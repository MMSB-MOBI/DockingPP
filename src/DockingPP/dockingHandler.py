from typeguard import typechecked
import DockingPP.typecheck as typecheck
import pyproteinsExt.structure.coordinates as PDB
from typing import Tuple, TypedDict, List, Optional
from DockingPP.pose import Pose
import DockingPP.error as error
import logging
import DockingPP.frequencies as frequencies

import threading
import ccmap

PARSER_PDB = PDB.Parser()

#Type defs 
PdbAtoms = TypedDict("PdbAtoms", {
    "x" : List[float],
    "y" : List[float],
    "z" : List[float],
    "seqRes" : List[str],
    "chainID" : List[str],
    "resName" : List[str],
     "name" : List[str]
})

@typechecked
class DockingHandler: 
    def __init__(self, grid_dimension, step, initial_euler, baryRec, baryLig):
       self.grid_dimension : int = grid_dimension
       self.step : float = step 
       self.initial_euler : Tuple[float, float, float] = initial_euler
       self.baryRec : Tuple[float, float, float] = baryRec
       self.baryLig : Tuple[float, float, float] = baryLig
       self.ligand : PdbAtoms = None 
       self.receptor : PdbAtoms = None 
       self.poses : List[Pose] = []
       self.offsetRec : Tuple[float, float, float] = tuple( [ -1 * bary for bary in self.baryRec]) 
       self.offsetLig : Tuple[float, float, float] = tuple( [ -1 * bary for bary in self.baryLig])
       self._raw_contact_map : List[List[int]] = None #Raw contact map from ccmap
       self._cmap_poses = None
       self.freq = None

    @property
    def cmap_poses(self):
        if self._cmap_poses == None : 
            self._cmap_poses = [p for p in self.poses if p.contact_computed]
        return self._cmap_poses

    def setLigand(self, ligand_pdb:str):
        self.ligand = PARSER_PDB.load(file = ligand_pdb)
    
    def setReceptor(self, receptor_pdb:str):
        self.receptor = PARSER_PDB.load(file = receptor_pdb)

    def computeContactMap(self, nb_threads:int, nb_poses:int, distance:float = 5):
        if nb_poses > len(self.poses):
            raise error.IncompatiblePoseNumber(f"You try to compute contact map on {nb_poses} and only {len(self.poses)} are loaded")
        logging.info(f"== Compute contact map ==\nnumber of threads : {nb_threads}\nnumber of poses : {nb_poses}\ndistance : {distance}")
        
        if not self.ligand:
            logging.error("Ligand is not set. Call setLigand first.")
            return
        
        if not self.receptor:
            logging.error("Receptor is not set. Call setReceptor first.")
            return
    
        output = [ None for i in range(nb_threads) ]
        threadPool = []
        for i, poses in self._split_poses(nb_poses, nb_threads):
            threadPool.append(threading.Thread(target = self._ccmap_thread, args = ([p.euler for p in poses], [p.translation for p in poses], i, output, distance)))

        for th in threadPool:
            th.start()

        for th in threadPool:
            th.join() 
        
        ccmap_result = [ pose for thread in output for pose in thread]
        self._decodeContactMap(ccmap_result)

    def _decodeContactMap(self, ccmap_result):
        self._raw_contact_map = []
        ligand_residue_number = self.ligand.residueNumber
        for pose_index in range(len(ccmap_result)): 
            pose_contact = ccmap_result[pose_index]
            pose_object = self.poses[pose_index]
            pose_object.contact_computed = True
            residues_index=[(int(i/ligand_residue_number), i % ligand_residue_number) for i in pose_contact]
            for index in residues_index:
                pose_object.addContact(index)
                pose_object.addResidueAtInferface("ligand", index[1])
                pose_object.addResidueAtInferface("receptor", index[0])
            self._raw_contact_map.append(residues_index)


    def computeFrequencies(self, nb_poses):
        if not self._raw_contact_map:
            logging.error("Contact map doesn't exist. Call computeContactMap first.")
            return 

        if nb_poses > len(self.cmap_poses):
            raise error.IncompatiblePoseNumber(f"You try to compute frequencies for {nb_poses} but only {len(self.cmap_poses)} have contact map.")

        self.freq = frequencies.Frequencies(self.cmap_poses[:nb_poses])

    def rescorePoses(self, type_score:str, nb_poses:int):
        if nb_poses > len(self.cmap_poses):
            raise error.IncompatiblePoseNumber(f"Impossible to rescore {nb_poses} poses, only {len(self.cmap_poses)} have contact map")
        if not self.freq:
            logging.error("Frequencies doesn't exist. Call computeFrequencies first.")
            return

        for pose in self.cmap_poses[:nb_poses]:
            pose.computeScore(type_score, self.freq)


    def _ccmap_thread(self, eulers, translations, thread_number:int, output:List[Optional[int]], distance:float):
        output[thread_number] = ccmap.lzmap(self.receptor.atomDictorize, self.ligand.atomDictorize, eulers, translations, d = distance, encode = True, offsetRec = self.offsetRec, offsetLig = self.offsetLig)
        return

    def _split_poses(self, nb_to_keep, nb_split):
        current_poses = self.poses[:nb_to_keep]
        assert(nb_split <= nb_to_keep)
        nWidth = int(nb_to_keep/nb_split)
        for i in range(nb_split):
            top = (i+1) * nWidth
            if i == (nb_split - 1):
                top += nb_to_keep%nb_split
            yield(i, current_poses[i*nWidth:top])


    def addPose(self, pose_index:int, euler:Tuple[float, float, float], translation:Tuple[float, float, float]):
        p = Pose(pose_index, euler, translation)
        self.poses.append(p)

    def __str__(self):
        return f"#DockingHandler object\nGrid dimension : {self.grid_dimension}\nStep : {self.step}\nInitial euler vector : {self.initial_euler}\nNumber of poses : {len(self.poses)}\nLigand offset : {self.offsetLig}\nReceptor offset : {self.offsetRec}"
