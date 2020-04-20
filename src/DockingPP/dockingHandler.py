from typeguard import typechecked
import DockingPP.typecheck as typecheck
import pyproteinsExt.structure.coordinates as PDB
from typing import Tuple, TypedDict, List, Optional
from DockingPP.pose import Pose
import DockingPP.error as error
import logging
from DockingPP.frequencies import Frequencies
import DockingPP.loader as loader

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

class DockingHandler: 
    def __init__(self, grid_dimension : int, step : float, initial_euler : Tuple[float, float, float], baryRec : Tuple[float, float, float], baryLig : Tuple[float, float, float]):
        """[summary]
        
        :param grid_dimension: [description]
        :type grid_dimension: int
        :param step: [description]
        :type step: float
        :param initial_euler: [description]
        :type initial_euler: Tuple[float, float, float]
        :param baryRec: [description]
        :type baryRec: Tuple[float, float, float]
        :param baryLig: [description]
        :type baryLig: Tuple[float, float, float]
        """
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
        self.freq : Frequencies = None
        self._nb_rescored_poses:int = 0
        self._nb_cmap_poses:int = 0

    @property
    def cmap_poses(self):
        """[summary]
        
        :return: [description]
        :rtype: [type]
        """
        return self.poses[:self._nb_cmap_poses]

    @property
    def rescored_poses(self):
        return self.poses[:self._nb_rescored_poses]

    def setLigand(self, ligand_pdb:str):
        """[summary]
        
        :param ligand_pdb: [description]
        :type ligand_pdb: str
        """
        self.ligand = PARSER_PDB.load(file = ligand_pdb)
    
    def setReceptor(self, receptor_pdb:str):
        """[summary]
        
        :param receptor_pdb: [description]
        :type receptor_pdb: str
        """
        self.receptor = PARSER_PDB.load(file = receptor_pdb)

    def computeContactMap(self, nb_threads:int, nb_poses:int, distance:float = 5):
        """[summary]
        
        :param nb_threads: [description]
        :type nb_threads: int
        :param nb_poses: [description]
        :type nb_poses: int
        :param distance: [description], defaults to 5
        :type distance: float, optional
        :raises error.IncompatiblePoseNumber: [description]
        """
        if nb_poses > len(self.poses):
            raise error.IncompatiblePoseNumber(f"You try to compute contact map on {nb_poses} and only {len(self.poses)} are loaded")
        logging.info(f"== Compute contact map ==\nnumber of threads : {nb_threads}\nnumber of poses : {nb_poses}\ndistance : {distance}")
        
        if not self.ligand:
            logging.error("Ligand is not set. Call setLigand first.")
            return
        
        if not self.receptor:
            logging.error("Receptor is not set. Call setReceptor first.")
            return
    
        self._nb_cmap_poses = nb_poses
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

    def _decodeContactMap(self, ccmap_result : List[List[int]]):
        """[summary]
        
        :param ccmap_result: [description]
        :type ccmap_result: List[List[int]]
        """
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


    def computeFrequencies(self, nb_poses:int):
        """[summary]
        
        :param nb_poses: [description]
        :type nb_poses: int
        :raises error.IncompatiblePoseNumber: [description]
        """
        logging.info(f"== Compute frequencies ==\nNumber of poses:{nb_poses}")

        if not self._raw_contact_map:
            logging.error("Contact map doesn't exist. Call computeContactMap first.")
            return 

        if nb_poses > self._nb_cmap_poses:
            raise error.IncompatiblePoseNumber(f"You try to compute frequencies for {nb_poses} but only {self._nb_cmap_poses} have contact map.")

        self.freq = Frequencies(self.cmap_poses[:nb_poses])

    def rescorePoses(self, type_score:str, nb_poses:int):
        """[summary]
        
        :param type_score: [description]
        :type type_score: str
        :param nb_poses: [description]
        :type nb_poses: int
        :raises error.IncompatiblePoseNumber: [description]
        """

        if nb_poses > self._nb_cmap_poses:
            raise error.IncompatiblePoseNumber(f"Impossible to rescore {nb_poses} poses, only {self._nb_cmap_poses} have contact map")
        if not self.freq:
            logging.error("Frequencies doesn't exist. Call computeFrequencies first.")
            return

        self._nb_rescored_poses = nb_poses
        for pose in self.rescored_poses:
            pose.computeScore(type_score, self.freq)


    def _ccmap_thread(self, eulers:List[Tuple[float, float, float]], translations: List[Tuple[float, float, float]], thread_number:int, output:List[Optional[int]], distance:float):
        """[summary]
        
        :param eulers: [description]
        :type eulers: List[Tuple[float, float, float]]
        :param translations: [description]
        :type translations: List[Tuple[float, float, float]]
        :param thread_number: [description]
        :type thread_number: int
        :param output: [description]
        :type output: List[Optional[int]]
        :param distance: [description]
        :type distance: float
        """
        output[thread_number] = ccmap.lzmap(self.receptor.atomDictorize, self.ligand.atomDictorize, eulers, translations, d = distance, encode = True, offsetRec = self.offsetRec, offsetLig = self.offsetLig)
        return

    def _split_poses(self, nb_to_keep:int, nb_split:int) -> Tuple[int, List['DockingPP.pose.Pose']]:
        """[summary]
        
        :param nb_to_keep: [description]
        :type nb_to_keep: int
        :param nb_split: [description]
        :type nb_split: int
        """
        current_poses = self.poses[:nb_to_keep]
        assert(nb_split <= nb_to_keep)
        nWidth = int(nb_to_keep/nb_split)
        for i in range(nb_split):
            top = (i+1) * nWidth
            if i == (nb_split - 1):
                top += nb_to_keep%nb_split
            yield(i, current_poses[i*nWidth:top])


    def addPose(self, pose_index:int, euler:Tuple[float, float, float], translation:Tuple[float, float, float]):
        """[summary]
        
        :param pose_index: [description]
        :type pose_index: int
        :param euler: [description]
        :type euler: Tuple[float, float, float]
        :param translation: [description]
        :type translation: Tuple[float, float, float]
        """
        p = Pose(pose_index, euler, translation)
        self.poses.append(p)

    def serializeRescoring(self, output_file:str, scores_to_write: List[str] = ["residues_sum", "residues_average", "residues_log_sum", "residues_square_sum", "contacts_sum", "contacts_average", "contacts_log_sum", "contacts_square_sum"]):
        """[summary]
        
        :param output_file: [description]
        :type output_file: str
        :param scores_to_write: [description], defaults to ["residues_sum", "residues_average", "residues_log_sum", "residues_square_sum", "contacts_sum", "contacts_average", "contacts_log_sum", "contacts_square_sum"]
        :type scores_to_write: List[str], optional
        """

        if not self._nb_rescored_poses : 
            logging.error("Poses has not been rescored, call rescorePoses")
            return

        o = open(output_file, "w")
        o.write(f"#Rescoring of {self._nb_rescored_poses} poses with frequencies computed on {self.freq.nb_poses_used} poses.\n")
        o.write("#Pose\t" + "\t".join(scores_to_write) + "\n")
        for p in self.rescored_poses:
            o.write(f"{p.index}\t{p.serializeScores(scores_to_write)}\n")
        o.close()
        logging.info(f"Scores writes to {output_file}")

    def getOriginalPoses(self, nb_poses):
        if nb_poses > len(self.poses):
           raise error.IncompatiblePoseNumber(f"Try to get {nb_poses} but only {len(self.poses)} have been load, including {self._nb_cmap_poses} with contact map.")

        return self.poses[:nb_poses]
        

    def getRankedPoses(self, score, nb_poses):
        if nb_poses > self._nb_rescored_poses:
            raise error.IncompatiblePoseNumber(f"Try to rank {nb_poses} but only {self._nb_rescored_poses} have been rescored.")

        sorted_poses = sorted(self.poses[:nb_poses], key=lambda pose:pose.getScore(score), reverse = True)
        return sorted_poses

    def storeRMSD(self, rmsd_file: str): 
        rmsds = loader.loadRMSD(rmsd_file, len(self.poses))
        for i in range(len(rmsds)):
            self.poses[i].setRMSD(rmsds[i])

    def getNativePoses(self, topX : int, rmsd_cutoff: float, ranked_by: str) -> List['DockingPP.pose.Pose']:
        if ranked_by == "original":
            if topX > len(self.poses):
                raise error.IncompatiblePoseNumber(f"Try to get native poses from top {topX} poses but only {len(self.poses)} have been stored.")
            poses_to_proceed = self.poses
        else:
            poses_to_proceed = self.getRankedPoses(ranked_by, self._nb_rescored_poses)

        return [p for p in poses_to_proceed[:topX] if p.isNative(rmsd_cutoff)]

    def __str__(self):
        return f"#DockingHandler object\nGrid dimension : {self.grid_dimension}\nStep : {self.step}\nInitial euler vector : {self.initial_euler}\nNumber of poses : {len(self.poses)}\nLigand offset : {self.offsetLig}\nReceptor offset : {self.offsetRec}"
