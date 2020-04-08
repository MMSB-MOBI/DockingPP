import DockingPP.types as types
from typeguard import typechecked
import DockingPP.typecheck as typecheck
import pyproteinsExt.structure.coordinates as PDB
from typing import Tuple, TypedDict, List
from DockingPP.pose import Pose
import DockingPP.error as error
import logging

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
    def __init__(self):
       self.grid_dimension : int = None
       self.step : float = None 
       self.initial_euler : Tuple[float, float, float] = None
       self.baryRec : Tuple[float, float, float] = None
       self.baryLig : Tuple[float, float, float] = None
       self.ligand_atoms : PdbAtoms = None 
       self.receptor_atoms : PdbAtoms = None 
       self.poses : List[Pose] = []

    def setLigand(self, ligand_pdb):
        self.ligand_atoms = PARSER_PDB.load(file = ligand_pdb).atomDictorize
    
    def setReceptor(self, receptor_pdb):
        self.receptor_atoms = PARSER_PDB.load(file = receptor_pdb).atomDictorize

    def computeContactMap(self, nb_threads:int, nb_poses:int):
        if nb_poses > len(self.poses):
            raise error.IncompatiblePoseNumber(f"You try to compute contact map on {nb_poses} and only {len(self.poses)} are loaded")
        logging.info(f"== Compute contact map ==\nnumber of threads : {nb_threads}\nnumber of poses : {nb_poses}")
        return

    def getContactMap(self):
        pass

    def addPose(self, pose_index:int, euler:Tuple[float, float, float], translation:Tuple[float, float, float]):
        p = Pose(pose_index, euler, translation)
        self.poses.append(p)