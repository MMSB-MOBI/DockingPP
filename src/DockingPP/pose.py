import DockingPP.error as error
from typing import Tuple, Set, List
import logging
import math

class Pose :
    """Object to handle single docking pose

    Attributes:
        index (int) : pose index
        euler (Tuple[float, float, float]) : pose euler angles
        translation (Tuple[float, float, float]) : pose translation vectors
        residues_interface (Dict[role, residues_index] where role is ligand|receptor and residues_index is Set[int]) : list of residues index at interface in ligand and receptor. Set with calls to addResidueAtInterface.
        contacts (Set[Tuple[int, int]]) : Set of (i,j) tuples with i index of receptor residue and j index of ligand residue and i is in contact with j. Set with calls to addContact
        rescoring (Dict[str, float]) : Dictionnary that store rescoring for the pose. Dictionnary with type of score as key and score value as value. Set with calls to computeScore.

    """
    def __init__(self, index:int, euler:Tuple[float, float, float], translation: Tuple[float, float, float]):
        self.index : int = index
        self.euler : Tuple[float, float, float] = euler
        self.translation : Tuple[float, float, float] = translation
        self.residues_interface : Dict[str,Set[int]] = {"ligand" : set(), "receptor" : set()}
        self.contacts : Set[Tuple[int, int]] = set() #tuple(receptor_index, ligand_index)
        self.rescoring : Dict [str, float] = {}

    def addContact(self, contact:Tuple[int, int]):
        """Add a contact to pose. 

        Args:
            contact (Tuple[int, int]): contact (i,j) between residue i of receptor and residue j of ligand
        """
        self.contacts.add(contact)

    def addResidueAtInferface(self, role:str, residue_i:int):
        """Add a residue at pose interface

        Args:
            role (str): To which residue belongs, ligand or receptor
            residue_i (int): Index of the residue

        Raises:
            error.InvalidRole: Raise if role is not valid, so is not ligand or receptor
        """

        try: 
            self.residues_interface[role].add(residue_i)
        except KeyError:
            raise error.InvalidRole(f"{role} is not a valid role.")

    def computeScore(self, type_score: str, name_score:str, score_fn):
        """Compute the given score with given function. Complete rescoring attribute.

        Args:
            type_score (str): Global type of the score, if it's derived from contacts frequencies or residues interface frequencies. Can be contacts|residues.
            name_score (str): Name of the score
            score_fn ([type]): Function to compute score

        Raises:
            error.InvalidScore: Raise if global type of the score is invalid
        """

        if type_score == "contacts":
            score = score_fn(self.contacts)
        elif type_score == "residues":
            score = score_fn(self.residues_interface)
        else:
            raise error.InvalidScore(f"{type_score} is invalid. Choose contacts or residues")

        self.rescoring[name_score] = score

    def serializeScores(self, scores:List[str]) -> str:
        """Str representation of rescoring

        Args:
            scores (List[str]): Scores for which you want the representation

        Raises:
            error.InvalidScore: Raises if a score is not computed or invalid.

        Returns:
            str: Str representation of scores, "\t" separated, in the order given in scores attributes.

        """
        serialized = ""
        for sc in scores: 
            try:
                serialized += f"\t{self.rescoring[sc]}"
            except KeyError:
                raise error.InvalidScore(f"{sc} score is invalid or not computed")
        
        return serialized

    def getScore(self, score) -> float:
        """Get a specific score

        Args:
            score ([type]): Score to get

        Raises:
            error.InvalidScore: Raise if score is not computed or invalid.

        Returns:
            float: score value
        """
        if not score in self.rescoring:
            raise error.InvalidScore(f"{score} score is invalid or not computed")

        return self.rescoring[score]

    def setRMSD(self, rmsd:float):
        self.rmsd = float(rmsd)

    def isNative(self, rmsd_cutoff:float) -> bool:
        if self.rmsd <= rmsd_cutoff : 
            return True
        return False
