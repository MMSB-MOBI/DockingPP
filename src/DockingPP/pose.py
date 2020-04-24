import DockingPP.error as error
from typing import Tuple, Set, List
import logging
import math

class Pose :
    """Object to handle single docking pose
    """
    def __init__(self, index:int, euler:Tuple[float, float, float], translation: Tuple[float, float, float]):
        self.index : int = index
        """pose index"""
        self.euler : Tuple[float, float, float] = euler
        """pose euler angles"""
        self.translation : Tuple[float, float, float] = translation
        """pose translation vector"""
        self.residues_interface : Dict[str,Set[int]] = {"ligand" : set(), "receptor" : set()}
        """Dictionary of residues index at the interface in ligand and receptor. It has molecule type as key and associated set of residues index as value. Set with calls to addResidueAtInterface."""
        self.contacts : Set[Tuple[int, int]] = set() #tuple(receptor_index, ligand_index)
        """Set of (i,j) tuples with i is the index of receptor residue and j is the index of ligand residue and i is in contact with j. Set with calls to addContact."""
        self.rescoring : Dict [str, float] = {}
        """Dictionary that stores rescoring for the pose. It has name of the score as key and score value as value. Set with calls to computeScore."""

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
