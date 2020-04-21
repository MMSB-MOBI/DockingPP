import DockingPP.error as error
from typing import Tuple, Set, List
import logging
import math

class Pose :
    def __init__(self, index:int, euler:Tuple[float, float, float], translation: Tuple[float, float, float]):
        self.index : int = index
        self.euler : Tuple[float, float, float] = euler
        self.translation : Tuple[float, float, float] = translation
        self.residues_interface : Dict[str,Set[int]] = {"ligand" : set(), "receptor" : set()}
        self.contacts : Set[Tuple[int, int]] = set() #tuple(receptor_index, ligand_index)
        self.rescoring : Dict [str, float] = {}
        self.contact_computed : bool = False
        self.rmsd: int = None

    @property
    def nb_residues(self):
        return len([residues for role in self.residues_interface for residues in self.residues_interface[role]])

    def addContact(self, contact:Tuple[int, int]):
        """[summary]
        
        :param contact: [description]
        :type contact: Tuple[int, int]
        """
        self.contacts.add(contact)

    def addResidueAtInferface(self, role:str, residue_i:int):
        """[summary]
        
        :param role: [description]
        :type role: str
        :param residue_i: [description]
        :type residue_i: int
        :raises error.InvalidRole: [description]
        """

        try: 
            self.residues_interface[role].add(residue_i)
        except KeyError:
            raise error.InvalidRole(f"{role} is not a valid role.")

    def computeScore(self, type_score: str, name_score:str, score_fn):
        """[summary]
        
        :param score_fn: [description]
        :type score_fn: 
        :raises error.InvalidScore: [description]
        """

        if type_score == "contacts":
            score = score_fn(self.contacts)
        elif type_score == "residues":
            score = score_fn(self.residues_interface)
        else:
            raise error.InvalidScore(f"{type_score} is not valid. Choose contacts or residues")

        self.rescoring[name_score] = score

    def serializeScores(self, scores:List[str]) -> str:
        """[summary]
        
        :param scores: [description]
        :type scores: List[str]
        :raises error.InvalidScore: [description]
        :return: [description]
        :rtype: [type]
        """
        serialized = ""
        for sc in scores: 
            try:
                serialized += f"\t{self.rescoring[sc]}"
            except KeyError:
                raise error.InvalidScore(f"{sc} is invalid or not computed")
        
        return serialized

    def getScore(self, score):
        if not score in self.rescoring:
            raise error.InvalidScore(f"{score} is invalid or not computed")

        return self.rescoring[score]

    def setRMSD(self, rmsd:float):
        self.rmsd = float(rmsd)

    def isNative(self, rmsd_cutoff:float) -> bool:
        if self.rmsd <= rmsd_cutoff : 
            return True
        return False
