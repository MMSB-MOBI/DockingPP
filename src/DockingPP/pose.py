import DockingPP.error as error
from typing import Tuple
import logging
import math

class Pose :
    def __init__(self, index:int, euler:Tuple[float, float, float], translation: Tuple[float, float, float]):
        self.index = index
        self.euler = euler
        self.translation = translation
        self.ligand_residues_interface = set()
        self.receptor_residues_interface = set()
        self.contacts = set() #tuple(receptor_index, ligand_index)
        self.rescoring = {}
        self.contact_computed = False

    def addContact(self, contact:Tuple[int, int]):
        self.contacts.add(contact)

    def addResidueAtInferface(self, role:str, residue_i:int):
        if role == "receptor":
            self.receptor_residues_interface.add(residue_i)
        elif role == "ligand":
            self.ligand_residues_interface.add(residue_i)
        else:
            raise error.InvalidRole("role must be ligand | receptor")

    def computeScore(self, type_score:str, freq):
        if type_score == "contacts_sum":
            score = sum([freq.rel_frequencies_contact.get(c, 1/freq.nb_poses_used) for c in self.contacts])

        elif type_score == "contacts_average":
            score = sum([freq.rel_frequencies_contact.get(c, 1/freq.nb_poses_used) for c in self.contacts]) / len(self.contacts)

        elif type_score == "contacts_log_sum":
            score = sum([math.log(freq.rel_frequencies_contact.get(c, 1/freq.nb_poses_used)) for c in self.contacts])

        elif type_score == "contacts_square_sum":
            score = sum([freq.rel_frequencies_contact.get(c, 1/freq.nb_poses_used)**2 for c in self.contacts])
        
        else:
            raise error.InvalidScore(f"{type_score} is invalid")

        self.rescoring[type_score] = score

    def getScore(self):
        return