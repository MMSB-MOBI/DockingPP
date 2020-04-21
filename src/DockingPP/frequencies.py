from typing import List, Dict, Tuple, Set
import math

class Frequencies:
    def __init__(self, poses: List['DockingPP.pose.Pose']):
        self.nb_poses_used : int = len(poses)
        self.count_contact : Dict[Tuple[int, int], int]= {}
        self.count_residue : Dict[str, Dict[int,int]]= {"ligand": {}, "receptor": {}}
        self._computeCount(poses)
        self.available_scores = { "contacts_sum" : ("contacts", self.getContactFrequenciesSum),\
                                  "contacts_average" : ("contacts", self.getContactFrequenciesAverage),\
                                  "contacts_log_sum" : ("contacts", self.getContactFrequenciesLogSum),\
                                  "contacts_square_sum" : ("contacts", self.getContactFrequenciesSquareSum),\
                                  "residues_sum" : ("residues", self.getResidueFrequenciesSum),
                                  "residues_average": ("residues", self.getResidueFrequenciesAverage),\
                                  "residues_log_sum": ("residues", self.getResidueFrequenciesLogSum),\
                                  "residues_square_sum": ("residues", self.getResidueFrequenciesSquareSum)}

    @property
    def rel_frequencies_contact(self) -> Dict[Tuple[int, int], float]:
        """[summary]
        
        :return: [description]
        :rtype: Dict[Tuple[int, int], float]
        """
        if not hasattr(self, "_rel_frequencies_contact"):
            self._rel_frequencies_contact = { key : val/self.nb_poses_used for key, val in self.count_contact.items()}
        return self._rel_frequencies_contact
    
    @property
    def rel_frequencies_residue(self) -> Dict[str, Dict[int, float]]:
        if not hasattr(self, "_rel_frequencies_residue"):
            self._rel_frequencies_residue = {"ligand":{}, "receptor": {}}
            for role in self.count_residue:
                for idx in self.count_residue[role]:
                    self._rel_frequencies_residue[role][idx] = self.count_residue[role][idx] / self.nb_poses_used
        return self._rel_frequencies_residue

    def _computeCount(self, poses : List['DockingPP.pose.Pose']):
        """[summary]
        
        :param poses: [description]
        :type poses: List[DockingPP.pose.Pose]
        """
        for p in poses: 
            for contact in p.contacts: 
                if not contact in self.count_contact:
                    self.count_contact[contact] = 0
                self.count_contact[contact] += 1
            
            for role in ["ligand", "receptor"]:
                for residue in p.residues_interface[role]:
                    if not residue in self.count_residue[role]:
                        self.count_residue[role][residue] = 0
                    self.count_residue[role][residue] += 1


    def getResidueFrequenciesSum(self, list_residue:Dict[str, Dict[int,Set[int]]]) -> float:
        """[summary]
        
        :param list_residue: [description]
        :type list_residue: Dict[int,Set[int]]
        :return: [description]
        :rtype: float
        """
        return sum([self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used) for role in list_residue for idx in list_residue[role]])
    
    def getResidueFrequenciesLogSum(self, list_residue:Dict[str, Dict[int,Set[int]]]) -> float:
        """[summary]
        
        :param list_residue: [description]
        :type list_residue: Dict[int,Set[int]]
        :return: [description]
        :rtype: float
        """
        return sum([math.log(self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used)) for role in list_residue for idx in list_residue[role]])

    def getResidueFrequenciesSquareSum(self, list_residue:Dict[str, Dict[int,Set[int]]]) -> float:
        """[summary]
        
        :param list_residue: [description]
        :type list_residue: Dict[int,Set[int]]
        :return: [description]
        :rtype: float
        """
        return sum([self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used)**2 for role in list_residue for idx in list_residue[role]])

    def getResidueFrequenciesAverage(self, list_residue:Dict[str, Dict[int,Set[int]]]) -> float:
        nb_residues = len([residues for role in list_residue for residues in list_residue[role]])
        return self.getResidueFrequenciesSum(list_residue) / nb_residues

    def getContactFrequenciesSum(self, list_contact:Dict[Tuple[int, int], int]) -> float: 
        """[summary]
        
        :param list_contact: [description]
        :type list_contact: Dict[Tuple[int, int], int]
        :return: [description]
        :rtype: float
        """
        return sum([self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used) for contact in list_contact])

    def getContactFrequenciesLogSum(self, list_contact:Dict[Tuple[int, int], int]) -> float: 
        """[summary]
        
        :param list_contact: [description]
        :type list_contact: Dict[Tuple[int, int], int]
        :return: [description]
        :rtype: float
        """
        return sum([math.log(self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used)) for contact in list_contact])

    def getContactFrequenciesSquareSum(self, list_contact:Dict[Tuple[int, int], int]) -> float: 
        """[summary]
        
        :param list_contact: [description]
        :type list_contact: Dict[Tuple[int, int], int]
        :return: [description]
        :rtype: float
        """
        return sum([self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used)**2 for contact in list_contact])
        
    def getContactFrequenciesAverage(self, list_contact:Dict[Tuple[int, int], int]) -> float:
        return self.getContactFrequenciesSum(list_contact) / len(list_contact)