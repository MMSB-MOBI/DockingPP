from typing import List, Dict, Tuple, Set
import math

class Frequencies:
    """An object that stores frequencies computed for a list of poses
    """

    def __init__(self, poses: List['DockingPP.pose.Pose']):
        self.nb_poses_used : int = len(poses)
        """The number of poses used for frequencies computation"""
        self.count_contact : Dict[Tuple[int, int], int]= {}
        """Dictionary that stores raw counts of contacts in poses. 
        """
        self.count_residue : Dict[str, Dict[int,int]]= {"ligand": {}, "receptor": {}}
        """Dictionary that stores raw counts of residues at interface in poses. 
        """
        self.available_scores = { "CONSRANK_U" : ("contacts", self.getContactFrequenciesSum),\
                                  "CONSRANK" : ("contacts", self.getContactFrequenciesAverage),\
                                  "contact_log_sum" : ("contacts", self.getContactFrequenciesLogSum),\
                                  "contact_square_sum" : ("contacts", self.getContactFrequenciesSquareSum),\
                                  "residue_sum" : ("residues", self.getResidueFrequenciesSum),
                                  "residue_average": ("residues", self.getResidueFrequenciesAverage),\
                                  "residue_log_sum": ("residues", self.getResidueFrequenciesLogSum),\
                                  "residue_square_sum": ("residues", self.getResidueFrequenciesSquareSum)}
        """A dictionary that stores available scores for rescoring. These scores are computed for a given list of contacts or a given list of residues at interface from contacts relative frequencies or relative frequencies of residues at interface  calculated here. The available scores are stored with their name, their type (if we need frequencies of contacts or residues at interface to compute them), and the function to compute them.
        """
        
        self._computeCount(poses)

    @property
    def rel_frequencies_contact(self) -> Dict[Tuple[int, int], float]:
        """The relative frequencies of contacts in poses.

        Returns:
            Dict[Tuple[int, int], float]: Dictionary with tuple (i,j) as key (where i is the index of residue in receptor and j is the index of residue in ligand) and score as value.
        """
        if not hasattr(self, "_rel_frequencies_contact"):
            self._rel_frequencies_contact = { key : val/self.nb_poses_used for key, val in self.count_contact.items()}
        return self._rel_frequencies_contact
    
    @property
    def rel_frequencies_residue(self) -> Dict[str, Dict[int, float]]:
        """The relative frequencies of residues at interface in poses.

        Returns:
            Dict[str, Dict[int, float]]: Dictionary with molecule role as key (ligand or receptor) and dictionary of frequencies for the molecule as value. Dictionary of frequencies has residue index as key and score as value.
        """
        if not hasattr(self, "_rel_frequencies_residue"):
            self._rel_frequencies_residue = {"ligand":{}, "receptor": {}}
            for role in self.count_residue:
                for idx in self.count_residue[role]:
                    self._rel_frequencies_residue[role][idx] = self.count_residue[role][idx] / self.nb_poses_used
        return self._rel_frequencies_residue

    def _computeCount(self, poses : List['DockingPP.pose.Pose']):
        """Compute raw count for contact and residues at interface in the given list of poses. Set count_contact and count_residue attributes.

        Args:
            poses (List[): List of poses to compute frequencies.
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


    #Maybe find a way to not have a single function for each score ???? 

    def getResidueFrequenciesSum(self, list_residue:Dict[str,Set[int]]) -> float:
        """For a given ensemble of residues, compute their relative frequencies sum. 

        Args:
            list_residue (Dict[str, Set[int]]): Ensemble of residues to compute frequencies sum. It's a dictionary with same format as DockingPP.pose.Pose residues_interface attribute (dictionary with molecule type as key and set of residues index as value). 

        Returns:
            float: The sum of residues relative frequencies.
        """
        return sum([self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used) for role in list_residue for idx in list_residue[role]])
    
    def getResidueFrequenciesLogSum(self, list_residue:Dict[str,Set[int]]) -> float:
        """For a given ensemble of residues, compute their relative frequencies log sum.

        Args:
            list_residue (Dict[str,Set[int]]): Ensemble of residues to compute frequencies log sum. It's a dictionary with same format as DockingPP.pose.Pose residues_interface attribute (dictionary with molecule type as key and set of residues index as value). 

        Returns:
            float: The log sum of residues relative frequencies.
        """
        return sum([math.log(self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used)) for role in list_residue for idx in list_residue[role]])

    def getResidueFrequenciesSquareSum(self, list_residue:Dict[str, Set[int]]) -> float:
        """For a given ensemble of residues, compute their relative frequencies square sum.

        Args:
            list_residue (Dict[str, Set[int]]): Ensemble of residues to compute frequencies square sum. It's a dictionary with same format as DockingPP.pose.Pose residues_interface attribute (dictionary with molecule type as key and set of residues index as value). 

        Returns:
            float: The square sum of residues relative frequencies.
        """
        return sum([self.rel_frequencies_residue[role].get(idx, 1 / self.nb_poses_used)**2 for role in list_residue for idx in list_residue[role]])

    def getResidueFrequenciesAverage(self, list_residue:Dict[str, Set[int]]) -> float:
        """For a given ensemble of residues, compute their relative frequencies sum normalized by the number of residues.

        Args:
            list_residue (Dict[str, Set[int]]): Ensemble of residues to compute frequencies square sum. It's a dictionary with same format as DockingPP.pose.Pose residues_interface attribute (dictionary with molecule type as key and set of residues index as value). 

        Returns:
            float: The normalized sum of residues relative frequencies.
        """
        nb_residues = len([residues for role in list_residue for residues in list_residue[role]])
        return self.getResidueFrequenciesSum(list_residue) / nb_residues

    def getContactFrequenciesSum(self, list_contact: Set[Tuple[int, int]]) -> float: 
        """For a given ensemble of contacts, compute their relative frequencies sum.

        Args:
            list_contact (Set[Tuple[int, int]]): Set of contacts to compute frequencies sum. It's a set of tuple(i,j) with i the index of residue in receptor and j the index of residue in ligand.

        Returns:
            float: The sum of contacts relative frequencies.
        """
        return sum([self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used) for contact in list_contact])

    def getContactFrequenciesLogSum(self, list_contact:Dict[Tuple[int, int], int]) -> float: 
        """For a given ensemble of contacts, compute their relative frequencies log sum.

        Args:
            list_contact (Set[Tuple[int, int]]): Set of contacts to compute frequencies sum. It's a set of tuple(i,j) with i the index of residue in receptor and j the index of residue in ligand.

        Returns:
            float: The log sum of contacts relative frequencies.
        """
        return sum([math.log(self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used)) for contact in list_contact])

    def getContactFrequenciesSquareSum(self, list_contact:Dict[Tuple[int, int], int]) -> float: 
        """For a given ensemble of contacts, compute their relative frequencies square sum.

        Args:
            list_contact (Set[Tuple[int, int]]): Set of contacts to compute frequencies sum. It's a set of tuple(i,j) with i the index of residue in receptor and j the index of residue in ligand.

        Returns:
            float: The square sum of contacts relative frequencies.
        """
        return sum([self.rel_frequencies_contact.get(contact, 1 / self.nb_poses_used)**2 for contact in list_contact])
        
    def getContactFrequenciesAverage(self, list_contact:Dict[Tuple[int, int], int]) -> float:
        """For a given ensemble of contacts, compute their relative frequencies sum normalized by the number of contacts.

        Args:
            list_contact (Set[Tuple[int, int]]): Set of contacts to compute frequencies sum. It's a set of tuple(i,j) with i the index of residue in receptor and j the index of residue in ligand.

        Returns:
            float: The normalized sum of contacts relative frequencies.
        """
        return self.getContactFrequenciesSum(list_contact) / len(list_contact)