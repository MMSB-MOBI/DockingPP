class Frequencies:
    def __init__(self, poses: List[DockingPP.pose.Pose]):
        self.nb_poses_used : int = len(poses)
        self.count_contact : Dict[Tuple[int, int], int]= {}
        self.count_ligand_residue : Dict[int, int]= {}
        self.count_receptor_residue : Dict[int, int] = {}
        self._computeCount(poses)

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
    def rel_frequencies_ligand_residue(self) -> Dict[int, float]:
        """[summary]
        
        :return: [description]
        :rtype: Dict[int, float]
        """
        if not hasattr(self, "_rel_frequencies_ligand_residue"):
            self._rel_frequencies_ligand_residue = { key : val/self.nb_poses_used for key, val in self.count_ligand_residue.items()}
        return self._rel_frequencies_ligand_residue

    @property
    def rel_frequencies_receptor_residue(self) -> Dict[int, float]:
        """[summary]
        
        :return: [description]
        :rtype: Dict[int, float]
        """
        if not hasattr(self, "_rel_frequencies_receptor_residue"):
            self._rel_frequencies_receptor_residue = { key : val/self.nb_poses_used for key, val in self.count_ligand_residue.items()}
        return self._rel_frequencies_receptor_residue


    def _computeCount(self, poses : List[DockingPP.pose.Pose]):
        """[summary]
        
        :param poses: [description]
        :type poses: List[DockingPP.pose.Pose]
        """
        for p in poses: 
            for contact in p.contacts: 
                if not contact in self.count_contact:
                    self.count_contact[contact] = 0
                self.count_contact[contact] += 1
            
            for ligand_residue in p.ligand_residues_interface:
                if not ligand_residue in self.count_ligand_residue:
                    self.count_ligand_residue[ligand_residue] = 0
                self.count_ligand_residue[ligand_residue] += 1

            for rec_residue in p.receptor_residues_interface:
                if not rec_residue in self.count_receptor_residue:
                    self.count_receptor_residue[rec_residue] = 0
                self.count_receptor_residue[rec_residue] += 1