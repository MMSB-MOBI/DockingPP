
from typeguard import typechecked
from DockingPP.dockingHandler import DockingHandler
import DockingPP.typecheck as typecheck
import pypstruct.coordinates as PDB
from typing import Tuple, TypedDict, List, Optional, Dict
from DockingPP.pose import Pose
import DockingPP.error as error
import logging
from DockingPP.frequencies import Frequencies
import DockingPP.loader as loader
import DockingPP.clustering as clustering
import threading
import ccmap

PARSER_PDB = PDB.Parser()


class DockingHandlerPDB(DockingHandler):
    """A class that handles docking results directly from PDB files
    """

    def __init__(self,  step: float):
        super().__init__(step)
        self.step: float = step
        self.ligand:  'List[PDB.Structure]' = []
        """list of ligand load with PDB.Parser()"""
        self.receptor: 'List[PDB.Structure]' = []
        """list of receptor load with PDB.Parser()"""
        self.ligandfilename:  'List[str]' = []
        """list of ligand filename"""
        self.receptorfilename: 'List[str]' = []
        """list of receptor filename"""
        self.poses: List[Pose] = []
        """List of docking poses. Set with calls to addPose method."""

        # Raw contact map from ccmap
        self._raw_contact_map: List[List[int]] = None
        self.freq: Frequencies = None
        """Object that stores frequencies computations. Set with computeFrequencies method."""
        self._nb_rescored_poses: int = 0
        self._nb_cmap_poses: int = 0
        self.clusters: Dict[str, Dict[Pose, List[Pose]]] = {}
        """Poses clusters. Dictionary with score as key and dictionnary as value with representative pose as key and other poses that belongs to cluster as value. Set with clusterPoses."""

    def __str__(self):
        return f"#DockingHandler from PDB files object \nNumber of poses : {len(self.poses)}\n"

    def addLigand(self, ligand_pdb: str):
        """Set ligand and ligandfilename attributes of class with ligand pdb.

        Args:
            ligand_pdb (str): Path to ligand pdb
        """
        logging.info(f"== Set ligand ==\nfile: {ligand_pdb}")
        #add filename 
        self.ligandfilename.append(ligand_pdb)
        #add loaded structure of this pdb file
        self.ligand.append(PARSER_PDB.load(file=ligand_pdb))

    def addReceptor(self, receptor_pdb: str):
        """Set receptor and receptorfilename attributes of class with receptor pdb.

        Args:
            receptor_pdb (str): Path to receptor pdb
        """
        logging.info(f"== Set receptor ==\nfile: {receptor_pdb}")
        #add filename 
        self.receptorfilename.append(receptor_pdb)
        #add loaded structure of this pdb file
        self.receptor.append(PARSER_PDB.load(file=receptor_pdb))

    def addPose(self, pose_index: int):
        """Add pose (number of ligand/receptor)

        Args:
            pose_index (int): Index of the pose

        """
        p = Pose(pose_index)
        self.poses.append(p)

    def computeContactMap(self, nb_threads: int, nb_poses: int = -1, distance: float = 5):
        """Function that compute contact map for given poses and distance. It uses ccmap module, decode and store its results.

        Args:
            nb_threads (int): Number of threads to compute contact map
            nb_poses (int): Number of poses to compute contact map
            distance (float, optional): Distance (in Angstrom) below which two residues are considered in contact. Defaults to 5.

        Raises:
            error.IncompatiblePoseNumber: Raise if you want to compute on more poses than loaded.
        """
        
        if nb_poses == -1:
            nb_poses = len(self.poses)
            print( "Number of poses : ",nb_poses)
            
            
        if nb_poses > len(self.poses):
            raise error.IncompatiblePoseNumber(
                f"You try to compute contact map on {nb_poses} and only {len(self.poses)} are loaded")
        logging.info(
            f"== Compute contact map ==\nnumber of threads : {nb_threads}\nnumber of poses : {nb_poses}\ndistance : {distance}")

        if not self.ligand:
            raise error.PdbNotSet("Ligand is not set. ")

        if not self.receptor:
            raise error.PdbNotSet(
                "Receptor is not set. ")

        self._nb_cmap_poses = nb_poses
        output = [None for i in range(nb_threads)]
        threadPool = []

        nb_to_keep: int = nb_poses
        nb_split: int = nb_threads
        nWidth = int(nb_to_keep/nb_split)

        #Cut the ligand and receptor structure into different part for threading 
        for i in range(nb_split):
            top = (i+1) * nWidth
            if i == (nb_split - 1):
                top += nb_to_keep % nb_split

            #Get subset of structure for each thread
            reclist = list(map(lambda  x  : x.atomDictorize, self.ligand[i*nWidth : top ]))
            liglist = list(map(lambda  x  : x.atomDictorize, self.receptor[i*nWidth : top ]))
            
            threadPool.append(threading.Thread(
                target=self._lcmap_thread,
                args=(i, reclist, liglist, output, distance)))

        for th in threadPool:
            th.start()

        for th in threadPool:
            th.join()

        ccmap_result = [pose for thread in output for pose in thread]
        self._decodeContactMap(ccmap_result)

    def _lcmap_thread(self, thread_number: int, receptor : List[PDB.Structure], ligand: List[PDB.Structure], output: List[Optional[int]], distance: float):
        """Prepare a thread for ccmap execution.
      Args:
            thread_number (int): Number of the tread
            receptor (List): PDB receptor structure for this thread
            ligand (List):  PDB ligand structure for this thread
            output (List[Optional[int]]): Where the output will be stored
            distance (float): Distance for contact computation 

        """
        output[thread_number] = ccmap.lcmap(
            receptor,
            ligand,
            d=distance,
            encode=True)

        return

    def _decodeContactMap(self, ccmap_result: List[List[int]]):
        """Decode ccmap int results. For each pose, decode int into index pairs, and add contact and residues to Pose object.

        Args:
            ccmap_result (List[List[int]]): ccmap results, list of list of int. Each element of the list is the list of encoded contact pairs for one pose.
        """
        self._raw_contact_map = []
        for pose_index in range(len(ccmap_result)):
            ligand_residue_number = len(self.ligand[pose_index].getResID)
            pose_contact = ccmap_result[pose_index]
            pose_object = self.poses[pose_index]
            pose_object.contact_computed = True

             
            residues_index = [(int(i/ligand_residue_number), i %
                               ligand_residue_number) for i in pose_contact]
            
            for index in residues_index:
                pose_object.addContact(index)
                pose_object.addResidueAtInferface("ligand", index[1])
                pose_object.addResidueAtInferface("receptor", index[0])
            self._raw_contact_map.append(residues_index)

    def computeDockingScore(self):
        """Compute a Docking score for each air of ligand/receptor.

        """
                
        #Build a dictionnary with each contact betwen residues for all the pdb 
        contact_scores={}
        #Number of poses
        poses_number = self._nb_cmap_poses
        #return a list of pared pdb and a score
        res = []
        p = 0
        while p < poses_number:
            indexes = self._raw_contact_map[p]            
            lig_residues = self.ligand[p].getResID 
            rec_residues = self.receptor[p].getResID 
           
            for pair in indexes:
                rec_residue = rec_residues[pair[0]]
                lig_residue = lig_residues[pair[1]]
                my_key = str(rec_residue)+" "+str(lig_residue)
                if my_key in contact_scores.keys():
                    contact_scores[my_key] += 1
                else:
                    contact_scores[my_key] = 1
            p+=1


        i = 0
        while i < poses_number:
            final_score = 0
            rec_name = self.receptorfilename[i].split("/")[-1]
            lig_name = self.ligandfilename[i].split("/")[-1]

            lig_residues = self.ligand[i].getResID 
            rec_residues = self.receptor[i].getResID 
            
            #Get decoded contact
            indexes = self._raw_contact_map[i]

            # énumère les contacts et recherche les résidus cibles dans le récepteur
            for pair in indexes:
                rec_residue = rec_residues[pair[0]]
                lig_residue = lig_residues[pair[1]]
                my_key = str(rec_residue)+" "+str(lig_residue)
                final_score+=contact_scores[my_key]
                 
            res.append(rec_name+" "+lig_name+" "+str(final_score))
            i += 1
        
        return res
    
    def getContactMap(self, index : int):
        """Get the contactcard for a index number (start at 1)
        
        Args:
            index (int):  index number.
        """
        return self.poses[index-1].contacts



