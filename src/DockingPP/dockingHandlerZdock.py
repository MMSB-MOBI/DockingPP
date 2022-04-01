import threading
from typeguard import typechecked
from DockingPP.dockingHandler import DockingHandler
from typing import Tuple, TypedDict, List, Optional, Dict
from DockingPP.pose import Pose
from DockingPP.frequencies import Frequencies
import logging
import ccmap
import pypstruct.coordinates as PDB

PARSER_PDB = PDB.Parser()

class DockingHandlerZdock(DockingHandler):
    """A class that handles docking results
    """

    def __init__(self, grid_dimension: int, step: float, initial_euler: Tuple[float, float, float], baryRec: Tuple[float, float, float], baryLig: Tuple[float, float, float]):
        self.grid_dimension: int = grid_dimension
        """size N of the NxNxN grid used in the docking"""
        self.step: float = step
        """spacing between grid cells"""
        self.initial_euler: Tuple[float, float, float] = initial_euler
        """initial rotation of the ligand in Euler angles"""
        self.baryRec: Tuple[float, float, float] = baryRec
        """receptor initial translation to center"""
        self.baryLig: Tuple[float, float, float] = baryLig
        """ligand initial translation to center"""
        self.ligand:  'pyproteinsExt.structure.coordinates.Structure' = None
        """ligand pdb parsed by pyproteinsExt. Set with setLigand method."""
        self.receptor: 'pyproteinsExt.structure.coordinates.Structure' = None
        """receptor pdb parsed by pyproteinsExt. Set with setReceptor method"""
        self.poses: List[Pose] = []
        """List of docking poses. Set with calls to addPose method."""
        self.offsetRec: Tuple[float, float, float] = tuple(
            [-1 * bary for bary in self.baryRec])
        """Corrected initial translation of receptor"""
        self.offsetLig: Tuple[float, float, float] = tuple(
            [-1 * bary for bary in self.baryLig])
        """Corrected initial translation of ligand"""
        self._raw_contact_map: List[List[int]
                                    ] = None  # Raw contact map from ccmap
        self.freq: Frequencies = None
        """Object that stores frequencies computations. Set with computeFrequencies method."""
        self._nb_rescored_poses: int = 0
        self._nb_cmap_poses: int = 0
        self.clusters: Dict[str, Dict[Pose, List[Pose]]] = {}
        """Poses clusters. Dictionary with score as key and dictionnary as value with representative pose as key and other poses that belongs to cluster as value. Set with clusterPoses."""

    def __str__(self):
        return f"#DockingHandler object\nGrid dimension : {self.grid_dimension}\nStep : {self.step}\nInitial euler vector : {self.initial_euler}\nNumber of poses : {len(self.poses)}\nLigand offset : {self.offsetLig}\nReceptor offset : {self.offsetRec}"

    def setLigand(self, ligand_pdb: str):
        """Set ligand attribute of class with ligand pdb. 

        Args:
            ligand_pdb (str): Path to ligand pdb
        """
        logging.info(f"== Set ligand ==\nfile: {ligand_pdb}")
        self.ligand = PARSER_PDB.load(file=ligand_pdb)

    def setReceptor(self, receptor_pdb: str):
        """Set receptor attribute of class with receptor pdb. 

        Args:
            receptor_pdb (str): Path to receptor pdb
        """
        logging.info(f"== Set receptor ==\nfile: {receptor_pdb}")
        self.receptor = PARSER_PDB.load(file=receptor_pdb)

    def addPose(self, pose_index: int, euler: Tuple[float, float, float], translation: Tuple[float, float, float]):
        """Add pose from zdock line informations. Used in loader.loadZdock

        Args:
            pose_index (int): Index of the pose
            euler (Tuple[float, float, float]): Euler vector of the pose
            translation (Tuple[float, float, float]): Translation vector of the pose

        """
        p = Pose(pose_index, euler, translation)
        self.poses.append(p)

    def computeContactMap(self, nb_threads: int, nb_poses: int, distance: float = 5):
        """Function that compute contact map for given poses and distance. It uses ccmap module, decode and store its results.  

        Args:
            nb_threads (int): Number of threads to compute contact map
            nb_poses (int): Number of poses to compute contact map
            distance (float, optional): Distance (in Angstrom) below which two residues are considered in contact. Defaults to 5.

        Raises:
            error.IncompatiblePoseNumber: Raise if you want to compute on more poses than loaded. 
        """
        if nb_poses > len(self.poses):
            raise error.IncompatiblePoseNumber(
                f"You try to compute contact map on {nb_poses} and only {len(self.poses)} are loaded")
        logging.info(
            f"== Compute contact map ==\nnumber of threads : {nb_threads}\nnumber of poses : {nb_poses}\ndistance : {distance}")

        if not self.ligand:
            raise error.PdbNotSet("Ligand is not set. Call setLigand first.")

        if not self.receptor:
            raise error.PdbNotSet(
                "Receptor is not set. Call setReceptor first.")

        self._nb_cmap_poses = nb_poses
        output = [None for i in range(nb_threads)]
        threadPool = []
        for i, poses in self._split_poses(nb_poses, nb_threads):
            threadPool.append(threading.Thread(target=self._ccmap_thread, args=(
                [p.euler for p in poses], [p.translation for p in poses], i, output, distance)))

        for th in threadPool:
            th.start()

        for th in threadPool:
            th.join()

        ccmap_result = [pose for thread in output for pose in thread]
        self._decodeContactMap(ccmap_result)

    def _ccmap_thread(self, eulers: List[Tuple[float, float, float]], translations: List[Tuple[float, float, float]], thread_number: int, output: List[Optional[int]], distance: float):
        """Prepare a thread for ccmap execution.

        Args:
            eulers (List[Tuple[float, float, float]]): List of euler vectors for the pose computed by the thread
            translations (List[Tuple[float, float, float]]): List of translation vectors for the pose computed by the thread
            thread_number (int): Number of the tread
            output (List[Optional[int]]): Where the output will be stored
            distance (float): Distance for contact computation 
        """
        output[thread_number] = ccmap.lzmap(self.receptor.atomDictorize, self.ligand.atomDictorize, eulers,
                                            translations, d=distance, encode=True, offsetRec=self.offsetRec, offsetLig=self.offsetLig)
        return
    
    def _decodeContactMap(self, ccmap_result : List[List[int]]):
        """Decode ccmap int results. For each pose, decode int into index pairs, and add contact and residues to Pose object. 

        Args:
            ccmap_result (List[List[int]]): ccmap results, list of list of int. Each element of the list is the list of encoded contact pairs for one pose.
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
    
    def getRankedClusters(self, ranked_by: str) -> List[Tuple[int, List['DockingPP.pose.Pose']]]:
        """Get clusters in decreasing order of given score

        Args:
            ranked_by (str): Score to rank by

        Raises:
            error.ClustersNotComputed: Raise if clusters are not computed
            error.InvalidScore : Raise if score is not computed or invalid

        Returns:
            List[Tuple[int, List[DockingPP.pose.Pose]]]: List of tuples where first element is cluster number and second is the list of poses inside the cluster.

        Examples:
            For 1BJ1, get CONSRANK_U clusters and display pose index.

            >>> clusters = DH.getRankedClusters("CONSRANK_U")
            >>> for clust in clusters[:2]:
            >>>     print("cluster :", clust[0], "poses :", [p.index for p in clust[1]])
            cluster : 0 poses : [16, 1, 5, ...]
            cluster : 1 poses : [74, 452, 1212, ...]

         """
        if not self.clusters:
            raise error.ClustersNotComputed(
                "Clusters have not been computed. Call clusterPoses first")

        if not ranked_by in self.clusters:
            raise error.InvalidScore(
                f"{ranked_by} cluster is invalid or not computed")

        nb_cluster = 0
        list_clusters = []
        for rep_pose in self.clusters[ranked_by]:
            list_clusters.append(
                (nb_cluster, [rep_pose] + self.clusters[ranked_by][rep_pose]))
            nb_cluster += 1

        return list_clusters
