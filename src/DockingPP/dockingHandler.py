from typeguard import typechecked
import DockingPP.typecheck as typecheck
import pypstruct.coordinates as PDB
from typing import Tuple, TypedDict, List, Optional, Dict
from DockingPP.pose import Pose
import DockingPP.error as error
import logging
from DockingPP.frequencies import Frequencies
import DockingPP.loader as loader
import DockingPP.clustering as clustering

PARSER_PDB = PDB.Parser()

# Type defs
PdbAtoms = TypedDict("PdbAtoms", {
    "x": List[float],
    "y": List[float],
    "z": List[float],
    "seqRes": List[str],
    "chainID": List[str],
    "resName": List[str],
    "name": List[str]
})


class DockingHandler:
    """A class that handles docking results
    """

    def __init__(self, step: float):

        self.step: float = step
        """spacing between grid cells"""

        self.ligand:  'pyproteinsExt.structure.coordinates.Structure' = None
        """ligand pdb parsed by pyproteinsExt. Set with setLigand method."""
        self.receptor: 'pyproteinsExt.structure.coordinates.Structure' = None
        """receptor pdb parsed by pyproteinsExt. Set with setReceptor method"""
        self.poses: List[Pose] = []
        """List of docking poses. Set with calls to addPose method."""

        """Corrected initial translation of ligand"""
        self._raw_contact_map: List[List[int]
                                    ] = None  # Raw contact map from ccmap
        self.freq: Frequencies = None
        """Object that stores frequencies computations. Set with computeFrequencies method."""
        self._nb_rescored_poses: int = 0
        self._nb_cmap_poses: int = 0
        self.clusters: Dict[str, Dict[Pose, List[Pose]]] = {}
        """Poses clusters. Dictionary with score as key and dictionnary as value with representative pose as key and other poses that belongs to cluster as value. Set with clusterPoses."""

    @property
    def cmap_poses(self):
        """Poses with contact map computed

        Returns:
            List[DockingPP.pose.Pose]: Poses with contact map computed
        """
        return self.poses[:self._nb_cmap_poses]

    @property
    def rescored_poses(self):
        """Rescored poses

        Returns:
            List[DockingPP.pose.Pose]: Poses with rescoring
        """
        return self.poses[:self._nb_rescored_poses]

    def computeFrequencies(self, nb_poses: int = -1):
        """Compute contact frequencies and residues at interface frequencies for given poses. Contact map has to be computed before.  

        Contact frequencies : For each contact between residue i of ligand and residue j of the receptor, count and relative frequency among poses are computed. 

        Frequencies of residues at interface : For each residue of the ligand and each residue of the receptor, count and relative frequency of the number of time the residue appears at interface among poses is computed. 

        Set freq attribute of the class. freq is a DockingPP.frequencies.Frequencies object. 

        Args:
            nb_poses (int): Number of poses to compute frequencies (optionnal)

        Raises:
            error.IncompatiblePoseNumber: Raise if you try to compute frequency with more poses than poses with contact map. 

        Examples: 
            Compute frequencies for 1BJ1 complex with 50 poses

            >>> DH.computeFrequencies(50)
            >>> DH.freq.rel_frequencies_contact
            {(312, 171): 0.5, (262, 170): 0.66, (313, 172): 0.5, (266, 5): 0.22, (320, 168): 0.64, (314, 161): 0.5, (259, 167): 0.34 ... }
            >>> DH.freq.rel_frequencies_residue
            {'ligand': {129: 0.74, 4: 0.7, 5: 0.68, 8: 0.74, 160: 0.5, 161: 0.5, 162: 0.74, 163: 0.54, 164: 0.72, 165: 0.7, 166: 0.66, 167: 0.72, 168: 0.7 ... }, 'receptor': {259: 0.48, 262: 1.0, 263: 0.66, 264: 1.0, 265: 0.7, 266: 1.0, 267: 0.66, 269: 0.62, 271: 0.94, 31: 0.9 ... }}


        """
        if nb_poses == -1:
            nb_poses = self._nb_cmap_poses

        logging.info(f"== Compute frequencies ==\nNumber of poses: {nb_poses}")

        if not self._raw_contact_map:
            raise error.ContactMapNotComputed(
                "Contact map doesn't exist. Call computeContactMap first.")

        if nb_poses > self._nb_cmap_poses:
            raise error.IncompatiblePoseNumber(
                f"You try to compute frequencies for {nb_poses} but only {self._nb_cmap_poses} have contact map.")

        self.freq = Frequencies(self.cmap_poses[:nb_poses])

    def rescorePoses(self, nb_poses: int, type_score: str):
        """Rescore N poses according to given type_score. A new score will be computed for each pose

        Args:
            nb_poses (int): Number of poses to rescore
            type_score (str): Score to use

                Available type score
                    * CONSRANK_U : sum of relative frequencies of each contact of the pose
                    * CONSRANK : CONSRANK_U normalised by number of contacts in the pose 
                    * contact_log_sum : sum of log of relative frequencies of each contact of the pose
                    * contact_square_sum : sum of square of relative frequencies of each contact of the pose
                    * residue_sum : sum of relative frequencies of each interface residue (ligand and receptor) of the pose
                    * residue_sum_ligand : sum of relative frequencies for ligand interface residues
                    * residue_sum_receptor : sum of relative frequencies for receptor interface residues
                    * residue_average : residues_sum normalised by number of interface residues in the pose
                    * residue_average_ligand : residues_sum normalised by number of interface residues in the pose, just for ligand residues
                    * residue_average_receptor : residues_sum normalised by number of interface residues in the pose, just for ligand receptor
                    * residue_log_sum : sum of log of relative frequencies of each interface residue of the pose
                    * residue_square_sum : sum of square of relative frequencies of each interface residue of the pose
                    * all : to compute all above scores

        Raises:
            error.IncompatiblePoseNumber: Raise if you try to rescore more poses than poses with contact map
            error.InvalidScore: Raise if you give an invalid type score. 

        Examples:
            For 1BJ1, rescore 2000 poses with all scores and display the first 2

            >>> DH.rescorePoses(2000, type_score = "all")
            >>> for p in DH.poses[:2]:
            >>>     print(p.index, p.rescoring)
            1 {'CONSRANK_U': 38.180000000000014, 'CONSRANK': 0.4061702127659576, 'contact_log_sum': -94.55436888131206, 'contact_square_sum': 17.90280000000001, 'residue_sum': 42.06, 'residue_average': 0.7378947368421053, 'residue_log_sum': -20.625622826512405, 'residues_square_sum': 33.8476}
            2 {'CONSRANK_U': 38.00000000000001, 'CONSRANK': 0.44186046511627913, 'contact_log_sum': -75.18567402702058, 'contact_square_sum': 18.285600000000002, 'residue_sum': 40.12, 'residue_average': 0.7569811320754717, 'residue_log_sum': -17.346104879674645, 'residues_square_sum': 32.776}

        """
        logging.info(
            f'== Rescore poses ==\nNumber of poses : {nb_poses}\n Scores : {type_score}')

        if nb_poses > self._nb_cmap_poses:
            raise error.IncompatiblePoseNumber(
                f"Impossible to rescore {nb_poses} poses, only {self._nb_cmap_poses} have contact map")

        if not self.freq:
            raise error.FrequenciesNotComputed(
                "Frequencies doesn't exist. Call computeFrequencies first.")

        self._nb_rescored_poses = nb_poses

        if type_score == "all":
            scores_to_compute = self.freq.available_scores

        else:
            if not type_score in self.freq.available_scores:
                raise error.InvalidScore(f"{type_score} score is not valid")
            scores_to_compute = {
                type_score: self.freq.available_scores[type_score]}

        for pose in self.rescored_poses:
            for score, score_info in scores_to_compute.items():
                pose.computeScore(score, *score_info)

    def _split_poses(self, nb_to_keep: int, nb_split: int) -> Tuple[int, List['DockingPP.pose.Pose']]:
        """Generator to split N first poses in X packets. 

        Args:
            nb_to_keep (int): N, the number of poses to split
            nb_split (int): X, the number of packets to create

        Yields:
            Tuple[int, List[DockingPP.pose.Pose]]: Tuple with index of the packet, and list of poses in the packet. 

        """
        current_poses = self.poses[:nb_to_keep]
        assert(nb_split <= nb_to_keep)
        nWidth = int(nb_to_keep/nb_split)
        for i in range(nb_split):
            top = (i+1) * nWidth
            if i == (nb_split - 1):
                top += nb_to_keep % nb_split
            yield(i, current_poses[i*nWidth:top])

    def serializeRescoring(self, output_file: str, type_score: List[str] = ["residue_sum", "residue_average", "residue_log_sum", "residue_square_sum", "CONSRANK_U", "CONSRANK", "contact_log_sum", "contact_square_sum"]):
        """Write rescoring results in a file. 

        Args:
            output_file (str): File to write results. The first line is a comment that resumes the number of poses used. The second line is the header with pose index followed by scores as given in scores_to_write list. The next lines are the scores for each pose.
            scores_to_write (List[str], optional): List of scores to write. The scores will be write in the given order. Defaults to ["residue_sum", "residue_average", "residue_log_sum", "residue_square_sum", "CONSRANK_U", "CONSRANK", "contact_log_sum", "contact_square_sum"].

        """
        if type_score == "all":
            scores_to_write = self.freq.available_scores

        else:
            if not type_score in self.freq.available_scores:
                raise error.InvalidScore(f"{type_score} score is not valid")
            scores_to_write = {
                type_score: self.freq.available_scores[type_score]}

        if not self._nb_rescored_poses:
            raise error.RescoringNotComputed(
                "Poses has not been rescored, call rescorePoses")

        o = open(output_file, "w")
        o.write(
            f"#Rescoring of {self._nb_rescored_poses} poses with frequencies computed on {self.freq.nb_poses_used} poses.\n")
        o.write("#Pose\t" + "\t".join(scores_to_write) + "\n")
        for p in self.rescored_poses:
            o.write(f"{p.index}{p.serializeScores(scores_to_write)}\n")
        o.close()
        logging.info(f"Scores writes to {output_file}")

    def getRankedPoses(self, score: str, nb_poses: int = -1) -> List['DockingPP.pose.Pose']:
        """Get poses ranked by given score

        Args:
            score (str): Score to rank by
            nb_poses (int): Number of poses to get

        Raises:
            error.IncompatiblePoseNumber: Raises if you try to rank more poses than rescored poses.

        Returns:
            List[DockingPP.pose.Pose]: List of poses in decreasing score order
        """

        if nb_poses == -1:
            nb_poses = self._nb_cmap_poses

        if score != "original" and nb_poses > self._nb_rescored_poses:
            raise error.IncompatiblePoseNumber(
                f"Try to rank {nb_poses} but only {self._nb_rescored_poses} have been rescored.")

        if score == "original":
            return self.poses[:nb_poses]

        sorted_poses = sorted(
            self.poses[:nb_poses], key=lambda pose: pose.getScore(score), reverse=True)
        return sorted_poses

    # Don't compute this here ??
    '''def storeRMSD(self, rmsd_file: str): 
        rmsds = loader.loadRMSD(rmsd_file, len(self.poses))
        for i in range(len(rmsds)):
            self.poses[i].setRMSD(rmsds[i])

    def getNativePoses(self, topX : int, rmsd_cutoff: float, ranked_by: str) -> List['DockingPP.pose.Pose']:
        if ranked_by == "original":
            if topX > len(self.poses):
                raise error.IncompatiblePoseNumber(f"Try to get native poses from top {topX} poses but only {len(self.poses)} have been stored.")
            poses_to_proceed = self.poses
        else:
            poses_to_proceed = self.getRankedPoses(ranked_by, self._nb_rescored_poses)

        return [p for p in poses_to_proceed[:topX] if p.isNative(rmsd_cutoff)]'''

    def clusterPoses(self, ranked_by: str, dist_cutoff: float, nb_poses: int):
        """Cluster the poses according to BSAS clustering. Clusters are computed for poses ranked by the given score in descending order. 

        BSAS clustering : Each pose is assigned to the first cluster when the distance between the pose and the representative of the cluster is less than given cutoff. If pose can't be assigned to any cluster, it creates a new one.  

        Set cluster attribute of the class. cluster is an ordered dictionary with representative Pose object as key and other Pose objects that belong to the cluster as values. 

        Args:
            ranked_by (str): Score for ranking
            dist_cutoff (float): Distance below which a pose is assigned to a cluster
            nb_poses (int): Number of poses to cluster

        """
        logging.info(f"== Cluster poses ==\nRanked by : {ranked_by}")
        poses_to_cluster = self.getRankedPoses(ranked_by, nb_poses)

        raw_clusters = clustering.BSAS([p.index for p in poses_to_cluster], [
                                       p.translation for p in poses_to_cluster], dist_cutoff)

        self.clusters[ranked_by] = {}
        for rep, others in raw_clusters.items():
            self.clusters[ranked_by][self.getPose(rep[0])] = [
                self.getPose(p[0]) for p in others]

    def getRankedClusterRepresentatives(self, ranked_by: str) -> List[Tuple[int, 'DockingPP.pose.Pose']]:
        """Get clusters representatives in decreasing order of given score

        Raises:
            error.ClustersNotComputed : Raise when clusters are not computed
            error.InvalidScore: Raise if score is not computed or invalid

        Returns:
            List[Tuple[int, DockingPP.pose.Pose]]: List of tuple where first element is the cluster number and second is the representative pose. 

        Examples:
            For 1BJ1, get CONSRANK_U clusters representatives and display the first 2 indexes 

            >>> representatives = DH.getRankedClusterRepresentatives("CONSRANK_U")
            >>> for clust in representatives[:2]:
            >>>       print("cluster", clust[0], "representative index", clust[1].index)
            cluster 0 representative index 16
            cluster 1 representative index 74

        """

        if not self.clusters:
            raise error.ClustersNotComputed(
                "Clusters have not been computed. Call clusterPoses first")

        if not ranked_by in self.clusters:
            raise error.InvalidScore(
                f"{ranked_by} cluster is invalid or not computed")

        cluster_nb = 0
        list_rep = []
        for rep_pose in self.clusters[ranked_by]:
            list_rep.append((cluster_nb, rep_pose))
            cluster_nb += 1

        return list_rep

    

    def loadScores(self, score_file: str):
        """Load scores from a score file serialized by serializeRescoring. 
        Load scores for each pose in rescoring Pose attribute, and store the number of rescored poses.

        Args:
            score_file (str): Path to score file

        Raises:
            error.IncompatiblePoseNumber: Raises if you try to load scores for more poses than loaded.
        """
        typecheck.validFile(score_file)
        with open(score_file) as f:
            f.readline()  # discard first lines with comments
            header = f.readline().rstrip()  # second line is header line
            scores = header.split("\t")[1:]  # get score names from header
            nb_poses = 0
            for pose_line in f:
                nb_poses += 1
                if nb_poses > len(self.poses):
                    raise error.IncompatiblePoseNumber(
                        f"There are more poses in scores file than in loaded zdock results.")

                pose_index = int(pose_line.split("\t")[0])
                pose_scores = pose_line.rstrip().split("\t")[1:]

                for i in range(len(scores)):
                    self.poses[pose_index - 1].rescoring[scores[i]
                                                         ] = float(pose_scores[i])
        self._nb_rescored_poses = nb_poses

    def getPose(self, pose_index: int):
        return self.poses[pose_index - 1]
