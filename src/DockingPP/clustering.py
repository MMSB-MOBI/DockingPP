from math import sqrt
from collections import OrderedDict 
from typing import List, Dict

def computeDistance(pose1: 'DockingPP.pose.Pose', pose2 : 'DockingPP.pose.Pose') -> float:
    """Compute distance between two poses. The distance is the distance between the center of mass of each pose, calculated from the translation vector.

    Args:
        pose1 ([DockingPP.pose.Pose])
        pose2 ([DockingPP.pose.Pose])

    Returns:
        float: The distance between the two poses
    """

    dist= sqrt(sum([(pose1.translation[i]-pose2.translation[i])**2 for i in range(3)]))
    return dist

def BSAS(poses : List['DockingPP.pose.Pose'], dist_cutoff: float) -> Dict['DockingP.pose.Pose', List['DockingPP.pose.Pose']]:
    """BSAS clustering from a list of poses. Browse the poses in the given order and add the pose to the first cluster where the distance between the pose and the representative pose of the cluster is below given distance cutoff. If the pose can't be added to any cluster, a new one is created.

    Args:
        poses (List[DockingPP.pose.Pose]): List of poses to cluster.
        dist_cutoff (float): Distance cutoff below which a pose is assigned to a cluster. 

    Returns:
        Dict[DockingPP.pose.Pose, List[DockingPP.pose.Pose]]: Dictionary that stores the clusters. It has representative pose object as key and list of other poses object in the cluster as value.
    """

    clusters = OrderedDict()
    for p in poses:
        in_cluster = False
        if not clusters: 
            clusters[p] = []
        else: 
            for representative in clusters:
                if computeDistance(representative, p) < dist_cutoff : 
                    in_cluster = True
                    clusters[representative].append(p)
                    break
            if not in_cluster:
                clusters[p] = []
            
    return clusters 


   