from math import sqrt
from collections import OrderedDict 
from typing import List, Dict, Tuple
import DockingPP.error as error

def computeDistance(transl1: Tuple[float, float, float], transl2 : Tuple[float, float, float]) -> float:
    """Compute distance between two poses. The distance is the distance between the center of mass of each pose, calculated from the translation vector.

    Args:
        transl1 (Tuple[float, float, float]) : The translation vector for pose 1
        transl2 (Tuple[float, float, float]) : The translation vector for pose 2

    Returns:
        float: The distance between the two poses
    """

    dist= sqrt(sum([(transl1[i]-transl2[i])**2 for i in range(3)]))
    return dist

def BSAS(poses_index:List[int], translations:List[Tuple[float, float, float]], dist_cutoff: float) -> Dict[Tuple[int, Tuple[float, float, float]], List[Tuple[int, Tuple[float, float, float]]]]:
    """
    BSAS clustering from a list of poses. Browse the poses in the given order and add the pose to the first cluster where the distance between the pose and the representative pose of the cluster is below given distance cutoff. If the pose can't be added to any cluster, a new one is created.

    Args:
        poses_index (List[int]): List of poses indexes
        poses_translation(List[Tuple[float, float, float]]) : List of poses translation vectors, in the order corresponding to poses_index
        dist_cutoff (float): Distance cutoff below which a pose is assigned to a cluster. 

    Returns:
        Dict[Tuple[int, Tuple[float, float, float]], List[Tuple[int, Tuple[float, float, float]]]]: Dictionary that stores the clusters. It has tuple with representative pose index and translation vector as key and list tuple with other poses index and translation vector that are present in the cluster as value.
        
    """

    if len(poses_index) != len(translations):
        raise error.InvalidArgument("poses_index and translations don't have the same size.")

    clusters = OrderedDict()

    for i in range(len(poses_index)):
        in_cluster = False
        idx = poses_index[i]
        transl = translations[i]
        if not clusters: 
            clusters[(idx, transl)] = []
        else: 
            for representative in clusters:
                if computeDistance(representative[1], transl) < dist_cutoff : 
                    in_cluster = True
                    clusters[representative].append((idx, transl))
                    break
            if not in_cluster:
                clusters[(idx, transl)] = []
            
    return clusters 


   