from math import sqrt
from collections import OrderedDict 
from typing import List, Dict

def computeDistance(pose1: 'DockingPP.pose.Pose', pose2 : 'DockingPP.pose.Pose') -> float:
    """Compute distance

    Args:
        pose1 ([type]): [description]
        pose2 ([type]): [description]

    Returns:
        float: [description]
    """

    dist= sqrt(sum([(pose1.translation[i]-pose2.translation[i])**2 for i in range(3)]))
    return dist

def BSAS(poses : List['DockingPP.pose.Pose'], dist_cutoff: float) -> Dict['DockingP.pose.Pose', List['DockingPP.pose.Pose']]:
    """[summary]

    Args:
        poses (List[DockingPP.pose.Pose]): [description]
        dist_cutoff (float): [description]

    Returns:
        Dict[DockingP.pose.Pose, List[DockingPP.pose.Pose]]: [description]
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


   