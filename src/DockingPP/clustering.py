from math import sqrt
from collections import OrderedDict 

def computeDistance(pose1, pose2):
    dist= sqrt(sum([(pose1.translation[i]-pose2.translation[i])**2 for i in range(3)]))
    return dist

def BSAS(poses, dist_cutoff) :
    """Can return list of pose's belonging if out is "list" or a dictionnary of clusters
    and the poses they contain if out is "dict"  """

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


   