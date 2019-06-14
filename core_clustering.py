#! /usr/bin/env python3
import sys ,pickle, time
from math import sqrt
from statistics import mean
import numpy as np

from sklearn.cluster import AgglomerativeClustering, Birch
from scipy.cluster.hierarchy import linkage, dendrogram , fcluster , ward
from scipy.spatial.distance import pdist

"""This script provides with clustering functions for zDock data. They use the coordinates of the poses'
barycenters in order to generate clusters with different algorithms. """

# from IPython import display

# path=["/Users/jprieto/Docking/modules/pyproteinsExt/src", "/Users/jprieto/Docking/scripts/megadock_module"]
# for way in path :
#     if way not in sys.path :
#         sys.path.append(way)
#
# from megadock import parse as megaparse
# from multiprocessing import Pool

# clustering validation
# from sklearn.metrics import silhouette_score

class Cluster(object):
    def __init__(self,poses):
        self.poses=poses

    @property
    def size(self):
        def testInsert(pair, new_val):
            if not pair :
                pair=[new_val,new_val]
                return pair
            else :
                if pair[0]>new_val:
                    pair[0]=new_val
                if pair[1]<new_val:
                    pair[1]=new_val
                return pair

        X=None
        Y=None
        Z=None
        for p in self.poses :
            x,y,z=p.translate
            X=testInsert(X,x)
            Y=testInsert(Y,y)
            Z=testInsert(Z,z)
        return (X,Y,Z)
    @property
    def repr(self):
        return self.poses[0]

# def reverseRankCluster(pList, ranked, maxd=5, out="list") :
#     """Returns list of pose's belonging if out is "list" or a dictionnary of clusters
#     and the poses they contain if out is "dict"  """
#     def calcul_dist(pose1,pose2):
#         dist= sqrt(sum([(pose1.translate[i]-pose2.translate[i])**2 for i in range(3)]))
#         return dist
#     assert len(ranked)==len(pList)
#     r_list=[pList[i] for i in ranked]
#     groups=[]
#     clusters_found=0
#     clusters={} #cluster_id : [poses]
#     for i in range(len(r_list)):
#     #     print(str(pose.id) + str(pose.translate))
#         in_cluster = False
#         for cluster_id in groups:
#             # For each cluster representative
#             representative = clusters[cluster_id][0]
#             if calcul_dist(r_list[i],representative) < maxd:
#                 clusters[cluster_id].append(r_list[i])
#                 in_cluster = True
#                 groups.insert(0,cluster_id)
#                 break
#         if not in_cluster:
#             clusters_found += 1
#             clusters[clusters_found] = [r_list[i]]
#             groups.insert(0,clusters_found)
#     if out=="list":
#         groups.reverse()
#         return groups
#     if out=="dict":
#         return clusters

def rankCluster(zD, ranked, maxd, out="list", stop=None) :
    """Can return list of pose's belonging if out is "list" or a dictionnary of clusters
    and the poses they contain if out is "dict"  """
    pList=zD.pList
    max=len(pList) if not stop else stop
    def calcul_dist(pose1,pose2):
        dist= sqrt(sum([(pose1.translate[i]-pose2.translate[i])**2 for i in range(3)]))
        return dist
    r_list=[pList[i] for i in ranked]
    groups=[]
    clusters_found=0
    clusters={} #cluster_id : [poses]
    for i in range(len(r_list)) :
        if i > max:
            break
    #     if i % 10000 ==0:
    #         print("Progress : " + str(i))
    # #     print(str(pose.id) + str(pose.translate))
        in_cluster = False
        for cluster_id in groups:
            # For each cluster representative
            representative = clusters[cluster_id][0]
            if calcul_dist(r_list[i],representative) < maxd:
                clusters[cluster_id].append(r_list[i])
                in_cluster = True
                groups.append(cluster_id)
                break
        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [r_list[i]]
            groups.append(clusters_found)

    if out=="list":
        return groups
    if out=="dict":
        return clusters

def cons_score(cluster,sc):
    rank=sc.rankPoses(element="con_fr_sum")
    #     print([(rank[p.id-1]) for p in cluster])
    return sum([(rank[p.id-1]) for p in cluster])/len(cluster)

def original_score(cluster,sc):
    return sum([(p.id) for p in cluster])/len(cluster)

def sortCluster(cluster,sc,fn="original_score"):
    _f={"original_score":original_score,"cons_score":cons_score}
    # sort clusters from bigger to smaller
    return sorted([cluster[c] for c in cluster], key=lambda o:_f[fn](o,sc))

def birchCluster(zD, maxd):
    data=zD.dictPos
    X =[[data['x'][i],data['y'][i],data['z'][i]] for i in range(7,len(zD.pList))]
    brc = Birch(branching_factor=50, n_clusters=None, threshold=maxd, compute_labels=True)
    #The radius of the subcluster obtained by merging a new sample and the closest subcluster should be lesser than the threshold.
    #Otherwise a new subcluster is started. Setting this value to be very low promotes splitting and vice-versa.
    brc.fit(X)
    # brc.set_params(n_clusters = 2500 )
    brc.partial_fit(np.matrix(X))
    groups=brc.predict(X)
    return groups

def wardCluster(zD, maxd):
    data=zD.dictPos
    Z = ward(pdist([[data['x'][i],data['y'][i],data['z'][i],data['a1'][i],data['a2'][i],data['a3'][i]] for i in range(7,len(zD.pList))]))
    groups = fcluster(Z,maxd, criterion = 'distance')
    return groups

def herarCluster(zD, maxc):
    cluster = AgglomerativeClustering(n_clusters=maxc, affinity='euclidean', linkage='complete')
    data=zD.dictPos
    cluster.fit_predict([[data['x'][i],data['y'][i],data['z'][i]] for i in range(7,len(zD.pList))])
    return cluster.labels_
