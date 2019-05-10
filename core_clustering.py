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


# f=sys.argv[1]
# zDock = pickle.load(open(f,'rb'))
# data= zDock.dictPos



def rankCluster(pList, maxd) :

    def calcul_dist(pose1,pose2):
        dist= sqrt(sum([(pose1.translate[i]-pose2.translate[i])**2 for i in range(3)]))
        return dist
    groups=[]
    clusters_found=0
    clusters={} #cluster_id : [poses]
    for i in range(7, len(pList)) :
    #     print(str(pose.id) + str(pose.translate))
        in_cluster = False
        for cluster_id in clusters.keys():
            # For each cluster representative
            representative = clusters[cluster_id][0]
            if calcul_dist(pList[i],representative) < maxd:
                clusters[cluster_id].append(pList[i])
                in_cluster = True
                groups.append(cluster_id)
                break
        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [pList[i]]
            groups.append(clusters_found)
    return clusters

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
