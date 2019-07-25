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

    def __str__(self):
        return str(self.poses) if self.poses else None

    def __repr__(self):
        return str([p.id for p in self.poses])

    def __getitem__(self, index):
        return self.poses[index]

    def __iter__(self):
        for i in self.poses :
            yield i


class ClusterColl(object):
    """ Must contain a dictionnary with integer keys """
    def __init__(self, clusters= None):
        if clusters:
            self.setClusters(clusters)
        else :
            self.clusters=None

    def addCluster(self, cluster):
        """ Takes a cluster object together in a cluster """
        self.clusters[max([key for key in self.clusters])+1]=cluster

    def setClusters(self,clusters):
        """ Build ClusterColl from a cluster's python dictionnary containing Docking poses
        and errase the previous clusters in collection"""
        self.clusters={ key : Cluster(clusters[key]) for key in clusters}

    def __str__(self):
        return str(self.clusters) if self.clusters else None

    def __repr__(self):
        return self.clusters

    def __getitem__(self, i ):
        return self.clusters[i]

    def __iter__(self):
        for key in self.clusters :
            yield key
    ###############################################
    ##                                           ##
    ##          Clustering Functions :           ##
    ##         - BSAS,                           ##
    ##         - Birch,                          ##
    ##         - ward,                           ##
    ##         - herarchical                     ##
    ##                                           ##
    ###############################################

def BSAS(rankedPoses, maxd, out="dict", start=0,stop=None) :
    """Can return list of pose's belonging if out is "list" or a dictionnary of clusters
    and the poses they contain if out is "dict"  """
    max=len(rankedPoses) if not stop else stop
    def calcul_dist(pose1,pose2):
        dist= sqrt(sum([(pose1.translate[i]-pose2.translate[i])**2 for i in range(3)]))
        return dist
    r_list=[p for p in rankedPoses]
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
        for cluster_id in clusters:
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
    elif out=="dict":
        return clusters

    else :
        raise Exception('out value can only be "list" or "dict"')

def birchCluster(zD, maxd, out='dict', N=None, start=0, stop=None ):
    #The radius of the subcluster obtained by merging a new sample and the closest subcluster should be lesser than the threshold.
    #Otherwise a new subcluster is started. Setting this value to be very low promotes splitting and vice-versa.
    data=zD.dictPos
    stop = len(zD.pList) if not stop else stop
    X =[[data['x'][i],data['y'][i],data['z'][i]] for i in range(start,stop)]
    brc = Birch(branching_factor=50, n_clusters=None, threshold=maxd, compute_labels=True)
    brc.fit(X)
    if N :
        brc.set_params(n_clusters = N )
    brc.partial_fit(np.matrix(X))
    groups=brc.predict(X)
    clusters=list2dict(zD,groups)
    if out=="list":
        return groups
    elif out=="dict":
        return clusters

def wardCluster(zD, maxd, start=0, stop=None ):
    stop = len(zD.pList) if not stop else stop
    data=zD.dictPos
    Z = ward(pdist([[data['x'][i],data['y'][i],data['z'][i],data['a1'][i],data['a2'][i],data['a3'][i]] for i in range(start,stop)]))
    groups = fcluster(Z,maxd, criterion = 'distance')
    clusters=list2dict(zD,groups)
    return clusters

def herarCluster(zD, maxc=None, linkage='complete', start=0, stop=None ):
    """linkage : {“ward”, “complete”, “average”, “single”} default is complete"""
    stop = len(zD.pList) if not stop else stop
    if maxc:
        cluster = AgglomerativeClustering(n_clusters=maxc, affinity='euclidean', linkage=linkage)
    else :
        cluster = AgglomerativeClustering(affinity='euclidean', linkage=linkage)
    data=zD.dictPos
    cluster.fit_predict([[data['x'][i],data['y'][i],data['z'][i]] for i in range(start,stop)])
    return list2dict(zD, cluster.labels_)




    ###############################################
    ##                                           ##
    ##        Clusters Utility Functions         ##
    ##                                           ##
    ###############################################


def clusScore(cluster,sc,rank='original_score'):
    if rank=='original_score':
        return sum([(p.id-1) for p in cluster])/len(cluster)
    else :
        try :
            rank=sc.ranks(element=rank)
            return sum([(rank[p.id-1]) for p in cluster])/len(cluster)
        except :
            raise Exception("Unknown rank type")


def sortCluster(cluster,sc,fn="original_score"):
    """Takes a clusters dictionary: {cluster1 : [ p1, p2, p3 ... ], cluster2 : [ p1, p2, p3 ... ]} """
    #Sort clusters using clusScore function. The lowest the score, the better the cluster.
    return sorted([cluster[c] for c in cluster], key=lambda o:clusScore(o,sc,fn))


def list2dict(zD,cluster):
    dictclus={}
    for u,v in enumerate(cluster):
        if v not in dictclus:
            dictclus[v]=[]
        dictclus[v].append(zD.pList[u])

    return dictclus
