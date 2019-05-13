#! /usr/bin/env python3

from dockingPP import parse,all_scores
from core_stats import writeScores

# set your data here :
pSize="size of packs for multiprocessing"
n="number of poses"
zF="DOCKING_FILE"
recPDB="RECEPTOR_PDB_FILE"
ligPDB="LIGAND_PDB_FILE"
resdir="results/"
name= "ComplexName"
c= "number of processes "

## Calculate scores and write Scores files
print(f"Parsing {zF} file")
zD= parse(zF, maxPose=n)
# print(n)
zD.setReceptor(recPDB)
zD.setLigand(ligPDB)
zD.setComplexName(name)
zD.ccmap(start=0,stop=int(n),dist=5, pSize=pSize, ncpu=c)

print("Calculating scores")
header = ["Surface size", "Residue freq sum", "Residue mean freq", "Residue log sum", "Residue square sum", "Number of contacts", "Contact freq sum", "Contact mean freq", "Contact log sum", "Contact square sum"]
fname=f"{resdir}/{zD.complexName}_{str(n)}.tsv"
print(fname)
writeScores(n-7, all_scores(zD, f"{resdir}/{zD.complexName}_{str(n)}" ), filename=fname, title=f"Scores for {zD.complexName}", header= header)

## VISUALISATION
from core_visuals import Scores

scores=Scores(fname)
scores.setPoses(zD.pList)
# generate a list of ranks for each pose : in position 0 we get pose[0]'s rank in rescoring
ranks=scores.rankPoses(element="con_log_sum", start=7)
# generate a ranked list of poses : in position 0 we get pose with best score
ranked=scores.rankedPoses(element="res_fr_sum", start=7, stop=8000)

# Generate 3D plot making visible rmsd and ranking
scores.plotly3D(ranks)

# Generate RMSD Graph
scores.rmsdGraphGenerator(ranked)

## Multiple Plots
wanted=["res_fr_sum","res_mean_fr","res_log_sum","con_fr_sum","con_mean_fr","con_log_sum" ] # list of wanted scores
complexes=["1BJ1","1KTZ","1GPW","1H9D"] # list of studied complexes
files= []"1BJ1_3000.tsv","1KTZ_3000.tsv","1GPW_10000.tsv","1H9D_8000.tsv"]

for i in range(len(complexes)):
    source=f"{complexes[i]}_r-{complexes[i]}_l.detail" #give access to zD file
    zD=parse(source,maxPose=3000)
    scores=Scores(files[i]) # give scores files generated with dockingPP
    scores.setPoses(zD.pList)
    plot=plt
    plt.title(complexes[i])
    axlist=multiplePlots(len(wanted)+1, size=(30,10))
    x=0
    for scor in wanted :
        ranked=scores.rankedPoses(element=scor, start=7,stop=2293) # get ranked poses
        scores.rmsdGraphGenerator(ranked, plot=axlist[x]) # generate rmsd plot
        axlist[x].set_title(scor) # set title
        x+=1
    scores.histRmsd(axlist[x])
    plt.title(complexes[i]+"_RMSDS")

## GENERATING CLUSTERS
from core_clustering import rankCluster, birchCluster, wardCluster,herarCluster
# Rank Clusters
# rankclusters can either be recovered in a list [clust1, cluster2, cluster1, cluster3, cluster4 ... ] or in a dictionnary {cluster1:[pose1, pose2], cluster2:[pose1,pose2] , ... }
# rankCluster can be performed over a chosen rank
ranked=scores.rankedPoses(element="con_mean_fr", start=7,stop=2293) # get ranked poses
clusters=rankCluster(zD.pList, ranked, 5, out="list")
# clusters=rankCluster(zD.pList, ranked, 5, out="dict")

# all other clusters methods return lists and take in DockingData objects:
groups2=birchCluster(zD,1.5)
groups3=wardCluster(zD,5)
groups4=herarCluster(zD,876)
