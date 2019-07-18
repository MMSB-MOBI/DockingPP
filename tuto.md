
# Import modules


```python
import sys
sys.path.append("/Users/jprieto/DockingPP")
from dockingPP import parse, zParse
from core_scores import Scores, countNative, eval_natives
from core_clustering import rankCluster as rC, sortCluster, birchCluster
%load_ext autoreload
```

## Create DockData object


```python
# Megadock version
DD=parse("/Users/jprieto/Docking/data/unbound-outputs/4CPA_r-4CPA_l.detail")
DD.setReceptor("/Users/jprieto/Docking/data/benchmark5/structures/4CPA_r_u.pdb")
DD.setLigand("/Users/jprieto/Docking/data/benchmark5/structures/4CPA_l_u.pdb")

# zDock version with pre-treated pdbs using 'cut -c1-54'
zD=zParse("/Users/jprieto/Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/results/4CPA.zd3.0.2.fg.fixed.out")
zD.setReceptor("/Users/jprieto/Docking/ZD_new_pdbs/4CPA_r_u.pdb.ms_2")
zD.setLigand("/Users/jprieto/Docking/ZD_new_pdbs/4CPA_l_u.pdb.ms_2")
```

## Calculate Contact maps


```python
DD.ccmap(start=0,stop=500,pSize=50)
```

    Created 10 data packets (50 zObjects each) for process pool
    unpacking


## Calculate all scores and write frequences and scores to files


```python
DD.write_all_scores(filename="/Users/jprieto/docking/Resultats/mycomplex")
```

    Warning : only 500 poses could be analysed
    Warning : File /Users/jprieto/docking/Resultats/mycomplex_resstats.tab already exists, do you wish to continue anyway and replace it ? (yes/no)yes
    Warning : File /Users/jprieto/docking/Resultats/mycomplex_constats.tab already exists, do you wish to continue anyway and replace it ? (yes/no)yes
    Warning : File /Users/jprieto/docking/Resultats/mycomplex.tsv already exists, do you wish to continue anyway and replace it ? (yes/no)yes





    '/Users/jprieto/docking/Resultats/mycomplex.tsv'



## Parse scores with the Scores Class


```python
SC=Scores(filename="/Users/jprieto/docking/Resultats/mycomplex.tsv")
# or 
DD.setScores(filename="/Users/jprieto/docking/Resultats/mycomplex.tsv")
```

## Make clusters with the BSAS algorithm
*different ranks can be used*


```python
natural_rank= [i for i in range(50)]

SC.setPoses(DD.pList)
res_fr_rank=SC.rankedPoses(element="res_fr_sum")
con_fr_rank=SC.rankedPoses(element="con_fr_sum")

c_clusters=rC(DD,con_fr_rank,5, out='dict', stop=500)

```

## You can also use birch Algorithm for instance


```python
b_clusters=birchCluster(DD, 5)
# print(b_clusters)
```

## Sort clusters using Ranks and get representatives


```python

sor_clus=sortCluster(c_clusters,SC, fn="cons_score")
sor_bclus=sortCluster(b_clusters,SC, fn="cons_score")

# These are the final candidate poses for prediction.
rep=[c[0] for c in sor_clus]

brep=[c[0] for c in sor_bclus]

```


```python
print([p.id for p in rep])
print([p.id for p in brep])

a=brep[:]
a.extend(rep)
a=list(set(a))
for p in a :
    if p in rep and p in brep:
        print(p)
```

    [175, 315, 166, 81, 360, 450, 394, 29, 295, 455, 28, 420, 240, 11, 113, 426, 203, 148, 94, 109, 283, 324, 75, 194, 22, 40, 99, 180, 247, 382, 105, 453, 397, 95, 70, 464, 425, 372, 344, 476, 254, 427, 216, 257, 458, 227, 380, 356, 87, 150, 359, 333, 483]
    [387, 293, 143, 367, 326, 352, 457, 278, 451, 222, 173, 337, 166, 378, 13, 91, 80, 38, 484, 23, 74, 21, 209, 100, 469, 155, 92, 33, 8, 5, 156, 130, 4, 17, 129, 1, 42, 104, 3, 12, 2, 317, 71, 45, 365, 220, 70, 443, 48, 79, 211, 107, 171, 373, 349, 448, 476]
    70) (0.42, 2.66, -2.19) (2.4, -21.599999999999998, -3.5999999999999996)
    166) (-0.42, 1.83, -0.47) (-4.8, -12.0, -14.399999999999999)
    476) (-2.83, 1.53, -0.3) (22.8, 4.8, 19.2)


## Analyse the performance of each method


```python
# The MEGADOCK PARSER already picks up the RMSD of each decoy in the results file
# For zDock, you will have to set it manually for each decoy using 'p.set_RMSD(rmsd)'
# You
rmsds=[p.rmsd for p in rep]
brmsds=[p.rmsd for p in brep]
print(countNative(rmsds))
print(countNative(brmsds))


```

    {5: 0, 10: 1, 20: 1, 100: 2, 200: 2, 'out': 0}
    {5: 0, 10: 0, 20: 2, 100: 4, 200: 4, 'out': 0}


## Now use it on a set of complexes and count how many you got right


```python
Natives={}
Natives["mycomplex"]=countNative(rmsds)
print(eval_natives(Natives, 10))
```

    (['mycomplex'], [])
