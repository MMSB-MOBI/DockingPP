
# Import modules 


```python
import sys
sys.path.append("/Users/jprieto/DockingPP")
from dockingPP import parse, zParse
from core_visuals import Scores, countNative
from core_clustering import rankCluster as rC, sortCluster, birchCluster
%load_ext autoreload
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload


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
%autoreload 2
DD.write_all_scores(filename="/Users/jprieto/docking/Resultats/mytest")
```

    Warning : only 500 poses could be analysed
    Warning : File /Users/jprieto/docking/Resultats/mytest_resstats.tab already exists, do you wish to continue anyway and replace it ? (yes/no)yes
    Warning : File /Users/jprieto/docking/Resultats/mytest_constats.tab already exists, do you wish to continue anyway and replace it ? (yes/no)yes
    Warning : File /Users/jprieto/docking/Resultats/mytest.tsv already exists, do you wish to continue anyway and replace it ? (yes/no)yes





    '/Users/jprieto/docking/Resultats/mytest'



## Parse scores with the Scores Class


```python
SC=Scores("/Users/jprieto/docking/Resultats/mytest.tsv")
```

## Make clusters with the BSAS algorithm 
*different ranks can be used*


```python
natural_rank= [i for i in range(50)]

SC.setPoses(DD.pList)
res_fr_rank=SC.rankPoses(element="res_fr_sum")
con_fr_rank=SC.rankPoses(element="con_fr_sum")

c_clusters=rC(DD,con_fr_rank,8, out='dict', stop=500)

```

## You can also use birch Algorithm for instance 


```python
b_clusters=birchCluster(DD, 8)
# print(b_clusters)
b_dict={}
for u,c in enumerate(b_clusters): 
    if u == 500: 
        break
    if c not in b_dict: 
        b_dict[int(c)]=[] 
    b_dict[int(c)].append(DD.pList[int(u)])
# print(b_dict)
# print(sum([len(b_dict[c]) for c in b_dict]))
# print(len(SC.rankPoses(element="res_log_sum")))
# print(max([p.id-1 for c in b_dict for p in b_dict[c]]))
```

## Sort clusters using Ranks and get representatives


```python

sor_clus=sortCluster(c_clusters,SC, fn="cons_score")
sor_bclus=sortCluster(b_dict,SC, fn="cons_score")

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

    [166, 81, 394, 29, 148, 295, 203, 109, 426, 94, 28, 75, 283, 99, 180, 95, 464, 49, 259, 372, 458, 344, 254, 483, 227, 356, 257, 333, 150, 359]
    [387, 143, 352, 451, 80, 390, 38, 156, 209, 21, 469, 42, 33, 5, 56, 48, 46, 78, 1, 12, 3, 365, 220, 349, 476]


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

    {5: 0, 10: 1, 20: 1, 100: 1, 200: 1, 'out': 0}
    {5: 0, 10: 1, 20: 1, 100: 1, 200: 1, 'out': 0}


## Now use it on set of complexes and count how many you got right 


```python
def eval_natives(natives,n):
    positives=0
    good=[]
    bad=[]
    for c in natives: 
    #     print(natives[c])
        if natives[c][n]>0:
            positives+=1
            good.append(c)
        if natives[c][200]==0 and natives[c]["out"]==0:
            bad.append(c)
    return good, bad
print(eval_natives({"4CPB":countNative(rmsds)}, 10))
```

    (['4CPB'], [])

 