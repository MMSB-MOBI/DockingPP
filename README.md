### DOCKING POST-PROCESSING Utility Class and functions


# Table of Contents
  * [Dependencies](#chapter-1)
  * [Quick Start](#chapter-2)
  * [The megadock format](#chapter-3)
  * [The ccmap format](#chapter-6)
  * [Getting started](#chapter-4)
  * [Available objects and functions](#chapter-5)

---------------------------------
<a id="chapter-1"></a>

</br> 

## Dependencies 
dockingPP has the following dependencies : 

* **Python 3**
* pyProteinsExt (<https://github.com/glaunay/pyproteinsExt>)
* pyProteins (<https://github.com/glaunay/pyproteins>) 
* ccmap (<https://github.com/glaunay/ccmap>) 
* NumPy (<http://www.numpy.org/>)
* sklearn (<https://scikit-learn.org/stable/>) 
* scipy (<https://www.scipy.org/>)
* matplotlib (<https://matplotlib.org/>) 
* plotly (<https://plot.ly/python/>) 
* mpl_toolkits (<https://matplotlib.org/1.4.3/mpl\_toolkits/index.html>)


Install dependencies using 

```sh
pip install -r LOCAL_PATH_TO_REPO/requirements.txt
``` 

<a id="chapter-2"></a>

</br> 
## Quick Start 


```sh
git clone https://github.com/juliaprieto/dockingPP.git  LOCAL_PATH_TO_REPO
``` 
Check for dependencies in LOCAL\_PATH\_TO\_REPO/requirments.txt

```python
sys.path.append(LOCAL_PATH_TO_REPO/src)
``` 


<a id="chapter-3"></a>

----------------------------

</br> 

## The megadock format 

A typical megadock results file is derived from the [ZDOCK format](http://zdock.umassmed.edu/zdock_output_format.php). 

- First line features the number of grid cells and the size of the cubic cell. 
- Line 2 indicates the initial random rotation of the **ligand**. 
- Line 3 indicates the initial position of the **receptor** barycenter within the original system of coordinates. 
- Line 4 indicates similar information for the **ligand**.
- From line 12 to 16, docking poses' information is provided.   

<span style="color:Crimson ;font-weight=500">In order to reconstruct a docking pose both **receptor** and **ligand** molecule barycenters have to be set to the origin.</span>  

- The first field is the **Euler angles (X,Z,X')** triplet to apply rotation to the **ligand** molecule to reconstruct its orientation in the pose.   
- The second triplet is the **translation vector** in mesh units to move the **ligand** molecule to its position in  the pose.

<pre style="font-family:Monaco;font-size:0.8em;padding:1em;border:gray 1px solid; background-color: Gainsboro;">
 1 126     1.20
 2         0.000000        0.000000        0.000000
 3 1a2k_u1.pdb     27.372499       8.099500        82.649002
 4 1a2k_u2.pdb     59.376999       0.796000        75.828003
 5 
 6 Elec score ratio:   1.00
 7 ACE  score ratio:   1.00
 8 
 9 *** docking results ***
10
11 Rank   Angle               Trans         total       rPSC     rPSC(+)  rPSC(-)  ELEC       RecACE     LCARMSD NearNat
12     1: (-2.72, 1.87, 0.62) ( 19, 15, 12)     4376.78     2977     7522    -4545     315.14    1084.64  56.351       0
13     2: (-2.72, 1.74, 0.66) ( 19, 14, 13)     4288.76     3066     7566    -4500     140.73    1082.03  56.310       0
14     3: (-2.93, 1.98,-3.13) ( 27,  8,  7)     4288.50     3394     7579    -4185     -27.59     922.10  64.465       0
15     4: ( 1.68, 1.19, 2.71) (121,104,106)     4270.28     3384     5994    -2610      94.39     791.89  52.556       0
16     5: ( 2.41, 1.63, 1.84) (  7, 96,117)     4269.77     3563     4688    -1125     337.08     369.69  61.792       0
</pre>
So for each pose, the ccmap libray has to:   
1. move ligand and receptor to center  
2. rotate ligand   
3. translate ligand  
4. compute interaction matrix an return contact map 
 
<a id="chapter-6"></a>

</br> 

## The ccmap format 
The ccmap Library (<https://github.com/glaunay/ccmap>) reconstructs docking poses using euler angles and translation values and analyses the formed complex to identify the contacts at the interface, and returns what is called the contact map. 
#### Docking data processing 
This library can be used with data from MEGADOCK or data from ZDOCK.
If more than one rotation is applied to the ligand, the dockingPP library computes a set of Euler angles equivalent to both rotations. </br>
The ligand translations are cropped in order to fit the docking grid. 
If translation value (t\_value) is larger than half grid size (nCells / 2) , translation value becomes :  **t\_value - nCells** </br>
Then the translation values in grid cells unit are converted into Amstrong unit : 
**t = (-1) * t * step**

#### Ccmap reconstruction 

In order to analyse poses, ccmap reconstructs the pose by computing the position of every atom. It applies a rotation of the ligand using euler (x-z-x) angles and then translates every atom ligand (x) from initial barycenter  of the ligand (l<sub>0</sub>) to center of the grid : x = x- l<sub>0</sub> then adds translation x = x + t. 

#### The contact map format 

For the identification of contacts in a pose, ccmap discretizes space again in order to reduce complexity. It returns a list of integers with maximal size of **(number of receptor residues * number of ligand residues)** the maximal number of contacts. </br>
Those integers represent the none-zero positions in a sparse matrix with indexes (i,j) </br>.
The integers are obtained by multiplying the number of residues in ligand(i_max) by j (receptor's index in matrix) and adding i (ligand's index in matrix). This way, each integer represents a contact between the ligand residue of index i and receptor residue of index j. 

----------------------
<a id="chapter-4"></a>
</br> 
## Getting started with Docking Post-Processing 
#### Docking prediction Output reading 

```python
from dockingPP import parse
DD= parse('PATH_TO_FILE/1BJ1_r-1BJ1_l.detail', maxPose =45) # maxpose default is 0
```
*Be careful, since docking outputs vary and may even provide near native solutions for the study.*   
*For now, the program expects 7 near native solutions in the beginning of the file. These are not taken into account for residue or contact statistics. Therefore if maxpose=4000 you will get 7 near native solutions and 3993 computed poses*


#### Add reference PDB for receptor and ligand

```python
DD.setReceptor('PATH/1BJ1_r_u.pdb')
DD.setLigand('PATH/1BJ1_l_u.pdb')
```
Make sure you use unbound molecules to do your predictions, since the bound version of the structures is supposedly unknown. 
#### Add name to Docking object

```python
DD.setComplexName('1BJ1')
```

#### Compute interaction interface
This action will be the most ram consuming. It may take a while. 

This function uses the ccmap module    
 1* Each pose keeps a dictionary of its amino-acids contacts 10000 poses are treated in 4mn, with a memory cost ~ 9286M   
 2* Releasing the GIL w/in C extension does reduce total computation cost to 3mn40s. Bottleneck is single structure computation in C 

```python
DD.ccmap(start=0,stop=n,dist=5, pSize=50)
# start default : 0 
# stop default : len of DD.pList 
# dist default : 4.5
# pSize default : 200
```
#### Calculate new scores for the Docking results

```python
DD.all_scores()
```
#### ... Or write them to a file 

```python
DD.write_all_scores(filename="1BJ1_new_scores", title="Rescoring 1BJ1")
```

#### ... Or import them from a file 

```python
DD.setScores(filename="1BJ1_new_scores.tsv" )
```

*See more in the tuto file : <https://github.com/juliaprieto/DockingPP/blob/master/tuto.md>*
</br>

<a id="chapter-5"></a>
## Available objects and functions

* <a href=#docking style="text-decoration: none;"> dockingPP </a>
	* <a href=#general style="text-decoration: none;"> DockData Object's Basics</a>
	* <a href=#rescoring style="text-decoration: none;"> Rescoring Utilities</a>
	* <a href=#clustering style="text-decoration: none;"> Clustering Utilities</a>
* <a href=#stats style="text-decoration: none;"> Statistics, src/core_stats.py</a>
* <a href=#clustering style="text-decoration: none;"> Clustering Object, src/core_clustering.py </a>


</br>

<a id="docking"></a>
### dockingPP

<span style="color:Crimson ;font-weight=500"> *class DockData* </span> _ *Docking results container*
we will consider as pose_index its original rank 

* Attributes
	- **complexName** : can be set using method setComplexName
	- **step** : size of the docking grid cells
	- **nCells** : size of the docking grid
	- **eulerREC** : euler angles of the receptors' reference position
	- **fileREC** : name of the receptor's original pdb file 
	- **fileLIG** : name of the ligand's original pdb file 
	- **baryREC** : receptor's barycenter position
	- **baryLIG** : ligand's reference's barycenter
	- **pList** : list of poses in the experiment
	- **pdbObjLigand** : structure object built from PDB file
	- **pdbObjReceptor** : structure object built from PDB file
	- **getStats** : returns resStats and contactStats objects into a list
			`resStats,conStats=DD.getStats`
	- **contactStats** : returns contactStats object
	- **resStats** : returns resStats object
* Available Methods<a id="general"></a>

	#### General Functions
	
	- **setComplexName**
	- **ccmap(start=0, stop=1000, dist=5, pSize=100, ncpu=10)** : calculate and store ccmap for each pose using multiprocessing </br>
*Use start and stop to choose the number of poses to be computed, dist to set the minimal distance to detect a contact in the interface of the two molecules, pSize (size of the paquets) and ncpu (number of workers) as multiprocessing parameters.* 
	- **push(pose\_id, pose\_euler, pose\_tr)** : add pose to pList 
	- **setReceptor(RecFile)** : add PDB reference
	- **setLigand(LigFile)** : add PDB reference
	- **loadRMSD(filename="")** : load RMSDs from zDock like rmsd file (tsv :  [ pose id ; RMSDS ] ) 
	- **dictPos()** : returns a dictionary useful for visualisation { 'x' : [ ] , 'y' : [ ] , 'z' : [ ] } 
	- **all\_scores()**
	- **write\_all\_scores()**<a id="rescoring"></a>


	#### Rescoring and visualization Methods
	Choose scoring scores from :'original\_rank', 'r\_size', 'res\_fr\_sum', 'res\_mean\_fr', 'res\_log\_sum', 'res\_sq\_sum', 'c\_size', 'con\_fr\_sum', 'con\_mean\_fr', 'con\_log\_sum', 'con\_sq\_sum'")
		
	- **.rankedPoses(element=**"original_rank", **start**=0, **stop**=None ): add PDB reference
	- **.rankedIDs(element**="original_rank", **start**=0, **stop**=None ) : add PDB reference
	- **.ranks(element**="original_rank") : returns rank values from 1 to n 
	- **.rmsds()** : add PDB reference
	- **.rankedRmsds(element**="original_rank", **start**=0, **stop**= None) : add PDB reference
	- **.rankedPoses(element**="original_rank", **cutoff**=5 ) : add PDB reference
	- **.countNatives(element**="original_rank", **cutoff**=5 ) : add PDB reference
	- **.plot3D(element**="original_rank", **name**="",  title="") : add PDB reference
	- **.mutliPlot3D( wanted_scores**,  **title**='Docking decoys', **size**=(600,400)) : add PDB reference
	- **.rmsdPlot(element**="original_rank", **start**=0, **stop**=None, **plot**=None, **title**=None ) : add PDB reference<a id="clustering"></a>

	#### Clustering Methods
		
	- **.BSAS(maxd , element**="original_rank"**, out**="dict"**, start**=0**,stop**=None**): cluster experiment's poses with Basic sequential algorithm scheme lead by element order using maximal distance maxd between barycenters as inclusion condition. 
	- **. birchCluster(maxd, out**='dict'**, N**=None): use birch clustering algorithm with maxd being maximal threshold for cluster's radius and N the number of clusters to build.
	- **.wardCluster(maxd,  start**=0**, stop**=None) : use ward linkage clustering with maxd the maximum variance of the clusters.
	- **.herarCluster(self, maxc**=None**, linkage**='complete'**, start**=0**, stop**=None ) : apply herarchichal clustering, linkage can be {“ward”, “complete”, “average”, “single”}. maxc is the number of clusters to be kept . 
	
</br>
<span style="color:Crimson ;font-weight=500"> *class Pose* </span> _ *Ligand Pose information*. 
  
* Attributes
	- **.id** : integer, similar to its position in the original ranking
	- **.euler : (psi, theta, phi)** rotation angles from the original pdb ligand to the pose
	- **.translate : (x,y,z)** translation distances from the original pdb ligand to the pose
	- **.ccmap**: json formatted contact map, dictionary of contacts between residues   </br> 
format : {'type': 'contactList', 'data': [{'root': {'resID': '  50 ', 'chainID': 'A'}, 'partners': [{'resID': ' 650 ', 'chainID': 'B'}]}, ... ]} 
	- **.belongsTo** : returns DockData object containing the pose
	- **.dictorizedReceptor** : positions for the receptors atoms
	- **.dictorizedLigand** : positions for the ligands atoms
	- **.rmsd** : if set, returns the rmsd between the pose and the native pose
	- **.resMapList** : returns the list of residues implied in the contact map for this pose in the form of the concatenation of the residue number, its chain and its original molecule </br> 
*format* : ['171ALig', '175ALig' , ... ]
	- **.resSize** : number of residues in pose contactmap
	- **.conSize** : number of contacts in pose contactmap 
	- **.scores** : returns a dictionnary of poses' scores in data set if scores have been calculated

* Available Methods
	- **.has_ccmap(error=False)** : Returns boolean, checks wether pose's ccmap has been calculated
	- **.has_rmsd()** : checks wether pose's ccmap has been calculated
	- **.ccmap(dist=5)** : calculate contact map
	- **.dump()** : returns a pdb-like content describing the pose, ligand and receptor , rotations have been applied if pose's ccmap has been calculated

</br>

<a id="stats"></a>
### src.core_stats
<span style="color:Crimson ;font-weight=500"> *class CmapRes* </span> _ *store residues from ccmap and their counts* 

* Attributes 
	- **.resID** : residue id
	- **.chainID** : residue chain in the form of the original PDB name
	- **.role** : Residue origin : receptor or ligand ('Rec' and 'Lig') 
	- **.count** : number of occurences in residues counts
	- **.cCount** : number of occurences in contacts counts
	- **.index** : unique index for each residue ( concatenation of the residue number, its chain and its original molecule ) 
* Available Methods
	- **.increase_count(count='plain')** : add 1 to the pose count number in residue counts if *count='plain'*, in contact counts if *count='pond'*
	- **.reset_all()** : resets all counts to 0

</br>
<span style="color:Crimson ;font-weight=500"> *class ResStats* </span> _ *storage and transformation of the statistics on the residues of the contact maps in a zDock or MegaDock experiment* 

* Attributes 
	- **.rDict** : dictionary allowing access to every residue from its index 
	- **.expSize** : Total number of decoys taken into account for statistical calculation
	- **.plainResDict** : counts dictionary {residue : count}
	- **.resFreq** : residue frequency dictionary {residue : frequency}
	- **.pondResDict** : counts in contacts dictionary 
* Available Methods
	- **.setSize(n)** : set a size used for normalisation in frequencies
	- **.addRes(residue)** : add a residue to the statistics 
	- **.write(filename)** : write the residue frequencies of the docking experiment to a file

</br>

<span style="color:Crimson ;font-weight=500"> *class ContactStats* </span> _ *storage and transformation of the statistics on the residues of the contact maps in a zDock or MegaDock experiment* </br>
*see also <https://github.com/glaunay/pyproteins/tree/master/src/pyproteins/container/Core.py> for inherited attributes and functions from MdTree*

* Attributes 
	- **.expSize** : Total number of decoys taken into account for statistical calculation
* Available Methods
	- **.incrMdTree(A,B)** : adds 1 to the count of contact ocurrences between A and B
	- **.setSize(n)** : set a size used for normalisation in frequencies
	- **.all()** : returns all contact counts in the dataset allowing to perform statistical analyses and overall view of the counts distribution
	- **.get(A,B)** : returns the number of contact ocurrences between A and B
	- **.render_table(n=None)** : returns a table of size n**2 with counts in the intersection of columns and lines, each contact has a starting count of 1/expSize
	- **.write(filename)** : write the residue frequencies of the docking experiment to a file

</br>
<a id="clustering"></a>
### src.core_clustering
<span style="color:Crimson ;font-weight=500"> *class Cluster* </span> _ *store clusters : group of poses* 


* Attributes 
	- **.size** : number of poses in cluster
	- **.poses**: list of poses objects in cluster
	- **.bounds** : Returns maximal and minimal value of points in cluster for each axis ([xmin,xmax],[ymin,ymax],[zmin,zmax])
	- **.representative** : returns 1st pose in cluster
* Available Methods
	- **.meanRank(ranks)** : Takes ranks of poses in the original order (by id) and returns the average rank of the cluster.

</br>

<span style="color:Crimson ;font-weight=500"> *class ClusterColl* </span> _ *store cluster :  group of clusters obtained from clustering* 

* Attributes 
	- **.size** : number of clusters in collection
* Available Methods
	- **.addCluster(cluster):** : Manually add a cluser object to the collection. It will be added last
	- **.setClusters(clusters):**  Build ClusterColl from a cluster's python dictionnary containing Docking poses and errase the previous clusters in collection
	- **.sorted(element='original_rank' ,min_size=None):** Returns a list of sorted cluster objects based on mean rank of the clusters 
	- **.representatives(element =None, min_size=None):** Returns cluster's representatives in a particular order (sorted with average rank of clusters). It also allows you to suppress from the list those clusters with sizes below threshold.
	- **.addCluster(cluster):** : Manually add a cluser object to the collection. It will be added last
countNatives(self, poseList, cutoff=5):
</br>

-------------------------------



#### 

