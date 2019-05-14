### DOCKING POST-PROCESSING Utility Class and functions


# Table of Contents
  * [Dependencies](#chapter-1)
  * [Quick Start](#chapter-2)
  * [The megadock format](#chapter-3)
  * [Getting started](#chapter-4)
  * [Available objects and functions](#chapter-5)

---------------------------------
<a id="chapter-1"></a>
## Dependencies 
dockingPP has the following dependencies:

* **Python 3**
* pyProteinsExt (<https://github.com/glaunay/pyproteinsExt>)
* pyProteins (<https://github.com/glaunay/pyproteins>) 
* ccmap (<https://github.com/glaunay/ccmap>) 
* NumPy (<http://www.numpy.org/>)

</br>


<a id="chapter-2"></a>
## Quick Start 


```sh
git clone https://github.com/juliaprieto/dockingPP.git  LOCAL_PATH_TO_REPO
``` 
Check for dependencies in LOCAL\_PATH\_TO\_REPO/setup.py

```python
sys.path.append(LOCAL_PATH_TO_REPO/src)
``` 
</br>


<a id="chapter-3"></a>
## The megadock/zDock format 

A typical megadock results file is derived from the [ZDOCK format](http://zdock.umassmed.edu/zdock_output_format.php). 

- First line features the number of grid cells and the size of the cubic cell. 
- Line 3 indicates the initial position of the **receptor** barycenter within the original system of coordinates. 
- Line 4 indicates similar information for the **ligand**.
- From line 12 to 16, docking poses' information is provided.   

<span style="color:Crimson ;font-weight=500">In order to reconstruct a docking pose both **receptor** and **ligand** molecule barycenters have to be set to the origin.</span>  

- The first field is the **Euler angles** triplet to apply rotation to the **ligand** molecule to reconstruct its orientation in the pose.   
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
 
 </br> 

----------------------
<a id="chapter-4"></a>
### Getting started with Docking Post-Processing 
#### Docking prediction Output reading 

```python
from dockingPP import parse
zD= parse('/PATH_TO_FILE/1BJ1_r-1BJ1_l.detail', maxPose =45) # maxpose default is 0
```
*Be careful, since docking outputs vary and may even provide near native solutions for the study.*   
*For now, the program expects 7 near native solutions in the beginning of the file. These are not taken into account for residue or contact statistics. Therefore if maxpose=4000 you will get 7 near native solutions and 3993 computed poses*


#### Add reference PDB for receptor and ligand

```python
zD.setReceptor('/PATH/1BJ1_r_u.pdb')
zD.setLigand('/PATH/1BJ1_l_u.pdb')
```
Make sure you use unbound molecules to do your predictions, since the bound version of the structures is supposedly unknown. 
#### Add name to Docking object

```python
zD.setComplexName('1BJ1')
```

#### Compute interaction interface
This action will be the most ram consuming. It may take a while. 

This function uses the ccmap module    
 1* Each pose keeps a dictionary of its amino-acids contacts 10000 poses are treated in 4mn, with a memory cost ~ 9286M   
 2* Releasing the GIL w/in C extension does reduce total computation cost to 3mn40s. Bottleneck is single structure computation in C 

```python
zD.ccmap(start=0,stop=n,dist=5, pSize=50)
# start default : 0 
# stop default : len of zD.zList 
# dist default : 4.5
# pSize default : 200
```

</br>

<a id="chapter-5"></a>
## Available objects and functions

* <a href=#docking style="text-decoration: none;"> dockingPP.docking </a>
* <a href=#stats style="text-decoration: none;">dockingPP.core_stats</a>
* <a href=#cluster style="text-decoration: none;">dockingPP.core_clustering</a>
* <a href=#visuals style="text-decoration: none;">dockingPP.core_visuals</a>

</br>

<a id="docking"></a>
### dockingPP.docking

<span style="color:Crimson ;font-weight=500"> *class DockData* </span> _ *Docking results container*
we will consider as pose_index its original rank 

* Attributes
	- **complexName** 
	- **step** 
	- **nCells**
	- **eulerREC** : euler angles of the receptors' reference position
	- **fileREC** 
	- **fileLIG** 
	- **baryREC** : receptor's barycenter position
	- **baryLIG** : ligand's reference's barycenter
	- **zList** : list of poses in the experiment
	- **pdbObjLigand** : structure object built from PDB file
	- **pdbObjReceptor** : structure object built from PDB file
	- **getStats** : returns resStats and contactStats objects into a list
			`resStats,conStats=zD.getStats`
	- **contactStats** : returns contactStats object
	- **resStats** : returns resStats object
* Available Methods
	- **setComplexName**
	- **ccmap** : calculate and store ccmap for each pose 
	- **push** : add pose to zList 
	- **setReceptor** : add PDB reference
	- **setLigand** : add PDB reference
	- **dictPos** : returns a dictionary useful for visualisation { 'x' : , 'y' : } 
	- **bestPoses** : returns a list of the n best poses according to 1 rescoring method
	- **poseScores** : returns a dictionary : { pose_index : score }

<span style="color:Crimson ;font-weight=500"> *class Pose* </span> _ *Ligand Pose information*. 
  
* Attributes
	- **id**
	- **euler**
	- **translate**
	- **ccmap**:json formatted contact map
	- **belongsTo** : DockData object
	- ligOffset 
	- recOffset
	- dictorizedReceptor : positions for the receptors atoms
	- dictorizedLigand : positions for the ligands atoms
	- resMapList : returns the list of residues implied in the contact map
	- resSize : number of residues in pose contactmap
	- conSize : number of contacts in pose contactmap
* Available Methods
	- **SumScore** :returns the sum of the residues statistics, method can take two values : 'plain' , 'freq'
	- **MeanScore** : returns the mean of the residues statistics, method can take two values : 'plain' , 'freq'
	- **SquareSumScore** : returns the square sum of the residues statistics, method can take three values : 'plain' , 'freq' or 'log'</br>
*And so for contacts*
	- **cmapSumScore**
	- **cmapSquareSumScore**
	- has_ccmap : checks wether pose's ccmap has been calculated
	- ccmap : calculate contact map

</br>

<a id="stats"></a>
### dockingPP.core_stats
<span style="color:Crimson ;font-weight=500"> *class CmapRes* </span> _ *store residues from ccmap and their counts* 

* Attributes 
	- chainID
	- resID 
	- role : Rec or Lig residue
	- **count** : number of occurences in residues
	- **cCount** : number of occurences in contacts
	- index
* Available Methods
	- increase_count
	- **reset_all** : resets all counts to 0

<span style="color:Crimson ;font-weight=500"> *class ResStats* </span> _ *storage and transformation of the statistics on the residues of the contact maps in a zDock or MegaDock experiment* 

* Attributes 
	- rDict
	- expSize : Total number of decoys taken into account for statistical calculation
	- expName
	- **plainResDict** : counts dictionary {residue : count}
	- **resFreq** : residue frequency dictionary {residue : frequency}
	- **pondResDict** : counts in contacts dictionary 
* Available Methods
	- setSize
	- addRes
	- **write** : write the statistics of the experiment to a file. 

</br>

<span style="color:Crimson ;font-weight=500"> *class ContactStats* </span> _ *storage and transformation or the statistics on the residues of the contact maps in a zDock or MegaDock experiment* 

* Attributes 
	- see also <https://github.com/glaunay/pyproteins/tree/master/src/pyproteins/container/Core.py> for inherited attributes
	- expSize
	- expName
	- plainResDict
	- resFreq
	- pondResDict
* Available Methods
	- incrMdTree
	- setSize
	- get
	- render_table
	- write

</br>

<a id="cluster"></a>
### dockingPP.core_clustering 

</br>

<a id="visuals"></a>
### dockingPP.core_visuals 

</br>

-------------------------------



