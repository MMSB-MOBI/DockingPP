#! /usr/bin/env python3
import sys, os
sys.path.append("DockingPP")
from dockingPP import zParse
from rotation_utils import where


if __name__ == '__main__' :

    # Ajouter documentation et parser d'arguments
    
    complex=sys.argv[1]
    DD=zParse(f"../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/results/{complex}.zd3.0.2.fg.fixed.out", maxPose = 500)
    DD.setReceptor(f"../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/{complex}_r_u.pdb.ms_2")
    DD.setLigand(f"../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/{complex}_l_u.pdb.ms_2")

    DD.ccmap(dist=5)
    DD.write_all_scores(filename=f"../Docking/Resultats/scores/zD_{complex}_500")
    C=DD.pList[0].resMapList
    print(C)
    print(where(C,'Rec'))
    print(where(C,'Lig'))
