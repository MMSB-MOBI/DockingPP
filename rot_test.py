#! /usr/bin/env python3
import sys, os
sys.path.append("DockingPP")
from dockingPP import zParse
from rotation_utils import where


if __name__ == '__main__' :
    DD=zParse("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/results/1HIA.zd3.0.2.fg.fixed.out", maxPose = 2)
    DD.setReceptor("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_r_u.pdb.ms_2")
    DD.setLigand("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_l_u.pdb.ms_2")

    DD.pList[0].ccmap(dist=5)
    C=DD.pList[0].resMapList
    print(C)
    print(where(C,'Rec'))
    print(where(C,'Lig'))
