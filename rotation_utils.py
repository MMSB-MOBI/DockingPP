#! /usr/bin/env python3
import sys, os, numpy as np
from math import cos, sin, acos, asin, atan2, sqrt



# print(pdbObj.atomDictorize['x'])
def where(O,astr):
    count=0
    for i in O:
        if astr in i:
            count+=1
    return count

def rotate(x,y,z,psi,theta,phi) :
    r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi)
    r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi)
    r31 = sin(theta)*sin(phi)

    r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi)
    r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi)
    r32 = sin(theta)*cos(phi)

    r13 = sin(psi)*sin(theta)
    r23 = -cos(psi)*sin(theta)
    r33 = cos(theta)

    new_x= r11 * x + r12 * y + r13 * z
    new_y= r21 * x + r22 * y + r23 * z
    new_z= r31 * x + r32 * y + r33 * z

    return (new_x, new_y, new_z)


def trans_matrix(psi,theta,phi):
    """Construction de la matrice de transition a partir des angles d'euler"""
    r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi)
    r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi)
    r31 = sin(theta)*sin(phi)

    r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi)
    r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi)
    r32 = sin(theta)*cos(phi)

    r13 = sin(psi)*sin(theta)
    r23 = -cos(psi)*sin(theta)
    r33 = cos(theta)

    tr_matrix = np.matrix([[r11, r12, r13 ],[r21, r22, r23],[r31, r32, r33]])

    return tr_matrix

def euleurFromMatrix(matrix):
    t1=acos(matrix[2,2])
    t2= -1 * t1
    if matrix[2,0] != 0 and matrix[2,1] != 0 and matrix[0,2] != 0 and matrix[1,2] != 0 :
        f=atan2(matrix[2,0], matrix[2,1])
        t=atan2(acos(matrix[2,2]),sqrt(matrix[2,0]**2+matrix[2,1]**2) )
        p=atan2(matrix[0,2]/sin(t1), matrix[1,2] / sin(t2))


    return (p,t1,f)

def vertical_matrix(P):
    ver_coordinates=np.matrix([[i] for i in P])
    return ver_coordinates


def translate(x,y,z,t1,t2,t3):

    return (x + t1, y + t2, z + t3)


def transpose(A,B):
    """A,B two numpy matrixes A a transition matrix, 'matrice de passage'
    B the points to be transferred """
    return A.dot(B)

def recenter(O,P):
    """O,P tuples with 3 values : the x, y and z coordinates of a rotation origin - O
    and of a rotated point or atom - P
    Returns a numpy array with dimensions (1,3) """

    # assert len(O) == len(P) == 3

    rel_coordinates=P+np.matrix([[o] for o in O])
    return rel_coordinates

if __name__ == '__main__':

    DD=zParse("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/results/1HIA.zd3.0.2.fg.fixed.out")
    DD.setReceptor("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_r_u.pdb.ms_2")
    DD.setLigand("../Docking/data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_l_u.pdb.ms_2")
    # DD=zParse("../data/decoys_bm4_zd3.0.2_6deg_fixed/results/1HIA.zd3.0.2.fg.fixed.out")
    # DD.setReceptor("../data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_r_u_2.pdb")
    # DD.setLigand("../data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_l_u_2.pdb")


    # print(DD.pdbObjLigand.atomDictorize['x'])

    # print(len(DD.pdbObjReceptor.model[0]))
    # print(DD.pdbObjReceptor.model[0][-1])
    # while i < len(DD.pdbObjReceptor.model[0]):
        # sys.stdout.write(str(i))
        # sys.stdout.flush()
        # atom=DD.pdbObjReceptor.model[0][i]
    # print(DD.pdbObjLigand)
    expected="ATOM      1  N   THR     5      13.140  23.779  -6.666"

    # print(f"baryLIG : {DD.baryLIG} \n random : {DD.eulerREC} \n pose rotation : {DD.pList[0].euler}\n pose translation : {DD.pList[0].translate} \n baryREC : {DD.baryREC}")

    # This works !!!!
    # for atom in DD.pdbObjLigand.model[0]:
    #     if atom.serial==1 :
    #         initial_coords=(atom.x, atom.y, atom.z)
    #         print("\nInitial state = " + str(initial_coords))
    #         neg_bary_lig=[-1*i for i in DD.baryLIG]
    #
    #         coords=translate(*initial_coords, *neg_bary_lig)
    #         print("Translation to origin = " + str(coords))
    #         #Apply Initial Random rotation
    #         coords2=rotate(*coords,*DD.eulerREC)
    #         #Apply pose Rotation
    #         coords3=rotate(*coords2,*DD.pList[0].euler)
    #         print("Rotated coordinates = " + str(coords3))
    #
    #         # trans=[i for i in DD.pList[0].translate]
    #         coords4=translate(*coords3, *DD.baryREC)
    #         coords5=translate(*coords4, *DD.pList[0].translate)
    #         # print(coords4)
    #
    #
    #         print("Translation to receptor then ligand = " + str(coords5))
    #         # print(atom.chainID)
    #         # if atom.resSeq==1 :
    #             # print("x = " + str(atom.x))
    #         atom.x=coords5[0]
    #         atom.y=coords5[1]
    #         atom.z=coords5[2]
    #
    #         break


    # Make rotation matrices
    # a=trans_matrix(*DD.eulerREC)
    # b=trans_matrix(*DD.pList[0].euler)
    #
    # # Combine into one matrix
    # double=b.dot(a)
    #
    # # Recover combined angles
    # new_euler=euleurFromMatrix(double)
    #
    # # Crash test
    # E1=rotate(*DD.baryLIG, *DD.eulerREC)
    # E2=rotate(*E1, *DD.pList[0].euler)
    #
    #
    #
    # DD.pList[0].euler=new_euler
    # # print(double.dot(vertical_matrix(DD.baryLIG)))
    # # print(E2)
    # # print(double)
    # # print(trans_matrix(*new_euler))
    # #
    # # print(trans_matrix(*new_euler).dot(vertical_matrix(DD.baryLIG)))
    #
    # # Apply negative function to baryLIG and baryREC for negative translation in ccmap module
    # DD.baryLIG=tuple([-i for i in DD.baryLIG])
    # DD.baryREC=tuple([-i for i in DD.baryREC])

    DD.pList[0].ccmap(dist=5)
    C=DD.pList[0].resMapList

    # print(where(C,'Rec'))
    # print(where(C,'Lig'))
