#! /usr/bin/env python3
import sys, os, math, numpy as np
import pyproteinsExt.structure.coordinates as PDB
sys.path.append("DockingPP")
from dockingPP import zParse

parser = PDB.Parser()
pdbObj = parser.load(file="data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_r_u_2.pdb")

phi=2
theta=5
psi=3

# print(pdbObj.atomDictorize['x'])

def trans_matrix(psi,theta,phi):
    """Construction de la matrice de transition a partir des angles d'euler"""
    cp=math.cos(psi)
    sp=math.sin(psi)
    ct=math.cos(theta)
    st=math.sin(theta)
    cf=math.cos(phi)
    sf=math.sin(phi)

    tr_matrix = np.matrix([[cp*cf-sp*ct*sf,sp*cf+cp*ct*sf,st*sf ],[-cp*sf-sp*ct*sf,-sp*sf+cp*ct*cf,st*cf],[sp*st,-cp*st,ct]])
    tr_matrix=tr_matrix.T
    return tr_matrix



def coord_matrix(P):
    """O,P tuples with 3 values : the x, y and z coordinates of a rotation origin - O
    and of a rotating point or atom - P
    Returns a numpy array with dimensions (1,3) """

    rel_coordinates=np.matrix([[P[i]] for i in range(3)])
    return rel_coordinates

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

# A revoir l'efficacité de cette fonction
# def is_orthogonal(M):
#     if M.T.all() == M.I.all() :
#         print("This is an orthogonal matrix")
#         return True
#     else :
#         return False

DD=zParse("data/decoys_bm4_zd3.0.2_6deg_fixed/results/1HIA.zd3.0.2.fg.fixed.out")
# print(DD.baryREC)
# print(DD.eulerREC)
DD.setReceptor("data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_r_u_2.pdb")
DD.setLigand("data/decoys_bm4_zd3.0.2_6deg_fixed/input_pdbs/1HIA_l_u_2.pdb")

# print(DD.pdbObjLigand.atomDictorize['x'])

# print(len(DD.pdbObjReceptor.model[0]))
# print(DD.pdbObjReceptor.model[0][-1])
# while i < len(DD.pdbObjReceptor.model[0]):
    # sys.stdout.write(str(i))
    # sys.stdout.flush()
    # atom=DD.pdbObjReceptor.model[0][i]
for atom in DD.pdbObjReceptor.model[0]:
    initial_coords=(atom.x, atom.y, atom.z)
    a=trans_matrix(*DD.eulerREC)
    b=coord_matrix(initial_coords)
    # print(b)
    coords=DD.baryREC,transpose(a,b)
    atom.x=coords[0,0]
    atom.y=coords[1,0]
    atom.z=coords[2,0]


# DD.pdbObjReceptor.atomdictorize['x']=new_coords['x']
# DD.pdbObjReceptor.atomdictorize['y']=new_coords['y']
# DD.pdbObjReceptor.atomdictorize['z']=new_coords['z']
# print(str(DD.pdbObjReceptor.model[0]))
DD.pList[0].ccmap(dist=5)
C=DD.pList[0].resMapList
# print(C)
# resS=DD.resStats
print(C)
# print(C['data'])

# print(DD.pdbObjReceptor)
# a=trans_matrix(*DD.eulerREC)
# b=coord_matrix(DD.baryREC,(50,48,32))
#
# print(transpose(a,b))
# print(recenter((2,3,2),transpose(a,b)))
