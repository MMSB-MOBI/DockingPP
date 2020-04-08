#! /usr/bin/env python3
import sys, os, numpy as np
from math import cos, sin, acos, asin, atan2, sqrt , pi

def is_zero(num):
    return ( num <= 1e-14 and num >= -1e-14)

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

def eulerFromMatrix(matrix):
    t=acos(matrix[2,2])
    t2=-1*t
    if not is_zero(sin(acos(matrix[2,2]))) :
        f=atan2(matrix[2,0]/t, matrix[2,1]/t)
        f2=atan2(matrix[2,0]/t2, matrix[2,1]/t2)
        # t=atan2(sqrt(matrix[2,0]**2+matrix[2,1]**2),matrix[2,2] )
        p=atan2(matrix[0,2]/t, matrix[1,2]/t2)
        p2=atan2(matrix[0,2]/t2, matrix[1,2]/t)
    else :
        if equals(matrix[2,2], -1) :
            p= 1-atan2(matrix[1,0], matrix[0,0])
            t=pi
            f=1
        if equals(matrix[2,2], 1) :
            f=atan2(matrix[1,0], matrix[0,0]) /2
            t=0
            p=atan2(matrix[1,0], matrix[0,0]) /2
        else :
            # print(matrix[2,2])
            return 4
    return (p,t,f)

