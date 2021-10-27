# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:22:11 2021

@author: Joshua Fung
"""
import numpy as np

def deformation(abd_inv, load):
    return abd_inv.dot(load)

def strain(laminate, mid_plane_deformation):
    mid_ply_z = []
    for k in range(len(laminate.z_s)):
        mid_ply_z.append((laminate.z_s[k][1] - laminate.z_s[k][0])/2)
    ply_deformation = np.empty([0,3])
    for z in mid_ply_z:
        ply_deformation = np.append(ply_deformation,
                                    [mid_plane_deformation[0:3] + z * mid_plane_deformation[3:6]],
                                    axis=0)
    return ply_deformation
