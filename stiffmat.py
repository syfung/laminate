# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 17:13:43 2021

Laminate stiffness matrix calculator 

Calculate the stiffness matrix of a given ply 

@author: Joshua Fung
"""
import numpy as np

def transformation_matrix(angle):
    angle = np.radians(angle)
    c = np.cos(angle)
    s = np.sin(angle)
    
    t = np.matrix([[c**2, s**2, 2*c*s],
                   [s**2, c**2, -2*c*s],
                   [-c*s, c*s, c**2 - s**2]])    
    return t

def stiffness_matrix(E_1, E_2, nu_12, G_12):
       """Stiffness Matrix"""
       nu_21 = nu_12 *E_2 / E_1
       q = np.matrix([[E_1/(1-nu_12*nu_21), (nu_12*E_2)/(1-nu_12*nu_21), 0],
                      [(nu_21*E_1)/(1-nu_12*nu_21), E_2/(1-nu_12*nu_21), 0],
                      [0, 0, G_12]])
       return q

def q_rotate(angle, q):
    t = transformation_matrix(angle)
    q_bar = np.linalg.inv(t) * q * np.linalg.inv(np.transpose(t))
    return q_bar, t

class Laminate:
    def __init__(self, angles, thickness, e_l, e_t, nu_lt, g_lt):
        self.plies = []
        self.thicknesses = []
        self.z_s = []
        self.mid_ply_zs = []
        
        for angle in angles:
            p = Ply(angle, thickness, e_l, e_t, nu_lt, g_lt)
            self.plies.append(p)
            self.thicknesses.append(thickness)
            
        self.h = sum(self.thicknesses) / 2
                
        for k in range(len(angles)):
            z_top = sum(self.thicknesses[:k+1]) - self.h
            z_bot = sum(self.thicknesses[:k]) - self.h
            
            self.z_s.append((z_bot, z_top))
            self.mid_ply_zs.append((z_top - z_bot)/2 + z_bot)
            
        self.abd = self._abd_matrix()
    
    @property
    def abd_inv(self):
        return np.linalg.inv(self.abd)
            
    def _abd_matrix(self):
        a = np.zeros([3,3])
        b = np.zeros([3,3])
        d = np.zeros([3,3])
        
        for k in range(len(self.thicknesses)):
            ai = self.plies[k].q_bar * (self.z_s[k][1] - self.z_s[k][0])
            bi = (1/2) * self.plies[k].q_bar * (self.z_s[k][1]**2 - self.z_s[k][0]**2)
            di = (1/3) * self.plies[k].q_bar * (self.z_s[k][1]**3 - self.z_s[k][0]**3)
            
            a = a + ai
            b = b + bi
            d = d + di
        
        ab = np.concatenate((a,b), axis=1)
        bd = np.concatenate((b, d), axis=1)
        abd = np.matrix(np.concatenate((ab, bd)))
        
        abd = np.matrix(np.where(np.abs(abd) < np.max(abd)*1e-6, 0, abd))
        
        return abd


class Ply:
    def __init__(self, angle, thickness, e_1, e_2, nu_12, g_12):
        """Ply class"""
        self.angle = angle
        self.thickness = thickness
        self.q = stiffness_matrix(e_1, e_2, nu_12, g_12)
        self.q_bar, self.t = q_rotate(self.angle, self.q)
        
        
if __name__ == "__main__":
    np.set_printoptions(precision=4)
    p = Ply(30, 0.125, 131, 9, 0.22, 6)
    print("Q_bar:\n", p.q_bar)
    
    l = Laminate([0,30,60,90], 0.125, 131, 9.8, 0.22, 5.8)
    print("ABD Matrix:\n", l.abd)
    

    

