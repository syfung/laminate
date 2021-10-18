# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 17:13:43 2021

Laminate stiffness matrix calculator 

Calculate the stiffness matrix of a given ply 

@author: Joshua Fung
"""
import numpy as np

class Laminate:
    def __init__(self, angles, thickness, e_l, e_t, nu_lt, g_lt):
        self.plies = []
        for angle in angles:
            p = Ply(angle, thickness, e_l, e_t, nu_lt, g_lt)
            self.plies.append(p)
class Ply:
    def __init__(self, angle, thickness, e_l, e_t, nu_lt, g_lt):
        """Ply class"""
        self.angle = angle
        self.thickness = thickness
        self.q = self._q_matrix(e_l, e_t, nu_lt, g_lt)
        self.t = self._t_matrix(self.angle)
        self.q_bar = self._q_rotate(self.angle, self.q, self.t)
        
    def _t_matrix(self, angle):
        angle = np.radians(angle)
        c = np.cos(angle)
        s = np.sin(angle)
        
        t = np.matrix([[c**2, s**2, 2*c*s],
                       [s**2, c**2, -2*c*s],
                       [-c*s, c*s, c**2 - s**2]])
        
        return t
        
    def _q_matrix(self, E_l, E_t, nu_lt, G_lt):
        """Stiffness Matrix"""
        nu_tl = nu_lt *E_t / E_l
        q = np.matrix([[E_l/(1-nu_lt*nu_tl), (nu_lt*E_t)/(1-nu_lt*nu_tl), 0],
                       [(nu_tl*E_l)/(1-nu_lt*nu_tl), E_t/(1-nu_lt*nu_tl), 0],
                       [0, 0, G_lt]])
        return q
    
    def _q_rotate(self, angle, q, t):
        q_bar = np.linalg.inv(t) * q * np.linalg.inv(np.transpose(t))
        return q_bar
        
if __name__ == "__main__":
    p = Ply(30, 0.125, 131, 9, 0.22, 6)
    print(p.q_bar)
    
    l = Laminate([0,30,60,90], 0.125, 131, 9, 0.22, 6)
    

    

