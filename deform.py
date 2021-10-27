# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:22:11 2021

@author: Joshua Fung
"""
import numpy as np

def deformation(abd_inv, load):
    return abd_inv.dot(load)
