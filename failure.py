# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 23:16:02 2021

@author: Joshua Fung
"""
def tsai_wu(s, sigma_lp, sigma_ln, sigma_tp, sigma_tn, tau_lt):
    s = s.A1
    f11 = 1/(sigma_lp*sigma_ln)
    f22 = 1/(sigma_tp*sigma_tn)
    f1 = 1/sigma_lp - 1/sigma_ln
    f2 = 1/sigma_tp - 1/sigma_tn
    f66 = 1/tau_lt**2
    f12 = -((f11*f22)**0.5)/2
    # print("F11: %g, F22: %g, F66:%g, F1: %g, F2: %g, F12: %g"  % (f11,f22,f66,f1,f2,f12))
    t = f11*s[0]**2+f22*s[1]**2+f66*s[2]**2 + f1*s[0] + f2*s[1] + 2*f12*s[0]*s[1] 
    if t >= 1:
        return True, t
    return False, t

def tsai_hill(s, sigma_lp, sigma_ln, sigma_tp, sigma_tn, tau_lt):
    s = s.A1
    # use +ve when sigma 11 is +vs, -ve when sigma 11 is -ve, mixed when 0
    if s[0] > 0:
        sigma_l = sigma_lp
    elif s[0] < 0:
        sigma_l = sigma_ln
    else:
        sigma_l = 0.5 * sigma_lp + 0.5 * sigma_ln
    
    if s[1] > 0:
        sigma_t = sigma_tp
    elif s[1] < 0:
        sigma_t = sigma_tn
    else:
        sigma_t = 0.5 * sigma_tp + 0.5 * sigma_tn

    t = (s[0]/sigma_l)**2 - ((s[0]*s[1])/(sigma_l**2)) + (s[1]/sigma_t)**2 + (s[2]/tau_lt)**2
    if t >= 1:
        return True, t
    return False, t