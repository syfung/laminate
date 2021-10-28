# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 17:22:11 2021

@author: Joshua Fung
"""
import numpy as np

def deformation(abd_inv, load):
    return abd_inv.dot(load)

def strain_mid_ply(laminate, mid_plane_deformation):
    ply_deformation = []
    for z in laminate.mid_ply_zs:
        strain = mid_plane_deformation[:3]
        curvature  = mid_plane_deformation[3:]
        ply_deformation.append(strain + z * curvature)
    return ply_deformation

def strain_top_bottom_ply(laminate, mid_plane_deformation):
    ply_deformation = []
    for z in laminate.z_s:
        strain = mid_plane_deformation[:3]
        curvature  = mid_plane_deformation[3:]
        ply_deformation.append((strain + z[0] * curvature, strain + z[1] * curvature))
    return ply_deformation

def stress_mid_ply(laminate, mid_plane_deformation):
    strain = strain_mid_ply(laminate, mid_plane_deformation)
    stress = []
    for k in range(len(laminate.plies)):
        stress.append(laminate.plies[k].q_bar.dot(strain[k]))
    return stress

def stress_top_bottom_ply(laminate, mid_plane_deformation):
    strain = strain_top_bottom_ply(laminate, mid_plane_deformation)
    stress = []
    for k in range(len(laminate.plies)):
        stress.append((laminate.plies[k].q_bar.dot(strain[k][0]), laminate.plies[k].q_bar.dot(strain[k][1])))
    return stress
        
if __name__ == "__main__":
    import stiffmat
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=6)

    # =============================================================================
    # Ply properties
    # =============================================================================
    # Elastic properties
    E1 = 125 * 1e3  # MPa
    E2 = 9.8 * 1e3  # MPa
    G12 = 5.5 * 1e3  # MPa
    nu12 = 0.24  # -
    
    # Failure properties
    sigma_lp = 900  # MPa
    sigma_ln = 800  # MPa
    sigma_tp = 55  # MPa
    sigma_tn = 170  # MPa
    tau_lt = 90  # MPa
    
    # Thickness
    t = 0.125 # mm
    
    # =============================================================================
    # Load Vector
    # =============================================================================
    load = np.matrix((240., 82., 4., -63., 0., 0.)).T
    print("Applied Load")
    print("({0[0]:2.2f}N/mm, {1[0]:2.2f}N/mm, {2[0]:2.2f}N/mm, {3[0]:2.2f}N, {4[0]:2.2f}N,{5[0]:2.2f}N).T\n".format(*np.array(load)))
    
    lam = stiffmat.Laminate([0, 60, 0, 30], t, E1, E2, nu12, G12)
    print("ABD Matrix:")
    print(lam.abd)
    print("Unit:")
    print(np.matrix([["N/mm (MPa.mm)", "N (MPa.mm2)"],["N (MPa.mm2)", "N.mm(MPa.mm3)"]]), "\n")
    
    mid_plane_deformation = deformation(lam.abd_inv, load)
    print("Mid plane deforamtion:")
    print("({0[0]:2.4f} {1[0]:2.4f} {2[0]:2.4f} {3[0]:2.4f}1/mm {4[0]:2.4f}1/mm {5[0]:2.4f}1/mm).T\n".format(*np.array(mid_plane_deformation)))
    m_deform = mid_plane_deformation.copy()
    m_deform[3:6] = m_deform[3:6]*1000
    print("({0[0]:2.4f} {1[0]:2.4f} {2[0]:2.4f} {3[0]:2.4f}1/m {4[0]:2.4f}1/m {5[0]:2.2f}1/m).T\n".format(*np.array(m_deform)))
    
    strain = strain_mid_ply(lam, mid_plane_deformation)
    print(strain)
    strain_top_bottom = strain_top_bottom_ply(lam, mid_plane_deformation)
    # print(strain_top_bottom)
    plt.figure()
    plt.plot(lam.mid_ply_zs, [s.item(0) for s in strain], "x")
    plt.plot(list(sum(lam.z_s, ())), [s.item(0) for s in (list(sum(strain_top_bottom, ())))])
    
    stress = stress_mid_ply(lam, mid_plane_deformation)
    stress_top_bottom = stress_top_bottom_ply(lam, mid_plane_deformation)
    
    plt.figure()
    plt.subplot(1, 3, 1)
    plt.plot([s.item(0) for s in (list(sum(stress_top_bottom, ())))], list(sum(lam.z_s, ())), "c-x")
    plt.plot([s.item(0) for s in stress], lam.mid_ply_zs, "kx")
    ax = plt.gca()
    ax.invert_yaxis()
    
    plt.subplot(1, 3, 2)
    plt.plot([s.item(1) for s in (list(sum(stress_top_bottom, ())))], list(sum(lam.z_s, ())), "c-x")
    plt.plot([s.item(1) for s in stress], lam.mid_ply_zs, "kx")
    ax = plt.gca()
    ax.invert_yaxis()
    
    plt.subplot(1, 3, 3)
    plt.plot([s.item(2) for s in (list(sum(stress_top_bottom, ())))], list(sum(lam.z_s, ())), "c-x")
    plt.plot([s.item(2) for s in stress], lam.mid_ply_zs, "kx")
    ax = plt.gca()
    ax.invert_yaxis()