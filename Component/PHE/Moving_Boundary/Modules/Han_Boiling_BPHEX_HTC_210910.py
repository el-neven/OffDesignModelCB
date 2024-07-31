# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 16:25:26 2021

@author: jvega

reference : https://ht.readthedocs.io/en/release/ht.boiling_plate.html
ht.boiling_plate.h_boiling_Han_Lee_Kim(m, x, Dh, rhol, rhog, mul, kl, Hvap, Cpl, q, A_channel_flow, wavelength, chevron_angle=45.0)

"""
from scipy.optimize import fsolve
#%%
def Han_Boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G, DT_log, Qdot, honv_h, Dh, theta, pitch_co):
    
    def iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        Bo_g = max(Bo_g, 1e-8)
        Nu = Ge1*Re_eq**Ge2*Bo_g**0.3*Pr_l**0.4
        h = Nu*k_l/Dh
        U = (1/h +  1/honv_h)**-1
        A_tp = AU_tp/U
        q = Qdot/A_tp
        Bo = q/(G_eq*i_fg)
        res_Bo = (Bo-Bo_g)/Bo_g
        return res_Bo, Nu, h, U, A_tp, q, Bo
    
    def res_iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh):
        res_Bo,_,_,_,_,_,_ = iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
        
        return res_Bo
    
    G_eq = G * ( (1 - x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    AU_tp = Qdot/DT_log
    Ge1 = 2.81*(pitch_co/Dh)**(-0.041)*(theta)**(-2.83)
    Ge2 = 0.746*(pitch_co/Dh)**(-0.082)*(theta)**(0.61)
    Bo_0 = 0.5
    f_Bo = lambda xx: res_iter_Han_boiling(xx, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    sol = fsolve(f_Bo, Bo_0,)
    Bo_sol = sol[0]
    _, Nu, h, _, _, _, _ = iter_Han_boiling(Bo_sol, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
    
    h_boiling = h
    
    return h_boiling, Nu