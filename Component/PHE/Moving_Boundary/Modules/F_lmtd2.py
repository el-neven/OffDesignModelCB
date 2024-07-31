# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 14:50:11 2021

@author: jvega
"""
import numpy as np
from scipy.optimize import fsolve
import math

#%%
def F_lmtd2(R,P):
    
    def res_NTU_crossflow(Ntu, epsilon, C_r):
        Ntu = max(1e-4, Ntu)
        epsilon_g = 1-np.exp((1/C_r)*Ntu**0.22*(np.exp(-C_r*Ntu**0.78)-1))
        res = (epsilon-epsilon_g)
        return res

    #R = (T_h_i - T_h_ex)/(T_c_ex - T_c_i) --> must be positive [0; infinte [
    #P = (T_c_ex - T_c_i)/(T_h_i - T_c_i) --> must be within 0 and 1
    #-----------------------------------------------------------------------------
    #R and P must be positive
    R = abs(R)
    P = abs(P)
    #-----------------------------------------------------------------------------
    R_min = 1e-4
    R_max = 100
    
    R = max(R_min,min(R_max,R))
    if R < 0.41:
        P_max_for_given_R = 0.99999;
    else:
        P_max_for_given_R = 0.995*((0.966475661367996)*R**3 + (-1.431274296912407)*R**2 + (0.247230065033875)*R + (0.309118607270897)) / (R**4 + (-1.766309686745371)*R**3 + (1.287179055148762)*R**2 + (-0.902512766020191)*R + (0.484880205333508))
    P = max(0,  min(P_max_for_given_R,P))
    
    if R <= 1: # Cdot_min = Cdot_h
        epsilon = P
        Cr = R
        Pbis = P
        Rbis = R
    else: # Cdot_max = Cdot_c
        epsilon = P*R
        Cr = 1/R
        Pbis = P*R
        Rbis = 1/R
    
    f = lambda x: res_NTU_crossflow(x, epsilon, Cr)
    
    out_fsolve = fsolve(f, 1, full_output= 1)#, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)    
    NTU, res_NTU,flag_ntu = float(out_fsolve[0]), float(out_fsolve[1].get('fvec')), out_fsolve[2]
    if abs(res_NTU)< 1e-2  and flag_ntu > 0:
        if  R==1:
            F=P/NTU/(1-P)
        elif P <= 1e-04:
            F = 1
        else:
            F = epsilon*math.log((Pbis-1)/(Rbis*Pbis -1))/NTU/Pbis/(Rbis-1);
        flag = flag_ntu
    else:
        #disp(['Prob compute NTU with F_lmtd2(' num2str(R) ',' num2str(P) ')'])
        F = 1;
        flag = flag_ntu
    return F, flag

#%%

# R = -118189046.70120256
# P = -4.5784805034544295e-10

    
# def res_NTU_crossflow(Ntu, epsilon, C_r):
#         Ntu = max(1e-4, Ntu)
#         epsilon_g = 1-np.exp((1/C_r)*Ntu**0.22*(np.exp(-C_r*Ntu**0.78)-1))
#         res = (epsilon-epsilon_g)
#         return res

#     #R = (T_h_i - T_h_ex)/(T_c_ex - T_c_i) --> must be positive [0; infinte [
#     #P = (T_c_ex - T_c_i)/(T_h_i - T_c_i) --> must be within 0 and 1
# #-----------------------------------------------------------------------------
# #The
# R = abs(R)
# P = abs(P)
# #-----------------------------------------------------------------------------
# R_min = 1e-4
# R_max = 100
# R = max(R_min,min(R_max,R))
# if R < 0.41:
#     P_max_for_given_R = 0.99999
# else:
#     P_max_for_given_R = 0.995*((0.966475661367996)*R**3 + (-1.431274296912407)*R**2 + (0.247230065033875)*R + (0.309118607270897)) / (R**4 + (-1.766309686745371)*R**3 + (1.287179055148762)*R**2 + (-0.902512766020191)*R + (0.484880205333508))

# P = max(0,  min(P_max_for_given_R,P))
    
# if R <= 1: # Cdot_min = Cdot_h
#         epsilon = P
#         Cr = R
#         Pbis = P
#         Rbis = R
# else: # Cdot_max = Cdot_c
#         epsilon = P*R
#         Cr = 1/R
#         Pbis = P*R
#         Rbis = 1/R
    
# f = lambda x: res_NTU_crossflow(x, epsilon, Cr)
    
# out_fsolve = fsolve(f, 1, full_output= 1)#, args=(), fprime=None, full_output=0, col_deriv=0, xtol=1.49012e-08, maxfev=0, band=None, epsfcn=None, factor=100, diag=None)    
# NTU, res_NTU,flag_ntu = float(out_fsolve[0]), float(out_fsolve[1].get('fvec')), out_fsolve[2]
# if abs(res_NTU)< 1e-2  and flag_ntu > 0:
#         if  R==1:
#             F=P/NTU/(1-P)
#         elif P <= 1e-06:
#             F = 1
#         else:
#             F = epsilon*math.log((Pbis-1)/(Rbis*Pbis -1))/NTU/Pbis/(Rbis-1)
#         flag = flag_ntu
# else:
#         #disp(['Prob compute NTU with F_lmtd2(' num2str(R) ',' num2str(P) ')'])
#         F = 1;
#         flag = flag_ntu
#     # return F, flag
