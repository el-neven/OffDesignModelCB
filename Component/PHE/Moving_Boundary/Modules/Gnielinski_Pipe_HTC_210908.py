# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:25:03 2021

Reference: VDI section G1, page 696

@author: jvega
"""

from math import log10, sqrt, inf
import warnings
#%%
def Gnielinski_Pipe_HTC(mu, Pr, k, G, Dh, L):
    # print(mu, Pr, k, G, Dh, L)
    #-------------------------------------------------------------------------
    def Gnielinski_Laminar(Re, Pr, Dh, L):
        Nu_1 = 4.364
        Nu_2 = 1.953*(Re*Pr*Dh/L)**(1/3)
        #print(Nu_2)
        Nu = (Nu_1**3 + 0.6**3 + (Nu_2 - 0.6)**3)**(1/3)
        return Nu
    def Gnielinski_Turbulent(Re, Pr):
        f = (1.8*log10(Re) - 1.5)**(-2)
        Nu = ((f/8)*(Re-1000)*Pr)/(1+12.7*sqrt(f/8)*(Pr**(2/3)-1))
        return Nu
    #-------------------------------------------------------------------------
    Re_min = 0
    Re_max = 1e06
    Re = G*Dh/mu
    #-------------------------------------------------------------------------
    if Re > 1e4: #fully turbulent
        Pr_min = 0.1
        Pr_max = 1000
        Nu = Gnielinski_Turbulent(Re, Pr)
    elif Re < 2300: #fully laminar
        Pr_min = 0.6
        Pr_max = inf
        Nu = Gnielinski_Laminar(Re, Pr, Dh, L)
    else: #transition zone
        Pr_min = 0.1
        Pr_max = 1000
        gamma = (Re - 2300)/(1e4 - 2300)
        Nu_lam2300 = Gnielinski_Laminar(2300, Pr, Dh, L)
        Nu_turb10000 = Gnielinski_Turbulent(1e4, Pr)
        Nu = (1-gamma)*Nu_lam2300 + gamma*Nu_turb10000
    #-------------------------------------------------------------------------
    hConv = Nu*k/Dh
    #-------------------------------------------------------------------------
    if Re >= Re_max or Pr <=Re_min:
        # warnings.warn('Gnielinski singe-phase: Out of validity range --> Re = ', Re, ' is out of [', Re_min, ' - ', Re_max, '] !!!')
        warnings.warn('Gnielinski singe-phase: Reynolds Out of validity range !!!')
    if Pr >= Pr_max or Pr <= Pr_min:
        # warnings.warn('Gnielinski singe-phase: Out of validity range --> Re = ', Pr, ' is out of [', Pr_min, ' - ', Pr_max, '] !!!')
        warnings.warn('Gnielinski singe-phase: Prandtl Out of validity range  !!!')
    #-------------------------------------------------------------------------
    return hConv, Nu