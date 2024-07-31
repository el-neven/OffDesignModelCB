# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:47:52 2023

@author: Basile
"""

from math import log10, sqrt, inf
from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve
import warnings
#%%
def water_Plate_HTC(mu, Pr, k, G, Dh):
    """
    Calibrated heat transfer coeffcient correlation for water side
    
    Inputs
    ----------
    mu : Viscosity [kg/(m*s)]
    
    Pr : Prandtl Number [/]
    
    k : thermal conductivity [W/(m*K)]
        
    G : Mass flux [kg/(m^2 * s)]
    
    Dh : Spacing between plates [m]

    Outputs
    -------
    h_conv : HTC in convection
    
    Reference
    -------
    Refrigerant R134a vaporisation heat transfer and pressure drop inside a small brazed plate heat exchanger
    G.A. Longo, A. Gasparella

    """
    # Bounds on Re (arbitrary) # !!! shall be checked
    Re_max = 1e6
    Re_min = 5
    
    # Reynolds number
    Re = G*Dh/mu
    
    if Re <= Re_max and Re >= Re_min:
        # Nusselt number
        Nu = 0.277*Re**(0.766)*Pr**(0.333)
    else: 
        print("Reynolds number for water out of bounds.")
        return 0
        
    # HTC
    h_conv = (k/Dh)*Nu

    return h_conv

def forced_conv_boil():
    """
    ---- Inputs : -------- 
    
    
    ---- Outputs : --------
    
    h_conv_boil : Heat transfer coefficient related to forced convection boiling [W/(m^2)]
    
    ---- Reference(s) : --------
    
    Foundations of heat transfer : Frank P. Incropera
    
    """
    
    # rho_l = PropsSI('D','T',T_sat,'Q',0,fluid) # density at saturated liquid in kg/m^3
    # rho_v = PropsSI('D','T',T_sat,'Q',1,fluid) # density at saturated vapor in kg/m^3
    # sigma = PropsSI('I','T',T_sat,'Q',0.5,fluid) # surface tension of liquid vapor equilibrium
    # Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific heat
    
    # Weber number
    # We_D = (rho_v*V_flow**2*D_out)/sigma
    # coef_low = (1/np.pi)*(1 + (4/We_D)**(1/3))
    
    # Lienhard and Eichhorn correlations
    
    # if coef_low >= (0.275/np.pi)*(rho_l/rho_v)**(1/2) + 1: # low velocity region
    #     q_max = rho_v*V_flow*Dh_evap*coef_low
    # else: # high velocity region
    #     term_1 = ((rho_l/rho_v)**(3/4))/(169*np.pi)
    #     term_2 = ((rho_l/rho_v)**(1/2))/(19.2*np.pi*We_D**(1/3))
    #     q_max = rho_v*V_flow*Dh_evap*(term_1 + term_2)
    
    # if T_tube - T_sat <= 0:
    #     h_conv_boil = 0
    # elif T_tube - T_sat <= 0.5:
    #     h_conv_boil = q_max/(T_tube - T_sat) * (T_tube - T_sat)/0.5
    # else:
        
    # h_conv_boil = q_max/(T_tube - T_sat)
    
    return 
