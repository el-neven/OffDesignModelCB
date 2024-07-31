# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:38:04 

Modification w/r to previous version:
    -T_mean is arbitrarely set to the critical temperature in order to avoid errors.

@author: jvega
"""

import CoolProp.CoolProp as CP
import numpy as np
from scipy.interpolate import interp1d
import math

def conducticity_R1233zd(T, P):
    """ Richard A. Perkins and Marcia L. Huber

        "Measurement and Correlation of the Thermal Conductivity of trans-1-Chloro-3,3,3-trifluoropropene (R1233zd(E))"
        J. Chem. Eng. Data 2017, 62, 2659-2665
        Note that the critical enhancement contribution is not implemented.

        Same source as EES
       """
    # Dilute gas thermal conductivity
    T_c = 382.52 #[K]
    A_0 = -0.0103589
    A_1 = 0.0308929
    A_2 = 0.000230348

    k_0 = A_0 + A_1*(T/T_c) + A_2*(T/T_c)**2

    # Residual gas thermal conductivity
    try:
        Rho = CP.PropsSI('D', 'T', T, 'P', P, 'R1233zd(E)')
    except:
        Rho = CP.PropsSI('D', 'T', T, 'Q', 0, 'R1233zd(E)')
    Rho_c = 489.24 #[kg/m^3]

    B_11 = -0.0428296
    B_12 = 0.0434288
    B_21 = 0.0927099
    B_22 = -0.0605844
    B_31 = -0.0702107
    B_32 = 0.0440187
    B_41 = 0.0249708
    B_42 = -0.0155082
    B_51 = -0.00301838
    B_52 = 0.0021019

    Delta_kr = (B_11 + B_12*(T/T_c))*Rho/Rho_c + (B_21 + B_22*(T/T_c))*(Rho/Rho_c)**2 + (B_31 + B_32*(T/T_c))*(Rho/Rho_c)**3 + (B_41 + B_42*(T/T_c))*(Rho/Rho_c)**4 + (B_51 + B_52*(T/T_c))*(Rho/Rho_c)**5

    k = k_0 + Delta_kr
    return k




def PropsFluid(T_mean, P_mean, T_wall, fluid, oil_name = ''):
    
    if fluid == 'Oil':
        #-----------------------------------------------------------------------
        # Oil (therminol66) properties with respect to temperature
                
        if oil_name == "Therminol":
            T_f_val =   np.array([100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370]) # °C
            rho_f_val = np.array([955,948,941,934,928,921,914,907,899,892,885,878,870,863,856,848,840,832,825,817,809,800,792,783,775,766,757,748]) # kg/m^3
            cp_f_val =  np.array([1.84,1.87,1.91,1.94,1.98,2.01,2.05,2.09,2.12,2.16,2.19,2.23,2.27,2.30,2.34,2.38,2.42,2.45,2.49,2.53,2.57,2.61,2.65,2.69,2.73,2.77,2.81,2.85])*1e3 # J/(kg*K)
            k_f_val =   np.array([0.1135,0.1128,0.1121,0.1114,0.1107,0.1099,0.1091,0.1083,0.1074,0.1065,0.1056,0.1046,0.1036,0.1026,0.1015,0.1004,0.0993,0.0982,0.0970,0.0958,0.0946,0.0933,0.0920,0.0906,0.0893,0.0879,0.0865,0.0850]) # W/(m*K)
            mu_f_val =  np.array([3.60,2.92,2.42,2.05,1.75,1.52,1.33,1.18,1.06,0.95,0.86,0.784,0.718,0.661,0.611,0.567,0.529,0.495,0.464,0.437,0.413,0.391,0.371, 0.353,0.336,0.321,0.308,0.295])*1e-3 # Pa*s
        
        elif oil_name == "Pirobloc":
            T_f_val =   np.array([50,100,150,200,250,300]) # °C
            rho_f_val = np.array([856.16,823.87,791.59,759,726.42,692.07]) # kg/m^3
            cp_f_val =  np.array([1.985,2.186,2.382,2.583,2.779,2.985])*1e3 # J/(kg*K)
            k_f_val =   np.array([0.129,0.1244,0.1208,0.1174,0.114,0.11]) # W/(m*K)
            mu_f_val =  np.array([0.0147,0.00365,0.00144,0.000767,0.000476,0.000338]) # Pa*s
        
        else:
            print("Oil name not known")
            return 0
        
        f_rho = interp1d(T_f_val, rho_f_val)
        f_cp = interp1d(T_f_val, cp_f_val)
        f_k = interp1d(T_f_val, k_f_val)
        f_mu = interp1d(T_f_val, mu_f_val)    
        
        # Oil supply properties
        rho = f_rho(T_mean - 273.15)
        cp = f_cp(T_mean-273.15)
        k = f_k(T_mean-273.15)
        mu = f_mu(T_mean-273.15)
        Pr = (mu*cp)/k
        mu_wall = f_mu(T_mean-273.15)
        mu_rat = mu/mu_wall
        
        return mu, Pr, k, mu_wall, mu_rat, cp, rho
    if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
        k = conducticity_R1233zd(T_mean, P_mean)
        mu = CP.PropsSI('V', 'T', T_mean, 'P', P_mean, fluid)
        cp = CP.PropsSI('C', 'T', T_mean, 'P', P_mean, fluid)
        Pr = (mu*cp)/k
        mu_wall = CP.PropsSI("V", "T", T_wall, "P", P_mean, fluid)
        mu_rat = mu/mu_wall
        return mu, Pr, k, mu_wall, mu_rat, 0, 0
    else:
        #Force T_wall to be under TCRIT
        T_crit = CP.PropsSI("TCRIT", fluid)
        # T_mean = min(T_mean, T_crit)
        # T_wall = min(T_wall, T_crit)
        #-----------------------------------------------------------------------
        mu = CP.PropsSI('V',        'T', T_mean, 'P', P_mean, fluid)
        Pr = CP.PropsSI('Prandtl',  'T', T_mean, 'P', P_mean, fluid)
        k = CP.PropsSI('L',        'T', T_mean, 'P', P_mean, fluid)
        mu_wall = CP.PropsSI("V", "T", T_wall, "P", P_mean, fluid)
        mu_rat = mu/mu_wall
        return mu, Pr, k, mu_wall, mu_rat, 0, 0


# T = 60 + 273.15
# P = 2e5
# T_wall = 60 + 273.15
# fluid = "R1233zd(E)"

# mu, Pr, k, mu_wall, mu_rat, cp, rho = PropsFluid(T, P, T_wall, fluid)
# print(mu, Pr, k, mu_wall, mu_rat, cp, rho)