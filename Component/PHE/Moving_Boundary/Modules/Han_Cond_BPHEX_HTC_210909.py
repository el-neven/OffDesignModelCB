# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:40:14 2021
% Condensation heat transfer correlation published by Han et. al in
% "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
@author: jvega
"""

import numpy as np

def Han_Cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, beta, L_v, N_cp, m_dot, D_p):
    """
    Inputs
    ------
    x : Vapor quality
    mu_l : Liquid viscosity
    k_l : Liquid thermal conductivity
    Pr_l : Liquid Prandtl Number
    rho_l : Liquid density
    rho_v : Vapor density
    G : Mass flux
    Dh : Plate spacing
    pitch_co : Corrugated pitch
    beta : Chevron angle
    L_v : Vertical length between fluid ports (268,2 mm) 
    N_cp : Number of canals
    m_dot : Flowrate
    
    Outputs
    -------
    
    
    Reference
    ---------
    "The caracteristics of condensation in brazed plate heat exchangers with different chevron angles", Journal of the Korean Physical Society, 2003
    Han et. al
    
    """
    
    # Preliminary calculations
    theta = np.pi/2 - beta
    g = 9.81 # gravity acceleration
    
    G_eq = G * ( (1-x) + x * (rho_l/rho_v)**0.5)
    Re_eq = G_eq*Dh/mu_l
    
    # Heat Transfer
    Ge1 = 11.22*(pitch_co/Dh)**(-2.83)*(theta)**(-4.5)
    Ge2 = 0.35*(pitch_co/Dh)**(0.23)*(theta)**(1.48)
    
    Nu = Ge1*Re_eq**Ge2*Pr_l**(1/3)
    h_cond = Nu*k_l/Dh
    
    # Pressure drop
    Ge3 = 3521.1*(pitch_co/Dh)**(4.17)*(theta)**(-7.75)
    Ge4 = -1.024*(pitch_co/Dh)**(0.0925)*(theta)**(-1.3)
    
    f = Ge3*Re_eq**Ge4
    
    # Two phase related pressure drop
    DP_tp = f*(L_v*N_cp/Dh)*G_eq**2*rho_l
    
    # Port pressure drop
    m_dot_eq = m_dot*(1 - x + x*(rho_l/rho_v)**0.5)
    G_p = 4*(m_dot_eq/(np.pi*D_p**2))
    rho_m = 1/( (x/rho_v) + (1 - x)/rho_l )

    DP_port = 1.4*G_p**2/(2*rho_m)

    # Static head loss
    DP_stat = -rho_m*g*L_v # negative because downward flow <-> Condenser

    # The acceleration pressure drop for condensation is expressed as : ??? 
    
    DP_tot = DP_tp + DP_port + DP_stat
    
    return h_cond, Nu, DP_tot