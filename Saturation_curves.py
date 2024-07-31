# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 14:35:37 2023

@author: Elise
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt


class TS_curve_generator:
    
    def __init__(self, Fluid):
        self.Fluid = Fluid
        self.TS_curve()
        
    def TS_curve(self):
        # Calculate critical pressure and minimum pressure
        P_crit = 0.99 * PropsSI('Pcrit', 'T', 0, 'P', 0, self.Fluid)
        P_min = PropsSI('P', 'T', 273.15 - 10, 'Q', 0, self.Fluid)

        # Create an array of pressures
        P_lin = np.linspace(P_min, P_crit, 100)
        
        # Initialize arrays to store data
        s_liq = np.zeros(100)
        T_liq = np.zeros(100)
        s_vap = np.zeros(100)
        T_vap = np.zeros(100)

        # Loop through the pressures and collect data
        for i in range(100):
            s_liq[i] = PropsSI('S', 'P', P_lin[i], 'Q', 0, self.Fluid)
            T_liq[i] = PropsSI('T', 'P', P_lin[i], 'Q', 0, self.Fluid)
            s_vap[i] = PropsSI('S', 'P', P_lin[i], 'Q', 1, self.Fluid)
            T_vap[i] = PropsSI('T', 'P', P_lin[i], 'Q', 1, self.Fluid)

        # Create the full TS curve by combining liquid and vapor data
        self.s_TS_curve = np.concatenate((s_liq, np.flip(s_vap)))
        self.T_TS_curve = np.concatenate((T_liq, np.flip(T_vap)))
        plt.plot(self.s_TS_curve, self.T_TS_curve, color='black', linestyle='--')
        plt.xlabel('Entropy (J/kg K)')
        plt.ylabel('Temperature (K)')
        # plt.show()

    def points(self, s_array, T_array):

        plt.plot(s_array, T_array, 'bo-', markersize=4)
        # plt.plot(self.s_TS_curve, self.T_TS_curve, color='blue')
        # plt.xlabel('Entropy (J/kg K)')
        # plt.ylabel('Temperature (K)')
        # plt.title('TS Curve')
        plt.show()






class PH_curve_generator:

    def __init__(self, Fluid):
        self.Fluid = Fluid
        self.PH_curve()

    def PH_curve(self):
        # Calculate critical pressure and minimum pressure
        P_crit = 0.99 * PropsSI('Pcrit', 'T', 0, 'P', 0, self.Fluid)
        P_min = PropsSI('P', 'T', 273.15 - 10, 'Q', 0, self.Fluid)

        # Create an array of pressures
        P_lin = np.linspace(P_min, P_crit, 100)
        
        # Initialize arrays to store data
        h_liq = np.zeros(100)
        h_vap = np.zeros(100)
        
        # Loop through the enthalpies and collect data
        for i in range(100):
            h_liq[i] = PropsSI('H', 'P', P_lin[i], 'Q', 0, self.Fluid)
            h_vap[i] = PropsSI('H', 'P', P_lin[i], 'Q', 1, self.Fluid)

        # Create the full PH curve by combining liquid and vapor data
        self.P_PH_curve = np.concatenate((P_lin, np.flip(P_lin)))
        self.h_PH_curve = np.concatenate((h_liq, np.flip(h_vap)))
        plt.plot(self.h_PH_curve, self.P_PH_curve, color='black', linestyle='--')
        plt.xlabel('Enthalpy (J/kg)')
        plt.ylabel('Pressure (Pa)')

    def Isotherm(self, T):
        # Calculate critical pressure and minimum pressure
        P_crit = 0.99 * PropsSI('Pcrit', 'T', 0, 'P', 0, self.Fluid)
        P_min = PropsSI('P', 'T', 273.15 - 10, 'Q', 0, self.Fluid)

        # Create an array of pressures
        P_lin = np.linspace(P_min, P_crit, 100)
        
        # Initialize arrays to store data
        h_liq = np.zeros(100)
        h_vap = np.zeros(100)
        
        # Loop through the enthalpies and collect data
        for i in range(100):
            h_liq[i] = PropsSI('H', 'P', P_lin[i], 'T', T, self.Fluid)
            h_vap[i] = PropsSI('H', 'P', P_lin[i], 'T', T, self.Fluid)

        # Create the full PH curve by combining liquid and vapor data
        #P_PH_curve = np.concatenate((P_lin, np.flip(P_lin)))
        self.h_T_curve = np.concatenate((h_liq, np.flip(h_vap)))
        

    def points(self, P_array, h_array):

        plt.plot(h_array, P_array, 'bo-', markersize=4)
        # plt.plot(self.s_TS_curve, self.T_TS_curve, color='blue')
        # plt.xlabel('Entropy (J/kg K)')
        # plt.ylabel('Temperature (K)')
        # plt.title('TS Curve')
        plt.show()

    









