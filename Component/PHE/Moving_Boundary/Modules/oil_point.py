#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class Oil_Point_on_cycle:
    def __init__(self, name):
        """
        Parameters
        ----------
        /

        """

        self.completely_known = False
        self.state_known = False
        
        self.name = name
        self.fluid = "Oil"
        self.m_dot = None           # Mass flow rate [kg/s]
        self.T = None               # Temperature [K]
        self.p = None               # Pressure [Pa]
        self.h = None               # Spec. enthalpy [J/kg]
        self.D = None               # Mass density [kg/m^3]
        self.x = 0                  # Quality [kg/kg]
        
        if self.name == 'Therminol':
        
            T_f_val =   np.array([100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370]) # °C
            rho_f_val = np.array([955,948,941,934,928,921,914,907,899,892,885,878,870,863,856,848,840,832,825,817,809,800,792,783,775,766,757,748]) # kg/m^3
            cp_f_val =  np.array([1.84,1.87,1.91,1.94,1.98,2.01,2.05,2.09,2.12,2.16,2.19,2.23,2.27,2.30,2.34,2.38,2.42,2.45,2.49,2.53,2.57,2.61,2.65,2.69,2.73,2.77,2.81,2.85])*1e3 # J/(kg*K)
                    
            self.f_rho = interp1d(T_f_val, rho_f_val)
            self.f_cp = interp1d(T_f_val, cp_f_val)
            
        elif self.name == 'Pirobloc': # Pirobloc HTF Basic : http://www.coolprop.org/_downloads/44c89ee563242455aa99c8c63aa3ce94/PBB_fitreport.pdf
                        
            T_f_val =   np.array([50,100,150,200,250,300]) # °C
            rho_f_val = np.array([847.2,814.3,782.3,748,715,682.6]) # kg/m^3
            cp_f_val =  np.array([1.985,2.186,2.382,2.583,2.779,2.985])*1e3 # J/(kg*K)
                    
            self.f_rho = interp1d(T_f_val, rho_f_val)
            self.f_cp = interp1d(T_f_val, cp_f_val)    
            
        else:
            print("Oil name not known")
            
        
    def check_completely_known(self):
        if self.fluid != None and self.T != None and self.p != None:
            self.state_known = True
        else:
            pass
        
        if self.m_dot != None and self.state_known and not self.completely_known:
            self.completely_known = True
            print("Point completely known")
        else:
            pass

        
    def set_m_dot(self, value):
        if self.m_dot != None:
            print("Error: Variable already defined")
        else:
            self.m_dot = value
            self.check_completely_known()
        
    def set_T(self, value):
        if self.T != None:
            print("Error: Variable already defined")
        else:
            self.T = value
            
            cp = self.f_cp(value-273.15) 
            self.h = cp*value
            self.D = self.f_rho(value-273.15) 
            self.check_completely_known()

    def set_p(self, value):
        if self.p != None:
            print("Error: Variable already defined")
        else:
            self.p = value
            self.check_completely_known()
            
    def print_resume(self, unit_T='K', unit_p='Pa'):
        """
        Parameters
        ----------
        unit_T = Temperature unit: 'K' or 'C'
        unit_p = Temperature unit: 'Pa' or 'bar'
        
        """
        
        print("Fluid: " + self.fluid + "")
        print("Mass flow rate: " + str(self.m_dot) + "[kg/s]")
        
        if unit_T == 'K':
            print("Temperature: " + str(self.T) + "[K]")
        elif unit_T == 'C':
            print("Temperature: " + str(self.T-273.15) + "[°C]")
        else:
            print("Error: Wrong argument unit_T in the method print_resume")

        if unit_p == 'Pa':
            print("Pressure: " + str(self.p) + "[Pa]")
        elif unit_p == 'bar':
            print("Pressure: " + str(self.p/1e5) + "[bar]")
        else:
            print("Error: Wrong argument unit_p in the method print_resume")
            
        
        print("Spec. enthalpy: " + str(self.h) + "[J/kg]")
        
        print("Spec. entropy: " + str(self.s) + "[J/kg/K]")
        
        print("Mass density: " + str(self.D) + "[kg/m^3]")
        print("Quality: " + str(self.x) + "[-]")
        
        
        
        
        
        
        
        
        

