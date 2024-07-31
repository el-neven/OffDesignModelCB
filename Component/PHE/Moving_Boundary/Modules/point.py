#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI

class Point_on_cycle:
    def __init__(self):
        """
        Parameters
        ----------
        /

        """

        self.completely_known = False
        self.state_known = False
        self.variables_input = []
        
        self.fluid = None
        self.m_dot = None           # Mass flow rate [kg/s]
        self.T = None               # Temperature [K]
        self.p = None               # Pressure [Pa]
        self.h = None               # Spec. enthalpy [J/kg]
        self.s = None               # Spec. entropy [J/kg/K]
        self.D = None               # Mass density [kg/m^3]
        self.x = None               # Quality [kg/kg]
        
        
    def check_completely_known(self):
        if self.fluid != None:
            if len(self.variables_input)>2:
                print("Error: Too many state variables")
            elif len(self.variables_input)<2:
                pass
            elif self.state_known:
                pass
            else:
                try:
                    self.T = PropsSI('T', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.p = PropsSI('P', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.h = PropsSI('H', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.s = PropsSI('S', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.D = PropsSI('D', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.x = PropsSI('Q', self.variables_input[0][0], self.variables_input[0][1], self.variables_input[1][0], self.variables_input[1][1], self.fluid)
                    self.state_known = True
                except:
                    print("Error: This pair of inputs is not yet supported")
        else:
            pass
        
        if self.m_dot != None and self.state_known and not self.completely_known:
            self.completely_known = True
            print("Point completely known")
        else:
            pass


    def set_fluid(self, value):
        if self.fluid != None:
            print("Error: Variable already defined")
        else:
            try:
                PropsSI('M', value)
                self.fluid = value
                self.check_completely_known()
            except:
                print("Error: Incorect fluid name")
        
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
            self.variables_input = self.variables_input+[['T',value]]
            self.check_completely_known()
        
    def set_p(self, value):
        if self.p != None:
            print("Error: Variable already defined")
        else:
            self.p = value
            self.variables_input = self.variables_input+[['P',value]]
            self.check_completely_known()
        
    def set_h(self, value):
        if self.h != None:
            print("Error: Variable already defined")
        else:
            self.h = value
            self.variables_input = self.variables_input+[['H',value]]
            self.check_completely_known()
        
    def set_s(self, value):
        if self.s != None:
            print("Error: Variable already defined")
        else:
            self.s = value
            self.variables_input = self.variables_input+[['S',value]]
            self.check_completely_known()
        
    def set_D(self, value):
        if self.D != None:
            print("Error: Variable already defined")
        else:
            self.D = value
            self.variables_input = self.variables_input+[['D',value]]
            self.check_completely_known()
            
    def set_x(self, value):
        if self.x != None:
            print("Error: Variable already defined")
        else:
            self.x = value
            self.variables_input = self.variables_input+[['Q',value]]
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
        
        
        
        
        
        
        
        
        

