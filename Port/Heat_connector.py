# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:30:09 2024

@author: Elise
"""


class Heat_connector:

    def __init__(self):
        """
        Parameters
        ----------
        /

        """

        self.completely_known = False
        self.state_known = False
        self.variables_input = []
        
        # self.fluid = None
        # self.m_dot = None           # Mass flow rate [kg/s]
        self.T_1 = None             # Temperature of the hot/cold fluid [K]
        self.T_2 = None             # Temperature of the cold/hot fluid [K]
        # self.h_1 = None             # Spec. enthalpy of the hot/cold fluid[J/kg]
        # self.h_2 = None             # Spec. enthalpy of the cold/hot fluid[J/kg]
        self.AU = None              # Heat transfer coefficient [W/K]
        self.Q_dot = None           # Heat transfer rate [W]

    def check_completely_known(self):

        if self.T_1 != None and self.T_2 != None and self.AU != None:
            self.calculate_heat_trasfer()
            self.completely_known = True
            print("Heat transfer calculated")
        else:
            pass

    
    def calculate_heat_trasfer(self):
        self.Q_dot = self.AU*(self.T_1 - self.T_2)


    def set_T1(self, value):
        self.T_1 = value
        self.check_completely_known()

    def set_T2(self, value):
        self.T_2 = value
        self.check_completely_known()

    def set_AU(self, value):
        self.AU = value
        self.check_completely_known()

    