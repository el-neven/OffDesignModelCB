# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 17:30:28 2024

@author: Elise
"""

import math
import numpy as np
from scipy.optimize import fsolve, least_squares
import scipy.special as bessel
import random
import csv
import matplotlib.pyplot as plt


class Dry_cooler:
    def __init__(self):

        self.calculable = False
        self.parametrized = False
        self.defined = False
        self.parameters = []

        "Initialize connectors"
        self.point_su1 = None
        self.point_su2 = None
        self.point_ex1 = None
        self.point_ex2 = None

        "Initialize calibration parameters"
        # The thermal resistances are based on theoretical models and need to be calibrated with experimental values
        self.C1 = None # Correction term for the resitance made by the convection in the tubes for the forced convection case (depends on the flow rate)
        self.C2 = None # Correction term for the resitance made by the conduction through the tube walls and the thermal resistance for the forced convection case
        self.C3 = None # Correction term for the resitance made by the convection in the tubes for the natural convection case (depends on the flow rate)
        self.C4 = None # Correction term for the resitance made by the conduction through the tube walls and the thermal resistance for the natural convection case

        "Initialize geometrical parameters"
        # Tubes data for one fictional exchanger
        self.nb_row = None # Number of row of tubes
        self.nb_column = None # Number of column of tubes
        self.D_ext = None # External diameter [m]
        self.th_wall = None # Wall thikness of the tubes [m]
        self.e_wall = None # Wall roughness [m] (hypothesis) 
        self.L_tube = None # Length of one tube [m]
        self.k_tube = None # Conductivity of steel [W/mK] (hypothesis on the material of the fins) 

        # Fins data for one fictional exchanger
        self.th_fin = None # Thickness of the fins [m] 
        self.p_fin = None # space between fins [m] (evaluated by counting the number of fins on 4cm and applying p_fin = (4-19*th_fin)/(nb_fins_4-1))
        self.h_fin = None # Height of the fins [m] (hypothesis, because two tubes of different temperature use the same fin, we divided the fin hight by two. 
        # It derives from the hypothesis that the fin is high enough so that the two heat transfers are decoupled)
        self.l_fin = None # Length of the fins [m]
        self.k_alu = None # Conductivity of aluminium [W/mK] (hypothesis on the material of the fins)

    def set_inputs(self, point_su1, point_su2, point_ex1, point_ex2, fan):
        self.point_su1 = point_su1
        self.point_su2 = point_su2
        self.point_ex1 = point_ex1
        self.point_ex2 = point_ex2

        self.water = point_su1.fluid
        self.Tw = point_su1.T # Water temperature [K]
        self.Qw = point_su1.V_dot # Water flow rate [m^3/s]

        self.air = point_su2.fluid
        self.Ta = point_su2.T # External air temperature [K]
        self.Qa = point_su2.V_dot # Air flow rate [m^3/s]

        self.fan = fan # List of boolean indicating if the fan is working or not

        self.check_calculable()

    def check_calculable(self):
        self.Required_inputs = [self.point_su1.T, self.point_su1.V_dot, self.point_su2.T, self.point_su2.V_dot]
        if all(Input is not None for Input in self.Required_inputs):
            self.calculable = True

    def check_parametrized(self):
        if self.nb_row is None and self.nb_column is None and self.D_ext is None and self.th_wall is None and self.e_wall is None and self.L_tube is None and self.k_tube is None and self.th_fin is None and self.p_fin is None and self.h_fin is None and self.l_fin is None and self.k_alu is None:
            self.parametrized = False
        else:
            self.parametrized = True

    def set_geometrical_parameters(self, **kwargs):
        """
        Set the geometrical parameters of the Dry cooler model.

        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                self.parameters.append((key, value))  # Append the parameter name and value as a tuple
            else:
                print(f"Warning: Parameter '{key}' not found in the parameters.")

        self.check_parametrized()
        if self.parametrized:
            self.calculate_geometrical_parameters()

        #return parameters  # Return the list of parameter values

    def calculate_geometrical_parameters(self):
        """
        Calculate the geometrical parameters of the Dry cooler model.

        """
        self.D_int = self.D_ext - 2*self.th_wall # Internal diameter [m]
        self.parameters.append(('D_int', self.D_int))
        self.R_int = self.D_int/2 # Internal radius [m]
        self.parameters.append(('R_int', self.R_int))
        self.L_tot = self.L_tube*self.nb_row*self.nb_column # Total cross sectional area of the tubes [m^2]
        self.parameters.append(('L_tot', self.L_tot))
        self.Ac_tot = self.nb_row*self.nb_column*math.pi*(self.D_ext/2-self.th_wall)**2 # Total cross sectional area of the tubes [m^2]
        self.parameters.append(('Ac_tot', self.Ac_tot))
        self.nb_fins = self.L_tube/(self.p_fin+self.th_fin) # Number of fins
        self.parameters.append(('nb_fins', self.nb_fins))
        self.A_fin = 2*self.l_fin*self.h_fin-self.nb_column*self.nb_row*math.pi*(self.D_ext/2)**2 # Area of a fin [m^2]
        self.parameters.append(('A_fin', self.A_fin))
        self.A_fins_tot = self.A_fin*self.nb_fins # Total area of the fins [m^2]
        self.parameters.append(('A_fins_tot', self.A_fins_tot))
        self.A_unfin = self.nb_row*self.nb_column*(self.L_tube-self.nb_fins*self.th_fin)*math.pi*self.D_ext # Unfined area of the pipe [m^2]
        self.parameters.append(('A_unfin', self.A_unfin))
        self.A_tot = self.A_fins_tot + self.A_unfin # Total area of the fins plus the unfined area of the pipe [m^2]
        self.parameters.append(('A_tot', self.A_tot))
        self.r_fin_eff = math.sqrt(self.A_fins_tot/(2*self.nb_fins*math.pi+(self.D_ext/2)**2))
        self.parameters.append(('r_fin_eff', self.r_fin_eff ))

    def set_calibration_parameters(self, **kwargs):
        """
        Set the calibration parameters of the Dry cooler model.

        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                self.parameters.append((key, value))

        if self.C1 is None:
            self.C1 = 1
            self.parameters.append(('C1', self.C1))

        if self.C2 is None:
            self.C2 = 0 
            self.parameters.append(('C2', self.C2))
        
        if self.C3 is None:
            self.C3 = 1
            self.parameters.append(('C3', self.C3))

        if self.C4 is None:
            self.C4 = 0
            self.parameters.append(('C4', self.C4))

    def solve(self):
        """
        Solve the Dry cooler model.
        
        """

        if self.calculable and self.parametrized:

            "Creation of the exchangers"
            Hx1 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx2 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx3 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx4 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx5 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx6 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx7 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)
            Hx8 = self.Exchanger(self.Tw, self.Ta, self.Qw, self.parameters)

            "Link between the exchangers (Water side)" # One water flow
            Hx1.set_next_exch_water(Hx2)
            Hx2.set_next_exch_water(Hx3)
            Hx3.set_next_exch_water(Hx4)
            Hx4.set_next_exch_water(Hx5)
            Hx5.set_next_exch_water(Hx6)
            Hx6.set_next_exch_water(Hx7)
            Hx7.set_next_exch_water(Hx8)

            "Link the heat exchangers to the fans (Air side)" # 4 Air flows
            Hx8.set_next_exch_air(Hx1, fan=self.fan[0]) #Air flow passing throught Hx8 passes through Hx1 and Fan1
            Hx7.set_next_exch_air(Hx2, fan=self.fan[1])
            Hx6.set_next_exch_air(Hx3, fan=self.fan[2])
            Hx5.set_next_exch_air(Hx4, fan=self.fan[3])

            self.Hxs = [Hx1, Hx2, Hx3, Hx4, Hx5, Hx6, Hx7, Hx8]

            "Evaluate the exchangers operating point"
            fsolve(self.equations, [self.Tw], args=(Hx1, Hx8)) #Guess: Tw_out = Tw_in

            self.point_ex1.set_fluid(self.point_su1.fluid)
            self.point_ex1.set_V_dot(self.point_su1.V_dot)
            self.point_ex1.set_T(Hx8.Tw_out)
            self.defined = True

        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")

    def equations(self, x, Hx1, Hx8):
        Hx1.compute_T_out() # Compute the temperature of the water at the outlet of the first exchanger (to strat the iteration)
        return [x[0]-Hx8.Tw_out] # Iterates until the Tw_out guess = Tw_calculated
    
    class Exchanger:

        def __init__(self, Tw_in, Ta, Qw, parameters):
            self.fan = None
            parameters = dict(parameters)
            self.parameters = parameters
            
            self.Tw_in = Tw_in
            self.Qw = Qw
            self.Tw_out = None
            
            self.Ta = Ta
            self.Ta_in = Ta
            self.Ta_out = None
            
            self.q_dot = None
            
            self.fin_link = None
            self.next_exch_water = None
            self.next_exch_air = None
            
            self.h_bar_out = None
            self.q_dot_one_fin_estimate = None
            self.m_dot_estimate = None
            
            self.R_natural_out = None

            "Fluid data (water)"
            self.rho_w = 1000  # Density [kg/m^3]
            self.u_w = 8.9 * 10 ** -4  # Dynamic viscosity at 25°C [Pa*s]
            self.k_w = 0.598  # Thermal conductivity [W/mK]
            self.cp_w = 4180  # Specific heat [J/kgK]

            "Fluid data (air)"
            self.rho_a = 1.225  # Density [kg/m^3]
            self.k_a = 25.5 * 10 ** -3  # Thermal conductivity of air at 15°C
            self.cp_a = 1007  # Specific heat of air at 15°C [J/kgK]

            "Flow data (water)"
            self.m_dot_w = Qw * self.rho_w  # Mass flow rate [kg/s]
            self.v_tube = Qw / (parameters['Ac_tot'])  # Velocity of the water in the tube [m/s]
            
            "Flow data (air)"
            self.Qa = 106500 / 3600 / 4  # Volume flow rate [m^3/s] FROM DATASHEET
            self.m_dot_a = self.Qa * self.rho_a  # Mass flow rate of the air through one fictional exchanger [kg/s]

            "Flow characteristics (water)"
            self.Re = self.rho_w * self.v_tube * parameters['D_int'] / self.u_w  # Reynold's number
            self.Pr = self.cp_w * self.u_w / self.k_w  # Prantl's number
            self.f = (-2 * math.log(2 * parameters['e_wall'] / (7.54 * parameters['D_int']) - 5.02 / self.Re * math.log(2 * parameters['e_wall'] / (7.54 * parameters['D_int']) + 13 / self.Re, 10), 10)) ** -2  # Friction factor, Zigrang and Sylvester formula (turbulent flow)
            if self.Pr > 0.5 and self.Pr < 6000 and self.Re > 2300 and self.Re < 5 * 10 ** 6:
                self.Nu = (self.f / 8 * (self.Re - 1000) * self.Pr) / (1 + 12.7 * (self.Pr ** (2 / 3) - 1) * math.sqrt(self.f / 8))  # Gnielinski, 0.5<Pr<6000 and 2300<Re<5*10**6
            elif self.Re < 2300:  # Flow is laminar
                Gz = parameters['D_int'] * self.Re * self.Pr / parameters['L_tube']  # Greatz number
                self.Nu = 3.66 + (0.049 + 0.02 / self.Pr) * Gz ** 1.12 / (1 + 0.065 * Gz ** 0.7)
            else:
                print("No correlation for the Nu of that flow")
            
            self.h_in = self.Nu * self.k_w / parameters['D_int']
            self.C_dot_w = self.cp_w * self.m_dot_w

        def set_next_exch_water(self, exch):
            """
            Sets the next heat exchanger for water flow.

            Parameters:
            exch (HeatExchanger): The next heat exchanger object for water flow.

            """
            self.next_exch_water = exch

        def set_next_exch_air(self, exch, fan=True): #Fan set on by default
            """
            Sets the next heat exchanger for air flow.

            """
            self.next_exch_air = exch # Next Heat exchanger the air will go through
            self.set_fin_link(exch, fan) # Air flow passes through the fan after the heat exchanger
    
        # Method useful when natural convection is on play.
        def set_fin_link(self, exch, fan):
            """
            Links the Heat exchangers to the fan.

            """
            self.fin_link = exch
            self.fan = fan
            if exch.fin_link == None:
                exch.set_fin_link(self, fan) # Set a reversible link between the two exchangers (Hx1 is connected to Hx8 and Hx8 is connected to Hx1)
            self.fin_link.fan = fan #Probably not necessary

        def set_Tw_in(self, Tw_in):
            self.Tw_in = Tw_in

        def set_Ta_in(self, Ta_in):
            self.Ta_in = Ta_in

        def compute_T_out(self):
            if not self.fan:
                self.natural_conv()
            else:
                self.forced_conv()
            self.q_dot = self.eps*self.C_dot_min*(self.Tw_in - self.Ta_in)
            self.Tw_out = self.Tw_in - self.q_dot/self.C_dot_w
            self.Ta_out = self.Ta_in + self.q_dot/self.C_dot_a
            
            if self.next_exch_air != None:
                self.next_exch_air.set_Ta_in(self.Ta_out)
                
            if self.next_exch_water != None:
                self.next_exch_water.set_Tw_in(self.Tw_out)
                self.next_exch_water.compute_T_out()

        def forced_conv(self):
            "Flow characteristics (air)"
            self.C_dot_a = self.cp_a * self.m_dot_a

            "Thermal resistance calculation for one fictional exchanger"
            self.R_in = 1 / (self.h_in * math.pi * self.parameters['D_int'] * self.parameters['L_tot'])  # Convection in the tubes
            self.R_cond = math.log(self.parameters['D_ext'] / self.parameters['D_int']) / (2 * self.parameters['k_tube'] * math.pi * self.parameters['L_tot'])  # Conduction through the tube walls
            self.R_out = 0.00166  # Thermal resistance from EES code, in forced condition, with full air draw.
            
            "Heat transfer coefficient calculation for forced convection"
            R_tot = self.parameters['C1'] * self.R_in + self.R_cond + self.R_out + self.parameters['C2']
            UA = 1 / R_tot

            "Evaluation of exchanger characteristics via the NTU-epsilon method for cross-flow exchanger, both fluid unmixed"
            self.C_dot_min = min(self.C_dot_w, self.C_dot_a)
            C_dot_max = max(self.C_dot_w, self.C_dot_a)
            self.C_R = self.C_dot_min / C_dot_max
            NTU = UA / self.C_dot_min
            self.eps = 1 - math.exp((math.exp(-NTU * self.C_R * NTU ** -0.22) - 1) / (self.C_R * NTU ** -0.22))

        def natural_conv(self): # based on the book of Sam
            self.Ts = (self.Tw_in + self.fin_link.Tw_in)/2 # Fin temperature. Big hypothesis: the fin temperature is the mean between the water inlet temperature bewteen the two heat exchangers (ie Hx1 and Hx 8)
            S = self.parameters['p_fin']
            g = 9.81 # Gravitational constant
            beta = 1/self.Ts # Thermal expension coefficient of a perfect gas
            mu_a = 1.81*10**-5 # Dynamic viscosity of air at 15°C
            v_a = mu_a/self.rho_a # Kinematic viscosity of air at 15°C
            alpha = 21.7*10**-6 # Thermal diffusivity of air at 20°C
            
            Ra = g*S**3*beta*(self.Ts-self.Ta)/(v_a*alpha)
            Nu = Ra/24*S/(self.parameters['h_fin']*2)*(1-math.exp(-35*self.parameters['h_fin']*2/(Ra*S)))**0.75
            self.h_bar_out = Nu*self.k_a/S
            
            eta_fin = self.eta_fin_annular_rect(self.parameters['th_fin'], self.parameters['D_int']/2, self.parameters['r_fin_eff'], self.h_bar_out, self.parameters['k_alu'])
            self.eta_o = 1-(self.parameters['A_fins_tot']/self.parameters['A_tot'])*(1-eta_fin)
            
            # Try to estimate C_dot_min
            Ta_out_estimate = (self.Tw_in+self.Ta_in)/2 # Big hypothesis which is prbably far from true
            self.q_dot_one_fin_estimate = self.eta_o*self.h_bar_out*(self.parameters['A_fin'] + (self.parameters['nb_column']+self.parameters['nb_row'])*self.parameters['D_ext']*math.pi)*(self.Tw_in-self.Ta)
            self.m_dot_estimate = self.q_dot_one_fin_estimate/(self.cp_a*(Ta_out_estimate-self.Ta))*self.parameters['nb_fins']
            self.C_dot_a = self.m_dot_estimate*self.cp_a
            self.C_dot_min = min(self.C_dot_a, self.C_dot_w)
            C_dot_max = max(self.C_dot_w, self.C_dot_a)
            self.C_R = self.C_dot_min / C_dot_max
            
            "Thermal resistance calculation for one fictional exchanger"
            self.R_in = 1 / (self.h_in * math.pi * self.parameters['D_int'] * self.parameters['L_tot'])  # Convection in the tubes
            self.R_cond = math.log(self.parameters['D_ext'] / self.parameters['D_int']) / (2 * self.parameters['k_tube'] * math.pi * self.parameters['L_tot'])  # Conduction through the tube walls
            self.R_natural_out = 1/(self.eta_o*self.h_bar_out*self.parameters['A_tot'])
            R_tot = self.parameters['C3']*self.R_in + self.R_cond + self.R_natural_out + self.parameters['C4']
            UA = 1/R_tot
            NTU = UA/self.C_dot_min
            self.eps = 1-math.exp((math.exp(-NTU*self.C_R*NTU**-0.22)-1)/(self.C_R*NTU**-0.22))


        def eta_fin_annular_rect(self, th_fin, r_in, r_out, h, k_alu):
            mro = r_out*math.sqrt(2*h/(k_alu*th_fin))
            rio = r_in/r_out
            eta_fin = 2*rio*(bessel.kv(1, mro*rio)*bessel.iv(1, mro)-bessel.iv(1, mro*rio)*
                            bessel.kv(1, mro))/(bessel.kv(0, mro*rio)*bessel.iv(1, mro)+
                                                bessel.iv(0, mro*rio)*bessel.kv(1, mro))/(mro*(1-(rio)**2))
            return eta_fin



