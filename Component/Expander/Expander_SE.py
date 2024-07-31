# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 11:03:15 2023

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

import numpy as np
from scipy.optimize import fsolve
import time

class Expander_Class:
    def __init__(self):

        self.calculable = False
        self.parametrized = False
        self.defined = False

        "Inputs"
        self.su = None
        self.ex = None

        self.N_rot = None
        self.T_amb = None

        "Parameters"
        self.AU_amb = None # global heat transfer coefficient for ambient heat losses      [W/K]
        self.AU_su_n = None # global heat transfer coefficient for supply heat transfer    [W/K]
        self.AU_ex_n = None # global heat transfer coefficient for exhaust heat transfer  	[W/K]
        self.d_su1 = None # pressure losses nozzle diameter                            	  [m]
        self.m_dot_n = None # nominal flow rate                                         	[kg/s]
        self.A_leak = None # nozzle cross section area for the leakage                     [m^2]
        self.W_dot_loss_0 = None # constant losses term                                    [W]
        self.alpha = None # proportional losses coefficient                                [-]
        self.C_loss = None # Couple losses                                                  [Nm]
        self.rv_in = None # built-in volumetric ratio                                  [-]
        self.V_s = None # machine displacement volume                                      [m^3]

    def inputs(self, su, ex, N_exp, T_amb):

        self.su = su
        self.ex = ex

        self.N_rot = N_exp
        self.T_amb = T_amb

        self.check_calculable()
        
    def check_calculable(self):

        self.Required_inputs = [self.su.p, self.su.h, self.ex.p, self.N_rot, self.T_amb]
        if all(Input is not None for Input in self.Required_inputs):
                self.calculable = True
        
    def check_parametrized(self):
        if self.AU_amb != None and self.AU_su_n != None and self.AU_ex_n != None and self.A_leak != None and self.d_su1 != None and self.W_dot_loss_0 != None and self.alpha != None and self.C_loss != None and self.rv_in != None and self.V_s != None and self.T_amb != None:
            self.parametrized = True

    def set_parameters(self, **parameters):
        for key in parameters:
            setattr(self, key, parameters[key])
        self.check_parametrized()
         
    #------------------------------------------------------------------------
    def System(self, x):
        
        # If the rotational speed is an input of the model while the mass flow rate is unknown
        self.m_dot, self.T_w = x
        self.N = self.N_rot/60
            
        #Boundary on the mass flow rate
        self.m_dot = max(self.m_dot, 1e-5)
        
        "1. Supply conditions: su"
        T_su = self.su.T
        P_su = self.su.p
        h_su = self.su.h
        s_su = self.su.s
        rho_su = self.su.D
        Fluid = self.su.fluid
        P_ex = self.ex.p
        self.P_ex = P_ex
        h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, Fluid)
        h_max = PropsSI('H', 'P', 4e6, 'T', 500, Fluid)
        T_sat_su = PropsSI('T', 'P', P_su, 'Q', 1, Fluid)

        if T_su<T_sat_su:
            print('----Warning the compressor inlet stream is not in vapor phase---')
        "2. Supply pressure drop: su->su1"
        h_su1 = h_su #Isenthalpic valve
        # Assumption that the density doesn't change
        
        A_su = np.pi*(self.d_su1/2)**2
        V_dot_su = self.m_dot/rho_su
        C_su = V_dot_su/A_su
        h_su_thr1 = h_su-(C_su**2)/2
        h_su_thr = max(h_su_thr1, h_ex_is)
        P_su_thr = PropsSI('P', 'S', s_su, 'H', h_su_thr, Fluid)
        
        P_su1 = max(P_su_thr, P_ex+1)
        self.DP_su = P_su-P_su1
        T_su1 = PropsSI('T', 'P', P_su1, 'H', h_su1, Fluid)
        
        "3. Cooling at the entrance: su1->su2"
        
        try:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cp_su1 = PropsSI('CPMASS', 'P', P_su1, 'Q', 0, Fluid)
        AU_su = self.AU_su_n*(self.m_dot/self.m_dot_n)**(0.8)
        C_dot_su = self.m_dot*cp_su1
        NTU_su = AU_su/C_dot_su
        epsilon_su1 = max(0,1-np.exp(-NTU_su))
        Q_dot_su = max(0, epsilon_su1*self.m_dot*cp_su1*(T_su1-self.T_w))
        
        h_su2 = min(h_max, max(max(h_ex_is, PropsSI('H','Q', 0.1, 'P', P_su1, Fluid)), h_su1 - Q_dot_su/self.m_dot))
        P_su2 = P_su1 #No pressure drop just heat transfer
        
        rho_su2 = PropsSI('D', 'P', P_su2, 'H', h_su2, Fluid)
        T_su2 = PropsSI('T', 'P', P_su2, 'H', h_su2, Fluid)
        s_su2 = PropsSI('S', 'P', P_su2, 'H', h_su2, Fluid)
        
        
        "4. Leakage"
        try:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'H', h_su1, Fluid)
        except:
            cv_su1 = PropsSI('CVMASS', 'P', P_su1, 'Q', 0, Fluid)
            
        #Essayer avec autre définition pour gamma!!!
        gamma = max(1e-2, cp_su1/cv_su1)
        P_crit = P_su2*(2/(gamma+1))**(gamma/(gamma-1))
        P_thr_leak = max(P_ex, P_crit)
        rho_thr_leak = PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)
        h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
        C_thr_leak = np.sqrt(2*(h_su2-h_thr_leak))
        m_dot_leak = self.A_leak*C_thr_leak*rho_thr_leak
        
        if self.su.m_dot == None:
            m_dot_in = self.N*self.V_s*rho_su2
            m_dot_leak_bis = self.m_dot-m_dot_in
        elif self.su.m_dot != None:
            m_dot_in = self.m_dot-m_dot_leak
            if self.N_rot == None:
                self.N = m_dot_in/(self.V_s*rho_su2)
        
        
        "5. Internal expansion"
        "Isentropic expansion to the internal pressure: su2->in"
        rho_in = rho_su2/self.rv_in
        #Trick to not have problems with CoolProp
        try:
            P_in = PropsSI('P', 'D', rho_in, 'S', s_su2, Fluid)
        except:
            delta = 0.0001
            P_in = 0.5*PropsSI('P', 'D', rho_in*(1+delta), 'S', s_su2, Fluid)+0.5*PropsSI('P', 'D', rho_in*(1-delta), 'S', s_su2, Fluid)
            
        h_in = PropsSI('H', 'D', rho_in, 'P', P_in, Fluid)
        w_in_s = h_su2-h_in
        
        "Expansion at constant volume: in->ex2"
        w_in_v = (P_in-P_ex)/rho_in
        h_ex2 = h_in-w_in_v
        
        "Total work"
        W_dot_in = m_dot_in*(w_in_s+w_in_v)
        
        
        "6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
        
        h_ex1 = max(min((m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot, h_su2), h_ex2)
        P_ex1 = P_ex
        T_ex1 = PropsSI('T', 'P', P_ex1, 'H', h_ex1, Fluid)
        
        try:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'H', h_ex1, Fluid)
        except:
            cp_ex2 = PropsSI('CVMASS', 'P', P_ex1, 'Q', 0, Fluid)
        
        AU_ex = self.AU_ex_n*(self.m_dot/self.m_dot_n)**(0.8)
        C_dot_ex = self.m_dot*cp_ex2
        NTU_ex=AU_ex/C_dot_ex
        epsilon_ex = max(0, 1-np.exp(-NTU_ex))
        Q_dot_ex = max(0, epsilon_ex*C_dot_ex*(self.T_w-T_ex1))
        
        self.h_ex = h_ex1+Q_dot_ex/self.m_dot
        
        
        "8. Energy balance"
        Q_dot_amb = self.AU_amb*(self.T_w-self.T_amb)
        
        W_dot_loss = self.alpha*W_dot_in + self.W_dot_loss_0 + self.C_loss*self.N*2*np.pi
        self.W_dot = W_dot_in - W_dot_loss

        "9. Performances"
        W_dot_s = self.m_dot*(h_su-h_ex_is)
        
        self.epsilon_is = self.W_dot/W_dot_s
        
        
        "10. Residuals"
        self.res_E = abs((Q_dot_su + W_dot_loss - Q_dot_ex - Q_dot_amb)/(Q_dot_su + W_dot_loss))
        self.res = self.res_E
        self.res_m_leak = abs((m_dot_leak_bis-m_dot_leak)/m_dot_leak)
        self.res = [self.res, self.res_m_leak]

        return self.res
           
                
    #------------------------------------------------------------------------
    def solve(self):
        
        if self.calculable and self.parametrized:
            
            start_time = time.time()
            
            ff_guess = [0.7, 1.2, 0.8, 1.3, 0.4, 1.7, 3] #guesses on the filling factor to provide suitable initial point for the iteration
            x_T_guess = [0.7, 0.95, 0.8, 0.99, 0.9] #For the iteration on the T_w
            stop = 0
        
            j = 0
            while not stop and j < len(x_T_guess):
                k = 0
                while not stop and k < len(ff_guess):
                    # Loop to permit multiple attempts to solve the implicit calculation 
                    m_dot_guess = ff_guess[k]*self.V_s*self.N_rot/60*PropsSI('D', 'P', self.su.p, 'H', self.su.h, self.su.fluid) #initial value for M_dot
                    T_w_guess = x_T_guess[j]*self.su.T+(1-x_T_guess[j])*self.T_amb #initial value for T_wall
                    #---------------------------------------------------------------------
                    args = ()
                    x = [m_dot_guess, T_w_guess]
                    #--------------------------------------------------------------------------
                    # ub = 2*x # upper bound for fsolve
                    try: 
                        fsolve(self.System, x, args = args)
                        res_norm = np.linalg.norm(self.res)
                    except:
                        res_norm = 1
                    
                    if res_norm < 1e-4: #flag > 0 
                        stop = 1
                    k = k + 1
                j = j + 1
            self.convergence = stop
            print(res_norm)
            elapsed_time = time.time() - start_time
            # print("Optimization time:", elapsed_time, "seconds")
            if self.convergence == 0:
                print("The system did not converge")
            if self.convergence == 1:
                # print("The system converged")
                self.ex.set_fluid(self.su.fluid)
                self.ex.set_m_dot(self.m_dot)
                self.ex.set_h(self.h_ex)
                self.ex.set_p(self.P_ex)
                self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")