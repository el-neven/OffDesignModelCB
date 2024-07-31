import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

class Pump_Wilo():
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
        self.A = None # ??
        self.B = None # ??
        self.Q_i_star = None # ??
        self.k_f = None # ??
        self.k_s = None # ??
        self.k_l = None # ??

        self.k_m= None # ??
        self.k_r = None # ??
        self.A_EC = None # ??
        self.B_EC = None # ??
        self.C_EC = None # ??
        self.D_EC = None # ??
        self.E_EC = None # ??
        self.F_EC = None # ??

        self.N_sim = None # Vitesse de similitude à laquelle les paramètres ont été calibré

    def inputs(self, su, N):
        self.su = su
        self.N_rot = N

        self.check_calculable()

    def check_calculable(self):
        self.Required_inputs = [self.su.m_dot, self.su.D, self.N_rot]
        if all(Input is not None for Input in self.Required_inputs):
                self.calculable = True

    def set_parameters(self, **kwargs):
            """
            Set parameters of the heat exchanger.
    
            Parameters
            ----------
            **kwargs : dict
                Key-value pairs representing parameters and their values.
                
                Example of call : heat_exchanger.set_parameters(D_o=0.05, t=0.002, H_core=2.0)
            """
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the parameters.")

            self.check_parametrized()

    def check_parametrized(self):
        self.Required_parameters = [self.A, self.B, self.Q_i_star, self.k_f, self.k_s, self.k_l, self.k_m, self.k_r, self.A_EC, self.B_EC, self.C_EC, self.D_EC, self.E_EC, self.F_EC]
        if all(Parameter is not None for Parameter in self.Required_parameters):
                self.parametrized = True

    def model_H_sim(self, Q):
        """
        Determines the head of the pump based on the volume flow rate.
        """
        Q_p = Q
        Q_l = self.k_l*Q_p #Loss in the pump
        self.Q_i = Q_p+Q_l #Total flow rate needed
        H_i = -self.A*self.Q_i+self.B
        epsilon_f = self.k_f*self.Q_i**2
        epsilon_s = self.k_s*((self.Q_i_star-self.Q_i)/self.Q_i_star)**2
        self.H_calc = H_i-epsilon_f-epsilon_s
        self.eta_hydro = (Q_p/self.Q_i)*(self.H_calc/H_i)
        self.P_hydro= 9810*Q_p*self.H_calc
        return (np.array([self.H_calc, self.Q_i, self.P_hydro, self.eta_hydro]))
    
    def model_shaft(self, P_hydro, eta_hydro, Q_i, Q_i_star_sh, N):
        """
        Model for the shaft power of the pump based on the volume flow rate and the frequency.
        """
        P_i = P_hydro/eta_hydro
        P_r = self.k_r*((Q_i_star_sh-Q_i)/Q_i_star_sh)**2 #Recirculation losses
        P_h = P_i+P_r
        P_m = self.k_m*(N/1000)**3 #disk friction losses
        self.P_shaft = P_h+P_m 
        self.eta_pump = P_hydro/self.P_shaft
        self.T_shaft = self.P_shaft/(N*(1/60)*2*np.pi)
        return (np.array([self.P_shaft]))
    
    def model_motor(self, P_shaft):
        """
        Model for the motor power of the pump based on the volume flow rate and the frequency.
        """
        "Zufen Wang"
        #print(P_shaft)
        self.eta_motor = self.A_EC*1e-5*(P_shaft**(1/3))**5 - self.B_EC*1e-4*(P_shaft**(1/3))**4 + self.C_EC*1e-3*(P_shaft**(1/3))**3 - self.D_EC*1e-2*(P_shaft**(1/3))**2 + self.E_EC*1e-1*(P_shaft**(1/3))**1 - self.F_EC*1e-1*(P_shaft**(1/3))**0
        self.P_motor = P_shaft/self.eta_motor
        return(np.array([self.P_motor, self.eta_motor]))

    def solve(self):
        m_dot = self.su.m_dot
        Q = m_dot/self.su.D #m^3/s	
        #print(Q)
        Q_sim = Q*(self.N_sim/self.N_rot)
        # print(Q_sim)
        H_sim = self.model_H_sim(Q_sim)
        self.H = H_sim[0]*(self.N_rot/self.N_sim)**2
        Q_i = H_sim[1]*(self.N_rot/self.N_sim)
        P_hydro = H_sim[2]*(self.N_rot/self.N_sim)**3
        eta_hydro = H_sim[3]
        Q_i_star = self.Q_i_star*(self.N_rot/self.N_sim)

        Calc_shaft = self.model_shaft(P_hydro, eta_hydro, Q_i, Q_i_star, self.N_rot)
        self.W_dot_sh = Calc_shaft[0]
        Calc_motor = self.model_motor(self.W_dot_sh)
        self.W_dot_motor = Calc_motor[0]
        #print(self.W_dot_motor)



         