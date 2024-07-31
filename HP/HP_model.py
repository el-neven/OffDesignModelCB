
from CoolProp.CoolProp import PropsSI
from scipy.optimize import least_squares
from scipy.optimize import root
from scipy.optimize import minimize
import numpy as np

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from Component.Compressor.Compressor_SE import Compressor_Class
from Component.PHE.Moving_Boundary.HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger
from Component.PHE.Simple_MB.HX_MovingBoundary import HeatExchanger
from Component.Circulateur.Simulation_model import Pump_Wilo
from Saturation_curves import TS_curve_generator
from Saturation_curves import PH_curve_generator

from Port.Mass_connector import Mass_connector
import time
import matplotlib.pyplot as plt

"Definir le cycle ici!!!"
class HP():
    def __init__(self):
        "Design parameters"
        self.T_sh = None # global heat transfer coefficient for supply heat transfer   [W/K]
        self.T_sc = None

    def inputs(self, N_comp, T_amb, wf, T_w_cd_su, glide_cd, T_w_ev_su, glide_ev, Nb_comp):
        self.N_comp = N_comp
        self.T_amb = T_amb
        self.wf = wf
        self.T_w_cd_su = T_w_cd_su
        self.glide_cd = glide_cd
        self.T_w_ev_su = T_w_ev_su
        self.glide_ev = glide_ev
        self.Nb_comp = Nb_comp

    def set_parameters(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Parameter '{key}' not found in the parameters.")
    
    def system(self, x):
        
        "GUESS"
        self.P_cd = x[0]
        self.P_ev = x[1]
        T_ev = PropsSI('T', 'P', self.P_ev, 'Q', 0.5, 'R1233ZDE')
        T_cd = PropsSI('T', 'P', self.P_cd, 'Q', 0.5, 'R1233ZDE')
        self.rp = self.P_cd/self.P_ev
        # print('P_ev', self.P_ev, 'P_cd', self.P_cd)

        "---------------------------------------1. Compressor---------------------------------------"
        "Define ports"
        self.su_comp = Mass_connector()
        self.ex_comp = Mass_connector()

        # Supply state (su_exp)
        T_ev = PropsSI('T', 'P', self.P_ev, 'Q', 0.5, 'R1233ZDE')
        T_cd = PropsSI('T', 'P', self.P_cd, 'Q', 0.5, 'R1233ZDE')
        self.su_comp.set_fluid('R1233ZDE')
        self.su_comp.set_T(T_ev+self.T_sh)
        self.su_comp.set_p(self.P_ev)

        # Exhaust state (ex_exp)
        self.ex_comp.set_p(self.P_cd)

        "Define class"
        COMP = Compressor_Class()
        COMP.inputs(self.su_comp, self.ex_comp, self.N_comp, self.T_amb)

        COMP.set_parameters(**{
            'AU_amb': 2.93416093e+00,
            'AU_su_n': 1.79782166e+01,
            'AU_ex_n': 9.56055293e+00,
            'd_ex': 9.79880377e-03,
            'm_dot_n': 0.1,
            'A_leak': 3.63305728e-07,
            'W_dot_loss_0': 4.29061955e-01,
            'alpha': 3.17325865e-03,
            'C_loss': 4.03523425e-01,
            'rv_in': 1.7,
            'V_s': 1.21e-4,
        })

        COMP.solve()
        self.W_dot_comp = self.Nb_comp*COMP.W_dot
        self.epsilon_is = COMP.epsilon_is
        # print('Convergence', COMP.convergence)

        "---------------------------------------2. Condenser---------------------------------------"

        "Define ports"
        self.ex_cd = Mass_connector()

        su_cd_w = Mass_connector()
        ex_cd_w = Mass_connector()

        # Supply state R1233ZDE (su_cd)
        self.su_cd = self.ex_comp # Copy of the outlet (need to find another way to do it)
        self.su_cd.set_m_dot(self.ex_comp.m_dot*self.Nb_comp)

        # Exhaust state (ex_cd)
        h_init_cd = PropsSI('H', 'T', T_cd-self.T_sc, 'P', self.P_cd, 'R1233ZDE')

        # Supply state (water)
        su_cd_w.set_fluid('Water')
        su_cd_w.set_T(self.T_w_cd_su)
        su_cd_w.set_p(2e5)

        #Exhaust state (water)
        ex_cd_w.set_fluid('Water')
        ex_cd_w.set_T(self.T_w_cd_su+self.glide_cd)
        ex_cd_w.set_p(2e5)

        "Define class"
        self.COND = HeatExchanger(self.su_cd, su_cd_w, self.ex_cd, ex_cd_w)

        self.COND.set_parameters(**{
            'HX_type': 'condenser',
            'HX_D': 0.06, # Inlet port diameter
            'HX_A': 17.8, # Exchange surface
            'min_pinch': 2,
            'dT_sub_or_sup': 5
        })

        self.COND.solve()
        self.Q_dot_cd = self.COND.Q
        self.m_dot_w_cd = ex_cd_w.m_dot

        "---------------------------------------3. Expansion valve---------------------------------------"
        "Define ports"
        self.ex_vlv = Mass_connector()

        # Supply state
        self.su_vlv = self.ex_cd

        # Exhaust state
        self.ex_vlv.set_fluid('R1233ZDE')
        self.ex_vlv.set_p(self.P_ev)
        self.ex_vlv.set_h(self.ex_cd.h)
        self.ex_vlv.set_m_dot(self.ex_cd.m_dot)
        # print(self.ex_vlv.h, 'h')

        # if self.ex_vlv.h > PropsSI('H', 'P', self.P_ev, 'Q', 0, self.wf):
        #     # 2P        
        #     self.ex_vlv.set_T(PropsSI('T','P',self.P_ev,'Q',0,self.wf))
        #     self.ex_vlv.set_x(PropsSI('Q', 'P', self.P_ev, 'H', self.ex_vlv.h, self.wf))
        # else:
        #     # 1P
        #     T_eva_in = refpropm('T','H',h_eva_in,'P',P_eva_guess,fluid_wf);
        #     x_eva_in = -1
        
        # print(self.ex_vlv.x, 'x')

        "---------------------------------------4. Evaporator---------------------------------------"

        "Define ports"
        self.ex_ev = Mass_connector()

        su_ev_w = Mass_connector()
        ex_ev_w = Mass_connector()

        # Supply state (R1233ZDE)
        self.su_ev = self.ex_vlv

        # Exhaust state (ex_ev)
        h_init_ev = PropsSI('H', 'T', T_ev+self.T_sh, 'P', self.P_ev, 'R1233ZDE')

        # Supply state (water)
        su_ev_w.set_fluid('Water')
        su_ev_w.set_T(self.T_w_ev_su)
        su_ev_w.set_p(2e5)

        #Exhaust state (water)
        ex_ev_w.set_fluid('Water')
        ex_ev_w.set_T(self.T_w_ev_su-self.glide_ev)
        ex_ev_w.set_p(2e5)
        # print(self.su_ev.T, 'T_in_ev')
        # print(self.su_ev.x, 'x_in_ev')
        "Define class"
        self.EVAP = HeatExchanger(self.su_ev, su_ev_w, self.ex_ev, ex_ev_w)

        self.EVAP.set_parameters(**{
            'HX_type': 'evaporator',
            'HX_D': 0.06, #Diamètre de port d'entré
            'HX_A': 17.8, #Surface d'échange
            'min_pinch': 2,
            'dT_sub_or_sup': 5
        })

        self.EVAP.solve()
        # print('ok')
        self.Q_dot_ev = self.EVAP.Q
        self.m_dot_w_ev = ex_ev_w.m_dot
        self.su_ev_w = su_ev_w
        self.ex_ev_w = ex_ev_w
        self.su_cd_w = su_cd_w
        self.ex_cd_w = ex_cd_w

        "---------------------------------------------------------------------------------------------------------------"
        "Water pump consumption"
        "Water pump for the tank"
        # Pump 1: water tank side
        # "Define ports"
        su_w_pp1 = ex_cd_w
        # ex = Mass_connector()

        "Work connector (pas encore coder!!)"
        N_pp = 3000 #50 Hz

        circu = "DN40"

        "Define class"
        PUMP1 = Pump_Wilo()
        PUMP1.inputs(su_w_pp1, N_pp)

        "Model: Wilo Stratos GIGA2.0-I 50/1-20/1,5"
        PUMP1.set_parameters(**{
            'A': 1.30030624e+03,
            'B': 2.55030878e+01,
            'Q_i_star': 3.41914286e-03,
            'k_f': 8.21293476e+02,
            'k_s': 2.65756790e+00,
            'k_l': 1.57368361e-02,
            'k_m': 8.87680941e+00,
            'k_r': 1.56474790e-30,
            'A_EC': 2.64061001e-12,
            'B_EC': 1.46222613e+00,
            'C_EC': 5.61373947e+00,
            'D_EC': 8.26672680e+00,
            'E_EC': 5.67766643e+00,
            'F_EC': 5.51663353e+00,
            'N_sim': 3450
        })

        PUMP1.solve()
        # print(PUMP1.W_dot_motor)
        W_dot_pp1 = PUMP1.W_dot_motor

        "Water pump for the water side (aerotherme)"
        # Pump 1: water tank side
        # "Define ports"
        su_w_pp2 = ex_ev_w
        # ex = Mass_connector()

        "Work connector (pas encore coder!!)"
        N_pp = 3000 #50 Hz

        circu = "DN40"

        "Define class"
        PUMP2 = Pump_Wilo()
        PUMP2.inputs(su_w_pp2, N_pp)

        "Model: Wilo Stratos GIGA2.0-I 50/1-20/1,5"
        PUMP2.set_parameters(**{
            'A': 1.30030624e+03,
            'B': 2.55030878e+01,
            'Q_i_star': 3.41914286e-03,
            'k_f': 8.21293476e+02,
            'k_s': 2.65756790e+00,
            'k_l': 1.57368361e-02,
            'k_m': 8.87680941e+00,
            'k_r': 1.56474790e-30,
            'A_EC': 2.64061001e-12,
            'B_EC': 1.46222613e+00,
            'C_EC': 5.61373947e+00,
            'D_EC': 8.26672680e+00,
            'E_EC': 5.67766643e+00,
            'F_EC': 5.51663353e+00,
            'N_sim': 3450
        })

        PUMP2.solve()
        # print(PUMP2.W_dot_motor)
        W_dot_pp2 = PUMP2.W_dot_motor

        "Performances"
        self.W_dot_tot = self.W_dot_comp+W_dot_pp1+W_dot_pp2
        self.COP_HP = self.Q_dot_cd/self.W_dot_comp
        self.COP_tot = self.Q_dot_cd/self.W_dot_tot
        self.m_dot_r = self.ex_cd.m_dot

        "Residuals"
        self.res_1 = abs(self.P_cd-self.ex_cd.p)
        self.res_2 = abs(self.P_ev-self.ex_ev.p)
        # self.res_1 = (h_init_cd-self.ex_cd.h)
        # self.res_2 = (h_init_ev-self.ex_ev.h)

        "For the graphs"
        if self.EVAP.flag == 1 or self.COND.flag == 1:
            self.res_1 = 10000
            self.res_2 = 10000
        
        
        self.res = [self.res_1, self.res_2]
        print(self.res)
        return self.res_1, self.res_2
    
    def plot_Ts(self):
        "Gaphs"
        #Create array with point of the cycle
        T_cd = PropsSI('T', 'P', self.P_cd, 'Q', 0.5, 'R1233ZDE')
        T_ev = PropsSI('T', 'P', self.P_ev, 'Q', 0.5, 'R1233ZDE')
        s_cd_1 = PropsSI('S', 'P', self.P_cd, 'Q', 1, 'R1233ZDE')
        s_cd_0 = PropsSI('S', 'P', self.P_cd, 'Q', 0, 'R1233ZDE')
        s_ev_1 = PropsSI('S', 'P', self.P_ev, 'Q', 1, 'R1233ZDE')
        s_ev_0 = PropsSI('S', 'P', self.P_ev, 'Q', 0, 'R1233ZDE')
        
        self.T_array = [self.su_comp.T, self.ex_comp.T, self.su_cd.T, T_cd, T_cd, self.ex_cd.T, self.su_ev.T, T_ev, self.su_comp.T]
        self.s_array = [self.su_comp.s, self.ex_comp.s, self.su_cd.s, s_cd_1, s_cd_0, self.ex_cd.s, self.su_ev.s, s_ev_1, self.su_comp.s]

        T_h = np.linspace(self.ex_cd_w.T, self.su_cd_w.T, 100)
        s_h = np.linspace(self.su_cd.s, self.ex_cd.s, 100)
        T_c = np.linspace(self.su_ev_w.T, self.ex_ev_w.T, 100)
        s_c = np.linspace(self.ex_ev.s, self.su_ev.s, 100)
        # print('TH', T_h, s_h)
        plt.plot(s_h, T_h, color='blue', linestyle='-', label='Heat sink')
        plt.plot(s_c, T_c, color='red', linestyle='-', label='Heat source')

        TS_curve = TS_curve_generator('R1233ZDE')
        TS_curve.points(self.s_array, self.T_array)

    def plot_Ph(self):
        self.P_array = [self.su_comp.p, self.ex_comp.p, self.su_cd.p, self.ex_cd.p, self.su_vlv.p, self.ex_vlv.p, self.su_ev.p, self.su_comp.p]
        self.h_array = [self.su_comp.h, self.ex_comp.h, self.su_cd.h, self.ex_cd.h, self.su_vlv.h, self.ex_vlv.h, self.su_ev.h, self.su_comp.h]

        Ph_curve = PH_curve_generator('R1233ZDE')
        Ph_curve.points(self.P_array, self.h_array)

    def solve(self):
            
            # x_pp_guess = [2, 3, 4, 5, 6, 7, 8, 9, 10]
            # stop = 0
            # j=0
            # while not stop and j < len(x_pp_guess):
            #     i=0
            #     while not stop and i < len(x_pp_guess):
            #         P_cd_guess = PropsSI('P', 'T', self.T_w_cd_su+x_pp_guess[j]+5, 'Q', 0, 'R1233ZDE')
            #         P_ev_guess = PropsSI('P', 'T', self.T_w_ev_su-x_pp_guess[i]-5, 'Q', 0, 'R1233ZDE')

            #         try:
            #             root(self.system, [P_cd_guess, P_ev_guess], tol=1e-2)
            #             res_norm = np.linalg.norm(self.res)
            #         except:
            #             res_norm = 1
            #             pass 
            #         if res_norm < 1:
            #             stop = 1
            #         i = i + 1
            #     j = j + 1

        P_cd_guess = PropsSI('P', 'T', self.T_w_cd_su+self.glide_cd+2+5, 'Q', 0, 'R1233ZDE') # 2: pp, 5: sc
        P_ev_guess = PropsSI('P', 'T', self.T_w_ev_su-self.glide_ev-2-5, 'Q', 0, 'R1233ZDE')
        # print(P_cd_guess, P_ev_guess)
        P_cd_res = 10000
        it_cd = 0
        # HP Iterative resolution
        while abs(P_cd_res)>100 and it_cd<20:
            # LP iterative resolution
            P_ev_res = 10000
            it_ev = 0
            while abs(P_ev_res)>100 and it_ev<20:
                # print('P_ev_guess', P_ev_guess, 'P_cd_guess', P_cd_guess)
                res = self.system([P_cd_guess, P_ev_guess])
                # print('res', res)   
                # P_ev_res = self.res_2
                #Problème ici: quand je commence avec ce guess, dans le modèle de l'évaporateur, il commence avec ub=lb et mets donc la les pression égales.
                if self.EVAP.flag == 1:
                    # print('coucou')
                    P_ev_guess = P_ev_guess - 2000
                    P_ev_res = 1000
                else:
                    P_ev_guess = self.ex_ev.p
                    P_ev_res = self.res_2
                it_ev = it_ev+1
                # print('it_ev', it_ev)
            P_cd_res = self.res_1
            # print('P_cd_guess')
            # print(P_cd_guess)
            P_cd_guess = self.ex_cd.p
            # print(P_cd_guess)
            it_cd = it_cd+1
            # print('it_cd', it_cd)
            print('it_cd', it_cd)
            if it_cd == 30:
                P_cd_guess = P_cd_guess + 20
        if it_cd > 100:
            print('No convergence')



if __name__ == '__main__':

    start_time = time.time()
    HP = HP()
    HP.inputs(6000, 293, 'R1233ZDE', 76+273.15, 14, 62+273.15, 14, 2)
    HP.set_parameters(**{
        'T_sh': 5,
        'T_sc' : 5,
    })
    HP.solve()
    print(HP.W_dot_tot)
    print(HP.W_dot_comp)
    print(HP.Q_dot_cd)
    print(HP.Q_dot_ev)
    print(HP.COP_HP)
    print(HP.m_dot_r)
    print(HP.P_cd)
    print(HP.P_ev)
    print(HP.epsilon_is)
    print(HP.m_dot_w_cd)
    print(HP.m_dot_w_ev)
    print(HP.res)

    HP.plot_Ts()




    # TS_curve.PH_curve()
