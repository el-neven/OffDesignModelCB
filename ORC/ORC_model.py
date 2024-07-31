from CoolProp.CoolProp import PropsSI
from scipy.optimize import least_squares
from scipy.optimize import root
import numpy as np

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')
# print(sys.path)
from Component.Expander.Expander_SE import Expander_Class
from Component.PHE.Moving_Boundary.HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger
from Component.Pump_ref.Model import PumpModel
from Component.PHE.Simple_MB.HX_MovingBoundary import HeatExchanger
from Saturation_curves import TS_curve_generator
from Saturation_curves import PH_curve_generator

from Port.Mass_connector import Mass_connector
import time
import matplotlib.pyplot as plt

class ORC():
    def __init__(self):
        self._model = None

        "Design parameters"
        self.T_sh = None # global heat transfer coefficient for supply heat transfer   [W/K]
        self.T_sc = None

    def inputs(self, N_exp, T_amb, wf, T_w_cd_su, glide_cd, T_w_ev_su, glide_ev, Nb_exp):
        self.N_exp = N_exp
        self.T_amb = T_amb
        self.wf = wf
        self.T_w_cd_su = T_w_cd_su
        self.glide_cd = glide_cd
        self.T_w_ev_su = T_w_ev_su
        self.glide_ev = glide_ev
        self.Nb_exp = Nb_exp

    def set_parameters(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: Parameter '{key}' not found in the parameters.")
    
    def system(self, x):
        
        "GUESS"
        self.P_ev = x[1]
        self.P_cd = x[0]
        # self.P_ev = x[0]
        # self.rp = x[1]
        # self.P_cd = self.P_ev/self.rp
        print('P_ev', self.P_ev, 'P_cd', self.P_cd)

        "---------------------------------------1. Expander---------------------------------------"
        "Define ports"
        self.su_exp = Mass_connector()
        self.ex_exp = Mass_connector()

        # Supply state (su_exp)
        T_ev = PropsSI('T', 'P', self.P_ev, 'Q', 1, 'R1233ZDE')
        self.su_exp.set_fluid('R1233ZDE')
        self.su_exp.set_T(T_ev+self.T_sh)
        self.su_exp.set_p(self.P_ev)

        # Exhaust state (ex_exp)
        self.ex_exp.set_p(self.P_cd)

        "Define class"
        EXP = Expander_Class()
        EXP.inputs(self.su_exp, self.ex_exp, self.N_exp, self.T_amb)

        EXP.set_parameters(**{
            'AU_amb': 8.33758799e+00,
            'AU_su_n': 6.67152053e-01,
            'AU_ex_n': 3.21181352e+01,
            'd_su1': 6.31789061e-03,
            'm_dot_n': 0.1,
            'A_leak': 1.00000000e-10,
            'W_dot_loss_0': 8.19123951e-01,
            'alpha': 7.79756524e-02,
            'rv_in': 1.7,
            'V_s': 0.0000712,
            'C_loss': 4.68294054e-01
        })

        EXP.solve()
        print(EXP.convergence)
        self.W_dot_exp = self.Nb_exp*EXP.W_dot

        "---------------------------------------2. Condenser---------------------------------------"

        "Define ports"
        self.su_cd = Mass_connector()
        self.ex_cd = Mass_connector()

        self.su_cd_w = Mass_connector()
        self.ex_cd_w = Mass_connector()

        # Supply state R1233ZDE (su_cd)
        self.su_cd = self.ex_exp
        self.su_cd.set_m_dot(self.ex_exp.m_dot*self.Nb_exp)

        # Supply state (water)
        self.su_cd_w.set_fluid('Water')
        self.su_cd_w.set_T(self.T_w_cd_su)
        self.su_cd_w.set_p(2e5)

        #Exhaust state (water)
        self.ex_cd_w.set_fluid('Water')
        self.ex_cd_w.set_T(self.T_w_cd_su+self.glide_cd)
        self.ex_cd_w.set_p(2e5)

        "Define class"
        COND = HeatExchanger(self.su_cd, self.su_cd_w, self.ex_cd, self.ex_cd_w)

        COND.set_parameters(**{
            'HX_type': 'condenser',
            'HX_D': 0.06, #Diamètre de port d'entré
            'HX_A': 17.8, #Surface d'échange
            'min_pinch': 2,
            'dT_sub_or_sup': 5
        })

        COND.solve()
        self.Q_dot_cd = COND.Q
        self.m_dot_w_cd = self.ex_cd_w.m_dot

        "---------------------------------------3. Pump---------------------------------------"
        "Define ports"
        self.su_pp = Mass_connector()
        self.ex_pp = Mass_connector()

        # Supply state
        self.su_pp = self.ex_cd
        # Exhaust state
        self.ex_pp.set_fluid('R1233ZDE')
        self.ex_pp.set_p(self.P_ev)
        self.ex_pp.set_m_dot(self.ex_cd.m_dot)

        Pump = PumpModel()
        Pump.set_inputs(self.su_pp, self.ex_pp)
        Pump.solve()
        if Pump.convergence == False:
            print('Error')
            self.ex_pp.set_T(self.su_pp.T+1)

        "---------------------------------------4. Evaporator---------------------------------------"

        "Define ports"
        self.su_ev = Mass_connector()
        self.ex_ev = Mass_connector()

        self.su_ev_w = Mass_connector()
        self.ex_ev_w = Mass_connector()

        # Supply state (R1233ZDE)
        self.su_ev = self.ex_pp

        # Supply state (water)
        self.su_ev_w.set_fluid('Water')
        self.su_ev_w.set_T(self.T_w_ev_su)
        self.su_ev_w.set_p(2e5)

        #Exhaust state (water)
        self.ex_ev_w.set_fluid('Water')
        self.ex_ev_w.set_T(self.T_w_ev_su-self.glide_ev)
        self.ex_ev_w.set_p(2e5)

        "Define class"
        EVAP = HeatExchanger(self.su_ev, self.su_ev_w, self.ex_ev, self.ex_ev_w)

        EVAP.set_parameters(**{
            'HX_type': 'evaporator',
            'HX_D': 0.06, #Diamètre de port d'entré
            'HX_A': 17.8, #Surface d'échange
            'min_pinch': 2,
            'dT_sub_or_sup': 5
        })

        EVAP.solve()
        self.Q_dot_ev = EVAP.Q
        self.m_dot_w_ev = self.ex_ev_w.m_dot

        "Performances"
        self.eta_ORC = self.W_dot_exp/EVAP.Q
        self.m_dot_r = self.ex_cd.m_dot

        "Residuals"
        self.res_1 = self.P_cd-self.ex_cd.p
        self.res_2 = self.P_ev-self.ex_ev.p
        print(self.res_1, self.res_2)

        return [self.res_1, self.res_2]
    
    def plot_Ts(self):
        "Gaphs"
        #Create array with point of the cycle
        T_cd = PropsSI('T', 'P', self.P_cd, 'Q', 0.5, 'R1233ZDE')
        T_ev = PropsSI('T', 'P', self.P_ev, 'Q', 0.5, 'R1233ZDE')
        s_cd_1 = PropsSI('S', 'P', self.P_cd, 'Q', 1, 'R1233ZDE')
        s_cd_0 = PropsSI('S', 'P', self.P_cd, 'Q', 0, 'R1233ZDE')
        s_ev_1 = PropsSI('S', 'P', self.P_ev, 'Q', 1, 'R1233ZDE')
        s_ev_0 = PropsSI('S', 'P', self.P_ev, 'Q', 0, 'R1233ZDE')
        
        self.T_array = [self.su_exp.T, self.su_cd.T, T_cd, T_cd, self.ex_cd.T, self.su_ev.T, T_ev, T_ev, self.ex_ev.T, self.su_exp.T]
        self.s_array = [self.su_exp.s, self.su_cd.s, s_cd_1, s_cd_0, self.ex_cd.s, self.su_ev.s, s_ev_0, s_ev_1,self.ex_ev.s, self.su_exp.s]
        print(self.T_array, self.s_array)
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
        self.P_array = [self.su_exp.p, self.ex_exp.p, self.su_cd.p, self.ex_cd.p, self.su_pp.p, self.ex_pp.p, self.su_ev.p, self.su_exp.p]
        self.h_array = [self.su_exp.h, self.ex_exp.h, self.su_cd.h, self.ex_cd.h, self.su_pp.h, self.ex_pp.h, self.su_ev.h, self.su_exp.h]
        print(self.P_array, self.h_array)
        Ph_curve = PH_curve_generator('R1233ZDE')
        Ph_curve.points(self.P_array, self.h_array)

    def solve(self):
        # # Automatic Boundary Conditions
        # P_cd_lb = max(PropsSI('P', 'Q', 0, 'T', self.T_w_cd_su-20, self.wf), PropsSI('P_min', 'Q', 0, 'T', 273.15, self.wf))
        # P_ev_ub = PropsSI('P', 'Q', 0, 'T', min(PropsSI('Tcrit', 'Q', 0, 'T',273, self.wf)-2, self.T_w_ev_su-1), self.wf)

        # rp_max = P_ev_ub/P_cd_lb
        # rp_min = min(1.01, rp_max)
        # P_ev_lb = rp_min*P_cd_lb
        # P_cd_ub = P_ev_ub/rp_min

        # bounds = [(P_ev_lb, rp_min), (P_ev_ub, rp_max)]
        # # print(bounds)

        # P_ev_guess = P_ev_lb+0.5*(P_ev_ub-P_ev_lb)

        # rp_guess = rp_min+0.5*(rp_max-rp_min)

        # least_squares(self.system, [P_ev_guess, rp_guess], xtol=1e-2, bounds=bounds)
        
        P_cd_guess = PropsSI('P', 'T', self.T_w_cd_su+self.glide_cd+2+5, 'Q', 0, 'R1233ZDE') # 2: pp, 5: sc
        P_ev_guess = PropsSI('P', 'T', self.T_w_ev_su-self.glide_ev-2-5, 'Q', 0, 'R1233ZDE')
        # print(P_cd_guess, P_ev_guess)
        P_cd_res = 10000
        it_cd = 0
        # HP Iterative resolution
        try: 
            while abs(P_cd_res)>100 and it_cd<20:
                # LP iterative resolution
                P_ev_res = 10000
                it_ev = 0
                while abs(P_ev_res)>100 and it_ev<20:
                    # print('P_ev_guess', P_ev_guess, 'P_cd_guess', P_cd_guess)
                    try:
                        res = self.system([P_cd_guess, P_ev_guess])
                        P_ev_guess = self.ex_ev.p
                        P_ev_res = self.res_2
                    except:
                        P_ev_guess = P_ev_guess + 2000
                        P_ev_res = 1000

                    # print('res', res)   
                    # P_ev_res = self.res_2
                    #Problème ici: quand je commence avec ce guess, dans le modèle de l'évaporateur, il commence avec ub=lb et mets donc la les pression égales.
                    # if self.EVAP.flag == 1:
                    #     # print('coucou')
                    #     P_ev_guess = P_ev_guess - 2000
                    #     P_ev_res = 1000
                    # else:
                    it_ev = it_ev+1
                    print('it_ev', it_ev)
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
        except:
            print('No convergence')
            self.eta_ORC = 0
            self.m_dot_r = 0
            self.W_dot_exp = 0
            self.P_ev = 0
            self.P_cd = 0
            self.Q_dot_cd = 0
            self.Q_dot_ev = 0
            self.m_dot_w_cd = 0
            self.m_dot_w_ev = 0
            






if __name__ == '__main__':

    start_time = time.time()

    "Model used"
    N_exp = 6000
    T_amb = 293
    T_tap_in = 30+273.15
    glide_cd = 10
    T_sto_HT = 60+273.15
    glide_ev = 10
    Nb_exp = 3

    ORC_instance = ORC()
    ORC_instance.inputs(N_exp, T_amb, 'R1233ZDE', T_tap_in, glide_cd, T_sto_HT, glide_ev, Nb_exp)
    ORC_instance.set_parameters(**{
        'T_sh': 5,
        'T_sc' : 5,
    })
    ORC_instance.solve()

    print("--- %s seconds ---" % (time.time() - start_time))
    print(ORC_instance.eta_ORC)
    print(ORC_instance.m_dot_r)
    print(ORC_instance.W_dot_exp)
    print(ORC_instance.P_ev)
    print(ORC_instance.P_cd)
    print(ORC_instance.Q_dot_cd)
    print(ORC_instance.Q_dot_ev)
    print(ORC_instance.m_dot_w_cd)
    print(ORC_instance.m_dot_w_ev)

    ORC_instance.plot_Ts()
    ORC_instance.plot_Ph()


    
