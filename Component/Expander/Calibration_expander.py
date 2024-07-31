# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 14:51:44 2023

@author: Elise
"""

import pandas as pd

from point import Point_on_cycle
from Expander_SE import Expander_Class

from scipy.optimize import least_squares


Datas_raw = pd.read_csv('Données_Exp_Camilla.csv')

# Extract the first 20 lines
Data_set = Datas_raw.head(150)


#-----------------------------------------------------------------------
"Inputs"
P_su = Data_set['P_su']
rp = Data_set['Press.Ratio']
N_exp_RPM = Data_set['RPM']
T_su = Data_set['T_su']

#------------------------------------------------------------------------
"Datas to calibrate"
W_dot_sh_meas = Data_set['W_shaft']
m_dot_meas = Data_set['m_dot_wf']
Torque_meas = Data_set['Torque']


# -------------------------------------------------------------

# Least square method to fit the parameters to the experimental data

def objective(params_opt):
    # Paramètres du modèle semi-empirique
    AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss = params_opt
    # print(AU_amb, AU_su_n, AU_ex_n, d_su, A_leak, W_dot_loss_0, alpha, C_loss)
    
    m_dot_calc = []
    W_dot_sh_calc = []
    diff = []
    
    for i in range(0, len(P_su)):
        
        EX = Expander_Class()
        
        # Fitting parameters
        EX.set_AU_amb(AU_amb)
        EX.set_AU_su_n(AU_su_n)
        EX.set_AU_ex_n(AU_ex_n)
        EX.set_d_su1(d_su)
        EX.set_A_leak(A_leak)
        EX.set_W_dot_loss_0(W_dot_loss_0)
        EX.set_alpha(alpha)
        EX.set_C_loss(C_loss)
        
        # Known parameters -> OK
        EX.set_rv_in(1.7)
        EX.set_V_s(0.0000712)
        EX.set_T_amb(293)
        EX.set_m_dot_n(0.1)
        
        point_su = Point_on_cycle()
        point_ex = Point_on_cycle()

        point_su.set_fluid('R245fa')
        
        # "Inputs"
        point_su.set_T(T_su[i])
        point_su.set_p(P_su[i])

        EX.set_N_rot(N_exp_RPM[i])
        EX.set_rp(rp[i])

        EX.set_su(point_su)
        EX.set_ex(point_ex)

        EX.solve()
        #print(EX.m_dot)
        m_dot_calc.append(EX.m_dot)
        W_dot_sh_calc.append(EX.W_dot)
        
        diff.append(((EX.W_dot - W_dot_sh_meas[i])/W_dot_sh_meas[i])*((EX.m_dot - m_dot_meas[i])/m_dot_meas[i]))
    
    return diff

# Initial guess on the parameters of the semi-empirical model
AU_amb_guess = 10
AU_su_n_guess = 5
AU_ex_n_guess = 5
d_su_guess = 6.5e-3
m_dot_n_guess = 0.1
A_leak_guess = 2e-7
W_dot_loss_guess = 1
alpha_guess = 0.05
C_loss_guess = 0.5

initial_params = [AU_amb_guess, AU_su_n_guess, AU_ex_n_guess, d_su_guess, A_leak_guess, W_dot_loss_guess, alpha_guess, C_loss_guess]

# Boundaries on the parameters
bounds = ([0.1, 0.1, 0.1, 0.1e-3, 1e-10, 0.001, 0.001, 0.001], [100, 100, 100, 20e-3, 1e-5, 1000, 0.999, 50])

# Minimization of the error by adjusting the parameters                                                                        
result = least_squares(objective, initial_params, bounds=bounds)

# Optimal parameters
optimal_params = result.x

print("Paramètres optimaux :", optimal_params)







