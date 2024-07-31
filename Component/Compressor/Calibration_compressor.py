# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 09:10:59 2023

@author: Elise
"""

import pandas as pd

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from Port.Mass_connector import Mass_connector
from Component.Compressor.Compressor_SE import Compressor_Class

from CoolProp.CoolProp import PropsSI

from scipy.optimize import least_squares

Datas_raw = pd.read_csv('Component\\Compressor\\Exp_data_cp.csv')

# Extract the first 55 lines
Data_set = Datas_raw.head(55)

#-----------------------------------------------------------------------
"Inputs"
P_su = Data_set['P_su']
P_ex = Data_set['P_ex']
rp = P_ex/P_su
N_cp_RPM = Data_set['RPM']
SH = Data_set['SH']

#------------------------------------------------------------------------
"Datas to calibrate"
W_dot_sh_meas = Data_set['W_dot_shaft']
m_dot_meas = Data_set['m_dot_wf']

#---------------------------------------------------------------------------
# Least square method to fit the parameters to the experimental data

def objective(params_opt):
    AU_amb, AU_su_n, AU_ex_n, d_ex, A_leak, W_dot_loss_0, alpha, C_loss = params_opt
    #print(AU_amb, AU_su_n, AU_ex_n, d_ex, A_leak, W_dot_loss_0, alpha, C_loss, V_s)
    
    total_error = 0
    m_dot_calc = []
    W_dot_sh_calc = []
    diff = []
    
    for i in range(0, len(P_su)):
        
        "Define port"
        su = Mass_connector()
        ex = Mass_connector()

        su.set_fluid('R134a')
        # su.set_T(338.40)
        su.set_T(PropsSI('T', 'P', P_su[i], 'Q', 0.5, 'R134a') + SH[i])
        su.set_p(P_su[i])
        # C'est comme si le problème venait de CoolProp, si on met une température d'entrée plus haute, PropsSI ne se trompe pas.$
        # Hors avec la P et la T qui aumente, S devrait augmenter non?

        ex.set_p(rp[i]*P_su[i])

        "Work connector (pas encore coder!!)"
        N_exp = 6000

        "Heat connector (pas encore coder!!)"
        T_amb = 293

        CP = Compressor_Class()
        CP.inputs(su, ex, N_cp_RPM[i], T_amb)

        CP.set_parameters(**{
            'AU_amb': AU_amb,
            'AU_su_n': AU_su_n,
            'AU_ex_n': AU_ex_n,
            'd_ex': d_ex,
            'm_dot_n': 0.1,
            'A_leak': A_leak,
            'W_dot_loss_0': W_dot_loss_0,
            'alpha': alpha,
            'C_loss': C_loss,
            'rv_in': 1.7,
            'V_s': 1.21e-4,
        })

        CP.solve()

    
        # CP.set_V_s(V_s)

        m_dot_calc.append(CP.m_dot)
        W_dot_sh_calc.append(CP.W_dot)
        
        diff.append(((CP.W_dot - W_dot_sh_meas[i])/W_dot_sh_meas[i])*((CP.m_dot - m_dot_meas[i])/m_dot_meas[i]))
    
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
C_loss_guess = 0.1
V_s_guess = 0.00012

initial_params = [AU_amb_guess, AU_su_n_guess, AU_ex_n_guess, d_su_guess, A_leak_guess, W_dot_loss_guess, alpha_guess, C_loss_guess]

# Boundaries on the parameters
bounds = ([0.1, 0.1, 0.1, 0.1e-3, 1e-10, 0.001, 0.001, 0.000001], [100, 100, 100, 20e-3, 1e-5, 1000, 0.999, 50])
                                                                          
# Minimization of the error by adjusting the parameters
result = least_squares(objective, initial_params, bounds=bounds)

# Optimal parameters
optimal_params = result.x

print("Paramètres optimaux :", optimal_params)


# Paramètres optimaux : [1.00110546e+01 4.99878877e+00 5.00028640e+00 1.00954403e-02 5.60763366e-10 1.28169361e+00 4.88233002e-02 1.02860767e-01]



