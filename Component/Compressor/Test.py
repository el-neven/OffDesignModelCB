#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 12:09:29 2023

@author: olivierthome
"""

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

# import pandas as pd

from Port.Mass_connector import Mass_connector
from Component.Compressor.Compressor_SE import Compressor_Class
import time

from CoolProp.CoolProp import PropsSI
# import matplotlib.pyplot as plt

"Define port"
su = Mass_connector()
ex = Mass_connector()

su.set_fluid('R1233ZDE')
su.set_T(331.033964665788)
su.set_p(319296.5575177148)


ex.set_p(606240.1433176235)

"Work connector (pas encore coder!!)"
N_exp = 6000

"Heat connector (pas encore coder!!)"
T_amb = 293

start_time = time.time()
CP = Compressor_Class()
CP.inputs(su, ex, N_exp, T_amb)

# CP.set_parameters(**{
#     'AU_amb': 1.16167919e+01,
#     'AU_su_n': 5.97522947e+00,
#     'AU_ex_n': 4.78129121e+00,
#     'd_ex': 1.85102372e-02,
#     'm_dot_n': 0.1,
#     'A_leak': 3.69895065e-07,
#     'W_dot_loss_0': 1.00951405e+00,
#     'alpha': 1.08884500e-04,
#     'C_loss': 2.27e-02,
#     'rv_in': 1.7,
#     'V_s': 1.17966761e-04,
# })

CP.set_parameters(**{
    'AU_amb': 9.96513290e+00,
    'AU_su_n': 1.02359773e+01,
    'AU_ex_n': 2.24133147e+00 ,
    'd_ex': 1.82304791e-02,
    'm_dot_n': 0.1,
    'A_leak': 3.66336680e-07,
    'W_dot_loss_0': 9.05482168e-01,
    'alpha': 3.22395090e-03,
    'C_loss': 1.11169710e-06,
    'rv_in': 1.7,
    'V_s': 1.17889079e-04  ,
})
# Paramètres optimaux : [1.00110546e+01 4.99878877e+00 5.00028640e+00 1.00954403e-02 5.60763366e-10 1.28169361e+00 4.88233002e-02 1.02860767e-01]





CP.solve()
print(CP.convergence)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.2f} [s]")
print(CP.W_dot)
print(CP.T_w-273.15)
print(su.s, ex.s)
print(CP.m_dot)
print(ex.T-273.15, su.T-273.15)
print(PropsSI('S', 'H', ex.h, 'P', ex.p, 'R1233ZDE'), PropsSI('S', 'T', su.T, 'P', su.p, 'R1233ZDE'))
print(PropsSI('Q', 'H', ex.h, 'P', ex.p, 'R1233ZDE'), PropsSI('Q', 'T', su.T, 'P', su.p, 'R1233ZDE'))


