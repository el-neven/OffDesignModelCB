# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 13:05:24 2023

@author: Elise
"""

import pandas as pd

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from Port.Mass_connector import Mass_connector
from Component.Compressor.Compressor_SE import Compressor_Class
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt

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

#------------------------------------------------------------------------
"Parameters"

m_dot_calc = []
W_dot_sh_calc = []

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
    
    "Found with optimization"
    CP.set_parameters(**{
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
    
    CP.solve()
    

    m_dot_calc.append(CP.m_dot)
    W_dot_sh_calc.append(CP.W_dot)
    




#----------------------------------------------------------------------------
# Parity plot of the mass flow rate

fig,ax = plt.subplots(figsize=(5.5,4.3),constrained_layout=True) 

plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas), max(m_dot_meas)], 'k', linewidth=1.7, label=r'45° line')
plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas*1.1), max(m_dot_meas)*1.1], '--', linewidth=1.4, color='gray', label=r'10% error')
plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas*0.9), max(m_dot_meas)*0.9], '--', linewidth=1.4, color='gray')


plt.scatter(m_dot_meas, m_dot_calc, marker='*', label='Data points')

plt.xlabel(r'$\dot{m}_{meas} [kg/s]$', fontsize=16)
plt.ylabel(r'$\dot{m}_{calc} [kg/s]$', fontsize=16)

plt.legend(fontsize=12)

plt.grid(True)

plt.savefig('Calibration_m_dot_cp.eps', format='eps', bbox_inches='tight')
plt.savefig('Calibration_m_dot_cp.svg', format='svg', bbox_inches='tight')
plt.show()


#----------------------------------------------------------------------------
# Parity plot of the work at the shaft

plt.plot([min(W_dot_sh_meas), max(W_dot_sh_meas)], [min(W_dot_sh_meas), max(W_dot_sh_meas)], 'k', linewidth=1.7, label=r'45° line')
plt.plot([min(W_dot_sh_meas), max(W_dot_sh_meas)], [min(W_dot_sh_meas)*1.1, max(W_dot_sh_meas)*1.1], '--', linewidth=1.4, color='gray', label=r'10% error')
plt.plot([min(W_dot_sh_meas), max(W_dot_sh_meas)], [min(W_dot_sh_meas)*0.9, max(W_dot_sh_meas)*0.9], '--', linewidth=1.4, color='gray')

plt.scatter(W_dot_sh_meas, W_dot_sh_calc, marker='*', label='Data points')

plt.xlabel(r'$\mathrm{\dot{W}_{sh, meas}}  [W]$', fontsize=16)
plt.ylabel(r'$\mathrm{\dot{W}_{sh, calc}} [W]$', fontsize=16)

plt.legend(fontsize=12)

plt.grid(True)

plt.savefig('Calibration_W_dot_cp.eps', format='eps', bbox_inches='tight')
plt.savefig('Calibration_W_dot_cp.svg', format='svg', bbox_inches='tight')
plt.show()



