# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 17:32:08 2023

@author: Elise
"""

import pandas as pd

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from Port.Mass_connector import Mass_connector
from Component.Expander.Expander_SE import Expander_Class
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt

Datas_raw = pd.read_csv('Component\\Expander\\Exp_Data_Camilla.csv')

import matplotlib.pyplot as plt

# Extract the first 20 lines
Data_set = Datas_raw.head(150)

#-----------------------------------------------------------------------
"Inputs"
P_su = Data_set['P_su']
rp = Data_set['Press.Ratio']
P_ex = Data_set['P_ex']
N_exp_RPM = Data_set['RPM']
T_su = Data_set['T_su']

#------------------------------------------------------------------------
"Datas to calibrate"
W_dot_sh_meas = Data_set['W_shaft']
m_dot_meas = Data_set['m_dot_wf']
Torque_meas = Data_set['Torque']

#------------------------------------------------------------------------
"Parameters"

m_dot_calc = []
W_dot_sh_calc = []
N_rot = []

print(P_su)
print(rp)
print(N_exp_RPM)
print(T_su)

for i in range(0, len(P_su)):
    
    "Define ports"
    su = Mass_connector()
    ex = Mass_connector()

    su.set_fluid('R245fa')
    su.set_T(T_su[i])
    su.set_p(P_su[i])
    print(su.T)
    print(su.p)
    ex.set_p(P_ex[i])
    print(ex.p)

    "Work connector (pas encore coder!!)"
    N_exp = N_exp_RPM[i]
    print(N_exp)

    "Heat connector (pas encore coder!!)"
    T_amb = 293

    "Define class"
    EX = Expander_Class()
    EX.inputs(su, ex, N_exp, T_amb)
    
    #Calibrated parameters
    EX.set_parameters(**{
        'AU_amb': 8.33758799e+00,
        'AU_su_n': 6.67152053e-01,
        'AU_ex_n': 3.21181352e+01,
        'd_su1': 6.31789061e-03,
        'm_dot_n': 0.1,
        'A_leak': 1.00000000e-10,
        'W_dot_loss_0': 8.19123951e-01,
        'alpha': 7.79756524e-02,
        'rv_in': 1.7,  #1.7,
        'V_s': 0.0000712,
        'C_loss': 4.68294054e-01
    })
    
    EX.solve()
    print(i)
    print(EX.convergence)
    m_dot_calc.append(EX.m_dot)
    W_dot_sh_calc.append(EX.W_dot)






#----------------------------------------------------------------------------
# Parity plot of the mass flow rate

fig,ax = plt.subplots(figsize=(5.5,4.3),constrained_layout=True) 

plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas), max(m_dot_meas)], 'k', linewidth=1.7, label=r'45° line')
plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas)*1.1, max(m_dot_meas)*1.1], '--', linewidth=1.4, color='gray', label=r'10% error')
plt.plot([min(m_dot_meas), max(m_dot_meas)], [min(m_dot_meas)*0.9, max(m_dot_meas)*0.9], '--', linewidth=1.4, color='gray')

plt.scatter(m_dot_meas, m_dot_calc, marker='*', label='Data points')


plt.xlabel(r'$\dot{m}_{meas} [kg/s]$', fontsize=16)
plt.ylabel(r'$\dot{m}_{calc} [kg/s]$', fontsize=16)

plt.legend(fontsize=12)

plt.grid(True)

plt.savefig('Calibration_m_dot.eps', format='eps', bbox_inches='tight')
plt.savefig('Calibration_m_dot.svg', format='svg', bbox_inches='tight')
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

plt.savefig('Calibration_W_dot.eps', format='eps', bbox_inches='tight')
plt.savefig('Calibration_W_dot.svg', format='svg', bbox_inches='tight')
plt.show()