# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:39:37 2023

@author: Elise
"""
import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

# import pandas as pd

from Port.Mass_connector import Mass_connector
from Expander_SE import Expander_Class

import numpy as np
import matplotlib.pyplot as plt

"Define ports"
su = Mass_connector()
ex = Mass_connector()

su.set_fluid('R245fa')
su.set_T(374.1)
su.set_p(955214.9)

ex.set_p(293940.1)

"Work connector (pas encore coder!!)"
N_exp = 1500

"Heat connector (pas encore coder!!)"
T_amb = 293

"Define class"
EX = Expander_Class()
EX.inputs(su, ex, N_exp, T_amb)

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
print(EX.convergence)
print(EX.W_dot)

# #---------------------------------------------------------------------------------
# # Graph varying the pressure ratio

# # Create lists to store results
# rp_values = np.linspace(1.3, 6, 25)
# epsilon_is = []

# # Loop through different values of CP.set_rp
# for rp in rp_values:
    
#     EX.set_rp(rp)
    
#     EX.solve()
#     epsilon_is.append(EX.epsilon_is)


# plt.plot(rp_values, epsilon_is, linewidth=1.9)

# plt.xlabel(r'$\mathrm{r_{p}}$ [-]', fontsize=16)
# plt.ylabel(r'$\mathrm{\epsilon_{is}}$ [-]', fontsize=16)

# plt.ylim([0, 1.05])
# plt.xlim([1.4, 6])

# plt.grid(True)
# plt.show()




