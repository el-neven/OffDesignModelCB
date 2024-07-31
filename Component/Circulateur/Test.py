import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

# import pandas as pd

from Port.Mass_connector import Mass_connector
from Component.Circulateur.Simulation_model import Pump_Wilo




import numpy as np
import matplotlib.pyplot as plt

"Define ports"
su = Mass_connector()
ex = Mass_connector()

su.set_fluid('Water')
su.set_D(1000)
su.set_m_dot(2)

"Work connector (pas encore coder!!)"
N_pp = 3000 #50 Hz

circu = "DN40"

"Define class"
PUMP = Pump_Wilo()
PUMP.inputs(su, N_pp)

"Model: Wilo Stratos GIGA2.0-I 50/1-20/1,5"

if circu == "DN40":
    PUMP.set_parameters(**{
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

else:
    PUMP.set_parameters(**{
        'A': 0,
        'B': 0,
        'Q_i_star': 0,
        'k_f': 0,
        'k_s': 0,
        'k_l': 0,
        'k_m': 0,
        'k_r': 0,
        'A_EC': 0,
        'B_EC': 0,
        'C_EC': 0,
        'D_EC': 0,
        'E_EC': 0,
        'F_EC': 0,
        'N_sim': 0
    })

PUMP.solve()
print(PUMP.W_dot_motor)

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




