import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from ORC.ORC_model import ORC

import numpy as np

import time

f = open("ORC\\Results\\tab_ORC_T_tap_30_glide_cd_10.csv", "w")
f.write("index, T_sto_HT [K], T_sto_LT [K], glide_ev [K], T_tap_in [K], T_tap_out [K], glide_cd [K], W_dot_exp [W], p_ev [Pa], p_cd [Pa], rp [-], Q_dot_cd [W], Q_dot_ev [W], eta_ORC [-], m_dot_w_cd [kg/s], m_dot_w_ev [kg/s], Nb_exp [-], m_dot_r [kg/s]")
to_csv = np.zeros(20, dtype=float)


start_time = time.time()
convergence = np.zeros(100, dtype=bool)
T_sto_step = 5 # [K]
T_sto_HT_min = 50 # Actually useless before T=45°C based on the simulations
T_sto_HT_max = 100
T_sto_HT_array = np.arange(T_sto_HT_min, T_sto_HT_max+1, T_sto_step, dtype=float)

T_tap_step = 5 # [K]
T_tap_su_min = 10
T_tap_su_max = 50
T_tap_su_array = np.arange(T_tap_su_min, T_tap_su_max+1, T_tap_step, dtype=float)

glide_ev_step = 2 # [K]
glide_ev_min = 5
glide_ev_max = 20
glide_ev_array = np.arange(glide_ev_min, glide_ev_max+1, glide_ev_step, dtype=float)

scenario = 0

scenario = 0

"Constantes"
N_exp = 6000 #RPM
T_amb = 293
Nb_exp = 2
# glide_ev = 15
glide_cd = 10
T_tap_su = 30

for i in range (0,len(T_sto_HT_array)):
    T_sto_HT = T_sto_HT_array[i] #T_ex_cd
    for j in range (0,len(glide_ev_array)):
        # T_tap_su = T_tap_su_array[j] #T_su_ev
        T_su_w_ev = T_sto_HT+273.15
        T_su_w_cd = T_tap_su+273.15
        glide_ev = glide_ev_array[j]

        ORC_instance = ORC()
        ORC_instance.inputs(N_exp, T_amb, 'R1233ZDE', T_su_w_cd, glide_cd, T_su_w_ev, glide_ev, Nb_exp)
        ORC_instance.set_parameters(**{
            'T_sh': 5,
            'T_sc' : 5,
        })
        ORC_instance.solve()

        to_csv[0] = T_sto_HT
        to_csv[1] = T_sto_HT-glide_cd
        to_csv[2] = glide_ev
        to_csv[3] = T_tap_su
        to_csv[4] = T_tap_su+glide_ev
        to_csv[5] = glide_cd
        to_csv[6] = ORC_instance.W_dot_exp
        to_csv[7] = ORC_instance.P_ev
        to_csv[8] = ORC_instance.P_cd
        to_csv[9] = 0 #Need to run the tests with rp
        to_csv[10] = ORC_instance.Q_dot_cd
        to_csv[11] = ORC_instance.Q_dot_ev
        to_csv[12] = ORC_instance.eta_ORC
        to_csv[13] = ORC_instance.m_dot_w_cd
        to_csv[14] = ORC_instance.m_dot_w_ev
        to_csv[15] = Nb_exp
        to_csv[16] = ORC_instance.m_dot_r

        
        scenario += 1
        f.write("\n"+str(scenario)+","+str(to_csv[0])+","+str(to_csv[1])+","+str(to_csv[2])+","+str(to_csv[3])+","+str(to_csv[4])+","+str(to_csv[5])+","+str(to_csv[6])+","+str(to_csv[7])+","+str(to_csv[8])+","+str(to_csv[9])+","+str(to_csv[10])+","+str(to_csv[11])+","+str(to_csv[12])+","+str(to_csv[15])+","+str(to_csv[16]))
        
f.close()

print("--- %s seconds ---" % (time.time() - start_time))
print(scenario, "scenarios simulated")
# print("Convergence [%]:", sum(convergence)/len(convergence)*100)
