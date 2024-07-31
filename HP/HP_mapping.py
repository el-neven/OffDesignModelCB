import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from HP.HP_model import HP

import numpy as np

import time

f = open("HP\\Results\\tab_HP_TDH_80_glide_ev_10.csv", "w")
f.write("index, T_sto_HT [K], T_sto_LT [K], glide_cd [K], T_DH+ [K], T_DH- [K], glide_ev [K], W_dot_cp [W], p_ev [Pa], p_cd [Pa], rp_cp [-], Q_dot_cd [W], Q_dot_ev [W], COP [-], COP_tot [-], m_dot_w_cd [kg/s], m_dot_w_ev [kg/s], Nb_comp [-], m_dot_r [kg/s]")
to_csv = np.zeros(20, dtype=float)

start_time = time.time()
convergence = np.zeros(100, dtype=bool)
T_sto_HT_step = 5 # [K]
T_sto_HT_min = 60 # Actually useless before T=45°C based on the simulations
T_sto_HT_max = 100
T_sto_HT_array = np.arange(T_sto_HT_min, T_sto_HT_max+1, T_sto_HT_step, dtype=float)

T_sto_LT_step = 5 # [K]
T_sto_LT_min = 10 # Actually useless before T=45°C based on the simulations
T_sto_LT_max = 90
T_sto_LT_array = np.arange(T_sto_LT_min, T_sto_LT_max+1, T_sto_LT_step, dtype=float)

T_DH_step = 5 # [K]
T_DH_su_min = 40
T_DH_su_max = 80
T_DH_su_array = np.arange(T_DH_su_min, T_DH_su_max+1, T_DH_step, dtype=float)

glide_cd_step = 1 # [K]
glide_cd_min = 5
glide_cd_max = 20
glide_cd_array = np.arange(glide_cd_min, glide_cd_max+1, glide_cd_step, dtype=float)

scenario = 0

"Constantes"
N_cp = 6000 #RPM
T_amb = 293
Nb_cp = 2
glide_ev = 10
# glide_cd = 15
T_DH_su = 80

for i in range (0,len(T_sto_LT_array)):
    T_sto_LT = T_sto_LT_array[i] #T_ex_cd
    for j in range (0,len(glide_cd_array)):
        # T_DH_su = T_DH_su_array[j] #T_su_ev
        glide_cd = glide_cd_array[j]
        T_su_w_ev = T_DH_su+273.15
        T_su_w_cd = T_sto_LT+273.15

        HP_instance = HP()
        HP_instance.inputs(N_cp, T_amb, 'R1233ZDE', T_su_w_cd, glide_cd, T_su_w_ev, glide_ev, Nb_cp)
        HP_instance.set_parameters(**{
            'T_sh': 5,
            'T_sc' : 5,
        })
        HP_instance.solve()

        to_csv[0] = T_sto_LT+glide_cd
        to_csv[1] = T_sto_LT
        to_csv[2] = glide_cd
        to_csv[3] = T_DH_su
        to_csv[4] = T_DH_su-glide_ev
        to_csv[5] = glide_ev
        to_csv[6] = HP_instance.W_dot_comp
        to_csv[7] = HP_instance.P_ev
        to_csv[8] = HP_instance.P_cd
        to_csv[9] = HP_instance.rp
        to_csv[10] = HP_instance.Q_dot_cd
        to_csv[11] = HP_instance.Q_dot_ev
        to_csv[12] = HP_instance.COP_HP
        to_csv[13] = HP_instance.COP_tot
        to_csv[14] = HP_instance.m_dot_w_cd
        to_csv[15] = HP_instance.m_dot_w_ev
        to_csv[16] = Nb_cp
        to_csv[17] = HP_instance.m_dot_r

        
        scenario += 1
        f.write("\n"+str(scenario)+","+str(to_csv[0])+","+str(to_csv[1])+","+str(to_csv[2])+","+str(to_csv[3])+","+str(to_csv[4])+","+str(to_csv[5])+","+str(to_csv[6])+","+str(to_csv[7])+","+str(to_csv[8])+","+str(to_csv[9])+","+str(to_csv[10])+","+str(to_csv[11])+","+str(to_csv[12])+","+str(to_csv[15])+","+str(to_csv[16])+","+str(to_csv[17]))
        
f.close()

print("--- %s seconds ---" % (time.time() - start_time))
print(scenario, "scenarios simulated")
# print("Convergence [%]:", sum(convergence)/len(convergence)*100)
