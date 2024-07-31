import numpy as np
from scipy.optimize import least_squares

from Simulation_model import Pump_Wilo

circu = "DN40"

# Valeurs max du circulateur
if circu == "DN40":
    n_max = 3450 #frequency
    Q_max = 27*(1/3600)
    H_max = 22
    P_max = 1400
    
elif circu == "DN50":
    n_max = 4020
    Q_max = 41*(1/3600)
    H_max = 19
    P_max = 1400
    
else:
    n_max = 0
    Q_max = 0
    H_max = 0
    P_max = 0

# Mesures depuis datasheet circulateur
if circu == "DN40":
    n = np.array([1000, 1500, 2000, 2000, 2668, 2668, 2668, 2668, 3000, 3000, 3000, 3450, 3450, 3450]) #[rpm]
    Q = np.array([4.5, 5.1, 5.9, 13.8, 2.1, 9.6, 11.5, 16.6, 4.7, 14.3, 20.8, 8.5, 20.2, 26.5])*(1/3600) #[m3/s]
    H = np.array([1.7, 4.1, 7.3, 5, 13.7, 12.5, 12, 9.8, 17, 14.5, 10.7, 22, 17, 12]) #[m_h2o]
    P = np.array([52, 113, 213, 307, 296, 513, 562, 661, 487, 815, 943, 883, 1346, 1467]) #[W]
    
elif circu == "DN50":
    n = np.array([0,0]) #[rpm]
    Q = np.array([0,0])*(1/3600) #[m3/h]
    H = np.array([0,0]) #[m_h2o]
    P = np.array([0,0]) #[W]
    
else:
    n = np.array([0,0]) #[rpm]
    Q = np.array([0,0]) #[m3/h]
    H = np.array([0,0]) #[m_h2o]
    P = np.array([0,0]) #[W]

n_sim = n_max*1
# Similarité (fait correspondre à chaque point de mesure son équivalent par similitude à n=n_sim)
# tout les Q_sim et H_sim sont déterminé par similitude à n=n_sim
print ("n de similarité: %.0f rpm\n" %(n_sim))
Q_sim = np.zeros(len(n), dtype=float)
H_sim = np.zeros(len(n), dtype=float)

for i in range (0,len(n)): #On reconstruit la courbe Q vs H pour n=nsim en utilisant la similitude. ca nous donne une courbe de "référence" oùsur la quelle on peut faire la similitude pour n'importe que N
    Q_sim[i] = Q[i]*(n_sim/n[i])
    H_sim[i] = H[i]*(n_sim/n[i])**2

Q_max_sim = Q_max*(n_sim/n_max)
H_max_sim = H_max*(n_sim/n_max)**2


import matplotlib.pyplot as plt

# plt.scatter(Q_sim, H_sim)
# plt.xlabel('Q_sim')
# plt.ylabel('H_sim')
# plt.title('Q_sim vs H_sim')
# plt.grid(True)
# plt.show()

"On va calibrer les paramètres basé sur la similarité"


#Calibration faite sur les points à N=N_sim
j = np.arange(len(n))
def calibration(parameters, Q_sim, j):
    circu = "DN40"
    #Paramètres à calibrer (pour N=N_sim)
    A = parameters[0]
    B = parameters[1]
    Q_i_star = parameters[2]
    k_f = parameters[3]
    k_s = parameters[4]
    k_l = parameters[5]
    
    k_m = parameters[6]
    k_r = parameters[7]
    
    A_EC = parameters[8]
    B_EC = parameters[9]
    C_EC = parameters[10]
    D_EC = parameters[11]
    E_EC = parameters[12]
    F_EC = parameters[13]
    
    # Calc_H -> det sur base de Q_sim et des paramètres
    # H_sim -> dét sur base de la similitude
    # C'est cette différence là qui va nous permettre de calibrer les paramètres
    Pump = Pump_Wilo()

    Pump.set_parameters(**{
    'A': A,
    'B': B,
    'Q_i_star': Q_i_star,
    'k_f': k_f,
    'k_s': k_s,
    'k_l': k_l,
    'k_m': k_m,
    'k_r': k_r,
    'A_EC': A_EC,
    'B_EC': B_EC,
    'C_EC': C_EC,
    'D_EC': D_EC,
    'E_EC': E_EC,
    'F_EC': F_EC,
    'circu': circu

    })

    Pump.model_H_sim(Q_sim)
    # sur le point de similitude
    # Calc_H = model_H_sim(Q_sim, A, B, Q_i_star, k_f, k_s, k_l)
    error_H = Pump.H_calc-H_sim[j]
    
    Pump.model_shaft(Pump.P_hydro*(n[j]/n_sim)**3, Pump.eta_hydro, Pump.Q_i*(n[j]/n_sim), Q_i_star*(n[j]/n_sim), n[j])
    # Sur le point de mesure
    # Calc_shaft = model_shaft(Calc_H[2]*(n[j]/n_sim)**3, Calc_H[3], Calc_H[1]*(n[j]/n_sim), Q_i_star*(n[j]/n_sim), n[j], k_m, k_r)
    Pump.model_motor(Pump.P_shaft, Pump.T_shaft, n[j])
    # Calc_motor = model_motor(model_EC, Calc_shaft[0], Calc_shaft[2], n[j], A_EC, B_EC, C_EC, D_EC, E_EC, F_EC)
    
    error_P = Pump.P_motor-P[j]
    #error_P = (9810*Q[j]*H[j])/P[j] - Calc_motor[0]*Calc_shaft[1]
    print(error_P, error_H)
    return(error_P*error_H)

# min | guess | max
x_0 = np.zeros(14, dtype=float)
bounds = (np.zeros(14, dtype=float), np.zeros(14, dtype=float))

bounds[0][0]=1          # A
x_0[0]      =8.89464288e+02
bounds[1][0]=10000

bounds[0][1]=0          # B
x_0[1]      =2.47558310e+01#H_max_sim
bounds[1][1]=H_max_sim*2

bounds[0][2]=0          # Q_i_star
x_0[2]      =3.38156129e-03#Q_max_sim/2
bounds[1][2]=Q_max_sim

bounds[0][3]=0          # k_f
x_0[3]      =3.08400247e+02
bounds[1][3]=np.inf

bounds[0][4]=0          # k_s
x_0[4]      =2.30939123e+00
bounds[1][4]=np.inf

bounds[0][5]=0          # k_l
x_0[5]      =1.39687488e-01#0.1
bounds[1][5]=0.5

bounds[0][6]=0          # k_m
x_0[6]      =0
bounds[1][6]=10

bounds[0][7]=0          # k_r
x_0[7]      =21
bounds[1][7]=1000

bounds[0][8:]=0          # eta_motor parameters
x_0[8:]      =5
bounds[1][8:]=100

print('coucou')
# Résolution
res = least_squares(calibration, x_0, bounds=bounds, args=(Q_sim, j))

# Résultats
parameters = res.x
A = parameters[0]
B = parameters[1]
Q_i_star = parameters[2]
k_f = parameters[3]
k_s = parameters[4]
k_l = parameters[5]

k_m = parameters[6]
k_r = parameters[7]
A_EC = parameters[8]
B_EC = parameters[9]
C_EC = parameters[10]
D_EC = parameters[11]
E_EC = parameters[12]
F_EC = parameters[13]


print("Parameters:")
print(parameters)

#PLOTS A VOIR!!!!!
############################################################################################################

# # Parity plot H
# erreur_parity_H_abs = np.zeros(len(n), dtype=float)
# erreur_parity_H_rel = np.zeros(len(n), dtype=float)
# plt.figure()
# #plt.title("Total head: parity plot")
# mediane = np.arange(25)
# plt.plot(mediane, mediane, color='C3')
# for i in range (0,len(n)):
#     H_calc = model_H_sim(Q_sim[i], A, B, Q_i_star, k_f, k_s, k_l)[0]*(n[i]/n_sim)**2
#     #print(model_H_sim(Q_sim[i], A, B, Q_i_star, k_f, k_s, k_l)[3])
#     plt.plot(H[i], H_calc, marker="o", markersize=4, markeredgecolor="C0", markerfacecolor="C0")
#     erreur_parity_H_abs[i] = abs(H[i]-H_calc)
#     erreur_parity_H_rel[i] = 100*abs(H[i]-H_calc)/H[i]

# plt.xlabel("Total head mesured [m_H2O]")
# plt.ylabel("Total head calculated [m_H2O]")
# plt.xlim(0,25)
# plt.ylim(0,25)
# plt.grid()
# plt.savefig("parity_circu_H_Wang.png", dpi=300)

# print ("\nErreur max sur H (abs): %.2f [m_H2o]" %(max(erreur_parity_H_abs)))
# print ("Erreur max sur H (rel): %.2f %%" %(max(erreur_parity_H_rel)))







# # Parity plot eta
# erreur_parity_P_abs = np.zeros(len(n), dtype=float)
# erreur_parity_P_rel = np.zeros(len(n), dtype=float)
# plt.figure()
# #plt.title("Total power consumed: parity plot")
# mediane = np.arange(1500)
# plt.plot(mediane, mediane, color="C3")
# for i in range (0,len(n)):
    
#     Calc_H = model_H_sim(Q_sim[i], A, B, Q_i_star, k_f, k_s, k_l)
#     Calc_shaft = model_shaft(Calc_H[2]*(n[i]/n_sim)**3, Calc_H[3], Calc_H[1]*(n[i]/n_sim), Q_i_star*(n[i]/n_sim), n[i], k_m, k_r)
#     #print(Calc_shaft[0])
#     #print(Calc_shaft[1])
#     Calc_motor = model_motor(model_EC, Calc_shaft[0], Calc_shaft[2], n[i], A_EC, B_EC, C_EC, D_EC, E_EC, F_EC)

#     P_calc = Calc_motor[0]

#     plt.plot(P[i], P_calc, marker="o", markersize=4, markeredgecolor="C0", markerfacecolor="C0")
#     erreur_parity_P_abs[i] = abs(P[i]-P_calc)
#     erreur_parity_P_rel[i] = 100*abs(P[i]-P_calc)/P[i]

# plt.xlim(0,1600)
# plt.ylim(0,1600)
# plt.xlabel("Power mesured [W]")
# plt.ylabel("Power calculated [W]")
# plt.grid()
# plt.savefig("parity_circu_P_Wang.png", dpi=300)

# print ("Erreur max sur P (abs): %.2f [W]" %(max(erreur_parity_P_abs)))
# print ("Erreur max sur P (rel): %.2f %%" %(max(erreur_parity_P_rel)))










############################################################################################################


#GRAPHE A REVOIR (PARITY PLOT), ... A voir!!!






















# # Graphe

# nbr_pts = round(Q_max_sim*3600)+1 # Nouveaux points simulés (pour n=n_sim pour que les paramètres correspondent)

# Q_calibration = np.zeros(nbr_pts, dtype=float)
# H_calibration = np.zeros(nbr_pts, dtype=float)

# for i in range (0,nbr_pts):
#     Q_calibration[i] = i/3600
    
#     Q_i = Q_calibration[i]+Q_l
#     H_i = A*Q_i+B
#     epsilon_f = k_f*Q_i**2
#     epsilon_s = k_s*((Q_i_star-Q_i)/Q_i_star)**2
#     H_calibration[i] = H_i-epsilon_f-epsilon_s
    
# plt.figure()
# plt.title("Total head")
# plt.plot(Q_calibration*3600, H_calibration)
# for i in range (0,nbr_pts):
#     plt.plot(Q_calibration[i]*(1000/n_sim)*3600, H_calibration[i]*(1000/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
#     plt.plot(Q_calibration[i]*(1500/n_sim)*3600, H_calibration[i]*(1500/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
#     plt.plot(Q_calibration[i]*(2000/n_sim)*3600, H_calibration[i]*(2000/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
#     plt.plot(Q_calibration[i]*(2668/n_sim)*3600, H_calibration[i]*(2668/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
#     plt.plot(Q_calibration[i]*(3000/n_sim)*3600, H_calibration[i]*(3000/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
#     plt.plot(Q_calibration[i]*(3450/n_sim)*3600, H_calibration[i]*(3450/n_sim)**2, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")
# for i in range (0,len(n)):
#     plt.plot(Q[i]*3600, H[i], marker="o", markersize=2, markeredgecolor="red", markerfacecolor="red")
#     plt.plot(Q_sim[i]*3600, H_sim[i], marker="o", markersize=1, markeredgecolor="green", markerfacecolor="green")

# plt.xlabel("Flow [m3/h]")
# plt.ylabel("Total head [m_H2O]")
# plt.xlim(0,30)
# plt.ylim(0,25)









# # Graphe Puissance global = f(Q) pour différents n + points de mesure + courbe de similarité
# P_calibration = np.zeros(nbr_pts, dtype=float)
# eta_pump_calibration = np.zeros(nbr_pts, dtype=float)
# for i in range (1,nbr_pts):
#     Q_i = Q_calibration[i]+Q_l
#     H_i = A*Q_i+B
    
#     # eta_m = k_m
#     P_i = 9810*Q_i*H_i
#     P_r = k_r*((Q_i_star-Q_i)/Q_i_star)**2
#     eta_r = P_i/(P_i+P_r)
#     # eta_v = (Q_i-Q_l)/Q_i
#     # eta_i = H_calibration[i]/H_i
#     # P_d = k_f*Q_i**2 # epsilon_f
#     # eta_d = P_i/(P_i+P_d)
    
#     # eta_pump_calibration[i] = eta_m*eta_r*eta_v*eta_i*eta_d # pour (Q_sim;n_sim), et donc (Q[j];n[j])
#     eta_pump_calibration[i] = (H_calibration[i]/H_i)*(Q_calibration[i]/Q_i)#*eta_r
    
#     P_hydro = 9810*Q_calibration[i]*H_calibration[i]
#     P_shaft = P_hydro/eta_pump_calibration[i]
#     T_shaft = P_shaft/(n_sim*(1/60)*2*np.pi) # pour Q_sim, n_sim
    
#     if model_EC == 1:
#         ratio_n = n_sim/B_EC
#         ratio_T = T_shaft/E_EC
#         eta_motor = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)
    
#     elif model_EC == 2:
#         ratio_n = n_sim/B_EC
#         ratio_T = T_shaft/E_EC
#         eta_motor = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)**2
    
#     elif model_EC == 3:
#         eta_motor = A_EC*1e-5*(P_shaft**(1/3))**5 - B_EC*1e-4*(P_shaft**(1/3))**4 + C_EC*1e-3*(P_shaft**(1/3))**3 - D_EC*1e-2*(P_shaft**(1/3))**2 + E_EC*1e-1*(P_shaft**(1/3))**1 - F_EC*1e-1*(P_shaft**(1/3))**0
        
#     elif model_EC == 4:
#         eta_motor = A_EC*n_sim + B_EC
#     else:
#         pass
    
#     P_calibration[i] = P_hydro/(eta_motor*eta_pump_calibration[i])



# plt.figure()
# plt.title("Total power consumed")
# plt.plot(Q_calibration[1:len(Q_calibration)], P_calibration[1:len(Q_calibration)])

# n_unique = np.unique(n)
# for i in range (0,len(n_unique)):
#     for ii in range (1,nbr_pts):
    
#         P_hydro = 9810*Q_calibration[ii]*H_calibration[ii]*(n_unique[i]/n_sim)**3
#         P_shaft = P_hydro/eta_pump_calibration[ii]
#         T_shaft = P_shaft/(n_unique[i]*(1/60)*2*np.pi) # pour Q_sim, n_sim
    
#         if model_EC == 1:
#             ratio_n = n_unique[i]/B_EC
#             ratio_T = T_shaft/E_EC
#             eta_motor = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)
    
#         elif model_EC == 2:
#             ratio_n = n_unique[i]/B_EC
#             ratio_T = T_shaft/E_EC
#             eta_motor = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)**2
    
#         elif model_EC == 3:
#             eta_motor = A_EC*1e-5*(P_shaft**(1/3))**5 - B_EC*1e-4*(P_shaft**(1/3))**4 + C_EC*1e-3*(P_shaft**(1/3))**3 - D_EC*1e-2*(P_shaft**(1/3))**2 + E_EC*1e-1*(P_shaft**(1/3))**1 - F_EC*1e-1*(P_shaft**(1/3))**0
        
#         elif model_EC == 4:
#             eta_motor = A_EC*n_unique[i] + B_EC
#         else:
#             pass

#         plt.plot(Q_calibration[ii]*(n_unique[i]/n_sim), P_hydro/(eta_motor*eta_pump_calibration[ii]), marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")

# for i in range (0, len(n)):
#     plt.plot(Q[i], P[i], marker="o", markersize=2, markeredgecolor="red", markerfacecolor="red")

# plt.xlabel("Flow [m3/h]")
# plt.ylabel("Power [W]")
# #plt.xlim(0,30)
# #plt.ylim(0,1600)







# # Graphe efficacity globale brut
# plt.figure()
# plt.title("Total efficiency mesured")
# for i in range (0,len(n)):
#     plt.plot(Q[i], 9810*Q[i]*(1/3600)*H[i]*(1/P[i])*100, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")

# plt.xlim(0,30)
# plt.ylim(0,100)
# plt.xlabel("Flow [m3/h]")
# plt.ylabel("Efficiency [%]")


# # Graphe efficacité global = f(Q) pour différents n + points de mesure + points de similarité + courbe de similarité
# plt.figure()
# plt.title("Total efficiency calculated")
# #plt.plot(Q_calibration[1:len(Q_calibration)], eta_pump_calibration[1:len(Q_calibration)])
# for i in range (0,len(n_unique)):
#     for ii in range (1,nbr_pts):
#         P_shaft_calibration = 9810*Q_calibration[ii]*(1/3600)*H_calibration[ii]*(n_unique[i]/n_sim)**3*(1/eta_pump_calibration[ii])
#         T_shaft_calibration = P_shaft_calibration/(n_unique[i]*(1/60)*2*np.pi)
        
#         if model_EC == 1:
#             ratio_n = n_unique[i]/B_EC
#             ratio_T = T_shaft_calibration/E_EC
#             eta_motor_calibration = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)
            
#         elif model_EC == 2:
#             ratio_n = n_unique[i]/B_EC
#             ratio_T = T_shaft_calibration/E_EC
#             eta_motor_calibration = A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)**2
        
#         elif model_EC == 3:
#             eta_motor_calibration = A_EC*1e-5*(P_shaft_calibration**(1/3))**5 - B_EC*1e-4*(P_shaft_calibration**(1/3))**4 + C_EC*1e-3*(P_shaft_calibration**(1/3))**3 - D_EC*1e-2*(P_shaft_calibration**(1/3))**2 + E_EC*1e-1*(P_shaft_calibration**(1/3))**1 - F_EC*1e-1*(P_shaft_calibration**(1/3))**0
            
#         elif model_EC == 4:
#             eta_motor_calibration = A_EC*n_unique[i] + B_EC
        
#         else:
#             pass
            
#         plt.plot(Q_calibration[ii]*(n_unique[i]/n_sim), eta_pump_calibration[ii]*eta_motor_calibration*100, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")


# for i in range (0,len(n)):
#     plt.plot(Q[i], 9810*Q[i]*(1/3600)*H[i]*(1/P[i])*100, marker="o", markersize=1, markeredgecolor="red", markerfacecolor="red")
    
# plt.xlim(0,30)
# plt.ylim(0,100)
# plt.xlabel("Flow [m3/h]")
# plt.ylabel("Efficiency [%]")
    
    


# # Graphe de eta pump = f(Q) pour différents n -> cohérence littérature
# plt.figure()
# plt.title("Pump efficiency")
# plt.plot(Q_calibration[1:len(Q_calibration)], 100*eta_pump_calibration[1:len(Q_calibration)])
# for i in range (0,len(n_unique)):
#     for ii in range (1,nbr_pts):
#         plt.plot(Q_calibration[ii]*(n_unique[i]/n_sim), eta_pump_calibration[ii]*100, marker="o", markersize=1, markeredgecolor="blue", markerfacecolor="blue")

# plt.xlim(0,30)
# plt.ylim(0,100)
# plt.xlabel("Flow [m3/h]")
# plt.ylabel("Efficiency [%]")




# # Graphe de eta motor = f(n,(T)) -> cohérence littérature

# if model_EC == 1 or model_EC == 2:
#     delta_mapping = 100
#     n_mapping = np.linspace(1000, 4000, delta_mapping)
#     T_mapping = np.linspace(0, 5, delta_mapping)
#     n_Mapping, T_Mapping = np.meshgrid(n_mapping, T_mapping)

#     eta_mapping = np.zeros((delta_mapping,delta_mapping), dtype=float)
#     for i in range(0,delta_mapping):
#         for ii in range(0,delta_mapping):
#             ratio_n = n_Mapping[i][ii]/B_EC
#             ratio_T = T_Mapping[i][ii]/E_EC
#             if model_EC == 1:
#                 eta_mapping[i][ii] =(A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1))*100
#             if model_EC == 2:
#                 eta_mapping[i][ii] =(A_EC*(ratio_n-1)**2+C_EC+D_EC*(ratio_T-1)**2)*100

#     lev_3  = [10,20,30,40,50,60,70,80,90,92,94,96]#,250]#400,500,600,700,800,900,1000]
#     #lev_3bis  = [90,91,92,93,94,95,96,97]#,250]#400,500,600,700,800,900,1000]

#     plt.figure()
#     plt.title("Motor efficiency")

#     CS2 = plt.contour(n_Mapping, T_Mapping, eta_mapping, cmap='binary_r', levels=lev_3)
#     plt.clabel(CS2, fontsize=8, fmt='%1.0f'+"%%")#, manual=manual_locations)
#     CS = plt.contourf(n_Mapping, T_Mapping, eta_mapping, cmap='bone_r', levels=lev_3, extend='max')

#     plt.xlabel("Speed [rpm]")
#     plt.ylabel("Torque [Nm]")
#     plt.grid(linewidth=0.5, which='both')

# elif model_EC == 3:
#     P_shaft_mapping = np.linspace(3.5, 11.5, 100) # P_shaft**(1/3)
#     eta_mapping = np.zeros(len(P_shaft_mapping), dtype=float)
#     for i in range (0, len(P_shaft_mapping)):
#         eta_mapping[i] = A_EC*1e-5*(P_shaft_mapping[i])**5 - B_EC*1e-4*(P_shaft_mapping[i])**4 + C_EC*1e-3*(P_shaft_mapping[i])**3 - D_EC*1e-2*(P_shaft_mapping[i])**2 + E_EC*1e-1*(P_shaft_mapping[i])**1 - F_EC*1e-1*(P_shaft_mapping[i])**0
    
#     plt.figure()
#     plt.title("Motor efficiency")
#     plt.plot(P_shaft_mapping, eta_mapping*100)
    
#     plt.xlim(3.5,12)
#     plt.ylim(0,100)
#     plt.xlabel("Shaft power**(1/3) [W**1/3]")
#     plt.ylabel("Efficiency [%]")

# else:
#     pass





###############################################################################




















# Q/P[bar] -> Power consumed
# refaire un parity plot du modèle

# refaire pts de mesure
# pompe 2


