##------------------------------ Compressor model -----------------------------------

1. Model description

Semi-empirical model developped by Vincent Lemort

- Inputs: T_su [K] or h_su, P_su [Pa], rp [-], N_rot [RPM] OR m_dot [kg/s]

 Rem: j'ai mit rp à la place de P_ex pour pouvoir tester le code avec plusieurs rp sans refaire une nouvelle classe à chaque fois mais c'est facilement interchangeable.

-Design parameters: AU_amb [W/K], AU_su_n [W/K], AU_ex_n [W/K], d_sex [m], m_dot_n [kg/s], A_leak [m^2], W_dot_loss_0 [W], alpha [-], C_loss [Nm], T_amb [K]
-Control parameters: rv_in [-], V_s [m^3]

-Ouputs: epsilon_is [-], W_dot_sh [W], T_ex [K], m_dot [kg/s] OR N_rot [RPM]


Iteration time: more or less 0.3 sec

4 different code:

	- Compressor_SE: Contains the Compressor class
	- Calibration_compressor: Contains the minimization code to determine the parameters based on the experimental data
	- Calibration_graphs: parity plots of the measured values vs the calculated ones
	- Compressor_test: example of how to use the semi-empirical model with a plot of epsilon_vs_rp

REM: J'ai dû adapter le swept volume car le débit ne suivait pas du tout aussinon. 


2. Exeprimental datas

Datas not precise at all, for expample just given a superheating constant for all the experiments, the pressure are fixed at round number (2 bar, 3 bar,...). Once again the mass flow rate fits very well while the work at the shaft not so well.


#---------------------------------------------------------------------------------------------------
TO DO:

# Adapater le calibration graphs et calibration compressor avec la structure updatée