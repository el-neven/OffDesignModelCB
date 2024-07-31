import numpy as np
import math
from scipy.constants import pi
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

def conducticity_R1233zd(T, P):
    """ Richard A. Perkins and Marcia L. Huber

        "Measurement and Correlation of the Thermal Conductivity of trans-1-Chloro-3,3,3-trifluoropropene (R1233zd(E))"
        J. Chem. Eng. Data 2017, 62, 2659-2665
        Note that the critical enhancement contribution is not implemented.

        Same source as EES
       """
    # Dilute gas thermal conductivity
    T_c = 382.52 #[K]
    A_0 = -0.0103589
    A_1 = 0.0308929
    A_2 = 0.000230348

    k_0 = A_0 + A_1*(T/T_c) + A_2*(T/T_c)**2

    # Residual gas thermal conductivity
    try:
        Rho = PropsSI('D', 'T', T, 'P', P, 'R1233zd(E)')
    except:
        Rho = PropsSI('D', 'T', T, 'Q', 0, 'R1233zd(E)')
    Rho_c = 489.24 #[kg/m^3]

    B_11 = -0.0428296
    B_12 = 0.0434288
    B_21 = 0.0927099
    B_22 = -0.0605844
    B_31 = -0.0702107
    B_32 = 0.0440187
    B_41 = 0.0249708
    B_42 = -0.0155082
    B_51 = -0.00301838
    B_52 = 0.0021019

    Delta_kr = (B_11 + B_12*(T/T_c))*Rho/Rho_c + (B_21 + B_22*(T/T_c))*(Rho/Rho_c)**2 + (B_31 + B_32*(T/T_c))*(Rho/Rho_c)**3 + (B_41 + B_42*(T/T_c))*(Rho/Rho_c)**4 + (B_51 + B_52*(T/T_c))*(Rho/Rho_c)**5

    k = k_0 + Delta_kr
    return k


def U_Gnielinski_calibrated(m, D, fluid, P):
    T = PropsSI('T', 'P', P, 'Q', 0, fluid)
    mu = PropsSI('V', 'Q', 0, 'P', P, fluid)
    cp = PropsSI('CPMASS', 'Q', 0, 'P', P, fluid)
    if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
        k = conducticity_R1233zd(T, P)
    else:
        k = PropsSI('L', 'Q', 0, 'P', P, fluid)
    Re = m * 4 / (pi * mu * D)
    Pr = mu * cp / k
    f = (0.79 * np.log(Re) - 1.64) ** (-2)
    c = 4.616163048309070 # calibrated for old model
    # c = 6
    Nu = c * (f / 8 * (Re - 1000) * Pr) / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    U = k * Nu / D
    return U

def U_DittusBoelter(m, heat_transfer_direction, D, fluid, P, T):

    if heat_transfer_direction == 'heated':
        a = 0.4
    elif heat_transfer_direction == 'cooled':
        a = 0.5

    if T is not None:
        mu = PropsSI('V', 'T', T, 'P', P, fluid)
        cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
        if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
            k = conducticity_R1233zd(T, P)
        else:
            k = PropsSI('L', 'T', T, 'P', P, fluid)
        Re = m * 4 / (pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re ** 0.8 * Pr ** a / D
    elif T is None:
        mu = (PropsSI('V', 'P', P, 'Q', 0, fluid) + PropsSI('V', 'P', P, 'Q', 1, fluid)) / 2
        cp = (PropsSI('C', 'P', P, 'Q', 0, fluid) + PropsSI('C', 'P', P, 'Q', 1, fluid)) / 2
        k = (PropsSI('L', 'P', P, 'Q', 0, fluid) + PropsSI('L', 'P', P, 'Q', 1, fluid)) / 2
        Re = m * 4 / (pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re ** 0.8 * Pr ** a / D
    return U

def U_Cooper_calibrater(Q, HX_A, P, fluid):
    pr = P/PropsSI('Pcrit', '', 0, ' ', 0, fluid)
    M = PropsSI('MOLARMASS', fluid)
    c = [0.978552748683140, 1.07700234277466] # calibrated for old HX model
    # c = [3, 1.5]
    U = c[0] * ((Q / HX_A)**(0.67 * c[1])) * 55 * (pr**(0.12 - 0.2 * np.log(pr))) * ((-np.log(pr))**(-0.55 * c[1])) * (M**(-0.5))
    return U

def U_Thonon(m, D, fluid, P, T):
    mu = PropsSI('V', 'T', T, 'P', P, fluid)
    cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
    if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
        k = conducticity_R1233zd(T, P)
    else:
        k = PropsSI('L', 'T', T, 'P', P, fluid)
    Re = m * 4 / (pi * mu * D)
    Pr = mu * cp / k
    U = 0.2946 * k * Re ** 0.7 * Pr ** (1 / 3) / D
    return U

class HeatExchanger:
    def __init__(self, su1, su2, ex1, ex2):

        # Inputs
        self.su1 = su1 # Working fluid
        self.su2 = su2 # Secondary fluid
        self.ex1 = ex1
        self.ex2 = ex2

        # Parameters
        self.HX_type = None
        self.HX_D = None
        self.HX_A = None
        self.min_pinch = None
        self.dT_sub_or_sup = None
        # self.TQ_flag = TQ_flag

    def set_parameters(self, **parameters):
        for key in parameters:
            setattr(self, key, parameters[key])

    def solve(self):
        # results: P_wf - T_wf_out - Q - m_dot_sf - Ai - Ui (i=sub-ssat-sup)

        ## Resolution of the problem - iterative variable P_wf

        # Get data from class
        # Working fluid
        T_wf_in = self.su1.T
        fluid_wf = self.su1.fluid
        m_dot_wf = self.su1.m_dot
        x_eva_in = self.su1.x

        # Secondary fluid
        T_sf_in = self.su2.T
        fluid_sf = self.su2.fluid
        P_sf = self.su2.p
        T_sf_out = self.ex2.T


        if self.HX_type=='evaporator':
            lb = PropsSI('P', 'T', T_wf_in, 'Q', 0, fluid_wf)
            # print(T_wf_in, 'T_wf_in')
            ub = PropsSI('P', 'T', T_sf_in-self.dT_sub_or_sup-self.min_pinch, 'Q', 0, fluid_wf)
            # ub = PropsSI('P', 'T', T_sf_in, 'Q', 0, fluid_wf)
            # print(T_wf_in, T_sf_in-self.dT_sub_or_sup-self.min_pinch)
            sf_plot_color = 'b'
        if self.HX_type=='condenser':
            lb = PropsSI('P', 'T', T_sf_in+self.dT_sub_or_sup+self.min_pinch, 'Q', 0, fluid_wf)
            # lb = PropsSI('P', 'T', T_sf_in, 'Q', 0, fluid_wf)
            try:
                ub = PropsSI('P', 'T', T_wf_in, 'Q', 0, fluid_wf)
            except:
                ub = PropsSI('P_max', 'T', T_wf_in, 'Q', 0, fluid_wf)
            sf_plot_color = 'r'

            "Ici on fait l'hypothèse que le pinch se trouve à la sortie du condenseur MAIS c'est pas toujours le cas! ATTENTION"
        self.flag = 0
        # Minimization problem
        if abs(ub - lb) < 10 or (ub < lb):
            print('The temperature of the working fluid is too high or too low to reach the secondary fluid temperature.')
            P_wf_f = (lb + ub)/2
            # print(lb, ub)
            # print('P_wf_f', P_wf_f, 'Type', self.HX_type)
            self.flag = 0
        else:
            P_wf0 = np.linspace(lb, ub, 200)
            res_vector = []
            for P_wf in P_wf0:
                try:
                    results_temp = MovingBoundariesHX(P_wf, self.HX_type, self.HX_D, self.HX_A,
                                                    m_dot_wf, T_wf_in,self. dT_sub_or_sup, fluid_wf,
                                                    P_sf, T_sf_in, T_sf_out, fluid_sf, x_eva_in)
                    res = results_temp.res
                    if math.isnan(res):
                        res = 1e6
                    else:
                        res = results_temp.res
                        # Append res to the list

                    res_vector.append(res)
                    
                except:
                    res_vector.append(1e6)

                min_indx = np.argmin(res_vector)
                self.res = res_vector[min_indx]
                P_wf_f = P_wf0[min_indx]

            # print(res_vector)
            # print(self.res)
            # print(min_indx)
        print('P_wf_f', P_wf_f, 'ub', ub, 'lb', lb, 'Type', self.HX_type)
        self.results = MovingBoundariesHX(P_wf_f, self.HX_type, self.HX_D, self.HX_A,
                                            m_dot_wf, T_wf_in, self.dT_sub_or_sup, fluid_wf,
                                            P_sf, T_sf_in, T_sf_out, fluid_sf, x_eva_in)

        self.ex1.set_T(self.results.T_wf_out)
        self.ex1.set_p(self.results.P_wf)
        self.ex1.set_m_dot(m_dot_wf)
        self.ex1.set_fluid(fluid_wf)
        self.Q = self.results.Q
        # self.res = res_vector[min_indx]

        self.ex2.set_T(self.results.T_sf_out)
        self.ex2.set_p(P_sf)
        self.ex2.set_m_dot(self.results.m_dot_sf)
        self.ex2.set_fluid(fluid_sf)

    def plot_heat_transfer_diagram(self):
        plt.figure(figsize=(6, 4.2))
        plt.ylabel('Temperature (°C)', fontweight='bold')
        plt.xlabel('Thermal power (kW)', fontweight='bold')
        plt.title(self.HX_type + ' heat transfer diagram')
        plt.box(True)
        plt.grid(True, linestyle=':', linewidth=0.5, color='black')
        plt.plot(self.results.Q_wf_plot, self.results.T_wf_plot, 'g', linewidth=2)
        plt.plot(self.results.Q_sf_plot,self. results.T_sf_plot, 'r', linewidth=2)
        plt.text(sum(self.results.Q_wf_plot)/2, self.results.T_wf_plot[1], str(round(self.results.P_wf/100)) + ' bar', backgroundcolor='white')
        plt.text(sum(self.results.Q_sf_plot)/2, sum(self.results.T_sf_plot)/2, str(round(self.results.m_dot_sf, 1)) + ' l/s', backgroundcolor='white')
        plt.show()

## Moving_boundaries function
class MovingBoundariesHX:
    def __init__(self, P_wf, HX_type, HX_D, HX_A, m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf, P_sf, T_sf_in, T_sf_out, fluid_sf, x_eva_in):
        self.P_wf = P_wf
        self.HX_type = HX_type
        self.HX_D = HX_D
        self.HX_A = HX_A
        self.m_dot_wf = m_dot_wf
        self.T_wf_in = T_wf_in
        self.dT_sub_or_sup = dT_sub_or_sup
        self.fluid_wf = fluid_wf
        self.P_sf = P_sf
        self.T_sf_in = T_sf_in
        self.T_sf_out = T_sf_out
        self.fluid_sf = fluid_sf
        self.x_eva_in = x_eva_in

        self.solve()

    def solve(self):

        # TEMPERATURES AND THERMAL POWERS
        # evaporating/consensing phase - wf side
        T_wf_sat = PropsSI('T','P',self.P_wf,'Q',0,self.fluid_wf)

        # case depending temperatures
        if self.HX_type == 'evaporator':
            sf_ht_dir = 'cooled'
            wf_ht_dir = 'heated'
            T_sf_h = self.T_sf_in
            T_sf_c = self.T_sf_out    #T_sf_h - self.glide_sf
            T_wf_h = T_wf_sat + self.dT_sub_or_sup
            T_wf_c = self.T_wf_in
            self.T_wf_out = T_wf_h
        elif self.HX_type == 'condenser':
            sf_ht_dir = 'heated'
            wf_ht_dir = 'cooled'
            T_sf_c = self.T_sf_in
            T_sf_h = self.T_sf_out   #T_sf_c + self.glide_sf
            T_wf_c = T_wf_sat - self.dT_sub_or_sup
            T_wf_h = self.T_wf_in
            self.T_wf_out = T_wf_c

        # Energy balance - wf side
        if self.HX_type == 'evaporator':
            if self.x_eva_in >= 0 and self.x_eva_in <= 1:
                h_wf_sat_liq = PropsSI('H','P', self.P_wf,'Q', self.x_eva_in, self.fluid_wf)
                h_wf_c = h_wf_sat_liq
            else:
                h_wf_sat_liq = PropsSI('H','P', self.P_wf,'Q', 0, self.fluid_wf)
                h_wf_c = PropsSI('H','P', self.P_wf,'T', T_wf_c, self.fluid_wf)
            h_wf_sat_vap = PropsSI('H','P', self.P_wf,'Q', 1, self.fluid_wf)
            h_wf_h = PropsSI('H','P', self.P_wf,'T', T_wf_h, self.fluid_wf)

        if self.HX_type == 'condenser':
            if self.x_eva_in >= 0 and self.x_eva_in <= 1:
                h_wf_sat_vap = PropsSI('H','P', self.P_wf,'Q', self.x_eva_in, self.fluid_wf)
                h_wf_h = h_wf_sat_vap
            else:
                h_wf_sat_vap = PropsSI('H','P', self.P_wf,'Q', 1, self.fluid_wf)
                h_wf_h = PropsSI('H','P', self.P_wf,'T', T_wf_h, self.fluid_wf)
            h_wf_sat_liq = PropsSI('H','P', self.P_wf,'Q', 0, self.fluid_wf)
            h_wf_c = PropsSI('H','P', self.P_wf,'T', T_wf_c, self.fluid_wf)

        Q_sub = max(0, self.m_dot_wf*(h_wf_sat_liq-h_wf_c))
        Q_sat = max(0, self.m_dot_wf*(h_wf_sat_vap-h_wf_sat_liq))
        Q_sup = max(0, self.m_dot_wf*(h_wf_h-h_wf_sat_vap))

        self.Q = Q_sub + Q_sat + Q_sup

        # Energy balance - sf side
        T_sf_mean = (T_sf_h + T_sf_c)/2
        cp_sf = PropsSI('C','T', T_sf_mean,'P', self.P_sf, self.fluid_sf)
        self.m_dot_sf = self.Q/(cp_sf*(T_sf_h - T_sf_c))
        T_sf_ch = T_sf_c + Q_sub/(self.m_dot_sf*cp_sf)
        T_sf_hc = T_sf_h - Q_sup/(self.m_dot_sf*cp_sf)

        # FLUID MASS DISTRIBUTION
        # Heat tranfer coefficients 1P - sf side: Dittus_Boelter
        U_sf = U_DittusBoelter(self.m_dot_sf, sf_ht_dir, self.HX_D, self.fluid_sf, self.P_sf, T_sf_mean)
        # Heat tranfer coefficients - wf side
        if self.HX_type == 'condenser':
            try:
                self.U_wf_h = U_Thonon(self.m_dot_wf, self.HX_D, self.fluid_wf, self.P_wf, (T_wf_h + T_wf_sat)/2)
            except:
                self.U_wf_h = 0 #On arrive déjà en 2 phases dans le condenseur
            self.U_wf_c = U_Thonon(self.m_dot_wf, self.HX_D, self.fluid_wf, self.P_wf, (T_wf_c + T_wf_sat)/2)
            self.U_wf_2P = U_Cooper_calibrater(self.Q, self.HX_A, self.P_wf, self.fluid_wf)

        elif self.HX_type == 'evaporator':
            self.U_wf_h = U_Thonon(self.m_dot_wf, self.HX_D, self.fluid_wf, self.P_wf, (T_wf_h + T_wf_sat)/2)
            try:
                self.U_wf_c = U_Thonon(self.m_dot_wf, self.HX_D, self.fluid_wf, self.P_wf, (T_wf_c + T_wf_sat)/2)
            except:
                self.U_wf_c = 0 #On arrive déjà en 2 phases dans l'évaporateur
            self.U_wf_2P = U_Gnielinski_calibrated(self.m_dot_wf, self.HX_D, self.fluid_wf, self.P_wf)

        try:
            self.U_sub = (1/U_sf + 1/self.U_wf_c)**(-1)
        except:
            self.U_sub = 0

        self.U_sat = (1/U_sf + 1/self.U_wf_2P)**(-1)

        try:
            self.U_sup = (1/U_sf + 1/self.U_wf_h)**(-1)
        except:
            self.U_sup = 0 #On arrive déjà en 2 phases

        # Mean logarithm temperature difference
        if self.HX_type == 'evaporator':
            tau_sub = [(T_sf_c-T_wf_c), (T_sf_ch-T_wf_sat)]
            tau_sat = [(T_sf_ch-T_wf_sat), (T_sf_hc-T_wf_sat)]
            tau_sup = [(T_sf_hc-T_wf_sat), (T_sf_hc-T_wf_h)]
        elif self.HX_type == 'condenser':
            tau_sub = [-(T_sf_c-T_wf_c), -(T_sf_ch-T_wf_sat)]
            tau_sat = [-(T_sf_ch-T_wf_sat), -(T_sf_hc-T_wf_sat)]
            tau_sup = [-(T_sf_hc-T_wf_sat), -(T_sf_h-T_wf_h)]

        if min(tau_sub) * max(tau_sub) < 0: #If are of opposite signs
            dT_ml_sub = 0
            self.flag_tau = 1
        else:
            dT_ml_sub = (max(tau_sub)-min(tau_sub))/np.log(max(tau_sub)/min(tau_sub))
            self.flag_tau = 0
        
        if min(tau_sat) * max(tau_sat) < 0: #If are of opposite signs
            dT_ml_sat = 0
            self.flag_tau = 1
        else:
            dT_ml_sat = (max(tau_sat)-min(tau_sat))/np.log(max(tau_sat)/min(tau_sat))
        
        if min(tau_sup) * max(tau_sup) < 0: #If are of opposite signs
            dT_ml_sup = 0
            self.flag_tau = 1
        else:
            dT_ml_sup = (max(tau_sup)-min(tau_sup))/np.log(max(tau_sup)/min(tau_sup))
            self.flag_tau = 0

        self.tau = [tau_sub, tau_sat, tau_sup]

        # Energy balance UA
        if Q_sub == 0:
            self.A_sub = 0
        else:
            self.A_sub = Q_sub/(self.U_sub*dT_ml_sub)
        if Q_sat == 0:
            self.A_sat = 0
        else:
            self.A_sat = Q_sat/(self.U_sat*dT_ml_sat)
        
        if Q_sup == 0: #Déja en 2 phases
            self.A_sup = 0
        else:
            self.A_sup = Q_sup/(self.U_sup*dT_ml_sup)

        self.A_sat = abs(self.HX_A - self.A_sub - self.A_sup)

        # Closure equations
        Q_sat_calc = self.A_sat*self.U_sat*dT_ml_sat
        # print('Q_sat_calc', Q_sat_calc, 'Q_sat', Q_sat)
        self.res = abs(Q_sat-Q_sat_calc)
        # print('res', self.res, 'type', self.HX_type)

        "Results for TQ plot"
        self.Q_wf_plot = [0, Q_sub, Q_sub+Q_sat, self.Q]
        self.Q_sf_plot = [0, self.Q]
        self.T_wf_plot = [T_wf_c, T_wf_sat, T_wf_sat, T_wf_h]
        self.T_sf_plot = [T_sf_c, T_sf_h]


