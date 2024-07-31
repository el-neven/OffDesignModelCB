from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

class ThermoclineStorage:
    def __init__(self, H, A, fill_material, z_init):
        self.H = H
        self.A = A
        self.fill_material = fill_material
        self.last_z_f = z_init  # initialize the level of the tank
        self.z_values = []
        self.U_values = []
        self.time_values = []
        self.total_time = 0

    def Mass_Energy_Balance(self, t, y):
        "Mass balance"
        dM_dt = self.m_dot

        return [dM_dt]

    def charge(self, m_dot_in, t_span, T_hot_in, T_cold_out):
        self.m_dot = m_dot_in
        self.T_hot_in = T_hot_in
        self.T_cold_out = T_cold_out

        "Initial state"
        self.z_i = self.last_z_f  # use the last z_f value as z_i
        Vol_cold_i = self.A * (self.z_i*self.H)  # Volume filled with cold water
        Vol_hot_i = self.A * self.H *(1- self.z_i)  # Volume filled with hot water
        # rho_hot = PropsSI('D', 'T', T_hot_in + 273.15, 'P', 1e5, 'water')
        # rho_cold = PropsSI('D', 'T', T_cold_out + 273.15, 'P', 1e5, 'water')
        rho_w = 1000 # Faux, il faut prendre le changemetn de densité en compte MAIS pas encore fait car aussinon le V_tot n'est pas conservé!!
        self.M_hot_i = Vol_hot_i * rho_w
        rho_cold = PropsSI('D', 'T', T_cold_out + 273.15, 'P', 1e5, 'water')
        self.M_cold_i = Vol_cold_i * rho_w

        u_cold = PropsSI('U', 'T', T_cold_out + 273.15, 'P', 1e5, 'water')
        u_hot = PropsSI('U', 'T', T_hot_in + 273.15, 'P', 1e5, 'water')
        self.U_tot_i = self.M_hot_i * u_hot + self.M_cold_i * u_cold

        "Integral"
        y0 = [0]  # No initial conditions
        t_eval = np.arange(t_span[0], t_span[1]+1, 1)  # specify time step of 1 second
        result = solve_ivp(self.Mass_Energy_Balance, t_span, y0, method='RK45', t_eval=t_eval)
        time_step = result.t[1] - result.t[0]
        time_step2 = result.t[-1] - result.t[-2]
        self.result = result

        DELTA_M = result.y[0]  # Total mass change at each time step

        "Final state"
        Vol_tot = self.A * self.H
        M_hot_max = Vol_tot * rho_w
        self.M_hot = self.M_hot_i + DELTA_M
        self.M_hot = np.where(self.M_hot >= M_hot_max, M_hot_max, self.M_hot)

        self.M_cold = self.M_cold_i - DELTA_M
        self.M_cold = np.where(self.M_cold <= 0, 0, self.M_cold)

        Vol_cold = np.where(self.M_cold == 0, 0, self.M_cold / rho_w)
        self.U_tot = self.M_hot * PropsSI('U', 'T', T_hot_in + 273.15, 'P', 1e5, 'water') + self.M_cold * PropsSI('U', 'T', T_cold_out + 273.15, 'P', 1e5, 'water')
        self.z = np.where(Vol_cold == 0, 0, Vol_cold / (self.A*H))
        # Vol_tot = self.A * self.H
        # Vol_hot = Vol_tot - Vol_cold
        self.last_z_f = self.z[-1]
        self.E = (self.U_tot - self.U_tot_i) * (result.t / 3600)
        self.E_f = self.E[-1]

        # Store results
        self.z_values.extend(self.z)
        self.U_values.extend(self.U_tot)
        self.time_values.extend(result.t + self.total_time)
        self.total_time += result.t[-1]

    def discharge(self, m_dot_out, t_span, T_hot_out, T_cold_in):
        self.m_dot = m_dot_out
        self.T_hot_out = T_hot_out
        self.T_cold_in = T_cold_in

        "Initial state"
        self.z_i = self.last_z_f  # use the last z_f value as z_i
        Vol_cold_i = self.A * (self.z_i*self.H)  # Volume filled with cold water
        Vol_hot_i = self.A * self.H *(1- self.z_i)  # Volume filled with hot water
        # rho_hot = PropsSI('D', 'T', T_hot_out + 273.15, 'P', 1e5, 'water')
        # rho_cold = PropsSI('D', 'T', T_cold_in + 273.15, 'P', 1e5, 'water')
        rho_w = 1000 #[kg/m3]
        self.M_hot_i = Vol_hot_i * rho_w
        self.M_cold_i = Vol_cold_i * rho_w

        u_cold = PropsSI('U', 'T', T_cold_in + 273.15, 'P', 1e5, 'water')
        u_hot = PropsSI('U', 'T', T_hot_out + 273.15, 'P', 1e5, 'water')
        self.U_tot_i = self.M_hot_i * u_hot + self.M_cold_i * u_cold

        "Integral"
        y0 = [0]  # No initial conditions
        t_eval = np.arange(t_span[0], t_span[1]+1, 1)  # specify time step of 1 second
        result = solve_ivp(self.Mass_Energy_Balance, t_span, y0, method='RK45', t_eval=t_eval)
        time_step = result.t[1] - result.t[0]
        time_step2 = result.t[-1] - result.t[-2]
        self.result = result

        DELTA_M = result.y[0]  # Total mass change at each time step

        "Final state"
        self.M_hot = self.M_hot_i - DELTA_M
        self.M_hot = np.where(self.M_hot <= 0, 0, self.M_hot)

        Vol_tot = self.A * self.H
        M_cold_max = Vol_tot * rho_w
        self.M_cold = self.M_cold_i + DELTA_M
        self.M_cold = np.where(self.M_cold >= M_cold_max, M_cold_max, self.M_cold)

        Vol_cold = np.where(self.M_cold == 0, 0, self.M_cold / rho_w)

        self.U_tot = self.M_hot * u_hot + self.M_cold * u_cold
        self.z = np.where(Vol_cold == 0, 0, Vol_cold / (self.A*H))
        self.last_z_f = self.z[-1]
        self.E = (self.U_tot - self.U_tot_i) / (3600000) #[kWh]
        self.E_f = self.E[-1]

        # Store results
        self.z_values.extend(self.z)
        self.U_values.extend(self.U_tot)
        self.time_values.extend(result.t + self.total_time)
        self.total_time += result.t[-1]

    def plot_results(self):
        # Plot z over time
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        plt.plot(self.time_values, self.z_values)
        plt.xlabel('Time (s)')
        plt.ylabel('z (m)')
        plt.title('Height of Cold Water Over Time')

        # Plot U over time
        plt.subplot(1, 2, 2)
        plt.plot(self.time_values, self.U_values)
        plt.xlabel('Time (s)')
        plt.ylabel('Internal Energy (J)')
        plt.title('Total Internal Energy Over Time')

        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    "Inputs"
    H = 3.48  # [m]
    # H = 3.1
    D = 1.8  # [m]
    Area = np.pi * D**2 / 4
    # Area = 2
    print("Vol", Area*H)
    m_dot_w = 2.5  # [kg/s]
    T_hot = 76  # [°C]
    T_cold = 62  # [°C]

    "Time"
    t_span1 = (0, 60 * 60)  # [s]
    t_span2 = (0, 15 * 60)

    TC = ThermoclineStorage(H, Area, 'water', z_init=1)

    TC.charge(m_dot_w, t_span1, T_hot, T_cold)
    # # print(TC.E_f)
    # # print(TC.z)
    # # print(TC.z[-1])
    # # print(TC.U_tot[-1])
    # # print(TC.M_cold_i)
    # # print(TC.M_hot_i)
    # # print(TC.M_cold[-1])
    # # print(TC.M_hot[-1])
    # TC.charge(3, t_span2, T_hot, T_cold)
    # print(TC.z[-1])
    # print(TC.U_tot[-1])
    # print(TC.M_cold_i)
    # print(TC.M_hot_i)
    # print(TC.M_cold[-1])
    # print(TC.M_hot[-1])
    # print(TC.E_f)
    # print(TC.z)

    # "Discharge"
    # T_hot_out = 76  # [°C]
    # T_cold_in = 62  # [°C]

    # t_span3 = (0, 15 * 60)  # [s]

    # TC.discharge(7, t_span3, T_hot_out, T_cold_in)
    # print(TC.E_f)
    # print(TC.E)
    # print(TC.z)
    # Pour la décharge thermique, il faut donc un beaucoup plus grand débit massique pour obtenir la même énergie que pour la charge thermique
    TC.plot_results()


