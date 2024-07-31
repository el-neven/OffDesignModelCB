
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

class HeatExchanger:
    def __init__(self, HX_type, HX_D, HX_A, min_pinch, TQ_flag):
        self.HX_type = HX_type
        self.HX_D = HX_D
        self.HX_A = HX_A
        self.min_pinch = min_pinch
        self.TQ_flag = TQ_flag

    def solve(self, m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf, P_sf, T_sf_in, glide_sf, fluid_sf, x_eva_in):
        # Define results dictionary
        results = {}

        # Resolution of the problem - iterative variable P_wf
        lb, ub = 0, 0
        if self.HX_type == 'evaporator':
            lb = PropsSI('P', 'T', T_wf_in, 'Q', 0, fluid_wf)
            ub = PropsSI('P', 'T', T_sf_in - dT_sub_or_sup - self.min_pinch, 'Q', 0, fluid_wf)
            sf_plot_color = 'b'
        elif self.HX_type == 'condenser':
            lb = PropsSI('P', 'T', T_sf_in + dT_sub_or_sup + self.min_pinch, 'Q', 0, fluid_wf)
            ub = PropsSI('P', 'T', T_wf_in, 'Q', 0, fluid_wf)
            sf_plot_color = 'r'
            print(lb, ub)

        # Minimization problem
        if abs(ub - lb) < 10 or ub < lb:
            P_wf0 = (lb + ub) / 2
        else:
            P_wf0 = np.arange(lb, ub + 1, 150)
        # print(P_wf0, ub, lb)
        res_temp = []
        flag_tau = []
        for P_wf in P_wf0:
            results_temp = MovingBoundariesHX(P_wf, self.HX_type, self.HX_D, self.HX_A,
                                              m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf,
                                              P_sf, T_sf_in, glide_sf, fluid_sf, x_eva_in)
            # print(P_wf)
            # flag_tau.append(results_temp.flag_tau)
            res_temp.append(results_temp.res)
            # P_wf_array = np.array(P_wf0)
        print(res_temp)
        # res_temp = np.array(res_temp)
        # res_temp = res_temp[flag_tau != 1]

        # if len(res_temp) == 0:  # Check compatibility of boundary conditions
        #     results['P_wf'] = np.nan
        #     results['flag_check'] = 0
        # else:
        #     P_wf0 = P_wf0[flag_tau != 1]
        #     min_indx = np.argmin(res_temp)
        #     results = MovingBoundariesHX(P_wf0[min_indx], self.HX_type, self.HX_D, self.HX_A,
        #                                  m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf,
        #                                  P_sf, T_sf_in, glide_sf, fluid_sf, x_eva_in)
        #     results['P_wf'] = P_wf0[min_indx]
        #     results['flag_check'] = 1

        # Heat transfer diagram
        # if self.TQ_flag == 1:
        #     fig, ax = plt.subplots()
        #     ax.set_ylabel('Temperature (°C)')
        #     ax.set_xlabel('Termal power (kW)')
        #     ax.set_title(self.HX_type + ' heat transfer diagram')
        #     ax.grid(True)
        #     ax.plot(results['plot']['Q_wf'] / 1e3, results['plot']['T_wf'] - 273.15, 'g', linewidth=2)
        #     ax.plot(results['plot']['Q_sf'] / 1e3, results['plot']['T_sf'] - 273.15, sf_plot_color, linewidth=2)
        #     ax.text(np.mean(results['plot']['Q_wf'] / 1e3), results['plot']['T_wf'][1] - 273.15,
        #             str(round(results['P_wf'] / 100)) + ' bar', backgroundcolor='w')
        #     ax.text(np.mean(results['plot']['Q_sf'] / 1e3), np.mean(results['plot']['T_sf']) - 273.15,
        #             str(round(results['m_dot_sf'], 1)) + ' l/s', backgroundcolor='w')

        # return results


class MovingBoundariesHX:
    def __init__(self, P_wf, HX_type, HX_D, HX_A, m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf, P_sf, T_sf_in,
                 glide_sf, fluid_sf, x_eva_in):
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
        self.glide_sf = glide_sf
        self.fluid_sf = fluid_sf
        self.x_eva_in = x_eva_in

        self.solve()

    def solve(self):
        # TEMPERATURES AND THERMAL POWERS
        # evaporating/condensing phase - wf side
        T_wf_sat = PropsSI('T', 'P', self.P_wf, 'Q', 0, self.fluid_wf)
        h_wf_sat_vap = PropsSI('H', 'P', self.P_wf, 'Q', 1, self.fluid_wf)

        if 0 <= self.x_eva_in < 1:
            h_wf_sat_liq = PropsSI('H', 'P', self.P_wf, 'Q', self.x_eva_in, self.fluid_wf)
            T_wf_in = T_wf_sat
        else:
            h_wf_sat_liq = PropsSI('H', 'P', self.P_wf, 'Q', 0, self.fluid_wf)

        if self.HX_type == 'evaporator':
            sf_ht_dir = 'cooled'
            wf_ht_dir = 'heated'

            T_sf_h = self.T_sf_in
            T_sf_c = T_sf_h - self.glide_sf
            T_wf_h = T_wf_sat + self.dT_sub_or_sup
            T_wf_c = self.T_wf_in
            self.T_wf_out = T_wf_h
        elif self.HX_type == 'condenser':
            sf_ht_dir = 'heated'
            wf_ht_dir = 'cooled'

            T_sf_c = self.T_sf_in
            T_sf_h = T_sf_c + self.glide_sf
            T_wf_c = T_wf_sat - self.dT_sub_or_sup
            T_wf_h = self.T_wf_in
            self.T_wf_out = T_wf_c

        # Energy balance - wf side
        h_wf_h = PropsSI('H', 'T', T_wf_h, 'P', self.P_wf, self.fluid_wf)
        

if __name__ == "__main__":
    # Define parameters
    param = {
        'HX_type': 'condenser',
        'HX_D': 0.05,
        'HX_A': 1.0,
        'min_pinch': 5,
        'TQ_flag': 1
    }

    # Define fluid properties
    fluid_wf = 'Water'
    fluid_sf = 'R134a'

    # Define other inputs
    m_dot_wf = 0.1  # kg/s
    T_wf_in = 50+273.15  # °C
    dT_sub_or_sup = 5  # °C
    P_sf = 2e5  # Pa
    T_sf_in = 20+273.15  # °C
    glide_sf = 20  # °C
    x_eva_in = 0  # -

    # Create HeatExchanger object
    hx = HeatExchanger(**param)

    # Solve and get results
    results = hx.solve(m_dot_wf, T_wf_in, dT_sub_or_sup, fluid_wf, P_sf, T_sf_in, glide_sf, fluid_sf, x_eva_in)

    # Print results
    print("Results:")
    print(results)
