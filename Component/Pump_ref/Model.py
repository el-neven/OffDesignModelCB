import numpy as np
from scipy.interpolate import RegularGridInterpolator

from CoolProp.CoolProp import PropsSI

class PumpModel:
    def __init__(self):

        self.calculable = False
        self.parametrized = False
        self.defined = False

        "Inputs"
        self.point_su = None
        self.point_ex = None

    def set_inputs(self, point_su, point_ex):
        self.point_su = point_su
        self.point_ex = point_ex

        self.check_calculable()

    def check_calculable(self):
        self.Required_inputs = [self.point_su.p, self.point_su.h, self.point_su.m_dot, self.point_ex.p]
        if all(Input is not None for Input in self.Required_inputs):
                self.calculable = True
    
    def frequency(self, Head_value, V_dot_value):
        # For the frequency w=60Hz, we can determine the Head corresponding to the volume flow rate based on the graph
        w1 = 60
        "Problème vient du V_dot_value bcp trop grand!"
        Head_w1 = -23.8 * V_dot_value + 94.76
        "Problème vient du Head_w1 négatif"
        # Using the similitudes, we can now determine the real frequency
        w = np.sqrt(Head_value/Head_w1) * w1
        return w
    
    def pump_consumption(self, w_value, V_dot_value):      
        # For the frequency w=60Hz, we can determine the Power corresponding to the volume flow rate
        w1 = 60
        W_dot_pp_w1 = -0.2 * V_dot_value + 1.8
        # Using the similitudes, we can now determine the real Power
        W_dot_pp_kW = W_dot_pp_w1 * (w_value/w1)**3
        W_dot_pp = W_dot_pp_kW*1000
        return W_dot_pp


    def solve(self):

        if self.calculable: 
            try:
                # Calculate the necessary parameters for the pump model
                m_dot = self.point_su.m_dot  # Mass flow rate at the suction point
                P_su = self.point_su.p  # Pressure at the suction point
                h_su = self.point_su.h  # Enthalpy at the suction point
                P_ex = self.point_ex.p  # Pressure at the discharge point
                rho_su = self.point_su.D  # Density at the suction point
                

                "Problème: rho_su bcp trop petit par rapport à l'eau et on finit donc avec des V_dot bcp trop grand"
                "Transform for water"
                rho_water = 1000  # Density of water [kg/m^3]
                m_dot_water = m_dot * rho_water / rho_su  # Mass flow rate of water [kg/s]
                V_dot_water = (m_dot_water/1000)*3600 # Volume flow rate of water [m^3/h]
                # print(V_dot_water)


                DeltaP_ref = P_ex - P_su  # Pressure difference across the pump
                DeltaP_water = DeltaP_ref * rho_water / rho_su # Pressure difference across the pump for water [Pa]
                g = 9.81  # Acceleration due to gravity
                # print('DeltaP_ref', DeltaP_ref)
                Head_water = DeltaP_water / (1000 * g)  # Calculate the head based on pressure difference and density
                # print(DeltaP_ref, DeltaP_water, Head_water)
                # Map the frequency and pump consumption based on the calculated head and mass flow rate
                self.w = self.frequency(Head_water, V_dot_water)  # Frequency based on head and volume flow rate [Hz]
                
                "Problème vient du self.w = nan"

                self.W_dot_pp = self.pump_consumption(self.w, V_dot_water)  # Pump consumption based on frequency and volume flow rate [kW]
                
                "Problème vient du self.W_dot_pp = nan"

                h_ex = h_su + self.W_dot_pp / m_dot  # Calculate the enthalpy at the discharge point
                #print(h_ex)
                # Set the calculated values for the discharge point
                self.point_ex.set_m_dot(m_dot)  # Set the mass flow rate at the discharge point
                self.point_ex.set_h(h_ex)  # Set the enthalpy at the discharge point
                if h_ex == None:
                    # print(h_ex)
                    self.convergence = False
                    # print(self.convergence)
                else:
                    self.convergence = True
                # self.convergence = True

            except:
                print("Error in the pump model")
                self.convergence = False

        else:
            if self.calculable == False:
                print("Input of the component not completely known")

    
    def frequency_mapping(self, V_dot_value, Head_value):
        "To plot the operating maps"
        # This function maps the frequency (w) based on the volume flow rate (V_dot) and head (Head) values.
        # Note: The mapping is approximate because the graph don't exactly respect the similitudes.

        # This function maps the frequency (w) based on the mass flow rate (m_dot) and head (Head) values.
        def frequency(V_dot_grid, Head_grid):
            # For the frequency w=60Hz, we can determine the Head corresponding to the mass flow rate based on the graph
            w1 = 60
            Head_w1 = -23.8 * V_dot_grid + 94.76
            # Using the similitudes, we can now determine the real frequency
            w = np.sqrt(Head_grid/Head_w1) * w1
            return w

        # Define the arrays for volume flow rate (V_dot) and head (Head)
        V_dot_arrays = np.linspace(0.01, 3.8, 100)
        Head_arrays = np.linspace(0.01, 100, 100)

        # Create a grid of volume flow rate (V_dot) and head (Head) values
        V_dot_grid, Head_grid = np.meshgrid(V_dot_arrays, Head_arrays, indexing='ij')

        # Calculate the frequency (w) values using the frequency mapping function
        w_grid = frequency(V_dot_grid, Head_grid)

        # Create an interpolator for the frequency (w) values
        w_interpolation = RegularGridInterpolator((V_dot_arrays, Head_arrays), w_grid,
                                                 bounds_error=False, fill_value=None)
        

    
    def W_dot_pp_mapping(self, w_value, V_dot_value):
        "To plot the operating maps"
        # This function maps the pump consumption (W_dot_pp) based on the volume flow rate (V_dot) and frequency (w) values.
        # Note: The mapping is approximate because the graph don't exactly respect the similitudes.

        def Pump_Consumption(V_dot_grid, w_grid):
            # For the frequency w=60Hz, we can determine the Power corresponding to the volume flow rate
            w1 = 60
            W_dot_pp_w1 = -0.2 * V_dot_grid + 1.8
            # Using the similitudes, we can now determine the real Power
            W_dot_pp = W_dot_pp_w1 * (w_grid/w1)**3
            return W_dot_pp
        
        # Define the arrays for volume flow rate (V_dot) and frequency (w)
        V_dot_arrays = np.linspace(0.01, 3.8, 100)
        w_arrays = np.linspace(0.01, 100, 100)

        # Create a grid of volume flow rate (V_dot) and frequency (w) values
        V_dot_grid, w_grid = np.meshgrid(V_dot_arrays, w_arrays, indexing='ij')
        # Calculate the Power (W_dot_pp) values using the frequency mapping function
        W_dot_pp_grid = Pump_Consumption(V_dot_grid, w_grid)

        # Create an interpolator for the Power (W_dot_pp) values
        W_dot_pp_interpolation = RegularGridInterpolator((V_dot_arrays, w_arrays), W_dot_pp_grid,
                                         bounds_error=False, fill_value=None)

        # Get the Power (W_dot_pp) value for a specific volume flow rate (m_dot) and frequency (w)
        point = np.array([V_dot_value, w_value])
        W_dot_pp_value = float(W_dot_pp_interpolation(point)) #otherwise it returns an array
        return W_dot_pp_value





