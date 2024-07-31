# from __future__ import division, print_function

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Off_design_Model_CB')

from Port.Mass_connector import Mass_connector
from Component.PHE.Moving_Boundary.HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger
from CoolProp.CoolProp import PropsSI 
import time


"--------- 1) Evaporator ------------------------------------------------------------------------------------------"

start_time = time.time()
"Refrigerant Su"
Evap_ref_su = Mass_connector()

# Set the fluid
Evap_ref_su.set_fluid('R1233zd(E)')

# For the test on the datasheet
Evap_ref_su.set_m_dot(0.33447089409117)  # Example mass flow rate [kg/s] (from SWEP)
Evap_ref_su.set_h(330790.6227404063) # Example temperature [K] (from SWEP)
P_sat = PropsSI('P','T', 60+273.15,'Q', 0.5, Evap_ref_su.fluid) # Example Pressure [Pa]
Evap_ref_su.set_p(550000)  # Example Pressure [Pa]

"Water Su"
Evap_w_su = Mass_connector()
Evap_w_su.set_fluid("water")

# Set other properties
Evap_w_su.set_m_dot(10)  # Example mass flow rate [kg/s]
Evap_w_su.set_T(90 + 273.15) # Example temperature [K]
Evap_w_su.set_p(2e5)  # Example Pressure [Pa]

#Pas de problème de "pair of inputs" until here
"Pressure Drop"

DP_H_ON = False
DP_C_ON = False

"Heat Exchanger parameters"

n_disc = 15 # number of discretization

calc = 1 # flag to compute the HX
plot = 1 # flag to plot the HX fluids temperature profiles
print_flag = 1 # flag to print the HX results

"Geometry"

HX_Evap_geom = Plate_HX_Geom_SWEP()
HX_Evap_geom.set_parameters_SWEP("P200THx140/1P_Evaporator")

"Flow and htc correlation types parameters"

htc_type = "Correlation" # not user-defined but the relation are implemented in the code htc = heat transfer coefficient
flow_type = "Counter_flow"

"Heat Exchanger initiation and computation"

print("\n")
HX_Evap =  Plate_HeatExchanger()
HX_Evap.inputs(Evap_w_su, Evap_ref_su, wf_T='cold')
print(HX_Evap.calculable)
HX_Evap.set_parameters(**{
    'geom': HX_Evap_geom,
    'n_disc': n_disc
})
HX_Evap.solve()
print(HX_Evap.C_ex.T-273.15)
print(HX_Evap.H_ex.T-273.15)
print(HX_Evap.Q)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.2f} [s]")
# # "--------- 3) Condenser ------------------------------------------------------------------------------------------"

# "Cyclopentane Su"
# Cond_ref_su = Mass_connector()

# # Set the fluid
# Cond_ref_su.set_fluid('R1233zd(E)')

# # Set other properties
# Cond_ref_su.set_m_dot(0.4582)  # Example mass flow rate [kg/s]
# Cond_ref_su.set_T(40 + 273.15) # Example temperature [K]
# P_su = PropsSI('P','T', 25+273.15,'Q', 0.5, Cond_ref_su.fluid) # Example Pressure [Pa]
# Cond_ref_su.set_p(P_su)  # Example Pressure [Pa]

# "Water Su"
# Cond_water_su = Mass_connector()

# # Set the fluid
# Cond_water_su.set_fluid("Water")

# # Set other properties
# Cond_water_su.set_m_dot(2.839)  # Example mass flow rate [kg/s]
# Cond_water_su.set_T(15 + 273.15) # Example temperature [K]
# Cond_water_su.set_p(2e5)  # Example Pressure [Pa]

# "Pressure Drop"

# DP_H_ON = True
# DP_C_ON = True

# "Heat Exchanger parameters"

# n_disc = 15 # number of discretization

# calc = 1 # flag to compute the HX
# plot = 0 # flag to plot the HX fluids temperature profiles
# print_flag = 1 # flag to print the HX results

# "Geometry"

# HX_Cond_geom = Plate_HX_Geom_SWEP()
# HX_Cond_geom.set_parameters_SWEP("P200THx140/1P_Condenser")

# "Flow and htc correlation types parameters"

# htc_type = "Correlation"
# flow_type = "Counter_flow"

# "Heat Exchanger initiation and computation"

# print("\n")
# HX_Cond =  Plate_HeatExchanger()
# HX_Cond.inputs(Cond_ref_su, Cond_water_su, wf_T='hot')
# print(HX_Cond.calculable)
# HX_Cond.set_parameters(**{
#     'geom': HX_Cond_geom,
#     'n_disc': n_disc
# })
# HX_Cond.solve()
# print(HX_Cond.C_ex.T-273.15)
# print(HX_Cond.H_ex.fluid)
# print(HX_Cond.H_ex.T-273.15)

# print(HX_Cond.Q)
# M_HX_Cond = sum(HX_Cond.Mvec_h)
    