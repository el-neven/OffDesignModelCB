import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/POC/Models/Off_design_Model_EN')

from Port.Mass_connector import Mass_connector
import time
from Component.PHE.Simple_MB.HX_MovingBoundary import HeatExchanger

#------------------------------------------------------------------------------
"Condenser test"

"Define ports"
su1 = Mass_connector()
su2 = Mass_connector()
ex1 = Mass_connector()
ex2 = Mass_connector()

# Refrigerant supply
su1.set_fluid('R1233zd(E)')
su1.set_T(40+273.15)
su1.set_m_dot(0.4582)
su1.set_x(1)
# su1.set_p(200000)

# Secondary fluid supply
su2.set_fluid('Water')
su2.set_T(15+273.15)
su2.set_p(2e5)
su2.set_m_dot(2.839)

# Secondary fluid exhaust
ex2.set_fluid('Water')
ex2.set_T(23+273.15+10) #10 glide
ex2.set_p(2e5)

"Define class"
HX = HeatExchanger(su1, su2, ex1, ex2)
# HX.inputs(su, ex, N_exp, T_amb)

HX.set_parameters(**{
    'HX_type': 'condenser',
    'HX_D': 0.06, #Diamètre de port d'entré
    'HX_A': 17.8, #Surface d'échange
    'min_pinch': 2,
    'dT_sub_or_sup': 5
})

HX.solve()
HX.plot_heat_transfer_diagram()
print(HX.Q)
print(HX.res)
print(HX.ex1.p, 'pressure')

#---------------------------------------------------------------------------------
#------------------------------------------------------------------------------
"Evaporator test"

"Define ports"
su1 = Mass_connector()
su2 = Mass_connector()
ex1 = Mass_connector()
ex2 = Mass_connector()

# Refrigerant supply
su1.set_fluid('R1233zd(E)')
su1.set_T(60+273.15)
su1.set_m_dot(0.4621)
su1.set_x(0)
# su1.set_p(400000)

# Secondary fluid supply
su2.set_fluid('Water')
su2.set_T(70+273.15)
su2.set_p(2e5)
su2.set_m_dot(2.425)

# Secondary fluid exhaust
ex2.set_fluid('Water')
ex2.set_T(62+273.15)
ex2.set_p(2e5)

"Define class"
HX = HeatExchanger(su1, su2, ex1, ex2)
# HX.inputs(su, ex, N_exp, T_amb)

HX.set_parameters(**{
    'HX_type': 'evaporator',
    'HX_D': 0.06, #Diamètre de port d'entré
    'HX_A': 17.8, #Surface d'échange
    'min_pinch': 2,
    'dT_sub_or_sup': 5
})

HX.solve()
HX.plot_heat_transfer_diagram()
print(HX.Q)
print(HX.res)
print(HX.ex1.p, 'pressure')
print(HX.ex1.T, 'temperature')