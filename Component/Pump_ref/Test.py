from Model import PumpModel
import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

from Port.Mass_connector import Mass_connector

"Define ports"
point_su = Mass_connector()
point_ex = Mass_connector()

point_su.set_fluid('R1233ZDE')
point_su.set_m_dot(0.13)
point_su.set_T(291.44817218160415)
point_su.set_p(200000)

point_ex.set_fluid('R1233ZDE')
point_ex.set_p(600000)


m_dot = 1.2
Head = 50
Pump = PumpModel()
Pump.set_inputs(point_su, point_ex)
Pump.solve()
print(Pump.w, Pump.W_dot_pp)
