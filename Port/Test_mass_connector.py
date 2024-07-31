# from Mass_connector import Mass_connector

# point = Mass_connector()

# point.set_fluid('water')
# point.set_p(101325)
# point.set_T(300)

# print(point.h)
# print(point.T)

# point.set_T(350)

# print(point.h)
# print(point.T)

# point.set_m_dot(10)
# print(point.m_dot)
# point.set_m_dot(20)
# print(point.m_dot)

from Heat_connector import Heat_connector

#Idée: par example échange de chaleur entre l'ambiance et l'expander.
# Faire un component Ambiance et connecter tout les composants à Ambiance.

Amb_Exp = Heat_connector()

Amb_Exp.set_T1(300)
Amb_Exp.set_T2(290)
Amb_Exp.set_AU(1000)

print(Amb_Exp.Q_dot)