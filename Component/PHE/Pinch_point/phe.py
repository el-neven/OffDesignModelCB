#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:54:19 2023

@author: olivierthome
"""

# Modèle d'échangeur SWEP

#    --------- 
#   |F1     F2|
#   |         |  
#   |         | 
#   |         | 
#   |         | 
#   |         | 
#   |F3     F4|
#    ---------

# Side 1 (S1): F1/F3
# Side 2 (S2): F2/F4

# Classical_config = True : F1 -> F3
# Classical_config = False: F3 -> F1

# Cross_flow = True:
#   if Classical_config:
#       F4 -> F2
#   else:
#       F2 -> F4

# Cross_flow = False:
#   if Classical_config:
#       F2 -> F4
#   else:
#       F4 -> F2


import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import math

class Pinch_method_PHE:
    def __init__(self, Classical_config, Cross_flow, Accuracy_order, Plot):
        """
        Parameters
        ----------
        Classical_config: Bool
        Cross_flow: Bool
        Accuracy_order = n: Power PHE defined with minimal error ± 10^n [W]
        Plot:   True pour sortir des graphes (à ne pas utiliser dans un process itératif)
                False pour ne pas sortir de graphes

        """
        
        # Pas de pertes de charge considérées
        # Pas de pertes à l'ambiance

        self.calculable = False
        self.parametrized = False
        self.defined = False
        
        self.classical_config = Classical_config
        self.cross_flow = Cross_flow
        self.order = 10**round(Accuracy_order)
        self.plot = Plot
        
        self.su_1 = None
        self.su_2 = None
        self.ex_1 = None
        self.ex_2 = None
        
        self.pinch = None
        
    def check_calculable(self):
        if self.su_1 != None and self.su_2 != None:
            if self.su_1.completely_known and self.su_2.completely_known:
                self.calculable = True
        
    def check_parametrized(self):
        if self.pinch != None:
            self.parametrized = True
            

    def set_su_1(self, value):
        self.su_1 = value
        self.check_calculable()
        
    def set_su_2(self, value):
        self.su_2 = value
        self.check_calculable()
        
    def set_ex_1(self, value):
        self.ex_1 = value
        self.check_calculable()
        
    def set_ex_2(self, value):
        self.ex_2 = value
        self.check_calculable()


    def set_pinch(self,value):
        self.pinch = value
        self.check_parametrized()
        
        
    def solve(self):
        if self.calculable and self.parametrized:
            
            if abs(self.su_1.T-self.su_2.T) <= self.pinch:
                self.ex_1.set_fluid(self.su_1.fluid)
                self.ex_2.set_fluid(self.su_2.fluid)
                self.ex_1.set_m_dot(self.su_1.m_dot)
                self.ex_2.set_m_dot(self.su_2.m_dot)
                self.ex_1.set_T(self.su_1.T)
                self.ex_2.set_T(self.su_2.T)
                self.ex_1.set_h(self.su_1.h)
                self.ex_2.set_h(self.su_2.h)
                print("dT input fluids < pinch")
                
            else:
                if self.cross_flow:
                    if self.su_1.T < self.su_2.T:
                        cold_fluid = self.su_1
                        hot_fluid = self.su_2
                    else:
                        cold_fluid = self.su_2
                        hot_fluid = self.su_1
                        
                    try:
                        dh1 = PropsSI('H', 'T', hot_fluid.T-self.pinch, 'P', cold_fluid.p, cold_fluid.fluid) - cold_fluid.h
                    except: # si point sous cloche de saturation
                        dh1 = PropsSI('H', 'T', hot_fluid.T-self.pinch, 'Q', 1, cold_fluid.fluid) - cold_fluid.h
                        # vérifie que la cause du expect était bien ça
                        if abs(PropsSI('P', 'T', hot_fluid.T-self.pinch, 'Q', 1, cold_fluid.fluid) - cold_fluid.p) > 1:
                            print('Error: input PHE not valid')
                    
                    try:
                        dh2 = hot_fluid.h - PropsSI('H', 'T', cold_fluid.T+self.pinch, 'P', hot_fluid.p, hot_fluid.fluid)
                    except: # si point sous cloche de saturation
                        dh2 = hot_fluid.h - PropsSI('H', 'T', cold_fluid.T+self.pinch, 'Q', 0, hot_fluid.fluid)
                        # vérifie que la cause du expect était bien ça
                        if abs(PropsSI('P', 'T', cold_fluid.T+self.pinch, 'Q', 0, hot_fluid.fluid) - hot_fluid.p) > 1:
                            print('Error: input PHE not valid')
                   
                    Q1 = dh1*cold_fluid.m_dot
                    Q2 = dh2*hot_fluid.m_dot
                    Q_max = round(min(Q1,Q2)/self.order)*self.order
                    
                    Q_list = np.zeros(int(Q_max/self.order)+1, dtype=int)
                    T_cold_fluid_list = np.zeros(int(Q_max/self.order)+1, dtype=float)
                    T_hot_fluid_list = np.zeros(int(Q_max/self.order)+1, dtype=float)
                    
                    for i in range (0,len(Q_list)):
                        Q_list[i] = i*self.order # [0,..,Q_max]
                        T_cold_fluid_list[i] = PropsSI('T', 'H', cold_fluid.h+(Q_list[i]/cold_fluid.m_dot), 'P', cold_fluid.p, cold_fluid.fluid) # [T_su,..,T_ex_max]
                        T_hot_fluid_list[i] = PropsSI('T', 'H', hot_fluid.h-(Q_list[i]/hot_fluid.m_dot), 'P', hot_fluid.p, hot_fluid.fluid) # [T_su,..,T_ex_max]
                        
                    solved = False
                    Q_i = -1 # Q_list[Q_i=-1] = Q_max
                    T_cold_i = 0
                    
                    # En partant de la température froide d'entrée, on essaye les points suivants (de T_cold) pour voir si le pinch est respecté avec Q_max
                    # Si un point non respecté est trouvé, on diminue la puissance de l'échangeur jusqu'à ce que le pinch soit respecté sur ce point
                    # Une fois le pinch respecté, on recommence l'opération à partit du point (T_cold) suivant
                    while not solved:
                        if T_hot_fluid_list[Q_i-T_cold_i]-T_cold_fluid_list[T_cold_i]<self.pinch: # point de pinch non respecté trouvé
                            Q_i-=1 # diminue puissance échangeur
                        else:
                            T_cold_i+=1 # prochain point à tester
                        
                        if Q_i-T_cold_i == -len(Q_list)-1: # prochain T_hot à tester hors de la liste, càd que tous les points ont été testés
                            solved = True
                        else:
                            pass
                        
                    self.ex_1.set_fluid(self.su_1.fluid)
                    self.ex_2.set_fluid(self.su_2.fluid)
                    self.ex_1.set_m_dot(self.su_1.m_dot)
                    self.ex_2.set_m_dot(self.su_2.m_dot)
                    self.ex_1.set_p(self.su_1.p)
                    self.ex_2.set_p(self.su_2.p)
                    
                    if self.su_1.T < self.su_2.T:
                        self.ex_1.set_h(self.su_1.h+(Q_list[Q_i]/self.su_1.m_dot))
                        self.ex_2.set_h(self.su_2.h-(Q_list[Q_i]/self.su_2.m_dot))
                    else:
                        self.ex_1.set_h(self.su_1.h-(Q_list[Q_i]/self.su_1.m_dot))
                        self.ex_2.set_h(self.su_2.h+(Q_list[Q_i]/self.su_2.m_dot))
                     
                    # Plot results
                    if self.plot:
                        if math.log10(Q_list[-1]) >= 3.7:
                            Q_list_kW = np.zeros(len(Q_list), dtype=float)
                            for i in range (0,len(Q_list)):
                                Q_list_kW[i] = Q_list[i]/1000
                            Q_list_plot = Q_list_kW
                            x_label = "$\.H$ [kW]"
                        else:
                            Q_list_plot = Q_list
                            x_label = "$\.H$ [W]"
                                
                        plt.plot(Q_list_plot, T_cold_fluid_list, linestyle='dotted', color='#1f77b4')
                        plt.plot(Q_list_plot, np.flip(T_hot_fluid_list), linestyle='dotted', color='#d62728')
                        if Q_i == -1:
                            plt.plot(Q_list_plot, T_cold_fluid_list, color='#1f77b4')
                            plt.plot(Q_list_plot, np.flip(T_hot_fluid_list), color='#d62728')
                        else:
                            plt.plot(Q_list_plot[0:Q_i+1], T_cold_fluid_list[0:Q_i+1], color='#1f77b4')
                            plt.plot(Q_list_plot[0:Q_i+1], np.flip(T_hot_fluid_list[0:Q_i+1]), color='#d62728')
                        plt.legend([cold_fluid.fluid+": pinch respected on su/ex", hot_fluid.fluid+": pinch respected on su/ex", cold_fluid.fluid+": pinch respected for all points", hot_fluid.fluid+": pinch respected for all points"], loc = "lower right", prop={'size': 6})
                        center = int((len(Q_list)+Q_i)/2)
                        print(center)
                        print(Q_list_plot[center])
                        plt.arrow(Q_list_plot[center], T_cold_fluid_list[center], Q_list_plot[center+1]-Q_list_plot[center], T_cold_fluid_list[center+1]-T_cold_fluid_list[center], width=0.5, color='#1f77b4')
                        plt.arrow(Q_list_plot[center], np.flip(T_hot_fluid_list[0:Q_i+1])[center], Q_list_plot[center-1]-Q_list_plot[center], np.flip(T_hot_fluid_list[0:Q_i+1])[center-1]-np.flip(T_hot_fluid_list[0:Q_i+1])[center], width=0.5, color='#d62728')
                        plt.xlim(left=0)
                        plt.ylabel("$T$ [K]")
                        plt.xlabel(x_label)
                        plt.grid()
                        plt.savefig("essai_phe_pinch.png", dpi=300)
                    else:
                        pass
                    
                    
                else: # parallele flow
                    pass
                
                self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")



















        