#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:53:29 2023

@author: olivierthome
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.interpolate import interp1d

class Plate_HX_Geom_SWEP(object):
    
    def __init__(self, **kwargs):

        "Heat transfer area"
        self.A_h = None # Hot side
        self.A_c = None # Cold side
        
        "HX Volume"
        self.H_V_tot = None # Hot side volume
        self.C_V_tot = None # Cold side volume
        
        "HX Dimensions"
        self.l = None # length
        self.w = None # width
        self.h = None # height
        self.l_v = None # length between ports
        self.casing_t = None # casing thickness
        
        "Plate related parameters"
        self.n_plates = None # number of plates
        self.t_plates = None # plate thickness
        self.plate_cond = None # Plate conduction
        self.fooling = None # Fooling Factor
        self.plate_pitch_co = None # corrugated pitch 
        self.chevron_angle = None # chevron angle

    def set_parameters_SWEP(self, name, **kwargs):
        
        if name == "B35TM0x10/1P":
            
            # Heat transfer area  
            self.A_h = 0.752 # m^2
            self.A_c = 0.752 # m^2
            
            # HTX dimensions
            self.l = 0.393 # [m] : length
            self.w = 0.243 # [m] : width
            self.h = 0.0446 # [m] : height
            self.l_v = 0.324 # [m] : length between ports
            self.casing_t = 0.005 # [m] : casing thickness # !!! arbitrary
            
            # Number and thickness of plates
            self.n_plates = 10
            self.t_plates = 0.0008 # [m] # !!! arbitrary
            
            # Fooling factor
            self.fooling = 98.51/100
        
            # Number of canals
            self.C_n_canals = 5
            self.H_n_canals = 4
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction
            self.plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
            self.chevron_angle = 20*np.pi/180 # !!! arbitrary
        
            # Total volume of each part
            self.H_V_tot = 0.9*1e-3 # [m^3]
            self.C_V_tot = 0.9*1e-3 # [m^3]
            
            # Canal thickness
            self.C_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.C_n_canals)
            self.H_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.H_n_canals)
        
            # Canal Surface
            self.C_CS = self.C_canal_t*(self.w-2*self.casing_t)
            self.H_CS = self.H_canal_t*(self.w-2*self.casing_t)
    
            # Dh : hydraulic diameter
            self.C_Dh = (4*self.C_canal_t*self.w)/(2*self.C_canal_t+2*self.w)
            self.H_Dh = (4*self.H_canal_t*self.w)/(2*self.H_canal_t+2*self.w)
        
        if name == "B20Hx24/1P":
            
            # Heat transfer area  
            self.A_h = 0.616 # m^2
            self.A_c = 0.616 # m^2
            
            # HTX dimensions
            self.l = 0.324 # [m] : length
            self.w = 0.094 # [m] : width
            self.h = 0.041 # [m] : height
            self.l_v = 0.2682 # [m] : length between ports
            self.casing_t = 0.005 # [m] : casing thickness # !!! arbitrary
            
            # Number and thickness of plates
            self.n_plates = 24
            self.t_plates = 0.0009 # [m] # !!! arbitrary
            
            # Fooling factor
            self.fooling = 0.705/100
        
            # Number of canals
            self.C_n_canals = 12
            self.H_n_canals = 11
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction
            self.plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
            self.chevron_angle = 20*np.pi/180 # !!! arbitrary
        
            # Total volume of each part
            self.H_V_tot = 0.34*1e-3 # [m^3]
            self.C_V_tot = 0.34*1e-3 # [m^3]        
            
            # Canal thickness
            self.C_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.C_n_canals)
            self.H_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.H_n_canals)
        
            # Canal Surface
            self.C_CS = self.C_canal_t*(self.w-2*self.casing_t)
            self.H_CS = self.H_canal_t*(self.w-2*self.casing_t)
    
            # Dh : hydraulic diameter
            self.C_Dh = (4*self.C_canal_t*self.w)/(2*self.C_canal_t+2*self.w)
            self.H_Dh = (4*self.H_canal_t*self.w)/(2*self.H_canal_t+2*self.w)  

        if name == "P200THx140/1P_Evaporator":
            
            # Heat transfer area  
            self.A_h = 17.8 # m^2 hot side OK
            self.A_c = 17.8 # m^2 cold side OK
            
            # HTX dimensions
            self.l = 0.525 # [m] :  OK
            self.w = 0.243 # [m] :  OK
            self.h = 0.3426 # [m] : height F??
            self.l_v = 0.450 # [m] : length between ports H??
            self.casing_t = 0.02 # [m] : casing thickness # !!! arbitrary
            
            # Number and thickness of plates
            self.n_plates = 140 # OK
            self.t_plates = 0.0019 # [m] # !!! arbitrary épaisseur plaque (certainement plus grand)
            
            # Fooling factor
            self.fooling = 0.111/100 # [m^2, °C/kW] OK
        
            # Number of canals per pass
            self.C_n_canals = 69 # Cold side OK
            self.H_n_canals = 70 # Hot side OK
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction acier inoxidable OK
            self.plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
            self.chevron_angle = 20*np.pi/180 # !!! arbitrary 20, 30 ou 45°
        
            # Total volume of each part
            self.H_V_tot = 16.63*1e-3 # [m^3] OK hot volume
            self.C_V_tot = 16.87*1e-3 # [m^3] OK cold volume
            
            # Canal thickness
            self.C_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.C_n_canals)
            self.H_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.H_n_canals)
        
            # Canal Surface
            self.C_CS = self.C_canal_t*(self.w-2*self.casing_t)
            self.H_CS = self.H_canal_t*(self.w-2*self.casing_t)
    
            # Dh : hydraulic diameter
            self.C_Dh = (4*self.C_canal_t*self.w)/(2*self.C_canal_t+2*self.w)
            self.H_Dh = (4*self.H_canal_t*self.w)/(2*self.H_canal_t+2*self.w)    
        
        if name == "P200THx140/1P_Condenser":
            
            # Heat transfer area  
            self.A_h = 17.8 # m^2 hot side OK
            self.A_c = 17.8 # m^2 cold side OK
            
            # HTX dimensions
            self.l = 0.525 # [m] :  OK
            self.w = 0.243 # [m] :  OK
            self.h = 0.346 # [m] : height F??
            self.l_v = 0.450 # [m] : length between ports H??
            self.casing_t = 0.06 # [m] : casing thickness # !!! arbitrary
            
            # Number and thickness of plates
            self.n_plates = 140 # OK
            self.t_plates = 0.0009 # [m] # !!! arbitrary épaisseur plaque (certainement plus grand)
            
            # Fooling factor
            self.fooling = 0.258/100 # [m^2, °C/kW] OK
        
            # Number of canals per pass
            self.C_n_canals = 70 # Cold side OK
            self.H_n_canals = 69 # Hot side OK
            
            # Plate values 
            self.plate_cond = 45 # [W/(m*K)] : plate conduction acier inoxidable OK
            self.plate_pitch_co = 0.0001 # 0.00745870973 # corrugated pitch # !!! arbitrary
            self.chevron_angle = 20*np.pi/180 # !!! arbitrary 20, 30 ou 45°
        
            # Total volume of each part
            self.H_V_tot = 16.63*1e-3 # [m^3] OK hot volume
            self.C_V_tot = 16.63*1e-3 # [m^3] OK cold volume
            
            # Canal thickness
            self.C_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.C_n_canals)
            self.H_canal_t = ((self.h-2*self.casing_t) - self.n_plates*self.t_plates)/(2*self.H_n_canals)
        
            # Canal Surface
            self.C_CS = self.C_canal_t*(self.w-2*self.casing_t)
            self.H_CS = self.H_canal_t*(self.w-2*self.casing_t)
    
            # Dh : hydraulic diameter
            self.C_Dh = (4*self.C_canal_t*self.w)/(2*self.C_canal_t+2*self.w)
            self.H_Dh = (4*self.H_canal_t*self.w)/(2*self.H_canal_t+2*self.w)    

        else: # manual setting
        
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger geometry.")
