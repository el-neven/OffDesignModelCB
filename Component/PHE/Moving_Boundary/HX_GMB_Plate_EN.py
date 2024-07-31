"""
Supplemental code for paper:
I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
"""

"""
Modification w/r to previous version:
    - Putting some order in the Objective Function "for" loops. Sparing some
    lines of code.
    - x_di_c correct calculation.
"""

# from __future__ import division, print_function

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot

import matplotlib.pyplot as plt
import numpy as np
from math import log, pi
import scipy.optimize

from Component.PHE.Moving_Boundary.Modules.F_lmtd2 import F_lmtd2
from Component.PHE.Moving_Boundary.Modules.PropsFluid_210907 import PropsFluid, conducticity_R1233zd
from Component.PHE.Moving_Boundary.Modules.Gnielinski_Pipe_HTC_210908 import Gnielinski_Pipe_HTC
from Component.PHE.Moving_Boundary.Modules.Han_Cond_BPHEX_HTC_210909 import Han_Cond_BPHEX_HTC
from Component.PHE.Moving_Boundary.Modules.Han_Boiling_BPHEX_HTC_210910 import Han_Boiling_BPHEX_HTC
from Component.PHE.Moving_Boundary.Modules.Kim_DryOutIncipience_210913 import Kim_DryOutIncipience
from Component.PHE.Moving_Boundary.Modules.find_2P_boundaries import find_2P_boundaries
from Component.PHE.Moving_Boundary.Modules.Plate_HTC import water_Plate_HTC

from Port.Mass_connector import Mass_connector
from Component.PHE.Moving_Boundary.Modules.Geometry_Plate_HX_SWEP import Plate_HX_Geom_SWEP

#%%
# Set to True to enable some debugging output to screen
debug = False

#%%

class Plate_HeatExchanger(): #useless to put (object) in the brackets
    
    class geom():
            pass      
    class H():
        pass
    class C():
        pass  
    
    def __init__(self, flow_type = 'Counter_flow', htc_type= 'Correlation', H_DP_ON = False, C_DP_ON = False, PropFluid = None):
        
        "Flow related parameters"
        self.flow_type = flow_type
        self.htc_type = htc_type
        
        "Pressure drop parameters"
        self.H_DP_ON = H_DP_ON
        self.C_DP_ON = C_DP_ON

        "Fluid properties"
        self.PropFluid = PropFluid
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.defined = False
        
        "Geometry sub-object"
        self.geom = None # parameters 
        
        "Supply"
        self.point_su = [None, None]
        
        "Exhaust"
        self.point_ex = [None, None]        
        
        "Is the working fluid the hot or cold fluid ?"
        self.wf_T = None
        
        "Hot Side"
        self.H_su = None
        self.H_ex = None
        
        "Cold side"
        self.C_su = None
        self.C_ex = None
        
        "Number of discretizations"
        self.n_disc = None
        
        "Type of flow"
        self.typeHEX = flow_type
                
#%%    
    def check_calculable(self):
        if self.point_su[0] != None and self.point_su[1] != None:
            if (self.point_su[0].completely_known) and (self.point_su[1].completely_known):
                self.calculable = True

    def inputs(self, H_su, C_su, wf_T):
        self.point_su = [H_su, C_su]
        self.H_su = H_su
        self.C_su = C_su
        self.wf_T = wf_T
        #self.point_ex = [H_ex, C_ex]
        self.check_calculable()

    def check_parametrized(self):
        if self.n_disc != None and self.geom != None:
            self.parametrized = True
           
#%%
    def set_parameters(self, **kwargs):
            """
            Set parameters of the heat exchanger.
    
            Parameters
            ----------
            **kwargs : dict
                Key-value pairs representing parameters and their values.
                
                Example of call : heat_exchanger.set_parameters(D_o=0.05, t=0.002, H_core=2.0)
            """
            for key, value in kwargs.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger.")
                    
                    "Type of HTC"
        
            if self.htc_type == "User-Defined":
                # !!! User-Defined Heat Transfer Coefficients (hot):
                self.H.h_liq = 100
                self.H.h_vap = 100
                self.H.h_twophase = 1000
                self.H.h_vapwet = 100
                self.H.h_tpdryout = 1000
                self.H.h_transcrit = 200
            
                # !!! User-Defined Heat Transfer Coefficients (cold):
                self.C.h_liq = 100
                self.C.h_vap = 100
                self.C.h_twophase = 1000
                self.C.h_vapwet = 100
                self.C.h_tpdryout = 10000
                self.C.h_transcrit = 200
            else:
                # Type 
                self.H.HeatExchange_Correlation = "Correlation"
                # params.H.HeatExchange_Correlation = "User-Defined"
                self.H.Correlation_1phase = "Gnielinski"
                self.H.Correlation_2phase = "Han_cond_BPHEX"
                
                self.C.HeatExchange_Correlation = "Correlation"
                # params.C.HeatExchange_Correlation = "User-Defined"
                self.C.Correlation_1phase = "Gnielinski"
                self.C.Correlation_2phase = "Han_Boiling_BPHEX_HTC"
    
            if self.H_DP_ON == True: # !!! arbitrarily set
                self.H.f_dp = {"K": 14.14, "B": 1.892}
            else:
                self.H.f_dp = {"K": 0, "B": 0}        
                
            if self.C_DP_ON == True:
                self.C.f_dp = {"K": 14.14, "B": 1.892}
            else:
                self.C.f_dp = {"K": 0, "B": 0}
                
            self.check_parametrized()


#%%
    def external_pinching(self): # To find boundaries for the iteration on Q_dot
        "Determine the maximum heat transfer rate based on the external pinching analysis"
        
        "1) Set of arbitrary bound values" # !!! Find out why      
        T_hmin = 218 # above CO2 freezing point (194.7 K)
            
        T_cmax = 273.15+250
        
        "2) Hot fluid side pinch"
        
        # Computation of outlet lowest possible enthalpy of the hot fluid using either the entrance cold fluid inlet temperature, or the arbitrary minimal
        #self.h_ho = CP.PropsSI('H','T', max(self.T_ci, T_hmin),'P',self.p_hi,self.H_su.fluid) # Equation 5 (Bell et al. 2015)
        #self.PropFluid.update(CP.PT_INPUTS, self.p_hi, max(self.T_ci, T_hmin))
        self.h_ho = self.PropFluid.update(CP.PT_INPUTS, self.p_hi, max(self.T_ci, T_hmin)).h()
        # Q_max computation
        Qmaxh = self.mdot_h*(self.h_hi-self.h_ho) # Equation 4 (Bell et al. 2015)
        
        if debug:
            print("Inlet Hot Enthalpy:", self.h_hi, "\n")
            print("Outlet Hot Enthalpy:", self.h_ho,"\n")
            print("mdot_h = ", self.mdot_h, "\n")
            print("Qmaxh =", Qmaxh, "\n")

        "3) Cold fluid side pinch"
                
        # Computation of highest possible outlet enthalpy of the hot fluid using either the entrance hot fluid inlet temperature, or the arbitrary maximal
        self.h_co = CP.PropsSI('H','T',min(self.T_hi, T_cmax),'P',self.p_ci,self.C_su.fluid) # Equation  (Bell et al. 2015)
        # Q_max computation
        Qmaxc = self.mdot_c*(self.h_co-self.h_ci) # Equation 6 (Bell et al. 2015)

        if debug:
            print("Inlet Cold Enthalpy:", self.h_ci, "\n")
            print("Outlet Cold Enthalpy:", self.h_co,"\n")
            print("mdot_c = ", self.mdot_c, "\n")
            print("Qmaxh =", Qmaxc, "\n")
            
        "4) Final pinch and cell boundaries computation"

        Qmax = min(Qmaxh, Qmaxc)
        
        if debug:
            print('Qmax (external pinching) is', Qmax)

        self.calculate_cell_boundaries(Qmax) # call calculate_cell_boundaries procedure

        return Qmax
    
#%%
    def calculate_cell_boundaries(self, Q):
        """ Calculate the cell boundaries for each fluid """
        
        "1) Re-calculate the outlet enthalpies of each stream"
                
        self.h_co = self.h_ci + Q/self.mdot_c
        self.h_ho = self.h_hi - Q/self.mdot_h

        "2) Calculate the dew and bubble pressures and enthalpies by accounting for pressure drops"
        
        "2.1) For the cold fluid"
        if not self.Transcritical_c: # if not in transcritical
            if (not (self.C.f_dp == {"K": 0, "B": 0}) or self.DP_c < 1e-2):
                # If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
                self.p_cdew = self.p_ci
                self.p_cbubble = self.p_co
                self.T_cbubble = self.T_cbubble_ideal
                self.T_cdew = self.T_cdew_ideal
                self.h_cbubble = self.h_cbubble_ideal
                self.h_cdew = self.h_cdew_ideal
            else:
                h_cbubble, h_cdew, p_cbubble, p_cdew, _, _, = find_2P_boundaries(self.C_su.fluid, self.h_ci, self.h_co, self.p_ci, self.p_co)
                self.p_cdew = p_cdew
                self.p_cbubble = p_cbubble
                self.h_cdew = h_cdew
                self.h_cbubble = h_cbubble
                self.T_cdew = CP.PropsSI("T", "P", self.p_cdew, "Q", 1, self.C_su.fluid)
                self.T_cbubble = CP.PropsSI("T", "P", self.p_cbubble, "Q", 0, self.C_su.fluid)
                
        "2.2) For the hot fluid"
        
        if not self.Transcritical_h:
            if (not (self.H.f_dp == {"K": 0, "B": 0}) or self.DP_h < 1e-2) and not self.h_incomp_flag:
                #If no pressure drop is calculated or if the pressure drop is neglectable, assign saturation conditions to the ideal case:
                self.p_hdew = self.p_hi
                self.p_hbubble = self.p_ho    
                self.T_hbubble = self.T_hbubble_ideal
                self.T_hdew = self.T_hdew_ideal
                self.h_hbubble = self.h_hbubble_ideal
                self.h_hdew = self.h_hdew_ideal
            elif not self.h_incomp_flag:
                h_hbubble, h_hdew, p_hdew, p_hbubble, _, _, = find_2P_boundaries(self.H_su.fluid, self.h_hi, self.h_ho, self.p_hi, self.p_ho)
                self.p_hdew = p_hdew
                self.p_hbubble = p_hbubble
                self.h_hdew = h_hdew
                self.h_hbubble = h_hbubble
                self.T_hdew = CP.PropsSI("T", "P", self.p_hdew, "Q", 1, self.H_su.fluid)
                self.T_hbubble = CP.PropsSI("T", "P", self.p_hbubble, "Q", 0, self.H_su.fluid)

        "3) Discretization of the heat exchanger, user defined"
        # The minimum number of cells desired is n_disc
        # n_disc is increased by 1. Now it is the minimum number of cells boundaries.
        n_disc = self.n_disc + 1
        
        # Force the minimum number of discretizations needed
        n_disc = max(n_disc, 2) #Three zone boundary model
        
        # The original Bell2015 code uses: "if h_cdew is not None". I don't see when that could happen without an error before.
        "3.1) Create a vector of enthalpies : initiates the cell boundaries. If phase change, dew and bubble point are added"
        
        # Cold side
        if not self.Transcritical_c:
            self.hvec_c = np.append(np.linspace(self.h_ci, self.h_co, n_disc), [self.h_cdew, self.h_cbubble])
        elif self.Transcritical_c:
            #In transcritical, just dont append the unexisting dew and bubble enthalpies
            self.hvec_c = np.linspace(self.h_ci, self.h_co, n_disc)
            
        self.hvec_c = np.sort(np.unique(self.hvec_c))
        # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
        self.hvec_c = [h for h in self.hvec_c if (h >= self.h_ci and h <= self.h_co)]
        
        # Hot side
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.hvec_h = np.append(np.linspace(self.h_hi, self.h_ho, n_disc), [self.h_hdew, self.h_hbubble])
        elif self.Transcritical_h or self.h_incomp_flag:
            #In transcritical, just dont append the unexisting dew and bubble enthalpies
            self.hvec_h = np.linspace(self.h_hi, self.h_ho, n_disc)
            
        self.hvec_h = np.sort(np.unique(self.hvec_h))
        # Ensure that the enthalpy boundaries are all enclosed by the inlet and outlet conditions
        self.hvec_h = [h for h in self.hvec_h if (h >= self.h_ho and h <= self.h_hi)]

        "4) Filling of the complementary cell boundaries"
        
        # Complementary cell boundaries serve to keep heat balances in check
        # Start at the first element in the vector
        k = 0
        while k < len(self.hvec_c)-1 or k < len(self.hvec_h)-1:
            if len(self.hvec_c) == 2 and len(self.hvec_h) == 2: # If there is only one cell
                break

            # Determine which stream is the limiting next cell boundary
            Qcell_hk = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])
            Qcell_ck = self.mdot_c*(self.hvec_c[k+1]-self.hvec_c[k])

            if debug:
              print("k+1 Hot Enthalpy:", self.hvec_h[k+1], "\n")
              print("k Hot Enthalpy:", self.hvec_h[k],"\n")
              print("mdot_h = ", self.mdot_h, "\n")
              print("Qcell_hk =", Qcell_hk, "\n-------\n")
              print("k+1 Cold Enthalpy:", self.hvec_c[k+1], "\n")
              print("k Cold Enthalpy:", self.hvec_c[k],"\n")
              print("mdot_c = ", self.mdot_c, "\n")
              print("Qcell_hk =", Qcell_ck, "\n-------\n")

            if round(Qcell_hk, 4) > round(Qcell_ck, 4):
                # Hot stream needs a complementary cell boundary as the available heat is higher than the cold cell heat absorption
                self.hvec_h.insert(k+1, self.hvec_h[k] + Qcell_ck/self.mdot_h)
                
            elif round(Qcell_hk, 4) < round(Qcell_ck, 4):
                # Cold stream needs a complementary cell boundary as the available heat absoprtion is higher than the hot cell heat 
                self.hvec_c.insert(k+1, self.hvec_c[k] + Qcell_hk/self.mdot_c)

            if debug:
                print(k,len(self.hvec_c),len(self.hvec_h),Qcell_hk, Qcell_ck)

            # Increment index
            k += 1

        if debug:
             print("Modified Length of hvec_c is ", len(self.hvec_c))
             print("Modified Length of hvec_h is ", len(self.hvec_h))

        assert(len(self.hvec_h) == len(self.hvec_c)) # Ensures that the two fluid streams have the same number of cell boundaries

        # Calculate the vectors of exchanged heat in each cell.
        self.Qvec_h = np.array([self.mdot_h*(self.hvec_h[i+1]-self.hvec_h[i]) for i in range(len(self.hvec_h)-1)])
        self.Qvec_c = np.array([self.mdot_c*(self.hvec_c[i+1]-self.hvec_c[i]) for i in range(len(self.hvec_c)-1)])

        if debug:
            if np.max(np.abs(self.Qvec_c/self.Qvec_h))<1e-5:
                print(self.Qvec_h, self.Qvec_c)
            
        "5) Computation of the thermal states at the cell boundaries"
        
        # Build the normalized enthalpy vectors (normalized by the total exchanged enthalpy)
        self.hnorm_h = (np.array(self.hvec_h)-self.hvec_h[0])/(self.hvec_h[-1]-self.hvec_h[0])
        self.hnorm_c = (np.array(self.hvec_c)-self.hvec_c[0])/(self.hvec_c[-1]-self.hvec_c[0])

        # Pressure distribution accross the HX. Method by RDickes (not validated).
        # Discretization of the pressure as a linear interpolation on enthalpy between in and out conditions # !!! Is this good ? 
        self.pvec_c = (1-self.hnorm_c)*self.p_ci + self.hnorm_c*self.p_co
        self.pvec_h = (1-self.hnorm_h)*self.p_hi + self.hnorm_h*self.p_ho

        # Calculate the temperature and entropy at each cell boundary
        self.Tvec_c = CP.PropsSI('T','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)
        
        self.Tvec_h = CP.PropsSI('T','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
        self.svec_h = CP.PropsSI('S','H',self.hvec_h,'P',self.pvec_h,self.H_su.fluid)
        
        self.svec_c = CP.PropsSI('S','H',self.hvec_c,'P',self.pvec_c,self.C_su.fluid)


        "6) Vapour quality calculation if in the two-phase zone. If not, force, the value to -2 or 2"
        # Take Transcritical into account. in that Case, force x = 3
        
        "6.1) Hot Fluid"
        self.x_vec_h = []
        for i, h_h in enumerate(self.hvec_h):
            if not self.Transcritical_h and not self.h_incomp_flag:
                if h_h < self.h_hbubble:
                    self.x_vec_h.append(-2)        
                elif h_h > self.h_hdew:
                    self.x_vec_h.append(2)
                else:
                    self.x_vec_h.append(min(1, max(0, (h_h - self.h_hbubble)/(self.h_hdew - self.h_hbubble))))
            else:
                self.x_vec_h.append(3)
                
        "6.2) Cold Fluid"
        self.x_vec_c = []
        for i, h_c in enumerate(self.hvec_c):
            if not self.Transcritical_c:
                if h_c < self.h_cbubble:
                    self.x_vec_c.append(-2)        
                elif h_c > self.h_cdew:
                    self.x_vec_c.append(2)
                else:
                    self.x_vec_c.append(min(1, max(0, (h_c - self.h_cbubble)/(self.h_cdew - self.h_cbubble))))
            else:
                self.x_vec_c.append(3)

        "7) Vector of saturation temperatures at quality = 0.5, considering the pressure drop at each cell"
        
        if not self.Transcritical_c:
            self.Tvec_sat_pure_c = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_c, self.C_su.fluid)
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.Tvec_sat_pure_h = CP.PropsSI("T", "Q", 0.5, "P", self.pvec_h, self.H_su.fluid)
            
        "8) Calculate pinch and heat exchanger border temperature deltas"
        
        self.DT_pinch = np.min((np.array(self.Tvec_h) - np.array(self.Tvec_c)))
        self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
        self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

#%%
    def internal_pinching(self, stream):
        """
        Determine the maximum heat transfer rate based on the internal pinching analysis
        NB : external pinching analysis has already been done
        """
        
        # Try to find the dew point enthalpy as one of the cell boundaries
        # that is not the inlet or outlet
        # Do this only if the fluid is sub-critical

        "1) Check for the hot stream"
        
        if stream == 'hot':
            if not self.Transcritical_h and not self.h_incomp_flag: # If the hot side is not transcritical
                for i in range(1,len(self.hvec_h)-1):
                    
                    # Check if enthalpy is equal to the dewpoint enthalpy of hot stream and hot stream is colder than cold stream (impossible)
                    if (abs(self.hvec_h[i] - self.h_hdew) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        
                        # Enthalpy of the cold stream at the pinch temperature
                        h_c_pinch = CP.PropsSI('H','T',self.T_hdew,'P',self.pvec_c[i], self.C_su.fluid) # Equation 10 (Bell et al. 2015)

                        # Heat transfer in the cell
                        Qright = self.mdot_h*(self.h_hi-self.h_hdew) # Equation 9 (Bell et al. 2015)
                        
                        # New value for the limiting heat transfer rate
                        Qmax = self.mdot_c*(h_c_pinch-self.h_ci) + Qright # Equation 12

                        # Recalculate the cell boundaries
                        self.calculate_cell_boundaries(Qmax)
                
                        return Qmax
                    
            elif self.Transcritical_h:
                #If the hot side is transcritical, do nothing
                pass
            
            "2) Check for the cold stream"
            
        elif stream == 'cold':
            if not self.Transcritical_c:
                # Check for the cold stream pinch point
                for i in range(1,len(self.hvec_c)-1):
                    
                    # Check if enthalpy is equal to the bubblepoint enthalpy of cold stream and hot stream is colder than cold stream (impossible)
                    if (abs(self.hvec_c[i] - self.h_cbubble) < 1e-6 and self.Tvec_c[i] > self.Tvec_h[i]):
                        
                        # Enthalpy of the cold stream at the pinch temperatur
                        h_h_pinch = CP.PropsSI('H','T',self.T_cbubble,'P',self.pvec_h[i], self.H_su.fluid) # Equation 14 (Bell et al. 2015)

                        # Heat transfer in the cell
                        Qleft = self.mdot_c*(self.h_cbubble-self.h_ci) # Equation 13 (Bell et al. 2015)

                        # New value for the limiting heat transfer rate
                        Qmax = Qleft + self.mdot_h*(self.h_hi-h_h_pinch) # Equation 16 (Bell et al. 2015)

                        # Recalculate the cell boundaries
                        self.calculate_cell_boundaries(Qmax)

                        return Qmax
                    
            elif self.Transcritical_c:
                #If the hot side is transcritical, do nothing
                pass
        else:
            raise ValueError
    
    #%%
    def solve(self, only_external = False, and_solve = True):
        
        """
        Parameters
        ----------
            - mdot_h : Hot fluid flowrate [kg/s]
            - p_hi : Hot fluid pressure [Pa]
            - h_hi : Hot fluid specific enthalpy [J/kg]
            - mdot_c : Cold fluid flowrate [kg/s]
            - p_ci : Cold fluid pressure [Pa]
            - h_ci : Cold fluid specific enthalpy [J/kg]
            - only_external (optional) : calls only exxternal_pinching function to determine Q_dot max if set to True (if there is no phase change)
            - and_solve (optional) : if set to true, solves the heat exchanger

        Returns
        -------
            - Q : Exchanged heat rate [W]

        """
        
        "1) Main Input variables"
        
        if self.wf_T == "C": # The working fluid is the cold fluid    
            self.H_su = self.point_su[1]
            self.C_su = self.point_su[0]
            
            self.H_ex = self.point_ex[1]
            self.C_ex = self.point_ex[0]   
            
        elif self.wf_T == "H": # The working fluid is the hot fluid    
            self.H_su = self.point_su[0]
            self.C_su = self.point_su[1]
            
            self.H_ex = self.point_ex[0]
            self.C_ex = self.point_ex[1]      

        self.h_incomp_flag = (self.H_su.fluid.find("INCOMP") != -1)
                
        # Hot fluid
        self.mdot_h = self.H_su.m_dot
        self.h_hi = self.H_su.h
        self.p_hi = self.H_su.p
        
        # Cold fluid 
        self.mdot_c = self.C_su.m_dot
        self.h_ci = self.C_su.h
        self.p_ci = self.C_su.p
                
        # Determine the inlet temperatures from the pressure/enthalpy pairs
        #Useless I think
        self.T_ci = CP.PropsSI('T', 'P', self.C_su.p, 'H', self.C_su.h, self.C_su.fluid)
        self.T_hi = CP.PropsSI('T', 'P', self.p_hi, 'H', self.h_hi, self.H_su.fluid)

        "2) Determine if the streams come in at a higher pressure than the transcritical pressure"
        
        # Initialization of the flags:    
        self.Transcritical_c = False
        self.Transcritical_h = False
        
        if not self.h_incomp_flag and (self.p_hi - CP.PropsSI("PCRIT", self.H_su.fluid)) >= 1e-06:
            self.Transcritical_h = True

        "3) Calculate the ideal bubble and dew temperatures/enthalpies for each stream IF the fluid is not transcritical"
        
        if not self.Transcritical_c:
            self.T_cbubble_ideal = CP.PropsSI('T', 'P', self.p_ci, 'Q', 0, self.C_su.fluid)
            self.T_cdew_ideal    = CP.PropsSI('T', 'P', self.p_ci, 'Q', 1, self.C_su.fluid)
            self.h_cbubble_ideal = CP.PropsSI('H', 'T', self.T_cbubble_ideal, 'Q', 0, self.C_su.fluid)
            self.h_cdew_ideal    = CP.PropsSI('H', 'T', self.T_cdew_ideal, 'Q', 1, self.C_su.fluid)  
            
        if not self.Transcritical_h and not self.h_incomp_flag:
            self.T_hbubble_ideal = CP.PropsSI('T', 'P', self.p_hi, 'Q', 0, self.H_su.fluid)
            self.T_hdew_ideal    = CP.PropsSI('T', 'P', self.p_hi, 'Q', 1, self.H_su.fluid)
            self.h_hbubble_ideal = CP.PropsSI('H', 'T', self.T_hbubble_ideal, 'Q', 0, self.H_su.fluid)
            self.h_hdew_ideal    = CP.PropsSI('H', 'T', self.T_hdew_ideal, 'Q', 1, self.H_su.fluid)
            
        "4) Calculate pressure drops"
        
        if self.H_DP_ON == True: # if the pressure drop are not neglected
            if self.H_su.fluid == 'Water':
                self.DP_h = self.H.f_dp["K"] * (self.mdot_h/self.H_su.D)**2
                self.p_ho = self.p_hi - self.DP_h
            else:
                self.DP_h = self.H.f_dp["K"] * self.mdot_h**(self.H.f_dp["B"]) # Empirical correlations : DP = K*m_dot**B
                self.p_ho = self.p_hi - self.DP_h
        else:
            self.DP_h = 0
            self.p_ho = self.p_hi
            
        if self.C.f_dp != {"K": 0, "B": 0}:
            if self.C_su.fluid == 'Water':
                self.DP_c = self.C.f_dp["K"] * (self.mdot_c/self.C_su.D)**2
                self.p_co = self.p_ci - self.DP_c
            else:
                self.DP_c = self.C.f_dp["K"] * self.mdot_c**(self.C.f_dp["B"])
                self.p_co = self.p_ci - self.DP_c
        else:
            self.DP_c = 0
            self.p_co = self.p_ci

        "5) Calculate maximum and actual heat rates"
        if (self.T_hi - self.T_ci) > 1e-2  and self.mdot_h  > 0 and self.mdot_c > 0: # Check that the operating conditions allow for heat transfer
            
            "5.1) Compute the external pinching & update cell boundaries"
            Qmax_ext = self.external_pinching() # Call to external-pinching procedure
            self.Qmax_ext = Qmax_ext
            Qmax = Qmax_ext
            
            if debug:
                print("External pinching calculation done. \n")
            
            "5.2) Compute the internal pinching & update cell boundaries"
            if not only_external: # If phase change is expected : Check the internal pinching
                for stream in ['hot','cold']:
                    Qmax_int = self.internal_pinching(stream) # Call to internal-pinching procedure
                    if Qmax_int is not None:
                        self.Qmax_int = Qmax_int
                        Qmax = Qmax_int
            
            # Maximum heat transfer rate determined by external or internal pinching
            self.Qmax = Qmax

            "5.3) Solve the heat exchanger to find the actual heat rate"
            if and_solve and not only_external:
                Q = self.solve_HX()

            self.epsilon_th = self.Q/self.Qmax # HTX efficiency
            self.residual = 1 - sum(self.w) # HTX residual # !!! (what is "w" ?)

            "5.4) Computation of the fluid mass inside the HTX"
            
            self.Vvec_h = self.geom.H_V_tot*np.array(self.w)
            self.Vvec_c = self.geom.C_V_tot*np.array(self.w)
            
            # Initiates the mass vectors # !!! (for each cell ?)
            self.Mvec_h = np.empty(len(self.hvec_h)-1)
            self.Mvec_c = np.empty(len(self.hvec_c)-1)
            
            # Mass computation # !!! understand + explain the formula
            for i in range(len(self.hvec_h)-1):
                
                self.Mvec_h[i] = self.Vvec_h[i]*0.5*(CP.PropsSI("D", "P", self.pvec_h[i], "H", self.hvec_h[i], self.H_su.fluid) + CP.PropsSI("D", "P", self.pvec_h[i+1], "H", self.hvec_h[i+1], self.H_su.fluid))
                self.Mvec_c[i] = self.Vvec_c[i]*0.5*(CP.PropsSI("D", "P", self.pvec_c[i], "H", self.hvec_c[i], self.C_su.fluid) + CP.PropsSI("D", "P", self.pvec_c[i+1], "H", self.hvec_c[i+1], self.C_su.fluid))

            if and_solve:
                # Hot side
                H_ex = Mass_connector()

                H_ex.set_fluid(self.H_su.fluid)
                H_ex.set_h(self.hvec_h[0]) #Vrai température déterminé sur base de Q
                H_ex.set_p(self.pvec_h[0]) 
                H_ex.set_m_dot(self.H_su.m_dot) 
                self.H_ex = H_ex

                if self.wf_T == "C": # The working fluid is the cold fluid    
                    self.point_ex[1] = H_ex
                    
                elif self.wf_T == "H": # The working fluid is the hot fluid    
                    self.point_ex[0] = H_ex
                    
                # Cold side
                C_ex = Mass_connector()
                
                C_ex.set_fluid(self.C_su.fluid)

                C_ex.set_h(self.hvec_c[-1]) # Example temperature [K] (last one in the vector Tvec_c)
                C_ex.set_p(self.pvec_c[-1]) # Example Pressure [Pa]
                C_ex.set_m_dot(self.C_su.m_dot) 
                # print(self.hvec_c) après vérification, l'erreur ne vient pas de là
                # print(self.h_ci + self.Q/self.mdot_c, self.hvec_c[-1], self.h_ci)
                # print(C_ex.T-273.15, self.Tvec_c[-1]-273.15)
                self.C_ex = C_ex
                
                if self.wf_T == "C": # The working fluid is the cold fluid    
                    self.point_ex[0] = C_ex
                    
                elif self.wf_T == "H": # The working fluid is the hot fluid    
                    self.point_ex[1] = C_ex
                                            
                self.defined = True
                return Q
        else: # Just a flag if the heat exchanger is not solved
            self.Q = 1
            
#%%    
    def objective_function(self, Q):
        "1) Initialize cell boundaries and results vectors"
        
        # Cell boundaries
        self.calculate_cell_boundaries(Q)

        # Initialize dry-out incipience identification
        self.x_di_c = 1
        self.dry_out_c = False
        self.ColdSide_Props_Calculated = False
        
        # Initialise results arrays
        w = [] # The wet-bulb potential temperature
        self.Avec_h = []
        self.alpha_c = []
        self.alpha_h = []
        self.phases_h = []
        self.phases_c = []

        "2) Iteration over all cells"
        
        # The length of the hvec vectors is the same. hvec_h is arbitrarily taken below
        for k in range(len(self.hvec_h)-1): # iterate over each cell
            
            "2.1) Diverse temperature parameter definition"
            
            Thi = self.Tvec_h[k+1] 
            Tci = self.Tvec_c[k]
            Tho = self.Tvec_h[k]
            Tco = self.Tvec_c[k+1]
            DTA = max(Thi - Tco, 1e-02)
            DTB = max(Tho - Tci, 1e-02)
            
            if DTA == DTB:
                LMTD = DTA
            else:
                try:
                    LMTD = (DTA-DTB)/log(abs(DTA/DTB))
                except ValueError as VE:
                    print(Q, DTA, DTB)
                    raise
            
            # Wall Temperature:
            T_wall = (Thi + Tho + Tci + Tco)/4
            
            # Cold side
            Tc_mean = 0.5*(Tci + Tco) # mean temperature over the cell
            if not self.Transcritical_c:
                Tc_sat_mean = 0.5*(self.Tvec_sat_pure_c[k] + self.Tvec_sat_pure_c[k+1]) # mean saturation temperature over the cell
            p_c_mean = 0.5*(self.pvec_c[k] + self.pvec_c[k+1]) # mean pressure over the cell
            
            # Hot side
            Th_mean = 0.5*(Thi + Tho) # mean temperature over the cell
            
            if not self.Transcritical_h and not self.h_incomp_flag:
                Th_sat_mean = 0.5*(self.Tvec_sat_pure_h[k] + self.Tvec_sat_pure_h[k+1]) # mean saturation temperature over the cell
            
            p_h_mean = 0.5*(self.pvec_h[k] + self.pvec_c[k+1]) # mean pressure over the cell

            "2.2) Hot side phase identification"
            
            # If not transcritical
            if not self.Transcritical_h and not self.h_incomp_flag:
                
                havg_h = (self.hvec_h[k] + self.hvec_h[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_h < self.h_hbubble: # below evaporation/bubble point
                    self.phases_h.append('liquid')
                    T_wall_h = T_wall
                    
                elif havg_h > self.h_hdew: # higher than condensation/dew point
                    #T_wall_h not used so far but it could appear in future Correlations implemented.
                    T_wall_h = max(T_wall, Th_sat_mean)
                    
                    if T_wall > Th_sat_mean:
                        self.phases_h.append('vapor')
                    else:
                        self.phases_h.append('vapor-wet')
                        
                else: # Between evaporation/bubble and condensation/dew point
                    self.phases_h.append('two-phase')
                    T_wall_h = T_wall
                    
            elif self.h_incomp_flag:
                self.phases_h.append('liquid')
                T_wall_h = T_wall
                
            elif self.Transcritical_h:
                self.phases_h.append("transcritical")
                T_wall_h = T_wall
            
            "2.3) Cold side phase identification"

            # If not transcritical
            if not self.Transcritical_c:
                
                havg_c = (self.hvec_c[k] + self.hvec_c[k+1])/2.0 # Average cell enthalpy over the cell
                
                if havg_c < self.h_cbubble: # below evaporation/bubble point
                    self.phases_c.append('liquid')
                    #The line comented below is the original line on the code. It is causing trouble.
                    #T_wall_c = min(T_wall, Tc_sat_mean)
                    T_wall_c = T_wall
                    
                elif havg_c > self.h_cdew: # higher than condensation/dew point
                    self.phases_c.append("vapor")
                    T_wall_c = T_wall
                    
                else: # Between evaporation/bubble and condensation/dew point
                    T_wall_c = T_wall
                    
                    if (0.5*self.x_vec_c[k] + 0.5*self.x_vec_c[k+1]) <= self.x_di_c:
                        self.phases_c.append('two-phase')  
                        
                    else:
                        self.phases_c.append("two-phase-dryout")
                        self.dry_out_c = True
                        
            elif self.Transcritical_c:
                self.phase_c.append('transcritical')
                
            
            G_c = (self.mdot_c/self.geom.C_n_canals)/self.geom.C_CS
            G_h = (self.mdot_h/self.geom.H_n_canals)/self.geom.H_CS
            
            # F correction factor for LMTD method:
            if self.typeHEX == "CrossFlow":
                self.R = (Thi - Tho)/(Tco - Tci)
                self.P = (Tco - Tci)/(Thi - Tci)
                self.F = F_lmtd2(self.R, self.P)[0]
            else:
                self.F = 1
            
            #UA_req including F correction factor in case of cross flow
            UA_req = self.mdot_h*(self.hvec_h[k+1]-self.hvec_h[k])/(self.F*LMTD)          
            
            "3) Cell heat transfer coefficients"
            
            "3.1) Hot side - User defined"
            
            if self.H.HeatExchange_Correlation == "User-Defined":
                if self.phases_h[k] == "liquid":
                    alpha_h = self.params.H.h_liq
                elif self.phases_h[k] == "vapor":
                    alpha_h = self.params.H.h_vap
                elif self.phases_h[k] == "two-phase":
                    alpha_h = self.params.H.h_twophase
                elif self.phases_h[k] == "vapor-wet":
                    alpha_h = self.params.H.h_vapwet
                elif self.phases_h[k] == "two-phase-dryout":
                    alpha_h = self.params.H.h_tpdryout
                elif self.phases_h[k] == "transcritical":
                    alpha_h = self.params.H.h_transcrit
                    
            elif self.H.HeatExchange_Correlation == "Correlation": 
                
                # Heat transfer coefficient calculated from Correlations.
                # Evaluate phases
                
                # 1 phase case
                
                if self.phases_h[k] == "liquid" or self.phases_h[k] == "vapor" or self.phases_h[k] == "vapor-wet" or self.phases_h[k] == "transcritical":
                    mu_h, Pr_h, k_h, mu_wall, mu_rat, _, _ = PropsFluid(Th_mean, p_h_mean, T_wall_h, self.H_su.fluid, self.h_incomp_flag)
                    
                    if self.H_su.fluid == 'water':
                        alpha_h = water_Plate_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh)
                    
                    if self.H.Correlation_1phase  == "Gnielinski":
                        alpha_h, _ = Gnielinski_Pipe_HTC(mu_h, Pr_h, k_h, G_h, self.geom.H_Dh, self.geom.l)
                        #print(alpha_h) #Whaaaaat c'est des nombres complexes!
                
                # 2 phases case
                elif self.phases_h[k] == "two-phase" or self.phases_h[k] == "vapor-wet":
                    if self.phases_h[k] == "two-phase":
                        x_h = min(1, max(0, 0.5*(self.x_vec_h[k] + self.x_vec_h[k])))
                    elif self.phases_h[k] == "vapor-wet":
                        x_h = 1
                    
                    # Thermodynamical variables
                    if self.H_su.fluid == 'R1233zd(E)' or self.H_su.fluid == 'R1233ZDE':
                        T_h_l = CP.PropsSI("T", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        k_h_l = conducticity_R1233zd(T_h_l, p_h_mean)

                        mu_h_l = CP.PropsSI("V", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        cp_h_l = CP.PropsSI("C", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        Pr_h_l = (mu_h_l*cp_h_l)/k_h_l
                        rho_h_l = CP.PropsSI("D", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        rho_h_v = CP.PropsSI("D", "Q", 1, "P", p_h_mean, self.H_su.fluid)
                        i_fg_h = CP.PropsSI('H', 'Q', 1, 'P', p_h_mean, self.H_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_h_mean, self.H_su.fluid)

                    else:    
                        mu_h_l = CP.PropsSI("V", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        k_h_l = CP.PropsSI("L", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        Pr_h_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        rho_h_l = CP.PropsSI("D", "Q", 0, "P", p_h_mean, self.H_su.fluid)
                        rho_h_v = CP.PropsSI("D", "Q", 1, "P", p_h_mean, self.H_su.fluid)
                        i_fg_h = CP.PropsSI('H', 'Q', 1, 'P', p_h_mean, self.H_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_h_mean, self.H_su.fluid)
                        
                    # i_fg_h = CP.PropsSI('H', 'Q', 1, 'P', p_h_mean, self.H_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_h_mean, self.H_su.fluid)
                    
                    # !!! Include different types of Correlation HERE
                    if self.H.Correlation_2phase == "Han_cond_BPHEX":
                        alpha_h_2phase, _, DP_H = Han_Cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, G_h, self.geom.H_Dh, self.geom.plate_pitch_co, self.geom.chevron_angle, self.geom.l_v, self.geom.H_n_canals, self.H_su.m_dot, self.geom.H_canal_t)

                    if self.phases_h[k] == "two-phase":
                        alpha_h = alpha_h_2phase
                        
                    elif self.phases_h[k] == "vapor-wet":
                        w_vap_wet = (Th_mean - Th_sat_mean)/(Th_mean - T_wall)
                        # The line below re-calculates alpha_h in case of having a vapor-wet condition
                        alpha_h = alpha_h_2phase - w_vap_wet*(alpha_h_2phase - alpha_h) # the last alpha_h in this equation is the 1 Phase calculation

            "3.2) Cold side - User defined"
            
            # User Defined
            if self.C.HeatExchange_Correlation == "User-Defined":
                if self.phases_c[k] == "liquid":
                    alpha_c = self.params.C.h_liq
                elif self.phases_c[k] == "vapor":
                    alpha_c = self.params.C.h_vap
                elif self.phases_c[k] == "two-phase":
                    alpha_c = self.params.C.h_twophase
                elif self.phases_c[k] == "vapor-wet":
                    alpha_c = self.params.C.h_vapwet
                elif self.phases_c[k] == "two-phase-dryout":
                    alpha_c = self.params.C.h_tpdryout
                elif self.phases_c[k] == "transcritical":
                    alpha_c = self.params.C.h_transcrit
                    
            elif self.C.HeatExchange_Correlation == "Correlation":

                # Heat transfer coefficient calculated from Correlations:
                # Evaluate phases
                
                # 1 phase case
                if self.phases_c[k] == "liquid" or self.phases_c[k] == "vapor" or self.phases_c[k] == "vapor-wet" or self.phases_c[k] == "transcritical":
                    mu_c, Pr_c, k_c, mu_wall, mu_rat, _, _ = PropsFluid(Tc_mean, p_c_mean, T_wall_c, self.C_su.fluid, False)
                    
                    if self.C_su.fluid == 'water':
                        alpha_c = water_Plate_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh)
                    
                    if self.C.Correlation_1phase  == "Gnielinski":
                        alpha_c, _ = Gnielinski_Pipe_HTC(mu_c, Pr_c, k_c, G_c, self.geom.H_Dh, self.geom.l)

                # 2 phases case
                elif self.phases_c[k] == "two-phase" or self.phases_c[k] == "vapor-wet":
                    if self.phases_c[k] == "two-phase":
                        x_c = min(1, max(0, 0.5*(self.x_vec_c[k] + self.x_vec_c[k])))
                    elif self.phases_c[k] == "vapor-wet":
                        x_c = 1
                        
                    # Thermodynamical variables
                    if self.C_su.fluid == 'R1233zd(E)' or self.C_su.fluid == 'R1233ZDE':
                        T_c_l = CP.PropsSI("T", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        k_c_l = conducticity_R1233zd(T_c_l, p_c_mean)

                        mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        cp_c_l = CP.PropsSI("C", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        Pr_c_l = (mu_c_l*cp_c_l)/k_c_l
                        rho_c_l = CP.PropsSI("D", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        rho_c_v = CP.PropsSI("D", "Q", 1, "P", p_c_mean, self.C_su.fluid)
                        i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)

                    else:    
                        mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        k_c_l = CP.PropsSI("L", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        Pr_c_l = CP.PropsSI("Prandtl", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        rho_c_l = CP.PropsSI("D", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        rho_c_v = CP.PropsSI("D", "Q", 1, "P", p_c_mean, self.C_su.fluid)
                        i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                        
                    # This boolean serves to spare to recalculate this properties in the dry-out analysis further ahead.
                    self.ColdSide_Props_Calculated = True
                    
                    # !!! Include different types of Correlation HERE
                    if self.C.Correlation_2phase == "Han_cond_BPHEX":
                        alpha_c_2phase, _ = Han_Cond_BPHEX_HTC(x_c, mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v, G_c, self.params.C.Dh, self.params.plate_pitch_co, self.params.chevron_angle)
                    
                    elif self.C.Correlation_2phase == "Han_Boiling_BPHEX_HTC":
                        alpha_c_2phase, _ = Han_Boiling_BPHEX_HTC(min(x_c, self.x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, LMTD*self.F, self.Qvec_c[k], alpha_h, self.geom.C_Dh, self.geom.chevron_angle, self.geom.plate_pitch_co)
                    
                    else:
                        raise ValueError("Correlation not found for Cold Side 2-Phase")
                            
                    if self.phases_c[k] == "two-phase":
                        alpha_c = alpha_c_2phase
                    elif self.phases_c[k] == "vapor-wet":
                        w_vap_wet = (Tc_mean - Tc_sat_mean)/(Tc_mean - T_wall)
                        #The line below re-calculates alpha_h in case of having a vapor-wet condition
                        alpha_c = alpha_c_2phase - w_vap_wet*(alpha_c_2phase - alpha_c) #the last alpha_c in this equation is the 1 Phase calculation
            
            "3.3) Store the heat transfer coefficients"
            
            self.alpha_c.append(alpha_c)
            self.alpha_h.append(alpha_h)
            
            "4) Recover some parameters and conductivity calculation"
            
            t = self.geom.t_plates
            # R_fooling_h = self.params.H.R_fooling 
            # R_fooling_c = self.params.C.R_fooling
            
            "5) Calculate w, the main objective of this function. This variable serves for residual minimization in the solver"
            
            # In the equation below, thickness resistance is given with respect to A_h arbitrarely
            UA_avail = 1/((1+self.geom.fooling)/(alpha_h*self.geom.A_h) + 1/(alpha_c*self.geom.A_c) + t/(self.geom.plate_cond)) 
            w.append(UA_req/UA_avail)
            
            self.w = w
            self.Avec_h.append(w[k]*self.geom.A_h)
            
            "5.1) Calculate dryout incipience only if subcritical"
            if not self.dry_out_c and self.phases_h[k] == "two-phase" and not self.Transcritical_c:
                    if not self.ColdSide_Props_Calculated:
                        #If these properties where not calculated before, do it:
                        mu_c_l = CP.PropsSI("V", "Q", 0, "P", p_c_mean, self.C_su.fluid)
                        rho_c_l = CP.PropsSI('D', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                        rho_c_v = CP.PropsSI('D', 'Q', 1, 'P', p_c_mean, self.H_su.fluid)
                        i_fg_c = CP.PropsSI('H', 'Q', 1, 'P', p_c_mean, self.C_su.fluid) - CP.PropsSI('H', 'Q', 0, 'P', p_c_mean, self.C_su.fluid)
                    
                    P_star_c = p_c_mean/CP.PropsSI("Pcrit", "Q", 1, "P", p_c_mean, self.C_su.fluid)
                    q_c = self.Qvec_c[k]/self.Avec_h[k]
                    sigma_c_l = CP.PropsSI("I", 'Q', 0, "P", p_c_mean, self.C_su.fluid)
                    self.x_di_c = Kim_DryOutIncipience(G_c, q_c, self.geom.C_Dh, P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c)
            
        if debug:
            print(Q, 1-sum(w))
        # print(Q, 1-sum(w))  
        return 1-sum(w)

#%%
    def solve_HX(self):
        """ 
        Solve the objective function using Brent's method and the maximum heat transfer 
        rate calculated from the pinching analysis
        """
        # print(self.objective_function(1e-5))
        # print(self.objective_function(self.Qmax))
        self.Q = scipy.optimize.brentq(self.objective_function, 1e-5, self.Qmax+2000, rtol = 1e-14, xtol = 1e-10)
        return self.Q
    
#%%
    def plot_objective_function(self, N = 100):
        """ Plot the objective function """
        # print(self.Qmax)
        Q = np.linspace(1e-5,self.Qmax,N)
        r = np.array([self.objective_function(_Q) for _Q in Q]) #Residuals!
        fig, ax = plt.subplots()
        ax.plot(Q, r)
        ax.grid(True)
        # plt.show()
        
#%%
    def plot_ph_pair(self):
        import warnings
        warnings.simplefilter("ignore")
        """ Plot p-h plots for the pair of working fluids """
        Ph_diagram_h = PropertyPlot(self.H_su.fluid, "PH", unit_system = "EUR")
        plt.plot(0.001*np.array(self.hvec_h), self.pvec_h*1e-05,'s-')
        Ph_diagram_h.calc_isolines()
        Ph_diagram_h.title("Hot side P-h diagram. Fluid: " + str(self.H_su.fluid))
        Ph_diagram_h.show()
        #----------------------------------------------------------------------
        Ph_diagram_c = PropertyPlot(self.C_su.fluid, "PH", unit_system = "EUR")
        plt.plot(0.001*np.array(self.hvec_c), self.pvec_c*1e-05,'s-')
        Ph_diagram_c.calc_isolines()
        Ph_diagram_c.title("Cold side P-h diagram. Fluid: " + str(self.C_su.fluid))
        Ph_diagram_c.show()
        warnings.simplefilter("default")

#%%
    def plot_Ts_pair(self):
        import warnings
        warnings.simplefilter("ignore")
        """ Plot a T-s plot for the pair of working fluids """
        Ph_diagram_h = PropertyPlot(self.H_su.fluid, "TS", unit_system = "EUR")
        plt.plot(0.001*self.svec_h,self.Tvec_h - 273.15,'s-')
        Ph_diagram_h.calc_isolines()
        Ph_diagram_h.title("Hot side T-s diagram. Fluid: " + str(self.H_su.fluid))
        Ph_diagram_h.show()
        # #----------------------------------------------------------------------
        Ph_diagram_c = PropertyPlot(self.C_su.fluid, "TS", unit_system = "EUR")
        plt.plot(0.001*self.svec_c,self.Tvec_c - 273.15,'s-')
        Ph_diagram_c.calc_isolines()
        Ph_diagram_c.title("Cold side T-s diagram. Fluid: " + str(self.C_su.fluid))
        Ph_diagram_c.show()
        warnings.simplefilter("default")

#%%
    def plot_cells(self, fName = '', dpi = 400):
        """ Plot the cells of the heat exchanger """
        plt.figure(figsize = (4,3))
        plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
        plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
        plt.xlim(0,1)
        plt.ylabel('T [K]') 
        plt.xlabel(r'$\hat h$ [-]')
        plt.grid(True)
        plt.tight_layout(pad = 0.2) 
        if fName != '':
            plt.savefig(fName, dpi = dpi)

#%%

def Plate_HX(wf_T, H_su, C_su, HX_geom, htc_type, flow_type, n_disc, DP_H_ON, DP_C_ON, calc, plot, print_flag):
    import time
    start_time = time.time()   

    # #Problem ici!!!    
    if H_su.T < C_su.T:
        temp = H_su
        H_su = C_su
        C_su = temp
        print("Hot and cold sides were inverted")
    
    HX = Plate_HeatExchanger() 
    
    if wf_T == "C": # The working fluid is the cold fluid    
        HX.set_parameters(point_su = [C_su,H_su], geom = HX_geom, n_disc = n_disc, wf_T = wf_T, flow_type = flow_type, htc_type = htc_type, H_DP_ON = DP_H_ON, C_DP_ON = DP_C_ON) 
    elif wf_T == "H": # The working fluid is the hot fluid    
        HX.set_parameters(point_su = [H_su,C_su], geom = HX_geom, n_disc = n_disc, wf_T = wf_T, flow_type = flow_type, htc_type = htc_type, H_DP_ON = DP_H_ON, C_DP_ON = DP_C_ON)
        
    
    if calc == 1:
        HX.check_calculable()
        HX.check_parametrized()
        if HX.calculable and HX.parametrized:
            #Actually run the HX code
            HX.solve()
        
        if print_flag == 1: 
            
            # print("Qmax =",HX.Qmax, "[W]")
            print("Q =", HX.Q, "[W]")
            print("-----------------------------------------------------------------------")
            print("Pinch point is:", HX.DT_pinch, "[K]")
            print("HX efficiency is:", HX.epsilon_th, "[-]")
            print("-----------------------------------------------------------------------")
            print("Phases in the hot side cells are: ", HX.phases_h)
            print("-----------------------------------------------------------------------")
            print("Phases in the cold side cells are: ", HX.phases_c)
            print("-----------------------------------------------------------------------")
            print("--- %s seconds ---" % (time.time() - start_time))
            print("\n")
            
        if plot == 1:
            HX.plot_cells('full.png')
        
        # Hot side
        H_ex = Mass_connector()
        
        H_ex.set_fluid(H_su.fluid)
        H_ex.set_T(HX.Tvec_h[-1])
        H_ex.set_p(HX.pvec_h[-1]) 
        H_ex.set_m_dot(HX.H_su.m_dot) 
        
        HX.H_ex = H_ex
        
        if wf_T == "C": # The working fluid is the cold fluid    
            HX.point_ex[1] = H_ex
            
        elif wf_T == "H": # The working fluid is the hot fluid    
            HX.point_ex[0] = H_ex
            
        # Cold side
        C_ex = Mass_connector()
        
        C_ex.set_fluid(C_su.fluid)
        C_ex.set_T(HX.Tvec_c[-1]) # Example temperature [K]
        C_ex.set_p(HX.pvec_c[-1]) # Example Pressure [Pa]
        C_ex.set_m_dot(HX.C_su.m_dot) 

        HX.C_ex = C_ex
        
        if wf_T == "C": # The working fluid is the cold fluid    
            HX.point_ex[0] = C_ex
            
        elif wf_T == "H": # The working fluid is the hot fluid    
            HX.point_ex[1] = C_ex
                    
        print("\n")
        
        HX.defined = True

        return HX, H_ex, C_ex
    
    print("\n")
    return HX, 0, 0

