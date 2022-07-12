# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:50:35 2020

@author: A. Praetorius and PradoDomercq
"""
import pandas as pd 
import sys


from objects.EnvCompartmentRiver import*
#class to generate environmental 

#import file storing required constants
from helpers.GlobalConstants import *


#import classes to generate objects
from objects.Particulates import Particulates #class to generate MP and SPM objects
from objects.ParticulatesBF import ParticulatesBF #class to generate MP and SPM objects
from objects.ParticulatesSPM import ParticulatesSPM #class to generate MP and SPM objects


#function to read lake data
def preProcessLayers(mode, compartments_prop, date, rs_indexs,riverSection):
    try:
        df= pd.DataFrame()
        if mode =="Standard":
            df = compartments_prop.loc[rs_indexs]
        elif mode == "Monthly":
            date =str(date.year)+"-"+str('%02d' % date.month)
            df = compartments_prop[compartments_prop.date == str(date)]#To be updated
                    
        surface = EnvCompartment(df, 0+4*int(riverSection))
        flowingWater = EnvCompartment(df, 1+4*int(riverSection))
        stagnantWater = EnvCompartment(df, 2+4*int(riverSection))
        sediment = EnvCompartment(df, 3+4*int(riverSection))
        
        surface.calc_dimensions() #Compartment 1
        flowingWater.calc_dimensions() #Compartment 2
        stagnantWater.calc_dimensions() #Compartment 3
        sediment.calc_dimensions() #Compartment 4
        
    except:    
        print("Error in data Preparation")
        sys.exit(1)  

    return surface, flowingWater, stagnantWater, sediment  

#fuction to genrate particles objects MPs and SPM
def preProcessElements(MP_prop, MP_index, SPM_index,compartments_prop, comp_index):
    #generate MicroPlastic object(s) --> A: pristine (free MP)
    
    #MP_index = 0 #currently only runnign for the first MP in the list

    MP1= Particulates(MP_prop, MP_index)
    MP1.calc_volume()
        
    #generate SPM object(s)
    #SPM_index = 11 #need to move SPM in own input file
    SPM1 = Particulates(MP_prop, SPM_index)
    SPM1.calc_volume()
    SPM1.calc_numConc(compartments_prop.SPM_mgL[comp_index], 0) #move this to river input file
    #SPM1.calc_settling(density_w_21C_kg_m3, mu_w_21C_kg_ms, g_m_s2, "Stokes")
    
    #B: heteroaggregated (MP attached to suspended particulate matter (SPM)) 
    MP1_SPM = ParticulatesSPM("MP1-SPM", MP1, SPM1) 
    MP1_SPM.calc_volume(MP1, SPM1)
    #MP1_SPM.calc_settling()
    
    #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
    MP1_BF = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
    MP1_BF.calc_volume()
    #MP1_BF.calc_settling()

    #D: biofilm-heteroaggregated (MP with BF layer attached to SPM)
    MP1_BF_SPM = ParticulatesSPM("MP1-BF-SPM", MP1_BF, SPM1) 
    MP1_BF_SPM.calc_volume(MP1_BF, SPM1)
    #MP1_BF_SPM.calc_settling()
   
    
    return MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM
