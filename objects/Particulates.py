# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:10:20 2019

@author: AntoniaPraetorius

#PradoDomercq added some modifications for Size bins

"""
import math
from helpers.GlobalConstants import *

#define Particulates class
class Particulates:
    "This is a class to create Particulate objects" 
    
    #class attribute
    species = "particulate"
    
    #constructor
    def __init__(self, plastic_prop, MP_index):
        self.name = plastic_prop.name.loc[MP_index]
        self.composition = plastic_prop.composition.loc[MP_index]
        self.density_kg_m3 = plastic_prop.density_kg_m3.loc[MP_index]
        self.shape = plastic_prop.MPshape.loc[MP_index]
        self.diameter_um = plastic_prop.diameter_um.loc[MP_index] #for spherical MPs and fibres. Should be 0 for all others.
        self.diameter_m = self.diameter_um*10**-6 #for spherical MPs and fibres. Should be 0 for all others.
        self.radius_m = self.diameter_um*10**-6/2
        self.length_a_um = plastic_prop.length_a_um.loc[MP_index] #longest length (for nonspherical MPs)
        self.length_a_m = self.length_a_um*10**-6 
        self.length_b_um = plastic_prop.length_b_um.loc[MP_index] #intermediate length (for nonspherical MPs)
        self.length_b_m = self.length_b_um*10**-6 
        self.length_c_um = plastic_prop.length_c_um.loc[MP_index] #shortest length (for nonspherical MPs)
        self.length_c_m = self.length_c_um*10**-6 
        
        
        
    #methods
    
    #volume calculation
    #different formulas for different particle shapes.
    #currently defined for spheres, fibres, cylinders, pellets and irregular fragments
    def calc_volume(self):
        
        if self.shape == "sphere":
            self.volume_m3 = 4/3*math.pi*(self.radius_m)**3
            #calculates volume (in m3) of spherical particles from MP radius  
            self.CSF = 1    
            #calculate corey shape factor (CSF) 
            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
            #particle number concentration calculation
         
        elif self.shape == "fibre" or self.shape == "fiber" or self.shape == "cylinder":
            self.volume_m3 = math.pi*(self.radius_m)**2*self.length_a_m
            #calculates volume (in m3) of fibres or cylinders from diameter and  
            #length assuming cylindrical shape
            self.CSF = self.radius_m/math.sqrt(self.length_a_m*self.radius_m)    
            #calculate corey shape factor (CSF) 
            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
            #particle number concentration calculation
           
        elif self.shape == "pellet" or self.shape == "fragment":
            self.volume_m3 = self.length_a_m*self.length_b_m*self.length_c_m
            #approximate volume calculation for irregular fragments
            #approximated as a cuboid using longest, intermediate and shortest length
            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
            #might adjust later!!
            self.CSF = self.length_c_m/math.sqrt(self.length_a_m*self.length_b_m)    
            #calculate corey shape factor (CSF) 
            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
            #particle number concentration calculation
            
        else:
            print("Error: unknown shape")
            #print error message for shapes other than spheres 
            #(to be removed when other volume calculations are implemented)
                   
        
        
        
    def calc_numConc(self, concMass_mg_L, concNum_part_L):
        
        if concNum_part_L == 0:
            self.concNum_part_m3 = concMass_mg_L/1000/self.density_kg_m3/self.volume_m3
            #if mass concentration is given, it is converted to number concentration
        else:
            self.concNum_part_m3 = concNum_part_L*1000
            #if number concentration is given, it is converted from part/L to part/m3
    
    def calc_projected_area(self):
        
        if self.shape == "sphere":
            self.projected_area_m2 = math.pi*(self.radius_m)**2
            #calculates projected area (m2) of spherical particles from MP radius
             
        elif self.shape == "fibre" or self.shape == "fiber" or self.shape == "cylinder":
            self.projected_area_m2 = math.pi*(self.radius_m)**2
            #calculates projected area (m2) of fibres or cylinders considering it as a circle
            #(reduces changes impact, as it looks more as a sphere)
            
        elif self.shape == "pellet" or self.shape == "fragment":
            self.projected_area_m2 = self.length_b_m*self.length_c_m
            #approximate projected area (m2) calculation for irregular fragments
            #approximated as a cuboid using intermediate and shortest length, making it closer to what the original model did 
            #(reduces changes impact)
            
        else:
            print("Error: unknown shape")
            #print error message for shapes other than spheres 
            #(to be removed when other volume calculations are implemented)
            
            
    def calc_sphericity(self):
        if self.shape == "sphere":
            self.sphericity = 1
        elif self.shape == "fibre" or self.shape == "fiber" or self.shape == "cylinder":
            self.sphericity == 0
        elif self.shape == "pellet" or self.shape == "fragment":
            self.sphericity ==
            
    def calc_reynolds_number(self):
        
        if self.shape == "sphere":
            #still need to implement speed_m2_s
            self.reynolds_num = (density_w_21C_kg_m3*diameter_m*speed_m2_s)/mu_w_21C_kg_ms
        
        else:
            #still need to implement speed_m2_s
            #using length_c for reduced impact in the results
            self.reynolds_num = (density_w_21C_kg_m3*length_c_m*speed_m2_s)/mu_w_21C_kg_ms
            
    def calc_drag_coef(self):
        
        if self.shape == "sphere":
            if self.reynolds_num < 1:
                self.drag_coef = 24/reynolds_num
            elif self.reynolds_num < 5:
                self.drag_coef = 24/reynolds_num*(1+(3/16)*reynolds_num)
            elif self.reynolds_num < 10**3:
                self.drag_coef = 24/reynolds_num*(1+0.15*reynolds_num**0.687)
            else:
                self.drag_coef = 24/reynolds_num*(0.0183*reynolds_num)
        else:
            print("Error: shape applicable for this model yet")
        
    
    #degradation estimations
    """ relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    
    

    