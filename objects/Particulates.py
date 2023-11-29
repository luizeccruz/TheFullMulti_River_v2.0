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
        self.length_a_um = plastic_prop.length_a_um.loc[MP_index] #longest length (for nonspherical MPs), 
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
            
        elif self.shape == "cube":
            self.volume_m3 = self.length_a_m**3
            #calculates volume (in m3) of cubes from its length
            self.CSF = 1
            
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
            
            
    #methods used for new settling model (and for further transport phenomena updates)
    def calc_area(self):
        
        if self.shape == "sphere":
            self.area = 4*math.pi*(self.diameter_m/2)**2
        elif self.shape == "fibre" or self.shape == "fiber" or self.shape == "cylinder":
            self.area = 2*math.pi*self.diameter_m/2*self.length_a_m + 2*math.pi*(self.diameter_m/2)**2
        elif  self.shape == "pellet" or self.shape == "fragment":
            self.area = 2*self.length_a_m*self.length_b_m + 2*self.length_a_m*self.length_c_m + 2*self.length_b_m*self.length_c_m
        elif self.shape == "cube":
            self.area = 6*self.length_a_m**2
        else:
            print(f"Error: shape {self.shape} not implemented")
    
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
            
        elif self.shape == "cube":
            self.projected_area_m2 = self.length_a_m**2
            
        else:
            print("Error: unknown shape")
            #print error message for shapes other than spheres 
            #(to be removed when other volume calculations are implemented)
            
            
            
    def calc_sphericity(self, sphericity):
        
        if sphericity == "calculated": 
            sphericity = (math.pi**(1/3)*(6*self.volume_m3)**(2/3))/self.area
            if (sphericity > 0.66):
                self.sphericity = sphericity
            else:
                print(f"sphericity {sphericity} not compatible with model")
        else:
            self.sphericity = float(sphericity)
            
            
    def calc_drag_coef(self):
        
        if self.shape != "sphere":
            self.sp_diameter_m = ((6*self.volume_m3)/math.pi)**(1/3)
        else:
            self.sp_diameter_m = self.diameter_m
                
        self.k1 = 0.843*math.log(self.sphericity/0.065, 10)
        self.k2 = 5.31 - 4.88*self.sphericity
        
        self.cd_re2 = (4*self.sp_diameter_m**3*density_w_21C_kg_m3*(self.density_kg_m3-density_w_21C_kg_m3)*g_m_s2)*(3*mu_w_21C_kg_ms**2)
        self.re = (((self.k1*self.cd_re2)/24)**(-1.2)+(self.cd_re2/self.k2)**(-0.6))**(-1/1.2)
        if (self.re<0.5):
            self.drag_coef = 24/(self.k1*self.re)
        elif (2*10**3<self.re<2*10**5):
            self.drag_coef = self.k2
        else:
            self.drag_coef = ((24/(self.k1*self.re))**(0.85)+self.k2**0.85)**(1/0.85)
            
            
    def calc_vSet(self):
        
        if self.shape != "sphere":
            self.sp_diameter_m = ((6*self.volume_m3)/math.pi)**(1/3)
        else:
            self.sp_diameter_m = self.diameter_m
        
        vset_m_s_0 =  (g_m_s2*(self.density_kg_m3-density_w_21C_kg_m3)*self.sp_diameter_m)/(18*mu_w_21C_kg_ms)
        Re_0 = (vset_m_s_0*self.density_kg_m3*self.sp_diameter_m)/mu_w_21C_kg_ms
        Re_corrected = self.CSF*Re_0
        
        Re = 0
        
           
        if Re_corrected < 1:
            self.re = Re_corrected
            self.vset = vset_m_s_0
        else:
            while Re_corrected<10**4 and Re_corrected>=1:
                if Re_corrected <= 25:
                    cd = 24/Re_corrected+3/(Re_corrected**(1/2))+0.34
                else:
                    cd = 24/Re_corrected + 3/(Re_corrected**(1/2))+0.92
                    
                vset_m_s = ((4*g_m_s2*(self.density_kg_m3-density_w_21C_kg_m3)*self.sp_diameter_m)/(3*cd*density_w_21C_kg_m3))**(1/2)
                Re_corrected = (vset_m_s*self.density_kg_m3*self.sp_diameter_m)/mu_w_21C_kg_ms
            self.vset = vset_m_s
        self.re = Re_corrected

        
    
    #degradation estimations
    """ relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    
    

    