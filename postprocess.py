# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 17:03:33 2022

@author: tytus
"""

import numpy as np
import math
from functions_FEM import *

global problem_dimension
problem_dimension = 2



def post_process (node_list, element_list, extended_node_list,scale):
    number_of_elements = np.size(element_list,0)
    nodes_per_element = 4
    
    
    disp,stress,strain = element_post_process(node_list, element_list, extended_node_list)
    
    stress_xx = np.zeros([nodes_per_element,number_of_elements])
    stress_xy = np.zeros([nodes_per_element,number_of_elements])
    stress_yx = np.zeros([nodes_per_element,number_of_elements])
    stress_yy = np.zeros([nodes_per_element,number_of_elements])
    
    strain_xx = np.zeros([nodes_per_element,number_of_elements])
    strain_xy = np.zeros([nodes_per_element,number_of_elements])
    strain_yx = np.zeros([nodes_per_element,number_of_elements])
    strain_yy = np.zeros([nodes_per_element,number_of_elements])
    
    disp_x = np.zeros([nodes_per_element,number_of_elements])
    disp_y = np.zeros([nodes_per_element,number_of_elements])
    
    X = np.zeros([nodes_per_element,number_of_elements])
    Y = np.zeros([nodes_per_element,number_of_elements])
    
    X = extended_node_list[element_list-1,0] + scale*extended_node_list[element_list-1,4*problem_dimension]
    Y = extended_node_list[element_list-1,1] + scale*extended_node_list[element_list-1,4*problem_dimension+1]
    
    X = X.T
    Y = Y.T
    
    stress_xx[:,:] = stress[:,:,0,0].T
    stress_xy[:,:] = stress[:,:,0,1].T
    stress_yx[:,:] = stress[:,:,1,0].T
    stress_yy[:,:] = stress[:,:,1,1].T
    
    strain_xx[:,:] = strain[:,:,0,0].T
    strain_xy[:,:] = strain[:,:,0,1].T
    strain_yx[:,:] = strain[:,:,1,0].T
    strain_yy[:,:] = strain[:,:,1,1].T
    
    disp_x = disp[:,:,0,0].T
    disp_y = disp[:,:,1,0].T
    
    disp_total = np.zeros([nodes_per_element, number_of_elements])

    for row in range(0,nodes_per_element):
        for column in range(0,number_of_elements):
            disp_total[row,column] = math.sqrt(disp_x[row,column]**2+disp_y[row,column]**2)
    
    return (stress_xx,stress_xy,stress_yx,stress_yy,disp_x,disp_y,X,Y,disp_total)
    
def element_post_process (node_list,element_list,extended_node_list):
    number_of_elements = np.size(element_list,0)
    nodes_per_element = 4
    
    gausianpoint_per_element = 4 
    
    disp = np.zeros([number_of_elements,nodes_per_element,problem_dimension,1])
    
    stress = np.zeros([number_of_elements,gausianpoint_per_element,problem_dimension,problem_dimension])
    strain = np.zeros([number_of_elements,gausianpoint_per_element,problem_dimension,problem_dimension])
    
    for e in range(1,number_of_elements+1):
        mini_node_list = element_list[e-1,0:nodes_per_element]
        
        #przypisanie przemieszczeń do poszczególnych nodów
        for i in range(1,nodes_per_element+1):
            for j in range(1,problem_dimension+1):
                disp[e-1,i-1,j-1,0] = extended_node_list[mini_node_list[i-1]-1,4*problem_dimension+j-1]
        
        x = np.zeros([nodes_per_element,problem_dimension])
        x[0:nodes_per_element,0:problem_dimension] = node_list[mini_node_list[0:nodes_per_element]-1,0:problem_dimension]
        
        u = np.zeros([problem_dimension,nodes_per_element])
        for i in range(1,nodes_per_element+1):
            for j in range(1,problem_dimension+1):
                u[j-1,i-1] = extended_node_list[mini_node_list[i-1]-1,4*problem_dimension+j-1]
        
        coor=x.T
        
        for gausianpoint_point in range(1,gausianpoint_per_element+1):
            
            epsilon = np.zeros([problem_dimension,problem_dimension])
            
            for i in range(1,nodes_per_element+1):
                
                jacobian = np.zeros([problem_dimension,problem_dimension])
                
                gradient = np.zeros([problem_dimension,nodes_per_element])
                
                (xi,eta,alpha) = GaussPoint(nodes_per_element,gausianpoint_per_element,gausianpoint_point)
                
                grad_nat = Grad_N_nat(nodes_per_element,xi,eta) 
                
                jacobian = coor @ grad_nat.T
                
                gradient = np.linalg.inv(jacobian).T @ grad_nat
                
                epsilon = epsilon + 1/2 *tensor_product(gradient[:,i-1],u[:,i-1]) + tensor_product(u[:,i-1],gradient[:,i-1])
                
            sigma = np.zeros([problem_dimension,problem_dimension])
            
            for a in range(1,problem_dimension+1):
                for b in range(1, problem_dimension+1):
                    for c in range(1,problem_dimension+1):
                        for d in range(1,problem_dimension+1):
                            sigma[a-b,b-1] = sigma[a-1,b-1]+ constitutive(a, b, c, d)*epsilon[c-d,d-1]
            
            
            for a in range(1,problem_dimension+1):
                for b in range(1,problem_dimension+1):
                    strain[e-1,gausianpoint_point-1,a-1,b-1] = epsilon[a-1,b-1]
                    stress[e-1,gausianpoint_point-1,a-1,b-1] = sigma[a-1,b-1]

    return disp,stress,strain
        
def tensor_product(u,v):
    u = u.reshape(len(v),1)
    v = v.reshape(len(v),1)
    
    A = u @ v.T
    
    return A
                
        