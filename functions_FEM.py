 # -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 18:38:46 2022

@author: Predator
"""

import numpy as np
import math

global problem_dimension
problem_dimension = 2

def assign_BCs (node_list, BC_flag, deformation_value):
    
    number_of_nodes = np.size(node_list,0)
    
    extended_node_list = np.zeros([number_of_nodes, 6 * problem_dimension])
    extended_node_list[:, 0:problem_dimension] = node_list
    
    if BC_flag == 'extension':
        
        for i in range(0, number_of_nodes):
            
            if  extended_node_list[i,0] == 0:
                
                #implementowanie warunku brzegowego Dirichleta
                extended_node_list[i,2] = -1 
                extended_node_list[i,3] = -1
                
                extended_node_list[i,8] = -deformation_value
                extended_node_list[i,9] = 0
                #8 oznacza kolumne zawierajaca odksztalcenia w kierunku x
                #9 oznacza kolumne zawuerająca odksztalcenia w kierunku y
            
            elif extended_node_list[i,0] == 1:
                
                extended_node_list[i,2] = -1
                extended_node_list[i,3] = -1
                
                extended_node_list[i,8] = deformation_value
                extended_node_list[i,9] = 0
                
            else:
                
                #warunek brzegowy Neumanna
                extended_node_list[i,2] = 1
                extended_node_list[i,3] = 1
                
                extended_node_list[i,10] = 0
                extended_node_list[i,11] = 0
                #siły są wyzerowane
            
    if BC_flag == 'expansion':
                
                for i in range(0, number_of_nodes):
                    
                    if  extended_node_list[i,0] == 0 or extended_node_list[i,0] == 1 or extended_node_list[i,1] == 0 or extended_node_list[i,1] == 1:
                        
                      
                        extended_node_list[i,2] = -1  
                        extended_node_list[i,3] = -1
                        
                        extended_node_list[i,8] = deformation_value*extended_node_list[i,0] 
                        extended_node_list[i,9] = deformation_value*extended_node_list[i,1]
                    else:
                        extended_node_list[i,2] = 1
                        extended_node_list[i,3] = 1
                        
                        extended_node_list[i,10] = 0
                        extended_node_list[i,11] = 0
                        
                        
    if BC_flag == 'shear':
                
                for i in range(0, number_of_nodes):
                    
                    if extended_node_list[i,1] == 0:
                        
                        extended_node_list[i,2] = -1
                        extended_node_list[i,3] = -1
                        
                        extended_node_list[i,8] = 0
                        extended_node_list[i,9] = 0
                    
                    elif extended_node_list[i,1] == 1:
                        
                        extended_node_list[i,2] = -1
                        extended_node_list[i,3] = -1
                        
                        extended_node_list[i,8] = deformation_value
                        extended_node_list[i,9] = 0
                    else:
                         
                        extended_node_list[i,2] = 1
                        extended_node_list[i,3] = 1
                        
                        extended_node_list[i,10] = 0
                        extended_node_list[i,11] = 0
                        
    DOFs = 0
    DOCs = 0
            
    for i in range(0,number_of_nodes):
        for j in range (0,problem_dimension):
            if extended_node_list[i,problem_dimension+j] == -1:
                        DOCs -= 1
                        extended_node_list[i,2*problem_dimension+j] = DOCs
            else:
                        DOFs += 1
                        extended_node_list[i,2*problem_dimension+j] = DOFs
                        
    for i in range(0,number_of_nodes):
        for j in range(0, problem_dimension):
            if extended_node_list[i,2*problem_dimension+j]<0:
                        extended_node_list[i,3*problem_dimension+j] = abs(extended_node_list[i,2*problem_dimension+j])+DOFs
            else:
                        extended_node_list[i,3*problem_dimension+j] = abs(extended_node_list[i,2*problem_dimension+j])
    DOCs = abs(DOCs)
    return(extended_node_list, DOFs, DOCs)
def assemble_displacments(extended_node_list, node_list):
    number_of_nodes = np.size(node_list,0)
    
    DOC = 0
    Up = []
    
    for i in range(0,number_of_nodes):
        for j in range(0,problem_dimension):
            if extended_node_list[i,problem_dimension+j] == -1:
                DOC +=1
                Up = np.append(Up,extended_node_list[i,4*problem_dimension+j])
            Up = np.vstack([Up]).reshape(-1,1)
    return Up

def assemble_forces(extended_node_list, node_list):
    number_of_nodes = np.size(node_list,0)
    
    degress_of_freedom = 0
    Fp = []
    
    for i in range(0,number_of_nodes):
        for j in range(0,problem_dimension):
            if extended_node_list[i, problem_dimension+j] == 1:
                degress_of_freedom +=1
                Fp = np.append(Fp,extended_node_list[i,5*problem_dimension+j])
                #Fp.append(extended_node_list[i,5*problem_dimension+j])
                Fp = np.vstack([Fp]).reshape(-1,1)
    return Fp

def update_nodes(extended_node_list, Unknown_displacments,Unknown_forces,node_list):
    number_of_nodes = np.size(node_list,0) 
    degress_of_freedom = 0
    DOCs = 0
    
    for i in range(0, number_of_nodes):
        for j in range(0, problem_dimension):
            if extended_node_list[i,problem_dimension+j] == 1:
                degress_of_freedom +=1
                extended_node_list[i,4*problem_dimension+j] = Unknown_displacments[degress_of_freedom-1]
            else:
                DOCs +=1
                extended_node_list[i,5*problem_dimension+j] = Unknown_forces[DOCs-1]
                
    return extended_node_list
    
def assemble_stiffness(extended_node_list, element_list, node_list):
    
    number_of_elements = np.size(element_list,0)
    nodes_per_element = np.size(element_list,1)
    number_of_nodes = np.size(node_list,0)
    
    stiffness_matrix = np.zeros([number_of_nodes*problem_dimension,number_of_nodes*problem_dimension,])
    
    for i in range(1, number_of_elements+1):
        nl = element_list[i-1,0:nodes_per_element]
        k = element_stiffness(nl,node_list)
         
        for r in range(0, nodes_per_element):
            for p in range(0, problem_dimension):
                for q in range(0, nodes_per_element):
                    for s in range(0, problem_dimension):
                        row = extended_node_list[nl[r]-1,p+3*problem_dimension]
                        column = extended_node_list[nl[q]-1,s+3*problem_dimension]
                        value = k[r*problem_dimension+p,q*problem_dimension+s]
                        stiffness_matrix[int(row)-1,int(column)-1] = stiffness_matrix[int(row)-1,int(column)-1] + value
    return stiffness_matrix
            
def element_stiffness(nl,node_list):
    
    nodes_per_element = np.size(nl,0)
    
    x = np.zeros([nodes_per_element, problem_dimension])
    x[0:nodes_per_element, 0:problem_dimension] = node_list[nl[0:nodes_per_element]-1,0:problem_dimension]
    
    stiffness_matrix = np.zeros([nodes_per_element*problem_dimension, nodes_per_element*problem_dimension])
    
    coordinates = x.T #tranzspozycja macierzy współrzędnych
    
    gausianpoint_per_element = 4 
    
    for i in range(1,nodes_per_element+1):
        for j in range(1, nodes_per_element+1):
               
            small_matrix = np.zeros([problem_dimension,problem_dimension])
            
            for gausian_point in range(1, gausianpoint_per_element+1):

                jacobian = np.zeros([problem_dimension,problem_dimension])
                  
                gradient = np.zeros([problem_dimension,nodes_per_element])
                
                (xi,eta,alpha) = GaussPoint(nodes_per_element, gausianpoint_per_element, gausian_point)
                
                grad_nat = Grad_N_nat(nodes_per_element,xi, eta)
                
                #mnozenie macierzowe macierzy koordynatow nodow i transponowanej funkjci ksztaltu
                jacobian = coordinates @ grad_nat.T
                
                gradient = np.linalg.inv(jacobian).T @ grad_nat
                
                for a in range(1,problem_dimension+1):
                    for c in range(1,problem_dimension+1):
                        for b in range(1,problem_dimension+1):
                            for d in range(1,problem_dimension+1):
                                
                                #linanlg.det() geeneruje wyzancznik macierzy
                                small_matrix[a-1,c-1] = small_matrix[a-1,c-1] + gradient[b-1,i-1] * constitutive(a, b, c, d) * gradient[d-1,j-1] * np.linalg.det(jacobian) * alpha 
                                                                                      
                                stiffness_matrix[((i-1)*problem_dimension+1)-1:i*problem_dimension , ((j-1)*problem_dimension+1)-1:j*problem_dimension] = small_matrix
    return stiffness_matrix

def GaussPoint (nodes_per_element, gausianpoint_per_element, gausian_point):
    
    if nodes_per_element == 4:
        
        if gausianpoint_per_element == 1:
            
            if gausian_point == 1:
                
                xi = 0
                eta = 0
                alpha = 4
        
        if gausianpoint_per_element == 4:
            
            if gausian_point == 1:
                
                xi = -1/math.sqrt(3)
                eta = -1/math.sqrt(3)
                alpha = 1
            
            if gausian_point == 2:
                
                xi = 1/math.sqrt(3)
                eta = -1/math.sqrt(3)
                alpha = 1
                
            if gausian_point == 3:
                
                xi = 1/math.sqrt(3)
                eta = 1/math.sqrt(3)
                alpha = 1
            
            if gausian_point == 4:
                
                xi = -1/math.sqrt(3)
                eta = 1/math.sqrt(3)
                alpha = 1
                
    return(xi,eta,alpha)

#funkcja odpowiedzialna za gradient funkcji kształtu
def Grad_N_nat(nodes_per_element,xi,eta):
    
    result = np.zeros([problem_dimension,nodes_per_element])
    

    result[0,0] = -1/4*(1-eta)
    result[0,1] = 1/4*(1-eta)
    result[0,2] = 1/4*(1+eta)
    result[0,3] = -1/4*(1+eta)
        
    result[1,0] = -1/4*(1-xi)
    result[1,1] = -1/4*(1+xi)
    result[1,2] = 1/4*(1+xi)
    result[1,3] = 1/4*(1-xi )
    
    return result

def constitutive(i,j,k,l):
    
    young_modulus = 8/3
    poisson_ration = 1/3
    
    c = (young_modulus/(2*(1+poisson_ration))) * (delta(i,l) * delta(j,k) + delta(i,k) * delta(j,l)) + (young_modulus*poisson_ration)/(1 - poisson_ration**2)  * delta(i,j) * delta(k,l)
                             
                            
    
    return c

def delta(i,j):
    if i == j:
        delta = 1
    else:
        delta = 0
    return delta
    
    
    
