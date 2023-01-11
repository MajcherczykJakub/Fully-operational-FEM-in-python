# -*- coding: utf-8 -*-
# cos jes skurwione jak chodzi o element list
"""
Created on Tue Nov  1 11:59:31 2022

@author: tytus
"""

import numpy as np
import math

global problem_dimension
problem_dimension = 2

def void_mesh (d1,d2,p,m,R):
    
    q =np.array([[0,0], [d1,0], [0,d2], [d1,d2]])
    
    
    number_of_nodes = 2*(p+1)*(m+1)+2*(p-1)*(m+1)
    number_of_elements = m*p*4
    nodes_per_element = 4
    
    ## Nodes
    
    node_list = np.zeros([number_of_nodes, problem_dimension])
    #0 to x, kierunek horyzontalny
    #1 to y, kierunek vertykalny  
    
    increment_in_horizontal_direction = (q[1,0]-q[0,0])/p
    increment_in_vertical_direction = (q[2,1]-q[0,1])/p
    
    #stworzenie listy nodów dla poszczegolnych regionow i potem ich sklejanie
    ##lista dla regionu 1
    coordinates11 = np.zeros([(p+1)*(m+1), problem_dimension])
    
    #dla nodów położonych na granicy plyty
    for column in range(1,p+2):
        coordinates11[column-1,0] = q[0,0] + (column-1)*increment_in_horizontal_direction
        coordinates11[column-1,1] = q[0,1]
    
    #dla nodów na granicy z otworem
    for column in range(1,p+2):
        coordinates11[m*(p+1)+column-1,0] = R*np.cos((5*math.pi/4) + (column-1)*((math.pi/2)/p)) +d1/2
        coordinates11[m*(p+1)+column-1,1] = R*np.sin((5*math.pi/4) + (column-1)*((math.pi/2)/p)) +d2/2
        
    #dla nodów pomiędzy
    for row in range(1,m):
        for column in range(1,p+2):
            
            increment_tau_horizontal = (coordinates11[m*(p+1)+column-1,0] - coordinates11[column-1,0] )/m
            increment_tau_vertical = (coordinates11[m*(p+1)+column-1,1] - coordinates11[column-1,1])/m
            
            coordinates11[row*(p+1)+column-1,0] = coordinates11[(row-1)*(p+1)+column-1,0] + increment_tau_horizontal
            coordinates11[row*(p+1)+column-1,1] = coordinates11[(row-1)*(p+1)+column-1,1] + increment_tau_vertical
    
    ## lista dla regionu 2
    
    coordinates22 = np.zeros([(p+1)*(m+1),problem_dimension])
    
    #2 w pierwszej wspolrzednej q to gorna granica plyty
    #dla granicy gornej plyty
    for column in range(1,p+2):
        coordinates22[column-1,0] = q[2,0] + (column-1)*increment_in_horizontal_direction
        coordinates22[column-1,1] = q[2,1]
    
    #dla nodów na granicy z otworem
    for column in range(1,p+2):
        coordinates22[m*(p+1)+column-1,0] = R*np.cos((3*math.pi/4) - (column-1)*((math.pi/2)/p)) +d1/2
        coordinates22[m*(p+1)+column-1,1] = R*np.sin((3*math.pi/4) - (column-1)*((math.pi/2)/p)) +d2/2
    
    #dla nodów pomiędzy
    for row in range(1,m):
        for column in range(1,p+2):
            
            increment_tau_horizontal = (coordinates22[m*(p+1)+column-1,0] - coordinates22[column-1,0] )/m
            increment_tau_vertical = (coordinates22[m*(p+1)+column-1,1] - coordinates22[column-1,1] )/m
            
            coordinates22[row*(p+1)+column-1,0] = coordinates22[(row-1)*(p+1)+column-1,0] + increment_tau_horizontal
            coordinates22[row*(p+1)+column-1,1] = coordinates22[(row-1)*(p+1)+column-1,1] + increment_tau_vertical
    
    ## lista dla regionu 3
    
    coordinates33 = np.zeros([(p-1)*(m+1),problem_dimension])
    
    #dla granicy gornej plyty
    for column in range(1,p):
        coordinates33[column-1,0] = q[0,0] 
        coordinates33[column-1,1] = q[0,1] + column*increment_in_vertical_direction
    
    for column in range(1,p):
        coordinates33[m*(p-1)+column-1,0] = R*np.cos((5*math.pi/4) - (column)*((math.pi/2)/p)) +d1/2
        coordinates33[m*(p-1)+column-1,1] = R*np.sin((5*math.pi/4) - (column)*((math.pi/2)/p)) +d2/2
    
    #dla nodów pomiędzy
    for row in range(1,m):
        for column in range(1,p):
            
            increment_tau_horizontal = (coordinates33[m*(p-1)+column-1,0] - coordinates33[column-1,0])/m
            increment_tau_vertical = (coordinates33[m*(p-1)+column-1,1] - coordinates33[column-1,1])/m
            
            coordinates33[row*(p-1)+column-1,0] = coordinates33[(row-1)*(p-1)+column-1,0] + increment_tau_horizontal
            coordinates33[row*(p-1)+column-1,1] = coordinates33[(row-1)*(p-1)+column-1,1] + increment_tau_vertical
    
    ## lista dla regionu 4
    
    coordinates44 = np.zeros([(p-1)*(m+1),problem_dimension])
    
    #dla granicy gornej plyty
    for column in range(1,p):
        coordinates44[column-1,0] = q[1,0] 
        coordinates44[column-1,1] = q[1,1] + column*increment_in_vertical_direction
    
     
    for column in range(1,p):
        coordinates44[m*(p-1)+column-1,0] = R*np.cos((7*math.pi/4) + (column)*((math.pi/2)/p)) +d1/2
        coordinates44[m*(p-1)+column-1,1] = R*np.sin((7*math.pi/4) + (column)*((math.pi/2)/p)) +d2/2
    
    for row in range(1,m):
        for column in range(1,p):
            
            increment_tau_horizontal = (coordinates44[m*(p-1)+column-1,0] - coordinates44[column-1,0])/m
            increment_tau_vertical = (coordinates44[m*(p-1)+column-1,1] - coordinates44[column-1,1] )/m
            
            coordinates44[row*(p-1)+column-1,0] = coordinates44[(row-1)*(p-1)+column-1,0] + increment_tau_horizontal
            coordinates44[row*(p-1)+column-1,1] = coordinates44[(row-1)*(p-1)+column-1,1] + increment_tau_vertical
    
    
    
    ##Reordering nodes
    
    for i in range(1,m+2):
        
        node_list[(i-1)*4*p:i*4*p,:] = np.vstack([coordinates11[(i-1)*(p+1):(i)*(p+1),:],
                                                  coordinates44[(i-1)*(p-1):(i)*(p-1),:],
                                                  np.flipud(coordinates22[(i-1)*(p+1):(i)*(p+1),:]),
                                                  np.flipud(coordinates33[(i-1)*(p-1):(i)*(p-1),:])])
    
    ##Generating elements
    
    elements_list = np.zeros([number_of_elements,nodes_per_element])
    
    for row in range(1,m+1):
        for j in range(1,4*p+1):
            
            if j==1:
                
                elements_list[(row-1)*(4*p)+j-1,0] = (row-1)*(4*p)+j
                elements_list[(row-1)*(4*p)+j-1,1] = elements_list[(row-1)*(4*p)+j-1,0] + 1
                elements_list[(row-1)*(4*p)+j-1,3] = elements_list[(row-1)*(4*p)+j-1,0] + 4*p
                elements_list[(row-1)*(4*p)+j-1,2] = elements_list[(row-1)*(4*p)+j-1,3] + 1
            
            elif j== 4*p:
                
                elements_list[(row-1)*(4*p)+j-1,0] = row*(4*p)
                elements_list[(row-1)*(4*p)+j-1,1] = (row-1)*(4*p) + 1
                elements_list[(row-1)*(4*p)+j-1,2] = elements_list[(row-1)*(4*p)+j-1,0] + 1
                elements_list[(row-1)*(4*p)+j-1,3] = elements_list[(row-1)*(4*p)+j-1,0] + 4*p
            
            else:
                
                elements_list[(row-1)*(4*p)+j-1,0] = elements_list[(row-1)*(4*p)+j-2,1]
                elements_list[(row-1)*(4*p)+j-1,3] = elements_list[(row-1)*(4*p)+j-2,2]
                elements_list[(row-1)*(4*p)+j-1,2] = elements_list[(row-1)*(4*p)+j-1,3] + 1
                elements_list[(row-1)*(4*p)+j-1,1] = elements_list[(row-1)*(4*p)+j-1,0] + 1
    
    elements_list=elements_list.astype(int)
    return(node_list,elements_list)
