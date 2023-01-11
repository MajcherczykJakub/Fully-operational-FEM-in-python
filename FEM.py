# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 18:26:33 2022

@author: Predator
"""

import numpy as np
from void_mesh import *
from matplotlib.figure import Figure
import math
import matplotlib.pyplot as plt
from functions_FEM import *
from postprocess import *
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)

from tkinter import *
from PIL import ImageTk, Image
from tkinter import messagebox
from functools import partial

root = Tk()
root.title("2d mechanics of materials FEM analysis")
root.iconbitmap("C:/Users/tytus/Desktop/Projekt_ZPO/FEM_icon.ico")

frame_mesh_picture = LabelFrame(root,text="Scheme of mesh")
frame_mesh_picture.grid(row=0,column=0,columnspan=3)

frame_mesh_settings = LabelFrame(root, text="mesh settings")
frame_mesh_settings.grid(row=1,column=0)

frame_mesh_gen = LabelFrame(root,text="Generated Mesh")
frame_mesh_gen.grid(row=1,column=1)

analysis_settings = LabelFrame(root,text="Analysis Settings")
analysis_settings.grid(row=1,column=2)

analysis_output = LabelFrame(root,text="Analysis Plots")
analysis_output.grid(row=0,column=3)

postprocess_settings = LabelFrame(root,text="Post process settings")
postprocess_settings.grid(row=1,column=3)


canvas = Canvas(frame_mesh_gen)
canvas.pack()

analysis_output = Canvas(analysis_output)
analysis_output.pack()

bc_flag_text = Label(analysis_settings,text="border conditions")
bc_flag_text.pack()

bc_modes = [
      ("shear","shear"),
      ("expansion","expansion"),
      ("extension","extension"),
]

bc = StringVar()
bc.set("extension")

for text, mode in bc_modes:
    Radiobutton(analysis_settings,text=text, variable=bc, value=mode).pack()
    
scale_text=Label(postprocess_settings ,text="scale in postprocess")
scale_text.pack()
scale_entry = Entry(postprocess_settings,width=49)
scale_entry.pack()

###########

deformation_text=Label(analysis_settings, text="deformation value")
deformation_text.pack()
deformation_entry = Entry(analysis_settings,width=49)
deformation_entry.pack()

display_mode_text = Label(postprocess_settings,text="analysis results")
display_mode_text.pack()

display_modes = [
      ("stress_XX","stress_XX"),
      ("stress_XY","stress_XY"),
      ("stress_YY","stress_YY"),
      ("stress_YX","stress_YX"),
      ("disp_total","disp_total"),
]

display = StringVar()
display.set('disp_total')

for text, mode in display_modes :
    Radiobutton(postprocess_settings,text=text, variable=display, value=mode).pack()

img = ImageTk.PhotoImage(Image.open("C:/Users/tytus/Desktop/Projekt_ZPO/Rysunek2.jpg"))
pictorial_drawing = Label(frame_mesh_picture, image = img)
pictorial_drawing.pack()#side = "top", fill = "both", expand = "yes"


d1_text = Label(frame_mesh_settings,text="d1:")
d1_text.pack()
d1_entry = Entry(frame_mesh_settings,width=49)
d1_entry.pack()

d2_text = Label(frame_mesh_settings,text="d2:")
d2_text.pack()
d2_entry = Entry(frame_mesh_settings,width=49)
d2_entry.pack()

m_text = Label(frame_mesh_settings,text="m:")
m_text.pack()
m_entry = Entry(frame_mesh_settings,width=49)
m_entry.pack()

p_text = Label(frame_mesh_settings,text="p:")
p_text.pack()
p_entry = Entry(frame_mesh_settings,width=49)
p_entry.pack()

R_text = Label(frame_mesh_settings,text="R:")
R_text.pack()
R_entry = Entry(frame_mesh_settings,width=49)
R_entry.pack()

def mesh_generate():
    
    global R
    R = float(R_entry.get())
    if R<0:
        messagebox.showerror('Incorect data','R value is to low')
        R_entry.delete(0,END)
        exit()
    global d1
    d1 = int(d1_entry.get())
    if d1<0 or 3*R>=d1:
        messagebox.showerror('Incorect data','d1 value is incorect')
        d1_entry.delete(0,END)
        exit()
    global d2
    d2 = int(d2_entry.get())
    if d2<0 or 3*R>=d2:
        messagebox.showerror('Incorect data','d2 value is incorect')
        d2_entry.delete(0,END)
        exit()
    global m
    m = int(m_entry.get())
    if m<0 :
        messagebox.showerror('Incorect data','m value is incorect')
        m_entry.delete(0,END)
        exit()
    global p
    p = int(p_entry.get())
    if p<0 :
        messagebox.showerror('Incorect data','p value is incorect')
        p_entry.delete(0,END)
        exit()
    node_list, element_list = void_mesh(d1,d2,p,m,R)

    number_of_nodes = np.size(node_list,0)
    number_of_elements = np.size(element_list,0)
        
    #plt.figure
    f = Figure(figsize=(3,3), dpi=100)
    plot1 = f.add_subplot(111)

    count = 1 # zlicza index nodow
    for i in range(0, number_of_nodes):
        plot1.annotate(count, xy = (node_list[i,0], node_list[i,1]))
        count +=1
    
        count2 =1
        for j in range(0, number_of_elements):
            
            plot1.annotate(count2, xy = ((node_list[element_list[j,0]-1,0]+node_list[element_list[j,1]-1,0]+node_list[element_list[j,2]-1,0]+node_list[element_list[j,3]-1,0])/4,
                                           (node_list[element_list[j,0]-1,1]+node_list[element_list[j,2]-1,1]+node_list[element_list[j,2]-1,1]+node_list[element_list[j,3]-1,1])/4),
                                           c = 'blue')
            count2+=1
            
        #plot lines
        x0,y0 = node_list[element_list[:,0]-1,0] , node_list[element_list[:,0]-1,1]
        x1,y1 = node_list[element_list[:,1]-1,0] , node_list[element_list[:,1]-1,1]
        x2,y2 = node_list[element_list[:,2]-1,0] , node_list[element_list[:,2]-1,1]
        x3,y3 = node_list[element_list[:,3]-1,0] , node_list[element_list[:,3]-1,1]
        
        
        plot1.plot(np.array([x0,x1]),np.array([y0,y1]), 'red', linewidth =3)
        plot1.plot(np.array([x1,x2]),np.array([y1,y2]), 'red', linewidth =3)
        plot1.plot(np.array([x2,x3]),np.array([y2,y3]), 'red', linewidth =3)
        plot1.plot(np.array([x3,x0]),np.array([y3,y0]), 'red', linewidth =3)
    canvas1 = FigureCanvasTkAgg(f, master = canvas)
    canvas1.draw()    
    canvas1.get_tk_widget().pack()

generate_mesh_button = Button(frame_mesh_settings,text='generate mesh',command= mesh_generate)
generate_mesh_button.pack(padx=10,pady=10)



def calculate():
    
    global R
    R = float(R_entry.get())
    if R<0:
        messagebox.showerror('Incorect data','R value is to low')
        R_entry.delete(0,END)
        exit()
    global d1
    d1 = int(d1_entry.get())
    if d1<0 or 3*R>=d1:
        messagebox.showerror('Incorect data','d1 value is incorect')
        d1_entry.delete(0,END)
        exit()
    global d2
    d2 = int(d2_entry.get())
    if d2<0 or 3*R>=d2:
        messagebox.showerror('Incorect data','d2 value is incorect')
        d2_entry.delete(0,END)
        exit()
    global m
    m = int(m_entry.get())
    if m<0 :
        messagebox.showerror('Incorect data','m value is incorect')
        m_entry.delete(0,END)
        exit()
    global p
    p = int(p_entry.get())
    if p<0 :
        messagebox.showerror('Incorect data','p value is incorect')
        p_entry.delete(0,END)
        exit()
    global scale
    scale = int(scale_entry.get())
    if scale<0 :
        messagebox.showerror('Incorect data','scale value is incorect')
        scale_entry.delete(0,END)
        exit()
    global deformation_value
    deformation_value = float(deformation_entry.get())
    if deformation_value<0.0 :
        messagebox.showerror('Incorect data','p value is incorect')
        deformation_value.delete(0,END)
        exit()
    global bc_flag
    bc_flag=bc.get()
    
    global display
    display = display.get()
    
    node_list, element_list = void_mesh(d1,d2,p,m,R)
    
    number_of_nodes = np.size(node_list,0)
    number_of_elements = np.size(element_list,0)
    
    
    (extended_node_list, degrees_of_freedom, DOCs) = assign_BCs(node_list,bc_flag,deformation_value)
    
    K = assemble_stiffness(extended_node_list, element_list, node_list)
    
    Fp = assemble_forces(extended_node_list, node_list)
    Up = assemble_displacments(extended_node_list, node_list)
    
    k_reduced = K[0:degrees_of_freedom, 0:degrees_of_freedom]
    k_Up = K[0:degrees_of_freedom, degrees_of_freedom:DOCs+degrees_of_freedom]
    k_Pu = K[degrees_of_freedom:DOCs+degrees_of_freedom,0 :degrees_of_freedom]
    k_Pp = K[degrees_of_freedom:DOCs+degrees_of_freedom, degrees_of_freedom: DOCs + degrees_of_freedom]
    
    F = Fp - (k_Up @ Up) 
    Unknown_displacments = np.linalg.solve(k_reduced, F)
    Unknown_forces = (k_Pu @ Unknown_displacments) + (k_Pp @ Up)
    
    extended_node_list = update_nodes(extended_node_list, Unknown_displacments,Unknown_forces,node_list)
    
    
    ##################################################################################
    
    
    (stress_xx, stress_xy, stress_yx, stress_yy, disp_x,disp_y,X,Y,disp_total) = post_process(node_list,element_list,extended_node_list,scale)
        
    stress_xxNormalized = (stress_xx - stress_xx.min())/(stress_xx.max()-stress_xx.min())
    disp_total_normalized = (disp_total - disp_total.min())/(disp_total.max()-disp_total.min())
    stress_xyNormalized = (stress_xy - stress_xy.min())/(stress_xy.max()-stress_xy.min())
    stress_yyNormalized = (stress_yy - stress_yy.min())/(stress_yy.max()-stress_yy.min())
    stress_yxNormalized = (stress_yx - stress_yx.min())/(stress_yx.max()-stress_yx.min())
    
    def truncate_colormap(cmap, minval = 0.0, maxval = 1.0, n =-1):
        if n== -1:
            n = cmap.N
        new_cmap = mcolors.LinearSegmentedColormap.from_list(
            'trunc({name},{a:2f},{b:2f})'.format(name = cmap.name, a = minval, b = maxval),
            cmap(np.linspace(minval,maxval, n)))
        
        return new_cmap
    
    if display=='stress_XX': 
        fig_2 = Figure(figsize=(3,3), dpi=100)
        fig_2 = plt.figure(2)
        plt.title('stress XX')
        axial_stress_xx = fig_2.add_subplot(111)
        
        for i in range(np.size(element_list,0)):
            x = X[:,i]
            y = Y[:,i]
            c = stress_xxNormalized[:,i]
            
            cmap = truncate_colormap(plt.get_cmap('jet'),c.min(), c.max())
            
            t = axial_stress_xx.tripcolor(x,y,c , cmap = cmap, shading = 'gouraud')
            
            p = axial_stress_xx.plot(x,y,'-k',linewidth = 0.5)
            

    elif display=='disp_total':
        fig_2 = Figure(figsize=(3,3), dpi=100)
        fig_2 = plt.figure(2)
        plt.title('total displacement')
        axial_stress_xx = fig_2.add_subplot(111)
        
        for i in range(np.size(element_list,0)):
            x = X[:,i]
            y = Y[:,i]
            c = disp_total_normalized[:,i]
            
            cmap = truncate_colormap(plt.get_cmap('jet'),c.min(), c.max())
            
            t = axial_stress_xx.tripcolor(x,y,c , cmap = cmap, shading = 'gouraud')
            
            p = axial_stress_xx.plot(x,y,'-k',linewidth = 0.5)
            
    
    elif display=='stress_XY':
        fig_2 = Figure(figsize=(3,3), dpi=100)
        fig_2 = plt.figure(2)
        plt.title('stress_XY')
        shear_stress_xy = fig_2.add_subplot(111)
        
        for i in range(np.size(element_list,0)):
            x = X[:,i]
            y = Y[:,i]
            c = stress_xyNormalized[:,i]
            
            cmap = truncate_colormap(plt.get_cmap('jet'),c.min(), c.max())
            
            t = shear_stress_xy.tripcolor(x,y,c , cmap = cmap, shading = 'gouraud')
            
            p = shear_stress_xy.plot(x,y,'-k',linewidth = 0.5)
            
    elif display=='stress_YY':
        fig_2 = Figure(figsize=(3,3), dpi=100)
        fig_2 = plt.figure(2)
        plt.title('stress_YY')
        axial_stress_yy = fig_2.add_subplot(111)
        
        for i in range(np.size(element_list,0)):
            x = X[:,i]
            y = Y[:,i]
            c = stress_yyNormalized[:,i]
            
            cmap = truncate_colormap(plt.get_cmap('jet'),c.min(), c.max())
            
            t = axial_stress_yy.tripcolor(x,y,c , cmap = cmap, shading = 'gouraud')
            
            p = axial_stress_yy.plot(x,y,'-k',linewidth = 0.5)
            
    elif display=='stress_YX':
        fig_2 = Figure(figsize=(3,3), dpi=100)
        fig_2 = plt.figure(2)
        plt.title('stress_YX')
        shear_stress_yx = fig_2.add_subplot(111)
        
        for i in range(np.size(element_list,0)):
            x = X[:,i]
            y = Y[:,i]
            c = stress_yxNormalized[:,i]
            
            cmap = truncate_colormap(plt.get_cmap('jet'),c.min(), c.max())
            
            t = shear_stress_yx.tripcolor(x,y,c , cmap = cmap, shading = 'gouraud')
            
            p = shear_stress_yx.plot(x,y,'-k',linewidth = 0.5)
            
    canvas1 = FigureCanvasTkAgg(fig_2, master = analysis_output)
    canvas1.draw()    
    canvas1.get_tk_widget().pack()

calculate_button = Button(analysis_settings,text='calculate',command= calculate)
calculate_button.pack() 



root.mainloop()


    