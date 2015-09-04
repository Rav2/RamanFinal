import Tkinter as tk
from Tkinter import *
import tkMessageBox as mb
import numpy as np
import aggregate_data as agg
from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib import cm
import fit_old as fo
import nanograin_size as ns
import load_files
import __builtin__

default_min_size = "1"
default_max_size = "10"
###############################################################################
#########################DEFINITIONS###########################################   
def map_maker():
    #LOADING DATA AND AGGREGATION
    data = load_files.load_mapping_file(file_name.get())
    aggregated_data = agg.agg_data(data, float(min_size.get()), float(max_size.get()), float(size_step.get()),
                                   float(omega_0.get()), float(gamma_0.get()), float(amplitude.get()), float(offset.get()))
    data_array = np.array(aggregated_data)
    a = f.add_subplot(1, 1, 1)
    min_value = min(__builtin__.map(min, data_array))
    max_value = max(__builtin__.map(max, data_array))
    if min_value == max_value:
        min_value -= 0.1
        max_value += 0.1
    #MAP MAKING
    b = a.pcolor(data_array, cmap='jet', vmin=min_value,
                 vmax=max_value)  #inne fajne cmap'y do wyboru: 'seismic' , 'coolwarm'
    f.colorbar(b)
    canvas = FigureCanvasTkAgg(f, master=GUI)
    canvas.get_tk_widget().place(x=0, y=0)
    
def single_fitting():
    #LOADING DATA AND AGGREGATION
    data = load_files.load_two_column_file(file_name.get())
    omega = data[2]
    intensity = data[3]
    result = ns.find_grain_diameter(omega, intensity, float(min_size.get()), float(max_size.get()), float(size_step.get()),
                                     float(omega_0.get()), float(gamma_0.get()), float(amplitude.get()), float(offset.get()))
    print type(result[1])
    g = result[1]
    canvas = FigureCanvasTkAgg(g, master=GUI)
    canvas.get_tk_widget().place(x=0, y=0)
    canvas.get_tk_widget().configure(border=0)    
    
def calibrate():
    #LOADING DATA AND AGGREGATION
    data = load_files.load_two_column_file(file_name.get())
    omega = data[2]
    intensity = data[3]
    result = fo.perform_fitting(omega, intensity, float(min_size.get()), float(max_size.get()))
    
    omega_0.delete(0, END)
    omega_0.insert(0,result[0][0])
    omega_0.config(foreground="green")
    gamma_0.delete(0, END)
    gamma_0.insert(0,result[0][1])
    gamma_0.config(foreground="green")
    amplitude.delete(0, END)
    amplitude.insert(0,result[0][2])
    amplitude.config(foreground="green")
    g = result[1]
    canvas = FigureCanvasTkAgg(g, master=GUI)
    canvas.get_tk_widget().place(x=0, y=0)
    canvas.get_tk_widget().configure(border=0) 
###############################################################################
########################GUI ELEMENTS###########################################   
GUI = tk.Tk()

GUI.geometry('900x505')
GUI.title("RamanMaps")

f = Figure(facecolor='white')
canvas = FigureCanvasTkAgg(f, master=GUI)
canvas.get_tk_widget().place(x=0, y=0)
canvas.get_tk_widget().configure(border=0)
notification_label = tk.Label()
notification_label.config(text="Please note, that mapping may take even a few hours, depending on number of measurements.")
notification_label.place(x=10,y=484)

file_name_label = tk.Label(text='File name:', font="Helvetica 10 bold")
file_name_label.place(x=680, y=10)
file_name = tk.Entry()
file_name.place(x=680, y=30)
file_name.insert(0,"map_file.txt")
#file_name.insert(0, "map_file.txt")

map_btn = tk.Button(text='Create map', command=map_maker, width=18)
map_btn.place(x=680, y=450)

def radio_button_select(event):
    if(event.widget["value"] ==1 ):
        map_btn.config(text="Create map!", command = map_maker)
        
        max_size_label.config(text="Max grain size to be checked [nm]:")
        max_size.delete(0, END)
        max_size.insert(0, 10)
        min_size.config(text="Min grain size to be checked [nm]:")
        min_size.delete(0, END)
        min_size.insert(0, 1)
        
        max_size.config(state=NORMAL)
        min_size.config(state=NORMAL)
        size_step.config(state=NORMAL)
        omega_0.config(state=NORMAL)
        gamma_0.config(state=NORMAL)
        amplitude.config(state=NORMAL)
        offset.config(state=NORMAL)
        notification_label.config(text="Please note, that mapping may take even a few hours, depending on number of measurements.")
    elif(event.widget["value"] == 2):
        map_btn.config(text="Perform fitting!", command = single_fitting)
        
        max_size_label.config(text="Max grain size to be checked [nm]:")
        max_size.delete(0, END)
        max_size.insert(0, 10)
        min_size.config(text="Min grain size to be checked [nm]:")
        min_size.delete(0, END)
        min_size.insert(0, 1)
        
        max_size.config(state=NORMAL)
        min_size.config(state=NORMAL)
        size_step.config(state=NORMAL)
        omega_0.config(state=NORMAL)
        gamma_0.config(state=NORMAL)
        amplitude.config(state=NORMAL)
        offset.config(state=NORMAL)
        notification_label.config(text="")
    elif(event.widget["value"] == 3):
        map_btn.config(text="Calibrate!", command = calibrate)
        
        max_size_label.config(text="Max Raman shift")
        max_size.delete(0, END)
        max_size.insert(0,"600")
        min_size_label.config(text="Min Raman shift")
        min_size.delete(0, END)
        min_size.insert(0,"450")
        
        size_step.config(state="readonly")
        omega_0.config(state="readonly")
        gamma_0.config(state="readonly")
        amplitude.config(state="readonly")
        offset.config(state="readonly")
        notification_label.config(text="")

v = tk.IntVar()
rbtn_one=tk.Radiobutton(GUI,
            text="Mapping file (four columns)",
            padx = 0, 
            variable=v, 
            value=1)
rbtn_one.place(x=670, y=60)
rbtn_one.select()
rbtn_one.bind("<Button-1>", radio_button_select)
rbtn_two=tk.Radiobutton(GUI,
            text="Single sample file (two columns)",
            padx = 0, 
            variable=v, 
            value=2)
rbtn_two.place(x=670, y=80)
rbtn_two.bind("<Button-1>", radio_button_select)
rbtn_three=tk.Radiobutton(GUI,
            text="Calibration file (two columns)",
            padx = 0, 
            variable=v,
            value=3)
rbtn_three.place(x=670, y=100)
rbtn_three.bind("<Button-1>", radio_button_select)


#initial_shift_label = tk.Label(text='Initial Raman Shift [cm^-1]:')
#initial_shift_label.place(x=680, y=260)
#initial_shift = tk.Entry()
#initial_shift.place(x=680, y=280)
#initial_shift.insert(0, "400")
#
#final_shift_label = tk.Label(text='Final Raman Shift [cm^-1]:')
#final_shift_label.place(x=680, y=210)
#final_shift = tk.Entry()
#final_shift.place(x=680, y=230)
#final_shift.insert(0, "600")

min_size_label = tk.Label(text="Min grain size to be checked [nm]:")
min_size_label.place(x=680, y=130)
min_size = tk.Entry()
min_size.place(x=680, y=145)
min_size.insert(0, 1)

max_size_label = tk.Label(text="Max grain size to be checked [nm]:")
max_size_label.place(x=680, y=170)
max_size = tk.Entry()
max_size.place(x=680, y=185)
max_size.insert(0, 10)

size_step_label = tk.Label(text="Step when checking size [nm]:")
size_step_label.place(x=680, y=210)
size_step = tk.Entry()
size_step.place(x=680, y=225)
size_step.insert(0, 0.5)

parameters_label = tk.Label(text="Parameters of main peak:", font= "Helvetica 10 bold")
parameters_label.place(x=680, y=255)

omega_0_label=tk.Label(text="omega [cm^-1]")
omega_0_label.place(x=680,y=275)
omega_0 = tk.Entry()
omega_0.place(x=680, y=290)
omega_0.insert(0, "521.604599812")

gamma_0_label=tk.Label(text="gamma [cm^-1]")
gamma_0_label.place(x=680,y=315)
gamma_0 = tk.Entry()
gamma_0.place(x=680, y=335)
gamma_0.insert(0, "4.42010161931")

amplitude_label=tk.Label(text="ampl [arb. u]")
amplitude_label.place(x=680,y=360)
amplitude = tk.Entry()
amplitude.place(x=680, y=375)
amplitude.insert(0, "71535.2586963")

offset_label=tk.Label(text="offset [arb. u]")
offset_label.place(x=680,y=400)
offset = tk.Entry()
offset.place(x=680, y=415)
offset.insert(0, "0.0")

GUI.mainloop()


