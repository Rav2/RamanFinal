import Tkinter as tk
from Tkinter import *
import aggregate_data as agg
from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import fit_old as fo
import nanograin_size as ns
import load_files
import __builtin__
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfilename
#import os

default_min_size = "1"
default_max_size = "10"
global first_mapping
first_mapping = 1
global datasets_for_map
datasets_for_map = []
global seleted_dataset
global map_styles
map_styles = ['Reds', 'YlOrRd', 'winter', 'rainbow']
global selected_style
global dataset_titles
dataset_titles = ['Grain diameter[nm]', 'Peak center shift[cm^-1]', 'Full width at half maximum[cm^-1]', 'Intensity[arb]']
#app_path = os.path.dirname(sys.argv[0])
###############################################################################
#########################DEFINITIONS###########################################
def on_closing():
    matplotlib.pyplot.close("all")
    root.quit()


def map_maker():
    # LOADING DATA AND AGGREGATION
    save_btn.place_forget()
    data = load_files.load_mapping_file(file_name.get())
    aggregated_data = agg.agg_data(data, float(min_size.get()), float(max_size.get())+float(size_step.get()), float(size_step.get()),
                                   float(omega_0.get()), float(gamma_0.get()), float(amplitude.get()),
                                   )

    global datasets_for_map
    for arg in datasets_for_map:
        datasets_for_map.remove(arg)
    for data_set in aggregated_data:
        data_array = np.array(data_set)
        min_value = min(__builtin__.map(min, data_set))
        max_value = max(__builtin__.map(max, data_set))
        if min_value == max_value:
            min_value -= 0.1 * min_value
            max_value += 0.1 * max_value
        datasets_for_map.append([data_array, min_value, max_value])
    global selected_dataset
    selected_dataset = 0
    
    # MAP MAKING
    global selected_style
    selected_style = 0
    global map_styles
    global dataset_titles
    global f
    f.clf()
    a = f.add_subplot(1, 1, 1)
    f.suptitle(dataset_titles[selected_dataset])
    b = a.pcolor(datasets_for_map[selected_dataset][0], cmap=map_styles[selected_style], vmin=datasets_for_map[selected_dataset][1],
                 vmax=datasets_for_map[selected_dataset][2])
    f.colorbar(b)
    global canvas
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().config(border=0)
    notification_label.config(text="")

    map_rbtn_one.place(x=10, y=484)
    map_rbtn_two.place(x=110, y=484)
    map_rbtn_three.place(x=210, y=484)
    map_rbtn_four.place(x=310, y=484)
    map_param_rbtn_one.place(x=10, y=5)
    map_param_rbtn_two.place(x=110, y=5)
    map_param_rbtn_three.place(x=210, y=5)
    map_param_rbtn_four.place(x=310, y=5)
    global first_mapping
    first_mapping = 0
    save_btn.place(x=545, y=1)


def single_fitting(are_four_col):
    # LOADING DATA AND AGGREGATION
    save_btn.place_forget()
    data = None
    if are_four_col :
        data = load_files.load_mapping_file(file_name.get())
    else:
        data = load_files.load_two_column_file(file_name.get())
    omega = data[2]
    intensity = data[3]
    result = ns.find_grain_diameter(omega, intensity, float(min_size.get()), float(max_size.get())+float(size_step.get()),
                                    float(size_step.get()),
                                    float(omega_0.get()), float(gamma_0.get()), float(amplitude.get()),
                                    True)
    global f
    f.clf()
    f = result[4]
    global canvas
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().configure(border=0, width=650, height=450)

    map_rbtn_one.place_forget()
    map_rbtn_two.place_forget()
    map_rbtn_three.place_forget()
    map_rbtn_four.place_forget()

    map_param_rbtn_one.place_forget()
    map_param_rbtn_two.place_forget()
    map_param_rbtn_three.place_forget()
    map_param_rbtn_four.place_forget()
    save_btn.place(x=545, y=1)


def calibrate():
    # LOADING DATA AND AGGREGATION
    save_btn.place_forget()
    data = load_files.load_two_column_file(file_name.get())
    omega = data[2]
    intensity = data[3]
    result = fo.perform_fitting(omega, intensity, float(min_size.get()), float(max_size.get())+float(size_step.get()))

    omega_0.delete(0, END)
    omega_0.insert(0, result[0][0])
    omega_0.config(foreground="green")
    gamma_0.delete(0, END)
    gamma_0.insert(0, result[0][1])
    gamma_0.config(foreground="green")
    amplitude.delete(0, END)
    amplitude.insert(0, result[0][2])
    amplitude.config(foreground="green")
    global f
    f.clf()
    f = result[1]
    f.suptitle("omega:{0:.4f}[cm^-1]   hwhm:{1:.4f}[cm^-1]   inten:{2:.4}[arb]".format(
        result[0][0], result[0][1], result[0][2]))
    global canvas
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().configure(border=0, width=650, height=450)

    map_rbtn_one.place_forget()
    map_rbtn_two.place_forget()
    map_rbtn_three.place_forget()
    map_rbtn_four.place_forget()

    map_param_rbtn_one.place_forget()
    map_param_rbtn_two.place_forget()
    map_param_rbtn_three.place_forget()
    map_param_rbtn_four.place_forget()
    save_btn.place(x=545, y=1)
###############################################################################
########################root ELEMENTS###########################################   
root = tk.Tk()

root.geometry('900x505')
root.title("RamanMaps")
global f
f = Figure(facecolor='white')
global canvas
canvas = FigureCanvasTkAgg(f, master=root)
canvas.get_tk_widget().place(x=0, y=30)
canvas.get_tk_widget().configure(border=0, width=650, height=450)
notification_label = tk.Label()
notification_label.config(
    text="Please note, that mapping may take several minutes, depending on number of measurements.")
notification_label.place(x=10, y=484)

# file_name_label = tk.Label(text='File name:', font="Helvetica 10 bold")
# file_name_label.place(x=680, y=10)
file_name = tk.Entry()
file_name.place(x=680, y=35)
file_name.config(width=26)

map_btn = tk.Button(text='Create map', command=map_maker, width=18)
map_btn.place(x=680, y=475)

def save_figure():
    file_opt = options = {}
    options['defaultextension'] = '.png'
    options['filetypes'] = [('vector graphics', '.png')]
    options['title'] = 'Save image'
    global f
    Tk().withdraw() # we don't want a full root, so keep the root window from appearing
    filepath = asksaveasfilename(**file_opt) # show an "Save" dialog box and return the path to the selected file
    f.savefig(filepath, bbox_inches='tight')


def open_file():
    file_opt = options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('text files', '.txt')]
    options['title'] = 'Open text file with data'
    Tk().withdraw()
    filepath = askopenfilename(**file_opt)
    file_name.delete(0, END)
    file_name.insert(0, filepath)
    file_name.xview_moveto(1.0)


save_btn = tk.Button(text='Save image', command=save_figure, width=10)
save_btn.place(x=546, y=1)
save_btn.place_forget()

open_btn = tk.Button(text='Open file', command=open_file, width=6)
open_btn.place(x=680, y=5)

def rbtn_select(event):
    if event.widget["value"] == 1:
        map_btn.config(text="Create map!", command=map_maker)

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
        if first_mapping == 1:
            notification_label.config(
                text="Please note, that mapping may take even a few hours, depending on number of measurements.")
    elif event.widget["value"] == 2:
        map_btn.config(text="Perform fitting!", command=lambda: single_fitting(True))
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
        if first_mapping == 1:
            notification_label.config(text="")
    elif event.widget["value"] == 3:
        map_btn.config(text="Perform fitting!", command=lambda: single_fitting(False))
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
        if first_mapping == 1:
            notification_label.config(text="")
    elif event.widget["value"] == 4:
        map_btn.config(text="Calibrate!", command=calibrate)

        max_size_label.config(text="Max Raman shift")
        max_size.delete(0, END)
        max_size.insert(0, "600")
        min_size_label.config(text="Min Raman shift")
        min_size.delete(0, END)
        min_size.insert(0, "450")

        size_step.config(state="readonly")
        omega_0.config(state="readonly")
        gamma_0.config(state="readonly")
        amplitude.config(state="readonly")
        if first_mapping == 1:
            notification_label.config(text="")


def map_rbtn_select(event):
    global canvas
    global selected_dataset
    global map_styles
    global datasets_for_map
    global dataset_titles
    global selected_style
    selected_style = event.widget["value"] - 1

    f.clf()
    a = f.add_subplot(1, 1, 1)
    f.suptitle(dataset_titles[selected_dataset])
    b = a.pcolor(datasets_for_map[selected_dataset][0], cmap=map_styles[selected_style], vmin=datasets_for_map[selected_dataset][1],
             vmax=datasets_for_map[selected_dataset][2])
    f.colorbar(b)
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().config(border=0, bg='yellow', width=650, height=450)


def map_dataset_rbtn_select(event):
    global canvas
    global dataset_titles
    global selected_dataset
    selected_dataset = event.widget["value"] - 1
    f.clf()
    a = f.add_subplot(1, 1, 1)
    f.suptitle(dataset_titles[selected_dataset])
    b = a.pcolor(datasets_for_map[selected_dataset][0], cmap='PuBuGn', vmin=datasets_for_map[selected_dataset][1],
             vmax=datasets_for_map[selected_dataset][2])  #inne fajne cmap'y do wyboru: 'seismic' , 'coolwarm'
    f.colorbar(b)
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().config(border=0, width=650, height=450)


v = tk.IntVar()
rbtn_one = tk.Radiobutton(root,
                          text="Mapping file (four columns)",
                          padx=0,
                          variable=v,
                          value=1)
rbtn_one.place(x=670, y=53)
rbtn_one.select()
rbtn_one.bind("<Button-1>", rbtn_select)
rbtn_two = tk.Radiobutton(root,
                          text="Single sample file (four columns)",
                          padx=0,
                          variable=v,
                          value=2)
rbtn_two.place(x=670, y=73)
rbtn_two.bind("<Button-1>", rbtn_select)
rbtn_three = tk.Radiobutton(root,
                          text="Single sample file (two columns)",
                          padx=0,
                          variable=v,
                          value=3)
rbtn_three.place(x=670, y=93)
rbtn_three.bind("<Button-1>", rbtn_select)
rbtn_four = tk.Radiobutton(root,
                            text="Calibration file (two columns)",
                            padx=0,
                            variable=v,
                            value=4)
rbtn_four.place(x=670, y=113)
rbtn_four.bind("<Button-1>", rbtn_select)



# initial_shift_label = tk.Label(text='Initial Raman Shift [cm^-1]:')
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

type_var = tk.IntVar()
type_rbtn_one = tk.Radiobutton(root,
                          text="policrystalline",
                          padx=0,
                          variable=type_var,
                          value=1)
type_rbtn_one.place(x=700, y=135)
type_rbtn_one.select()
# type_rbtn_one.bind("<Button-1>", type_rbtn_select)
type_rbtn_two = tk.Radiobutton(root,
                          text="semiamorphous",
                          padx=0,
                          variable=type_var,
                          value=2)
type_rbtn_two.place(x=700, y=155)
# type_rbtn_two.bind("<Button-1>", type_rbtn_select)
type_rbtn_three = tk.Radiobutton(root,
                          text="amorphous",
                          padx=0,
                          variable=type_var,
                          value=3)
type_rbtn_three.place(x=700, y=175)
# type_rbtn_three.bind("<Button-1>", type_rbtn_select)

min_size_label = tk.Label(text="Min grain size to be checked [nm]:")
min_size_label.place(x=680, y=200)
min_size = tk.Entry()
min_size.place(x=680, y=220)
min_size.insert(0, 1)

max_size_label = tk.Label(text="Max grain size to be checked [nm]:")
max_size_label.place(x=680, y=240)
max_size = tk.Entry()
max_size.place(x=680, y=260)
max_size.insert(0, 10)

size_step_label = tk.Label(text="Step when checking size [nm]:")
size_step_label.place(x=680, y=280)
size_step = tk.Entry()
size_step.place(x=680, y=300)
size_step.insert(0, 0.5)

parameters_label = tk.Label(text="Initial parameters of main peak:", font="Helvetica 10 bold")
parameters_label.place(x=680, y=325)

omega_0_label = tk.Label(text="omega [cm^-1]")
omega_0_label.place(x=680, y=345)
omega_0 = tk.Entry()
omega_0.place(x=680, y=365)
omega_0.insert(0, "521.604599812")

gamma_0_label = tk.Label(text="gamma [cm^-1]")
gamma_0_label.place(x=680, y=385)
gamma_0 = tk.Entry()
gamma_0.place(x=680, y=405)
gamma_0.insert(0, "4.42010161931")

amplitude_label = tk.Label(text="ampl [arb. u]")
amplitude_label.place(x=680, y=425)
amplitude = tk.Entry()
amplitude.place(x=680, y=445)
amplitude.insert(0, "71535.2586963")

map_var = tk.IntVar()
map_rbtn_one = tk.Radiobutton(root,
                          text="mode1",
                          padx=0,
                          variable=map_var,
                          value=1)
map_rbtn_one.select()
map_rbtn_one.bind("<Button-1>", map_rbtn_select)

map_rbtn_two = tk.Radiobutton(root,
                          text="mode2",
                          padx=0,
                          variable=map_var,
                          value=2)
map_rbtn_two.bind("<Button-1>", map_rbtn_select)
map_rbtn_three = tk.Radiobutton(root,
                            text="mode3",
                            padx=0,
                            variable=map_var,
                            value=3)
map_rbtn_three.bind("<Button-1>", map_rbtn_select)
map_rbtn_four = tk.Radiobutton(root,
                            text="mode4",
                            padx=0,
                            variable=map_var,
                            value=4)
map_rbtn_four.bind("<Button-1>", map_rbtn_select)
map_rbtn_one.place_forget()
map_rbtn_two.place_forget()
map_rbtn_three.place_forget()
map_rbtn_four.place_forget()

map_param_var = tk.IntVar()
map_param_rbtn_one = tk.Radiobutton(root,
                          text="diameter",
                          padx=0,
                          variable=map_param_var,
                          value=1)
map_param_rbtn_one.select()
map_param_rbtn_one.bind("<Button-1>", map_dataset_rbtn_select)

map_param_rbtn_two = tk.Radiobutton(root,
                          text="omega_0",
                          padx=0,
                          variable=map_param_var,
                          value=2)
map_param_rbtn_two.bind("<Button-1>", map_dataset_rbtn_select)
map_param_rbtn_three = tk.Radiobutton(root,
                            text="hwhm",
                            padx=0,
                            variable=map_param_var,
                            value=3)
map_param_rbtn_three.bind("<Button-1>", map_dataset_rbtn_select)
map_param_rbtn_four = tk.Radiobutton(root,
                            text="intensity",
                            padx=0,
                            variable=map_param_var,
                            value=4)
map_param_rbtn_four.bind("<Button-1>", map_dataset_rbtn_select)
map_param_rbtn_one.place_forget()
map_param_rbtn_two.place_forget()
map_param_rbtn_three.place_forget()
map_param_rbtn_four.place_forget()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()


