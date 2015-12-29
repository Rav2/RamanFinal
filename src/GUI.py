try:
    # for Python2
    import Tkinter as tk
    from tkFileDialog import askopenfilename, asksaveasfilename
except ImportError:
    # for Python3
    import tkinter as tk
    from tkinter.filedialog import askopenfilename, asksaveasfilename
# import Tkinter as tk
# from Tkinter import *
import aggregate_data as agg
from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import fit_old as fo
import nanograin_size as ns
import load_files

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
global number_of_peaks_to_fit
number_of_peaks_to_fit = 3
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
                                   float(omega_0_1.get()), float(gamma_0_1.get()), float(amplitude_1.get()), number_of_peaks_to_fit
                                   )

    global datasets_for_map
    del datasets_for_map[:]
    # for arg in datasets_for_map:
    #     datasets_for_map.remove(arg)
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
    a.set_xlabel("x ["+u"\u00B5"+"m]")
    a.set_ylabel("y ["+u"\u00B5"+"m]")
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
    """
    Gets parameters from window, executes fitting for the chosen file (four or two columns depending on argument's value)
    displays fitting curves and data points from file
    :param are_four_col:
    :return:
    """
    # LOADING DATA AND AGGREGATION
    save_btn.place_forget()
    data = None
    if are_four_col :
        data = load_files.load_mapping_file(file_name.get())
    else:
        data = load_files.load_two_column_file(file_name.get())
    omega = data[2]
    intensity = data[3]

    peaks_params = []
    if active_var_1.get():
        peak1 = []
        peak1.append(float(omega_0_1.get()))
        peak1.append(float(gamma_0_1.get()))
        peak1.append(float(amplitude_1.get()))
        peak1.append(float(place_from_1.get()))
        peak1.append(float(place_to_1.get()))
        peak1.append(visible_var_1.get())
        peak1.append(is_lorentz_var_1.get())
        peaks_params.append(peak1)
    if active_var_2.get():
        peak2 = []
        peak2.append(float(omega_0_2.get()))
        peak2.append(float(gamma_0_2.get()))
        peak2.append(float(amplitude_2.get()))
        peak2.append(float(place_from_2.get()))
        peak2.append(float(place_to_2.get()))
        peak2.append(visible_var_2.get())
        peak2.append(is_lorentz_var_2.get())
        peaks_params.append(peak2)
    if active_var_3.get():
        peak3 = []
        peak3.append(float(omega_0_3.get()))
        peak3.append(float(gamma_0_3.get()))
        peak3.append(float(amplitude_3.get()))
        peak3.append(float(place_from_3.get()))
        peak3.append(float(place_to_3.get()))
        peak3.append(visible_var_3.get())
        peak3.append(is_lorentz_var_3.get())
        peaks_params.append(peak3)


    result = ns.find_grain_diameter(omega, intensity, float(min_size.get()), float(max_size.get())+float(size_step.get()),
                                    float(size_step.get()),
                                    float(omega_0_1.get()), float(gamma_0_1.get()), float(amplitude_1.get()),
                                    True, number_of_peaks_to_fit, peaks_params)
    global f
    f.clf()
    f = result[4]
    f.suptitle(file_name.get().split('/')[-1])
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

    omega_0_1.delete(0, tk.END)
    omega_0_1.insert(0, result[0][0])
    omega_0_1.config(foreground="green")
    gamma_0_1.delete(0, tk.END)
    gamma_0_1.insert(0, result[0][1])
    gamma_0_1.config(foreground="green")
    amplitude_1.delete(0, tk.END)
    amplitude_1.insert(0, result[0][2])
    amplitude_1.config(foreground="green")
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

root.geometry('1300x505')
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
    tk.Tk().wm_withdraw() # we don't want a full root, so keep the root window from appearing
    filepath = asksaveasfilename(**file_opt) # show an "Save" dialog box and return the path to the selected file
    f.savefig(filepath, bbox_inches='tight')


def open_file():
    file_opt = options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('text files', '.txt')]
    options['title'] = 'Open text file with data'
    tk.Tk().wm_withdraw()
    filepath = askopenfilename(**file_opt)
    file_name.delete(0, tk.END)
    file_name.insert(0, filepath)
    file_name.xview_moveto(1.0)


save_btn = tk.Button(text='Save image', command=save_figure, width=10)
save_btn.place(x=546, y=1)
save_btn.place_forget()

open_btn = tk.Button(text='Open file', command=open_file, width=6)
open_btn.place(x=680, y=5)


def rbtn_select(event):
    """
    This function changes programmes appearence depending on the selected mode
    :param event:
    :return:
    """
    if event.widget["value"] == 1:
        map_btn.config(text="Create map!", command=map_maker)

        max_size_label.config(text="Max grain size to be checked [nm]:")
        max_size.delete(0, tk.END)
        max_size.insert(0, 10)
        min_size_label.config(text="Min grain size to be checked [nm]:")
        min_size.delete(0, tk.END)
        min_size.insert(0, 1)

        peak_active_1.config(state="disabled")
        peak_active_2.config(state="normal")
        peak_active_3.config(state="normal")
        peak_visible_1.config(state="normal")
        peak_visible_2.config(state="normal")
        peak_visible_3.config(state="normal")
        is_lorentz_1.config(state="normal")
        is_lorentz_2.config(state="normal")
        is_lorentz_3.config(state="normal")
        max_size.config(state=tk.NORMAL)
        min_size.config(state=tk.NORMAL)
        size_step.config(state=tk.NORMAL)

        place_from_1.config(state="normal")
        place_to_1.config(state="normal")
        omega_0_1.config(state=tk.NORMAL)
        gamma_0_1.config(state=tk.NORMAL)
        amplitude_1.config(state="disabled")

        place_from_2.config(state="normal")
        place_to_2.config(state="normal")
        omega_0_2.config(state=tk.NORMAL)
        gamma_0_2.config(state=tk.NORMAL)
        amplitude_2.config(state="disabled")

        place_from_3.config(state="normal")
        place_to_3.config(state="normal")
        omega_0_3.config(state=tk.NORMAL)
        gamma_0_3.config(state=tk.NORMAL)
        amplitude_3.config(state="disabled")
        if first_mapping == 1:
            notification_label.config(
                text="Please note, that mapping may take some time, depending on number of measurements.")
    elif event.widget["value"] == 2:
        map_btn.config(text="Perform fitting!", command=lambda: single_fitting(True))
        max_size_label.config(text="Max grain size to be checked [nm]:")
        max_size.delete(0, tk.END)
        max_size.insert(0, 10)
        min_size_label.config(text="Min grain size to be checked [nm]:")
        min_size.delete(0, tk.END)
        min_size.insert(0, 1)

        peak_active_1.config(state="disabled")
        peak_active_2.config(state="normal")
        peak_active_3.config(state="normal")
        peak_visible_1.config(state="normal")
        peak_visible_2.config(state="normal")
        peak_visible_3.config(state="normal")
        is_lorentz_1.config(state="normal")
        is_lorentz_2.config(state="normal")
        is_lorentz_3.config(state="normal")
        max_size.config(state=tk.NORMAL)
        min_size.config(state=tk.NORMAL)
        size_step.config(state=tk.NORMAL)

        place_from_1.config(state="normal")
        place_to_1.config(state="normal")
        omega_0_1.config(state=tk.NORMAL)
        gamma_0_1.config(state=tk.NORMAL)
        amplitude_1.config(state="disabled")

        place_from_2.config(state="normal")
        place_to_2.config(state="normal")
        omega_0_2.config(state=tk.NORMAL)
        gamma_0_2.config(state=tk.NORMAL)
        amplitude_2.config(state="disabled")

        place_from_3.config(state="normal")
        place_to_3.config(state="normal")
        omega_0_3.config(state=tk.NORMAL)
        gamma_0_3.config(state=tk.NORMAL)
        amplitude_3.config(state="disabled")
        if first_mapping == 1:
            notification_label.config(text="")
    elif event.widget["value"] == 3:
        map_btn.config(text="Perform fitting!", command=lambda: single_fitting(False))
        max_size_label.config(text="Max grain size to be checked [nm]:")
        max_size.delete(0, tk.END)
        max_size.insert(0, 10)
        min_size_label.config(text="Min grain size to be checked [nm]:")
        min_size.delete(0, tk.END)
        min_size.insert(0, 1)

        peak_active_1.config(state="disabled")
        peak_active_2.config(state="normal")
        peak_active_3.config(state="normal")
        peak_visible_1.config(state="normal")
        peak_visible_2.config(state="normal")
        peak_visible_3.config(state="normal")
        is_lorentz_1.config(state="normal")
        is_lorentz_2.config(state="normal")
        is_lorentz_3.config(state="normal")

        max_size.config(state=tk.NORMAL)
        min_size.config(state=tk.NORMAL)
        size_step.config(state=tk.NORMAL)

        place_from_1.config(state="normal")
        place_to_1.config(state="normal")
        omega_0_1.config(state=tk.NORMAL)
        gamma_0_1.config(state=tk.NORMAL)
        amplitude_1.config(state="normal")

        place_from_2.config(state="normal")
        place_to_2.config(state="normal")
        omega_0_2.config(state=tk.NORMAL)
        gamma_0_2.config(state=tk.NORMAL)
        amplitude_2.config(state="normal")

        place_from_3.config(state="normal")
        place_to_3.config(state="normal")
        omega_0_3.config(state=tk.NORMAL)
        gamma_0_3.config(state=tk.NORMAL)
        amplitude_3.config(state="normal")
        if first_mapping == 1:
            notification_label.config(text="")
    elif event.widget["value"] == 4:
        map_btn.config(text="Calibrate!", command=calibrate)

        max_size_label.config(text="Max Raman shift")
        max_size.delete(0, tk.END)
        max_size.insert(0, 600)
        min_size_label.config(text="Min Raman shift")
        min_size.delete(0, tk.END)
        min_size.insert(0, 450)

        size_step.config(state="readonly")
        place_from_1.config(state="readonly")
        place_to_1.config(state="readonly")
        omega_0_1.config(state="readonly")
        gamma_0_1.config(state="readonly")
        amplitude_1.config(state="readonly")
        peak_active_1.config(state="disabled")
        peak_visible_1.config(state="disabled")
        is_lorentz_1.config(state="disabled")

        place_from_2.config(state="readonly")
        place_to_2.config(state="readonly")
        omega_0_2.config(state="readonly")
        gamma_0_2.config(state="readonly")
        amplitude_2.config(state="readonly")
        peak_active_2.config(state="disabled")
        peak_visible_2.config(state="disabled")
        is_lorentz_2.config(state="disabled")

        place_from_3.config(state="readonly")
        place_to_3.config(state="readonly")
        omega_0_3.config(state="readonly")
        gamma_0_3.config(state="readonly")
        amplitude_3.config(state="readonly")
        peak_active_3.config(state="disabled")
        peak_visible_3.config(state="disabled")
        is_lorentz_3.config(state="disabled")
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
    a.set_xlabel("x ["+u"\u00B5"+"m]")
    a.set_ylabel("y ["+u"\u00B5"+"m]")
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
    a.set_xlabel("x ["+u"\u00B5"+"m]")
    a.set_ylabel("y ["+u"\u00B5"+"m]")
    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.get_tk_widget().place(x=0, y=30)
    canvas.get_tk_widget().config(border=0, width=650, height=450)


def type_rbtn_select(event):
    global number_of_peaks_to_fit
    if event.widget["value"] == 1:
        number_of_peaks_to_fit = 3
    elif event.widget["value"] == 2:
        number_of_peaks_to_fit = 2
    elif event.widget["value"] == 3:
        number_of_peaks_to_fit = 2

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
type_rbtn_one.bind("<Button-1>", type_rbtn_select)
type_rbtn_two = tk.Radiobutton(root,
                          text="semiamorphous",
                          padx=0,
                          variable=type_var,
                          value=2)
type_rbtn_two.place(x=700, y=155)
type_rbtn_two.bind("<Button-1>", type_rbtn_select)
type_rbtn_three = tk.Radiobutton(root,
                          text="amorphous",
                          padx=0,
                          variable=type_var,
                          value=3)
type_rbtn_three.place(x=700, y=175)
type_rbtn_three.bind("<Button-1>", type_rbtn_select)

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
############################################################################
########################  PARAMETRY PIKOW  #################################
shift_from_1 = tk.Label()


parameters_label = tk.Label(text="Initial parameters of peaks:", font="Helvetica 10 bold")
parameters_label.place(x=930, y=25)


place_from_label_1 =tk.Label(text="place from [cm^-1]")
place_from_label_1.place(x=1130, y=45)
place_from_1 = tk.Entry()
place_from_1.place(x=1130, y=65)
place_from_1.insert(0, "510.0")

place_to_label_1 =tk.Label(text="place to [cm^-1]")
place_to_label_1.place(x=1130, y=85)
place_to_1 = tk.Entry()
place_to_1.place(x=1130, y=105)
place_to_1.insert(0, "530.0")

omega_0_1_label_1 = tk.Label(text="omega [cm^-1]")
omega_0_1_label_1.place(x=930, y=45)
omega_0_1 = tk.Entry()
omega_0_1.place(x=930, y=65)
omega_0_1.insert(0, "521.604599812")

gamma_0_1_label_1 = tk.Label(text="gamma [cm^-1]")
gamma_0_1_label_1.place(x=930, y=85)
gamma_0_1 = tk.Entry()
gamma_0_1.place(x=930, y=105)
gamma_0_1.insert(0, "4.42010161931")

amplitude_label_1 = tk.Label(text="ampl [arb. u]")
amplitude_label_1.place(x=930, y=125)
amplitude_1 = tk.Entry()
amplitude_1.place(x=930, y=145)
amplitude_1.insert(0, "14000")
# amplitude_1.config(state="disabled")


active_var_1 = tk.BooleanVar()
active_var_1.set(True)
peak_active_1 = tk.Checkbutton(text="active", variable=active_var_1, state="disabled")
peak_active_1.place(x=930, y=170)
visible_var_1 = tk.BooleanVar()
visible_var_1.set(True)
peak_visible_1 = tk.Checkbutton(text="visible", variable=visible_var_1)
peak_visible_1.place(x=1020, y=170)
is_lorentz_var_1 = tk.BooleanVar()
is_lorentz_1 = tk.Checkbutton(text="lorentz", variable=is_lorentz_var_1)
is_lorentz_1.place(x=1120, y=170)
is_lorentz_var_1.set(True)

############################################################
place_from_label_2 =tk.Label(text="place from [cm^-1]")
place_from_label_2.place(x=1130, y=195)
place_from_2 = tk.Entry()
place_from_2.place(x=1130, y=215)
place_from_2.insert(0, "470.0")

place_to_label_2 =tk.Label(text="place to [cm^-1]")
place_to_label_2.place(x=1130, y=235)
place_to_2 = tk.Entry()
place_to_2.place(x=1130, y=255)
place_to_2.insert(0, "490.0")

omega_0_1_label_2 = tk.Label(text="omega [cm^-1]")
omega_0_1_label_2.place(x=930, y=195)
omega_0_2 = tk.Entry()
omega_0_2.place(x=930, y=215)
omega_0_2.insert(0, "480.0")

gamma_0_1_label_2 = tk.Label(text="gamma [cm^-1]")
gamma_0_1_label_2.place(x=930, y=235)
gamma_0_2 = tk.Entry()
gamma_0_2.place(x=930, y=255)
gamma_0_2.insert(0, "17")

amplitude_label_2 = tk.Label(text="ampl [arb. u]")
amplitude_label_2.place(x=930, y=275)
amplitude_2 = tk.Entry()
amplitude_2.place(x=930, y=295)
amplitude_2.insert(0, "10000")
# amplitude_2.config(state="disabled")


active_var_2 = tk.BooleanVar()
active_var_2.set(True)
peak_active_2 = tk.Checkbutton(text="active", variable=active_var_2)
peak_active_2.place(x=930, y=315)
visible_var_2 = tk.BooleanVar()
visible_var_2.set(True)
peak_visible_2 = tk.Checkbutton(text="visible", variable=visible_var_2)
peak_visible_2.place(x=1020, y=315)
is_lorentz_var_2 = tk.BooleanVar()
is_lorentz_2 = tk.Checkbutton(text="lorentz", variable=is_lorentz_var_2)
is_lorentz_2.place(x=1120, y=315)
is_lorentz_var_2.set(True)

############################################################
place_from_label_3 =tk.Label(text="place from [cm^-1]")
place_from_label_3.place(x=1130, y=340)
place_from_3 = tk.Entry()
place_from_3.place(x=1130, y=360)
place_from_3.insert(0, "485.0")

place_to_label_3 =tk.Label(text="place to [cm^-1]")
place_to_label_3.place(x=1130, y=380)
place_to_3 = tk.Entry()
place_to_3.place(x=1130, y=400)
place_to_3.insert(0, "500.0")

omega_0_1_label_3 = tk.Label(text="omega [cm^-1]")
omega_0_1_label_3.place(x=930, y=340)
omega_0_3 = tk.Entry()
omega_0_3.place(x=930, y=360)
omega_0_3.insert(0, "495.0")

gamma_0_1_label_3 = tk.Label(text="gamma [cm^-1]")
gamma_0_1_label_3.place(x=930, y=380)
gamma_0_3 = tk.Entry()
gamma_0_3.place(x=930, y=400)
gamma_0_3.insert(0, "20")

amplitude_label_3 = tk.Label(text="ampl [arb. u]")
amplitude_label_3.place(x=930, y=420)
amplitude_3 = tk.Entry()
amplitude_3.place(x=930, y=440)
amplitude_3.insert(0, "10000")
# amplitude_3.config(state="disabled")

active_var_3 = tk.BooleanVar()
active_var_3.set(True)
active_3 = tk.BooleanVar()
peak_active_3 = tk.Checkbutton(text="active", variable=active_var_3)
peak_active_3.place(x=930, y=460)
visible_var_3 =tk.BooleanVar()
visible_var_3.set(True)
peak_visible_3 = tk.Checkbutton(text="visible", variable=visible_var_3)
peak_visible_3.place(x=1020, y=460)
is_lorentz_var_3= tk.BooleanVar()
is_lorentz_3 = tk.Checkbutton(text="lorentz", variable=is_lorentz_var_3)
is_lorentz_3.place(x=1120, y=460)
is_lorentz_var_3.set(True)
############################################################################
############################################################################


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
                          text="omega_0_1",
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


