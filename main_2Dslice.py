import matplotlib.pyplot as plt
import numpy as np
from element import *
from thermal_model import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import math
            
def h_to_k(h, axis, x, y, z):

    """
    approximate a convection coefficient h by thru-half-element conduction

    args:
        h       desired convection coefficient, W/m^2K
        axis    must be 'x', 'y', 'z', or 'all'. This describes the intended convection faces of the element.
        x       element size
        y       element size
        z       element size

    returns:
        k       a conductivity value in W/mK, that if used in the prescribed element dimensions and convected on the prescribed face,
                will produce a conduction Rth equivalent to the desired convection Rth
    """

    assert isinstance(axis, str), "my_str is not of string type"

    k = None
    if axis == "x":
        k = x*h/2 
    elif axis == "y":
        k = y*h/2 
    elif axis == "z":
        k = z*h/2 
    elif axis == "all":
        if x == y and y == z:
            k = x*h/2 
        else:
            raise ValueError("element must be cubic (x=y=z)")
    else:
        raise ValueError("invalid input, axis must be 'x', 'y', 'z', or 'all'")

    return k
# =====================================================================================================
# LABELS AND MATERIALS --------------------------------------------------------------------------------
# =====================================================================================================

null =                      {'label':   'null',       # this might change if we want locally powered or observed chunks
                             'ref':     '.',
                             'eltype':  'null',     
                             'k':       1,
                             'cp':      1,    
                             'rho':     1}

reservoir =                 {'label':   'infinite heatsink',       # this might change if we want locally powered or observed chunks
                             'ref':     '~',
                             'eltype':  'reservoir',     
                             'k':       h_to_k(h=10026, axis="z", x=0.0001, y=0.001, z=0.001),
                             'cp':      1,    
                             'rho':     1}

aluminum =                  {'label':   'aluminum',
                             'ref':     'A',
                             'eltype':  'solid',
                             'k':       200,
                             'cp':      0.9,    
                             'rho':     2700}

tim =                       {'label':   'tim',
                             'ref':     't',
                             'eltype':  'solid',
                             'k':       3.5,
                             'cp':      0.9,    
                             'rho':     3000}

glass =                     {'label':   'glass',
                             'ref':     'g',
                             'eltype':  'solid',
                             'k':       1,
                             'cp':      0.83,    
                             'rho':     2250}

resistor =                 {'label':   'resistor',
                             'ref':     'r',
                             'eltype':  'solid',
                             'k':       3,
                             'cp':      0.7,    
                             'rho':     3200}

alumina =                   {'label':   'alumina',
                             'ref':     'a',
                             'eltype':  'solid',
                             'k':       30,
                             'cp':      0.9,    
                             'rho':     4000}

solder =                    {'label':   'solder',
                             'ref':     'L',
                             'eltype':  'solid',
                             'k':       59,
                             'cp':      0.232,    
                             'rho':     7370} 

mask =                      {'label':   'mask',
                             'ref':     'm',
                             'eltype':  'solid',
                             'k':       1.75,
                             'cp':      1.3,    
                             'rho':     1300}

copper =                    {'label':   'copper',
                             'ref':     'C',
                             'eltype':  'solid',
                             'k':       400,
                             'cp':      0.385,    
                             'rho':     8000}

thermistor =                {'label':   'thermistor',
                             'ref':     'T',
                             'eltype':  'solid',
                             'k':       30,
                             'cp':      0.9,    
                             'rho':     4000}

# =====================================================================================================
# 3D MODELING -----------------------------------------------------------------------------------------
# =====================================================================================================
frequency = 60  # Frequency in Hz
dt = 1 / (frequency * 10)  # Time step, smaller than 1/60 s for smooth curve
P_of_t = lambda t: 100 if t <= 0.25 else 0
#P_of_t = lambda t: math.sin(2*math.pi*frequency*t)
T_of_t = lambda t: math.sin(2*math.pi*frequency/10*t)

void_sizes = []
T_avg_list = []
T_max_list = []

my_model = ThermalModel()
# copper
my_model.create_volume(x_vol=0.003,y_vol=0.001,z_vol=0.000035,x_el=0.0001,y_el=0.001,z_el=0.000035,properties=null,T=0,P=0)                                      
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=30,y_ext=1,z_ext=1,properties=copper,T=0,P=0)             
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=6,y_ext=1,z_ext=1,properties=null,T=0,P=0)             

# solder
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.0001,thickness_element=0.0001,properties=null)            
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=30,y_ext=1,z_ext=1,properties=mask,T=0,P=0)             
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=24,y_ext=1,z_ext=1,properties=solder,T=0,P=0)             
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=6,y_ext=1,z_ext=1,properties=null,T=0,P=0)             

# Alumina
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.00066,thickness_element=0.00011,properties=tim)            
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=14,y_ext=1,z_ext=6,properties=alumina,T=0,P=0)             

# resistor
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.00002,thickness_element=0.00002,properties=tim)            
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=10,y_ext=1,z_ext=1,properties=resistor,T=0,P=1)             

# resistor
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.00002,thickness_element=0.00002,properties=tim)            
my_model.modify_volume_centered_slice(axis="z",polarity=-1,index_start_extrusion=-1,x_ext=10,y_ext=1,z_ext=1,properties=glass,T=0,P=0)             

# tim, remainder
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.000309,thickness_element=0.000103,properties=tim)            

# aluminum
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.002,thickness_element=0.001,properties=aluminum)            

# convection
my_model.extend_volume(axis="z",polarity=1,size_extrusion=0.001,thickness_element=0.001,properties=reservoir)            

#nullify half
#my_model.modify_volume_centered_slice(axis="y",polarity=1,index_start_extrusion=0,x_ext=30,y_ext=10,z_ext=16,properties=null,T=0,P=0)             

#my_model.visualize()
#exit()

my_model.generate_nodes()
my_model.converge_this_label_T(resistor)
my_model.observe_these_labels_T_avg([resistor])
my_model.observe_these_labels_T_max([resistor, solder])
#my_model.observe_these_labels_P([silicon])

df = my_model.solve(dt=0.00003,dt_sampling=0.001,t_max=0.01,dT_dt_converge=0.1)
print(df)
    

#enclose all of the above in a loop
#add a mod_vol_center to the solder TIM layer, introducing a null space of variable size
#grab the last value of the second column
#



# ==================================================================================================
# PLOTTING -----------------------------------------------------------------------------------------
# ==================================================================================================

if(1):

    import matplotlib.pyplot as plt

    # Create a figure and axis for the temperature plot
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # TEMPERATURES -------------------------------------------------------------------------------------
    # Plot temperature traces on the left y-axis
    for column in df.columns:
        if column.startswith("T_"):  # Check for temperature columns
            ax1.plot(df['time'], df[column], label=column)

    # Label the left y-axis
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature and Power Traces Over Time')

    # Create a second y-axis for the power traces
    ax2 = ax1.twinx()

    # POWERS -------------------------------------------------------------------------------------
    # Plot power traces on the right y-axis
    for column in df.columns:
        if column.startswith("P_"):  # Check for power columns
            ax2.plot(df['time'], df[column], linestyle='--', label=column)

    # Label the right y-axis
    ax2.set_ylabel('Power (W)')

    # COMBINE AND PLOT -------------------------------------------------------------------------------------
    # Combine legends from both axes
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    legend = ax1.legend(lines_1 + lines_2, labels_1 + labels_2, title='Traces', framealpha=1.0)  # Opaque legend

    # Add grid - snap to temperature (primary y-axis) and ensure both horizontal and vertical lines
    ax1.grid(True, which='both', axis='both')  # Apply grid to both axes
    ax2.grid(False)  # Disable grid for the secondary y-axis to avoid double lines

    # Show plot
    plt.show()


# DATA VISUALIZATION -----------------------------------------------------------------------------------------

label_to_color = {
    aluminum['label']:   "turquoise",
    tim['label']:        "pink",
    glass['label']:      "darkgrey",
    resistor['label']:   "purple",
    alumina['label']:    "white",
    solder['label']:     "grey",
    mask['label']:       "green",
    copper['label']:     "orange",
    reservoir['label']:  "blue",
    null['label']:       "none",
    thermistor['label']: "red",
}
if(1):

    # Initialize min and max values for x, y, and z
    x_min = float('inf')
    y_min = float('inf')
    z_min = float('inf')
    x_max = float('-inf')
    y_max = float('-inf')
    z_max = float('-inf')

    # Iterate over each element to find the min and max for each axis
    for x in range(my_model.x_el_qty_old):
        for y in range(my_model.y_el_qty_old):
            for z in range(my_model.z_el_qty_old):
                element = my_model.elements[x][y][z]
                x_min = min(x_min, element.x_corner)
                y_min = min(y_min, element.y_corner)
                z_min = min(z_min, element.z_corner)
                x_max = max(x_max, element.x_corner + element.x)
                y_max = max(y_max, element.y_corner + element.y)
                z_max = max(z_max, element.z_corner + element.z)

    # Find the overall range for each axis
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)

    # Set up the figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot each element
    for x in range(my_model.x_el_qty_old):
        for y in range(my_model.y_el_qty_old):
            for z in range(my_model.z_el_qty_old):
                element = my_model.elements[x][y][z]
                color = label_to_color[element.label]
                if color != "none":  # Skip plotting if the color is "none"
                    ax.bar3d(element.x_corner, element.y_corner, element.z_corner, element.x, element.y, element.z, color=color, alpha=0.8, edgecolor="k")

    # Set equal scaling on all axes
    ax.set_xlim([x_min, x_min + max_range])
    ax.set_ylim([y_min, y_min + max_range])
    ax.set_zlim([z_min, z_min + max_range])

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Grid Visualization')

    # Show the plot
    plt.show()


if(0):
    # Initialize min and max values for x, y, and z
    x_min = float('inf')
    y_min = float('inf')
    z_min = float('inf')
    x_max = float('-inf')
    y_max = float('-inf')
    z_max = float('-inf')

    # Iterate over each element to find the min and max for each axis
    for x in range(my_model.x_el_qty_old):
        for y in range(my_model.y_el_qty_old):
            for z in range(my_model.z_el_qty_old):
                element = my_model.elements[x][y][z]
                x_min = min(x_min, element.x_corner)
                y_min = min(y_min, element.y_corner)
                z_min = min(z_min, element.z_corner)
                x_max = max(x_max, element.x_corner + element.x)
                y_max = max(y_max, element.y_corner + element.y)
                z_max = max(z_max, element.z_corner + element.z)

    # Find the overall range for each axis
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)

    # Determine min and max temperature for color scaling
    T_min = float('inf')
    T_max = float('-inf')

    for x in range(my_model.x_el_qty_old):
        for y in range(my_model.y_el_qty_old):
            for z in range(my_model.z_el_qty_old):
                element = my_model.elements[x][y][z]
                T_min = min(T_min, element.T)
                T_max = max(T_max, element.T)

    # Normalize the temperature values to [0, 1] for color mapping
    norm = plt.Normalize(T_min, T_max)
    cmap = cm.get_cmap('jet')  # You can choose other colormaps as well

    # Set up the figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot each element with color based on temperature
    for x in range(my_model.x_el_qty_old):
        for y in range(my_model.y_el_qty_old):
            for z in range(my_model.z_el_qty_old):
                element = my_model.elements[x][y][z]
                color = label_to_color[element.label]
                if color != "none":  # Skip plotting if the color is "none"
                    color = cmap(norm(element.T))  # Map temperature to a color
                    ax.bar3d(element.x_corner, element.y_corner, element.z_corner, element.x, element.y, element.z, color=color, alpha=0.8, edgecolor="k")

    # Set equal scaling on all axes
    ax.set_xlim([x_min, x_min + max_range])
    ax.set_ylim([y_min, y_min + max_range])
    ax.set_zlim([z_min, z_min + max_range])

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Grid Visualization with Temperature Color Scale')

    # Add a color bar to indicate the temperature scale
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    fig.colorbar(mappable, ax=ax, label='Temperature (K)')

    # Show the plot
    plt.show()