# change notes
#   changed decimal points on T printout to 3 places


from element import *
from node import *
import pandas as pd
import math

#TODO rename mesh.py to thermal_model.py

class ThermalModel:
    def __init__(self):

        self.elements = None                # 3D list of elements.
        self.nodes = []                     # 1D list of nodes (element pairs). 
        self.elements_to_power = []         # 1D list of element references. Usually a minority of elements generate heat. Run math on only those!
        self.elements_to_observe_T_avg = [] # 2D list of element references. We monitor the T_avg of these elements for data collection purposes.
        self.labels_to_observe_T_avg = []   # list of str,  corresponds to above
        self.elements_to_observe_T_max = [] # 2D list of element references. We monitor the T_max of these elements for data collection purposes.
        self.labels_to_observe_T_max = []   # list of str,  corresponds to above
        self.elements_to_observe_P = []     # 2D list of element references. We monitor the P of these elements for data collection purposes.
        self.labels_to_observe_P = []       # list of str,  corresponds to above
        self.elements_to_converge = []      # 1D list of element references. We monitor the dT/dt of these elements for covergence purposes.
        self.label_to_converge = None       # str,          corresponds to above
        self.elements_to_CoE = []           # 1D list of element references. These are non-null, non-reservoir elements that we care to compute conservation of energy for.
        self.elements_to_prescribe_T = []   # 1D list of element references. These are elements that were assigned a non-constant T(t)
        self.x_el_qty_old = 0               # Rectangular maximum envelope size of mesh. QTY elements in X
        self.y_el_qty_old = 0               # "                                                          Y
        self.z_el_qty_old = 0               # "                                                          Z
        self.df = None                      # pandas df, stores T vs t history
        self.x_vol_max = None               # maximum size in meters of the volume, x axis
        self.x_vol_max = None               # "                                     y axis
        self.x_vol_max = None               # "                                     z axis

    # ========================================================================================
    # CAD 
    # ========================================================================================

    def create_volume(self, 
            x_vol,          #   x dimension of rectangular volume,  units: m 
            y_vol,          #   y dimension of rectangular volume,  units: m 
            z_vol,          #   x dimension of rectangular volume,  units: m 
            x_el,           #   x dimension of single element,      units: m        
            y_el,           #   y dimension of single element,      units: m 
            z_el,           #   x dimension of single element,      units: m 
            properties,     #   properties dictionary
            #eltype,         #   element type, "null" or "solid"             
            #label=".",      #   string, describes the element. ex. Silicon.
            #ref=".",        #   char, a reference symbol used in the text-based visualization
            #k=1,            #   volume thermal conductivity         units: W/mK
            T=0,            #   volume temperature,                 units: K
            P=0             #   volume heating power generation     units: W
            #cp=10,          #   material specific heat              units: J/gK 
            #rho=1000       #   material density                    units: kg/m^3 
            ):
        
        """
        This is the first 3D modeling function that must be used.
        Establishes a new rectangular volume of material, with independent element sizes in x y and z

        Arguments: See above

        Return: None
        """

        if self.elements != None:
            raise ValueError("create_first_volume method called when mesh already exists")

        # get element count on each axis
        #TODO: raise value error if QTYs do not reflect intergers
        #TODO: modulo check
        x_el_qty = int(x_vol/x_el)
        y_el_qty = int(y_vol/y_el)
        z_el_qty = int(z_vol/z_el)

        # get watts per element
        total_el_qty = x_el_qty * y_el_qty * z_el_qty
        P_el = self.compute_element_wise_P(P, total_el_qty)


        # create 3d list
        # when accessing index, follow convention elements[x][y][z]
        self.elements = [[[Element(x_el, y_el, z_el, properties['eltype'], properties['label'], properties['ref'], properties['k'], T, P_el, properties['cp'], properties['rho']) for _ in range(z_el_qty)] for _ in range(y_el_qty)] for _ in range(x_el_qty)]

        # update the volume max sizes
        self.x_el_qty_old = x_el_qty
        self.y_el_qty_old = y_el_qty
        self.z_el_qty_old = z_el_qty

        
    def extend_volume(self,
            axis,               #   str, "x", "y", or "z"
            polarity,           #   int, 1, -1
            size_extrusion,     #   float, size of the extrusion off the main volume                                           units: m
            thickness_element,  #   float, thickness of a single element in the axis of choice (as specified by "direction")   units: m
            properties,                  #   properties dictionary
            #eltype,             #   str, element type, "null" or "solid"             
            #label=".",          #   str, describes the element. ex. Silicon.
            #ref=".",            #   char, a reference symbol used in the text-based visualization
            #k=1,                #   volume thermal conductivity         units: W/mK
            T=0,                #   volume temperature,                 units: K
            P=0,                #   volume heating power generation     units: W
            #cp=10,              #   material specific heat              units: J/gK 
            #rho=1000            #   material density                    units: kg/m^3 
            ):
        """
        This is a 3D modeling operation, where we extrude off of 1 face of the existing rectangular volume/envelope,
        whether it is homogenous or composite of many modeling operations.
        This add-on volume inherits element dimensions orthoganal to extrusion axis.
        A new element thickness can be assigned in the axis of this extrusion.

        Arguments: see above

        Returns: none
        """
        
        #TODO: modulo check the extruded element count, must be interger, nonzero
        #TODO: polarity check, +1 or -1

        # define volumes bounds (old, new, and dividing line)
        x_el_qty_new = self.x_el_qty_old 
        y_el_qty_new = self.y_el_qty_old 
        z_el_qty_new = self.z_el_qty_old 
        divider = None      # this is the index at which the next zone starts (new volume for polarity +1, old volume for polarity -1)
        if axis == "x":
            x_el_qty_new += int(size_extrusion/thickness_element)
            divider = self.x_el_qty_old if polarity == 1 else x_el_qty_new
        elif axis == "y":
            y_el_qty_new += int(size_extrusion/thickness_element)
            divider = self.y_el_qty_old if polarity == 1 else y_el_qty_new
        elif axis == "z":
            z_el_qty_new += int(size_extrusion/thickness_element)
            divider = self.z_el_qty_old if polarity == 1 else z_el_qty_new
        else:
            raise ValueError("invalid argument, axis")

        # initialize new 3d list
        elements_new = [[[None for _ in range(z_el_qty_new)] for _ in range(y_el_qty_new)] for _ in range(x_el_qty_new)]

        # calculate P per element for new volume
        P_el = self.compute_element_wise_P(P, (x_el_qty_new*y_el_qty_new*z_el_qty_new))

        # recreate mesh, duplicating old mesh in old volume, creating new mesh in new volume
        for x in range(int(x_el_qty_new)):
            for y in range(int(y_el_qty_new)):
                for z in range(int(z_el_qty_new)):
                    if axis == "x":
                        # Determine if this index is in the old or new volume
                        if (polarity == 1 and x < divider) or (polarity == -1 and x >= divider):
                            # Duplicate old volume element into this index if this is old volume
                            old_x = x if polarity == 1 else x - divider
                            elements_new[x][y][z] = self.elements[old_x][y][z]    
                        else:
                            # Initialize new element if this is new volume
                            # new element sizes orthoganal to extrusion/extension axis are inherited from old volume elements
                            this_y_el = self.elements[0][y][z].get_size("y")
                            this_z_el = self.elements[0][y][z].get_size("z")
                            elements_new[x][y][z] = Element(thickness_element, this_y_el, this_z_el, properties['eltype'], properties['label'], properties['ref'], properties['k'], T, P_el, properties['cp'], properties['rho'])

                    elif axis == "y":
                        # Determine if this index is in the old or new volume
                        if (polarity == 1 and y < divider) or (polarity == -1 and y >= divider):
                            # Duplicate old volume element into this index if this is old volume
                            old_y = y if polarity == 1 else y - divider
                            elements_new[x][y][z] = self.elements[x][old_y][z]  
                        else:
                            # Initialize new element if this is new volume
                            # new element sizes orthoganal to extrusion/extension axis are inherited from old volume elements
                            this_x_el = self.elements[x][0][z].get_size("x")
                            this_z_el = self.elements[x][0][z].get_size("z")
                            elements_new[x][y][z] = Element(this_x_el, thickness_element, this_z_el, properties['eltype'], properties['label'], properties['ref'], properties['k'], T, P_el, properties['cp'], properties['rho'])

                    elif axis == "z":
                        # Determine if this index is in the old or new volume
                        if (polarity == 1 and z < divider) or (polarity == -1 and z >= divider):
                            # Duplicate old volume element into this index if this is old volume
                            old_z = z if polarity == 1 else z - divider
                            elements_new[x][y][z] = self.elements[x][y][old_z]     
                        else:
                            # Initialize new element if this is new volume
                            # new element sizes orthoganal to extrusion/extension axis are inherited from old volume elements
                            this_x_el = self.elements[x][y][0].get_size("x")
                            this_y_el = self.elements[x][y][0].get_size("y")
                            elements_new[x][y][z] = Element(this_x_el, this_y_el, thickness_element, properties['eltype'], properties['label'], properties['ref'], properties['k'], T, P_el, properties['cp'], properties['rho'])

        # update total vol bounds, update element list
        self.x_el_qty_old = x_el_qty_new
        self.y_el_qty_old = y_el_qty_new
        self.z_el_qty_old = z_el_qty_new
        self.elements = elements_new

        return


    def modify_volume(self,
                      x_coord,      #   int, starting element index coordinate (lowest magnitude) for the corner of this rectangle volume scope
                      y_coord,      #   "
                      z_coord,      #   "
                      x_ext,        #   int, number of element indices to modify, from start coordinate moving positive           
                      y_ext,        #   " 
                      z_ext,        #   " 
                      properties,            #   properties dictionary
                      #eltype,       #   str, element type, "null" or "solid"             
                      #label=".",    #   str, describes the element. ex. Silicon.
                      #ref=".",      #   char, a reference symbol used in the text-based visualization
                      #k=1,          #   volume thermal conductivity         units: W/mK
                      T=0,          #   volume temperature,                 units: K
                      P=0,          #   volume heating power generation     units: W
                      #cp=10,        #   material specific heat              units: J/gK 
                      #rho=1000      #   material density                    units: kg/m^3 
                      ):
        """
        This is a mesh modification operation.
        The scope of the modification is a rectangular volume whose size is specified by element counts, not meters.
        The original element dimensions are unchanged, but their physical properties, element-type designation, and char/str labels are changed.
        This operation starts a corner of a rectangular modification scope at some coordinate, and extrudes that out in 3 dimensions of size.
        
        Arguments: see above

        Returns: none
        """

        #TODO: modification scope cannot exceed existing volume

        indeces_x = []
        indeces_y = []
        indeces_z = []
        for i in range(x_ext):
            indeces_x.append(x_coord + i)
        for i in range(y_ext):
            indeces_y.append(y_coord + i)
        for i in range(z_ext):
            indeces_z.append(z_coord + i)

        P_el = self.compute_element_wise_P(P, (x_ext*y_ext*z_ext))

        for x in indeces_x:
            for y in indeces_y:
                for z in indeces_z:
                    self.elements[x][y][z].modify(properties['eltype'], properties['label'], properties['ref'], properties['k'], T, P_el, properties['cp'], properties['rho'])
                    #TODO: check what happens if i create a new element and assign that 3d list index to that new element. what happens to the old element


    def modify_volume_centered_slice(self,
        axis,                       #   str, "x", "y", or "z"
        polarity,                   #   int (1 or -1). +1 starts the first extrusion index from axis zero. -1 starts it from opposite face.
        index_start_extrusion,      #   int, 
        x_ext,                      #   int, number of element indices to modify, orthoganal sizes are centered, in-axis size starts from "index_start_extrusion"          
        y_ext,                      #   " 
        z_ext,                      #   " 
        properties,                          #   properties dictionary
        #eltype,                     #   str, element type, "null" or "solid"             
        #label=".",                  #   str, describes the element. ex. Silicon.
        #ref=".",                    #   char, a reference symbol used in the text-based visualization
        #k=1,                        #   volume thermal conductivity         units: W/mK
        T=0,                        #   volume temperature,                 units: K
        P=0,                        #   volume heating power generation     units: W
        #cp=10,                      #   material specific heat              units: J/gK 
        #rho=1000                    #   material density                    units: kg/m^3 
        ):
        """
        This is a mesh modification operation.
        The scope of the modification is a rectangular volume whose size is specified by element counts, not meters.
        The original element dimensions are unchanged, but their physical properties, element-type designation, and char/str labels are changed.
        This operation modifies a centered rectangular volume/subvolume in a select axis and layer "slice" of the main volume.

        Arguments: see above

        Returns: none
        """


        max_index_modified_slice_thickness = None    
        min_index_modified_slice_thickness = None    


        # starting index check
        # convert index to positive, if negative value was input
        if axis == "x":
            if index_start_extrusion < 0:
                index_start_extrusion = self.x_el_qty_old + index_start_extrusion
            if index_start_extrusion > self.x_el_qty_old-1:
                raise ValueError("x index cannot start outside of existing volume")
        elif axis == "y":
            if index_start_extrusion < 0:
                index_start_extrusion = self.y_el_qty_old + index_start_extrusion
            if index_start_extrusion > self.y_el_qty_old-1:
                raise ValueError("y index cannot start outside of existing volume")
        elif axis == "z":
            if index_start_extrusion < 0:
                index_start_extrusion = self.z_el_qty_old + index_start_extrusion
            if index_start_extrusion > self.z_el_qty_old-1:
                raise ValueError("z index cannot start outside of existing volume")

        # ending index check
        if polarity == 1:
            if axis == "x":
                max_index_modified_slice_thickness = index_start_extrusion + x_ext - 1
                if max_index_modified_slice_thickness > self.x_el_qty_old - 1:
                    raise ValueError("x index cannot end outside of existing volume")
            elif axis == "y":
                max_index_modified_slice_thickness = index_start_extrusion + y_ext - 1
                if max_index_modified_slice_thickness > self.y_el_qty_old - 1:
                    raise ValueError("y index cannot end outside of existing volume")
            elif axis == "z":
                max_index_modified_slice_thickness = index_start_extrusion + z_ext - 1
                if max_index_modified_slice_thickness > self.z_el_qty_old - 1:
                    raise ValueError("z index cannot end outside of existing volume")
        elif polarity == -1:
            if axis == "x":
                min_index_modified_slice_thickness = index_start_extrusion - x_ext + 1
                if min_index_modified_slice_thickness > self.x_el_qty_old - 1:
                    raise ValueError("x index cannot end outside of existing volume")
            if axis == "y":
                min_index_modified_slice_thickness = index_start_extrusion - y_ext + 1
                if min_index_modified_slice_thickness > self.y_el_qty_old - 1:
                    raise ValueError("y index cannot end outside of existing volume")
            if axis == "z":
                min_index_modified_slice_thickness = index_start_extrusion - z_ext + 1
                if min_index_modified_slice_thickness > self.z_el_qty_old - 1:
                    raise ValueError("z index cannot end outside of existing volume")


        if axis == "x":
            if y_ext > self.y_el_qty_old: raise ValueError("y size of modification scope cannot exceed original volume")
            if z_ext > self.z_el_qty_old: raise ValueError("z size of modification scope cannot exceed original volume")
            if polarity ==  1: x_coord = index_start_extrusion
            if polarity == -1: x_coord = min_index_modified_slice_thickness
            y_coord = math.floor((self.y_el_qty_old-y_ext)/2)
            z_coord = math.floor((self.z_el_qty_old-z_ext)/2)
        elif axis == "y":
            if x_ext > self.x_el_qty_old: raise ValueError("x size of modification scope cannot exceed original volume")
            if z_ext > self.z_el_qty_old: raise ValueError("z size of modification scope cannot exceed original volume")
            x_coord = math.floor((self.x_el_qty_old-x_ext)/2)
            if polarity ==  1: y_coord = index_start_extrusion
            if polarity == -1: y_coord = min_index_modified_slice_thickness
            z_coord = math.floor((self.z_el_qty_old-z_ext)/2)
        elif axis == "z":
            if x_ext > self.x_el_qty_old: raise ValueError("x size of modification scope cannot exceed original volume")
            if y_ext > self.y_el_qty_old: raise ValueError("y size of modification scope cannot exceed original volume")
            x_coord = math.floor((self.x_el_qty_old-x_ext)/2)
            y_coord = math.floor((self.y_el_qty_old-y_ext)/2)
            if polarity ==  1: z_coord = index_start_extrusion
            if polarity == -1: z_coord = min_index_modified_slice_thickness

        self.modify_volume(x_coord=x_coord, y_coord=y_coord, z_coord=z_coord, x_ext=x_ext, y_ext=y_ext, z_ext=z_ext, properties=properties, T=T, P=P)


    def compute_element_wise_P(self, P, total_el_qty):
        """
        process element-wise power, whether power is constant or variable P(t)

        args: power (function or float)
              total element qty assigned to this volume of material

        returns: power (function or float), divided between element qty
        """
        if isinstance(P, (float, int)):
            return P / total_el_qty
        elif callable(P):
            return lambda t: P(t) / total_el_qty
        else:
            raise ValueError("P must be valid input")


    # ========================================================================================
    # SOLVER
    # ========================================================================================

    def generate_nodes(self):
        """
        Preprocessing function to be executed before solving.
        Generates a 1D list of nodes that connects active elements only, excluding inactive and null elements.
        This optimized node list facilitates quicker computation during the solving phase by reducing the complexity 
        associated with the original 3D mesh structure.

        Args: None

        Returns: None
        
        NOTE:
        there are no redundant linkings of elements
        there are no "existing neighbors" check when establishing neighbor status.
        as of writing, the neighbor list of every element is not utilized
        """


        # Step through all elements in 3D list
        for x in range(self.x_el_qty_old):
            for y in range(self.y_el_qty_old):
                for z in range(self.z_el_qty_old):
                    element = self.elements[x][y][z]

                    if element.get_eltype() == "null":
                        continue  # Skip if the current element is a null type

                    # Since we start in a corner of the mesh, we only need to check in 3 directions: X+ Y+ Z+
                    # Nodes are NOT created if:
                    #   the element is "null" type
                    #   the element is "reservoir" type and the neighbor is also "reservoir" type (pointless to heat exchange inside reservoire)

                    
                    # Check for neighbors in the x+ direction
                    if x + 1 < self.x_el_qty_old:
                        neighbor_x = self.elements[x + 1][y][z]
                        if neighbor_x.get_eltype() != "null" or (element.get_eltype() == "reservoir" and neighbor_x.get_eltype() != "reservoir"):
                            element.neighbors.append(neighbor_x)
                            neighbor_x.neighbors.append(element)
                            node = Node(element, neighbor_x, "x")
                            self.nodes.append(node)

                    # Check for neighbors in the y+ direction
                    if y + 1 < self.y_el_qty_old:
                        neighbor_y = self.elements[x][y + 1][z]
                        if neighbor_y.get_eltype() != "null" or (element.get_eltype() == "reservoir" and neighbor_x.get_eltype() != "reservoir"):
                            element.neighbors.append(neighbor_y)
                            neighbor_y.neighbors.append(element)
                            node = Node(element, neighbor_y, "y")
                            self.nodes.append(node)

                    # Check for neighbors in the z+ direction
                    if z + 1 < self.z_el_qty_old:
                        neighbor_z = self.elements[x][y][z + 1]
                        if neighbor_z.get_eltype() != "null" or (element.get_eltype() == "reservoir" and neighbor_x.get_eltype() != "reservoir"):
                            element.neighbors.append(neighbor_z)
                            neighbor_z.neighbors.append(element)
                            node = Node(element, neighbor_z, "z")
                            self.nodes.append(node)

        self.establish_corner_coordinates()
        

    def establish_corner_coordinates(self):
        """
        Assigns corner coordinates to every element in the 3D list self.elements.
        Each element has attributes x, y, z for size in meters.
        Also records the cumulative dimensions at the last index of x, y, and z.
        """

        # Initialize cumulative positions
        z_cumulative = 0  # Start at the origin for z-axis

        for z_index in range(self.z_el_qty_old):
            y_cumulative = 0  # Start at the origin for y-axis

            for y_index in range(self.y_el_qty_old):
                x_cumulative = 0  # Start at the origin for x-axis

                for x_index in range(self.x_el_qty_old):

                    # Retrieve the current element
                    element = self.elements[x_index][y_index][z_index]

                    # Assign corner coordinates based on cumulative positions
                    element.assign_corner_coordinates(x_cumulative, y_cumulative, z_cumulative)

                    # Update the cumulative x position for the next element in the x-direction
                    x_cumulative += element.x

                # Update the cumulative y position for the next element in the y-direction
                y_cumulative += self.elements[x_index][y_index][z_index].y

            # Update the cumulative z position for the next element in the z-direction
            z_cumulative += self.elements[x_index][y_index][z_index].z

        # Record the cumulative dimensions at the last indices
        self.x_vol_max = x_cumulative  # This will be the cumulative length in x-direction after the last element
        self.y_vol_max = y_cumulative  # This will be the cumulative length in y-direction after the last element
        self.z_vol_max = z_cumulative  # This will be the cumulative height in z-direction after the last element


    def observe_these_labels_T_avg(self, list_of_labels_T):
        """
        Mark these labels for AVERAGE TEMPERATURE observation
        Identifies and stores references to elements based on a list of labels. 
        Each label corresponds to a list of elements, which are added to self.elements_to_observe_T_avg.

        Args:
            list_of_labels (list of dicts): A list of labels to observe. Each label corresponds to a group of elements to be monitored.

        Returns:
            None

        TODO: announce unused/unmatched labels. raise error.
        TODO: assert input type
        """

        if not isinstance(list_of_labels_T, list):
            raise ValueError("input must be a list")


        if len(self.labels_to_observe_T_avg) == 0: 
            self.labels_to_observe_T_avg = list_of_labels_T
        else:
            raise ValueError("cannot input observed labels twice")

        for label in list_of_labels_T:                                                
            elements_current_label = []                                             
            for x in range(self.x_el_qty_old):  
                for y in range(self.y_el_qty_old): 
                    for z in range(self.z_el_qty_old): 
                        if self.elements[x][y][z].label == label['label']:
                            elements_current_label.append(self.elements[x][y][z])

            self.elements_to_observe_T_avg.append(elements_current_label)  # For this one label, add one column of element references 
                                                                     # After loop, we will have n columns of element references,
                                                                     # one for each label

        return
    

    def observe_these_labels_T_max(self, list_of_labels_T):
        """
        Mark these labels for MAXIMUM TEMPERATURE observation
        Identifies and stores references to elements based on a list of labels. 
        Each label corresponds to a list of elements, which are added to self.elements_to_observe_T_max.

        Args:
            list_of_labels (list of dicts): A list of labels to observe. Each label corresponds to a group of elements to be monitored.

        Returns:
            None

        TODO: announce unused/unmatched labels. raise error.
        TODO: assert input type
        """

        if not isinstance(list_of_labels_T, list):
            raise ValueError("input must be a list")


        if len(self.labels_to_observe_T_max) == 0: 
            self.labels_to_observe_T_max = list_of_labels_T
        else:
            raise ValueError("cannot input observed labels twice")

        for label in list_of_labels_T:                                                
            elements_current_label = []                                             
            for x in range(self.x_el_qty_old):  
                for y in range(self.y_el_qty_old): 
                    for z in range(self.z_el_qty_old): 
                        if self.elements[x][y][z].label == label['label']:
                            elements_current_label.append(self.elements[x][y][z])

            self.elements_to_observe_T_max.append(elements_current_label)  # For this one label, add one column of element references 
                                                                     # After loop, we will have n columns of element references,
                                                                     # one for each label

        return


    def observe_these_labels_P(self, list_of_labels_P):
        """
        Mark these label for POWER observation
        Identifies and stores references to elements based on a list of labels. 
        Each label corresponds to a list of elements, which are added to self.elements_to_observe_P.

        Args:
            list_of_labels (list of dicts): A list of labels to observe. Each label corresponds to a group of elements to be monitored.

        Returns:
            None

        TODO: announce unused/unmatched labels. raise error.
        TODO: assert input type
        """

        if not isinstance(list_of_labels_P, list):
            raise ValueError("input must be a list")

        if len(self.labels_to_observe_P) == 0: 
            self.labels_to_observe_P = list_of_labels_P
        else:
            raise ValueError("cannot input observed labels twice")

        for label in list_of_labels_P:                                                
            elements_current_label = []                                             
            for x in range(self.x_el_qty_old):  
                for y in range(self.y_el_qty_old): 
                    for z in range(self.z_el_qty_old): 
                        if self.elements[x][y][z].label == label['label']:
                            elements_current_label.append(self.elements[x][y][z])

            self.elements_to_observe_P.append(elements_current_label)  # For this one label, add one column of element references 
                                                                     # After loop, we will have n columns of element references,
                                                                     # one for each label

        return


    def converge_this_label_T(self, my_label):
        """
        Identifies and stores references to elements based on a label. 
        The label corresponds to a list of elements, which are added to self.elements_to_converge.

        Args:
            my_label (dict): Label to observe for solver convergence criteria. This label corresponds to a group of elements.

        Returns:
            None
        """

        if not isinstance(my_label, dict):
            raise ValueError("input must be a dict")

        if self.label_to_converge == None: 
            self.label_to_converge = my_label['label']
        else:
            raise ValueError("cannot input convergence label twice")

        for x in range(self.x_el_qty_old):  
            for y in range(self.y_el_qty_old): 
                for z in range(self.z_el_qty_old): 
                    if self.elements[x][y][z].label == my_label['label']:
                        self.elements_to_converge.append(self.elements[x][y][z])

        return


    def observe_T_avg(self, elements_observed):
        """
        computes average temperature of a group of elements

        Args: list of elements

        Returns: T_average
        """
        T_accumulated = 0
        for element in elements_observed:
            T_accumulated += element.T

        return T_accumulated / len(elements_observed)


    def observe_T_max(self, elements_observed):
        """
        returns the maximum observed temperature in a group of elements

        Args: list of elements

        Returns: T_max
        """
        T_max = 0
        for element in elements_observed:
            if element.T > T_max:
                T_max = element.T

        return T_max


    def observe_P(self, elements_observed):
        """
        computes total power of a group of elements

        Args: list of elements

        Returns: P_total
        """
        P_accumulated = 0
        for element in elements_observed:
            P_accumulated += element.P

        return P_accumulated


    def solve(self, dt, dt_sampling, t_max, dT_dt_converge):
        """
        Solves the thermal simulation by iterating through time steps until either the maximum time or 
        the temperature convergence criteria is met.

        Args:
            dt (float):                             Time step for the simulation in seconds.
            dt_sampling (float):                    add dt_sampling interval, which must be a multiple of dt (seconds)
            t_max (float):                          The maximum time in seconds at which the simulation will halt, regardless of convergence.
            dT_dt_converge (float):                 The rate of temperature change (in K/s) that is considered converged.

        Returns:
            pandas dataframe of T vs t histories of all specified labeled element groups
        """

        # create 1D line-up for:
        #   elements that are powered
        #   elements that we need to compute conservation of energy for (non-null, non-reservoir)
        for x in range(self.x_el_qty_old):  
            for y in range(self.y_el_qty_old): 
                for z in range(self.z_el_qty_old): 
                    this_element = self.elements[x][y][z]
                    # add elements with constant, non-zero power OR prescribed P(t)
                    if (this_element.P != 0 and this_element.P_is_constant) or this_element.P_is_constant is False:
                        self.elements_to_power.append(this_element)
                    # add elements whose temperatures change throughout simulation
                    if this_element.get_eltype() not in ["null", "reservoir"]:
                        self.elements_to_CoE.append(this_element)
                    # add elements that follow a prescribed T(t)
                    if this_element.T_is_prescribed:
                        self.elements_to_prescribe_T.append(this_element)

        iteration = int(0)
        iteration_interval_to_sample = int(dt_sampling/dt)
        t_cur = 0
        T_to_converge_prev = None
        T_to_converge_cur = None
        dT_dt_cur = None

        # Initialize DataFrame for T vs t
        prefixed_labels_T_avg = ["T_avg_" + label['label'] for label in self.labels_to_observe_T_avg]
        prefixed_labels_T_max = ["T_max_" + label['label'] for label in self.labels_to_observe_T_max]
        suffixed_labels_P = ["P_" + label['label'] for label in self.labels_to_observe_P]
        column_names = ["time"] + prefixed_labels_T_avg + prefixed_labels_T_max + suffixed_labels_P


        self.df = pd.DataFrame(columns=column_names)
        decimal_points = 3

        # iterate and solve until convergence or timeout
        while t_cur < t_max:

            if iteration % iteration_interval_to_sample == 0:

                # observe temperature of converging elements
                T_to_converge_prev = T_to_converge_cur
                T_to_converge_cur = self.observe_T_avg(self.elements_to_converge)

                # calculate dT/dt of converging elements
                if T_to_converge_cur != None and T_to_converge_prev != None:
                    dT_dt_cur = (T_to_converge_cur-T_to_converge_prev)/dt_sampling

                # check for convergence
                if dT_dt_cur != None and abs(dT_dt_cur) < abs(dT_dt_converge):
                    print("dT/dt convergence reached, terminating solver")
                    break

                # Gather this row of data
                row_data = [t_cur]
                for element_group in self.elements_to_observe_T_avg:
                    temperature = self.observe_T_avg(element_group)
                    row_data.append(round(temperature, decimal_points))
                for element_group in self.elements_to_observe_T_max:
                    temperature = self.observe_T_max(element_group)
                    row_data.append(round(temperature, decimal_points))
                for element_group in self.elements_to_observe_P:
                    power = self.observe_P(element_group)
                    row_data.append(round(power, decimal_points))

                # Get this row into df
                self.df.loc[len(self.df)] = row_data

                # print progress
                if t_cur != None and T_to_converge_cur != None and dT_dt_cur != None:
                    print(f"t:\t{t_cur:.6f}\tT:\t{T_to_converge_cur:.3f}\tdT/dt:\t{dT_dt_cur:.4f}")

            # calculate dE from interelement heat transfer
            for node in self.nodes:
                node.calculate_heattransfer(dt)

            # calculate dE from element self heating 
            for element in self.elements_to_power:
                element.self_heat(dt,t_cur)

            # apply conservation of energy
            for element in self.elements_to_CoE:
                element.apply_energy_buffer()

            # reset temperatures only for elements with prescribed/forced T(t)
            for element in self.elements_to_prescribe_T:
                element.prescribe_temperature(t_cur)

            iteration += 1
            t_cur += dt

        if(t_cur) >= t_max:
            print("t_max reached, terminating solver")
    
        return self.df



    # ========================================================================================
    # DATA ANALYSIS
    # ========================================================================================

    #TODO
    #text-based visualization
    #3D visualization
    #
    
    def visualize(self, mode="ref"):
        
        # Iterate over each z layer
        for z in range(len(self.elements[0][0])):  
            for y in range(len(self.elements[0])):  
                for x in range(len(self.elements)): 
                    if mode == "ref":
                        print(self.elements[x][y][z].ref, end=" ")  # Print value at (x, y, z)
                    if mode == "T":
                        print(self.elements[x][y][z].T, end="\t")  # Print value at (x, y, z)
                        print(f"{self.elements[x][y][z].T:.1f}", end="\t")  # Print value at (x, y, z) with 1 decimal place

                print()  # New line after each row
            print()  # Blank line after each z layer

        return