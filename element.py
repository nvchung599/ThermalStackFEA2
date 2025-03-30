class Element:
    """
    The Element is the fundamental unit of a thermal model mesh, representing a rectangular volume of homogenous material.
    Each element defines material properties, such as thermal conductivity, specific heat, and density, which affect heat
    transfer and temperature changes. Element sizes in the X, Y, and Z dimensions are independent, allowing for high aspect
    ratios (useful for modeling thin substrates or anisotropic materials as isotropic).
    
    Attributes:
    -----------
    x : float
        Size of the element in the x-axis, units: meters.
    y : float
        Size of the element in the y-axis, units: meters.
    z : float
        Size of the element in the z-axis, units: meters.
    eltype : str
        Type of the element (e.g., "null", "solid", "reservoir").
    label : str
        A descriptive label for the material of the element (e.g., "Silicon").
    ref : str
        A reference symbol used in text-based visualizations for this element.
    k : float
        Thermal conductivity of the material, units: W/mK.
    T : float or function
        Initial temperature of the element or a function prescribing temperature over time, units: Kelvin.
    P : float or function
        Heating power generation within the element, units: Watts. Can be constant or a function of time.
    cp : float
        Specific heat capacity of the material, units: J/gK.
    rho : float
        Density of the material, units: kg/m^3.
    rth_x_half : float
        Thermal resistance in the x-direction for half the element, units: K/W.
    rth_y_half : float
        Thermal resistance in the y-direction for half the element, units: K/W.
    rth_z_half : float
        Thermal resistance in the z-direction for half the element, units: K/W.
    capacitance : float
        Heat capacity of the element, units: J/K.
    energy_buffer : float
        Accumulates energy changes during the simulation, applied all at once at the end of each timestep.
    neighbors : list
        Neighboring elements that this element can interact with thermally.
    """

    def __init__(self, x, y, z, eltype, label, ref, k, T, P, cp, rho):
        """
        Initializes an element with its material properties and dimensions.
        
        The thermal resistance values for each axis are calculated based on the element's size and material properties.
        Heating power (P) and temperature (T) can either be constants or time-varying functions.

        Args:
        -----
        x : float
            Size of the element in the x-axis, units: meters.
        y : float
            Size of the element in the y-axis, units: meters.
        z : float
            Size of the element in the z-axis, units: meters.
        eltype : str
            Type of the element (e.g., "null", "solid", "reservoir").
        label : str
            A descriptive label for the material of the element (e.g., "Silicon").
        ref : str
            A reference symbol used in visualizations.
        k : float
            Thermal conductivity of the material, units: W/mK.
        T : float or function
            Initial temperature or a time-varying function for the element, units: Kelvin.
        P : float or function
            Heating power or a time-varying function, units: Watts.
        cp : float
            Specific heat capacity, units: J/gK.
        rho : float
            Density, units: kg/m^3.

        Raises:
        -------
        ValueError:
            If P or T are not numbers or callable functions.
        """
        self.x = x      
        self.y = y 
        self.z = z 
        self.eltype = eltype
        self.label = label
        self.ref = ref

        # Compute thermal resistances based on element size and conductivity in each axis
        self.rth_x_half = (x / (y * z * k)) * 0.5
        self.rth_y_half = (y / (x * z * k)) * 0.5
        self.rth_z_half = (z / (y * x * k)) * 0.5

        # Volumetric heat capacity calculation
        cpv = cp * rho * 1000  # Units: J/m^3K
        self.capacitance = x * y * z * cpv  # Total heat capacity, J/K

        # Initialize heating power (P), which may be constant or a function of time
        self.P_is_constant = isinstance(P, (float, int))
        if self.P_is_constant:
            self.P = P
        elif callable(P):
            self.P_is_constant = False
            self.P_of_t = P
            self.P = self.P_of_t(0)
        else:
            raise ValueError("P must be a number or a function of time")

        # Initialize temperature (T), which may be constant or a function of time
        self.T_is_prescribed = callable(T)
        if self.T_is_prescribed:
            self.T_of_t = T
            self.T = self.T_of_t(0)
        elif isinstance(T, (float, int)):
            self.T = T
        else:
            raise ValueError("T must be a number or a function of time")

        # Buffer to accumulate energy changes
        self.energy_buffer = 0
        self.neighbors = []

        # Coordinates for rendering purposes
        self.x_corner = None
        self.y_corner = None
        self.z_corner = None

    def modify(self, eltype, label, ref, k, T, P, cp, rho):
        """
        Modifies the material properties of the element without changing its dimensions.

        Args:
        -----
        eltype : str
            New element type.
        label : str
            New descriptive label.
        ref : str
            New reference symbol for visualizations.
        k : float
            New thermal conductivity, units: W/mK.
        T : float or function
            New temperature or temperature function, units: Kelvin.
        P : float or function
            New heating power or power function, units: Watts.
        cp : float
            New specific heat capacity, units: J/gK.
        rho : float
            New density, units: kg/m^3.
        
        Returns:
        --------
        None
        """
        self.__init__(self.x, self.y, self.z, eltype, label, ref, k, T, P, cp, rho)

    def apply_energy_buffer(self):
        """
        Applies the accumulated energy buffer to update the element's temperature.

        After all energy transfers are calculated, this method is used to update the temperature
        based on the total accumulated energy and the element's heat capacity.

        Args: None
        Returns: None
        """
        self.T += self.energy_buffer / self.capacitance
        self.energy_buffer = 0

    def accumulate_energy(self, energy):
        """
        Accumulates energy (dE) in the buffer within a timestep.

        Args:
        -----
        energy : float
            Energy to accumulate, units: Joules.
        
        Returns:
        --------
        None
        """
        self.energy_buffer += energy

    def self_heat(self, dt, t=None):
        """
        Accumulates energy from self-heating based on internal power generation over a timestep.

        Args:
        -----
        dt : float
            Timestep, units: seconds.
        t : float, optional
            Current time, used if power is a time-varying function.
        
        Returns:
        --------
        None
        """
        if self.P_is_constant:
            self.accumulate_energy(self.P * dt)
        elif not self.P_is_constant and t is not None:
            self.P = self.P_of_t(t)
            self.accumulate_energy(self.P * dt)
        else:
            raise ValueError("If P is variable, time (t) must be provided.")

    def prescribe_temperature(self, t):
        """
        Overwrites the element's temperature with a prescribed temperature value at the given time.

        Args:
        -----
        t : float
            Current time, units: seconds.
        
        Returns:
        --------
        None
        """
        if self.T_is_prescribed:
            self.T = self.T_of_t(t)

    def get_eltype(self):
        """
        Returns the type of the element.

        Args: None
        Returns:
        --------
        str
            Element type (e.g., "null", "solid", "reservoir").
        """
        return self.eltype

    def get_size(self, axis):
        """
        Returns the size of the element along a specified axis.

        Args:
        -----
        axis : str
            Axis along which to get the size ('x', 'y', or 'z').
        
        Returns:
        --------
        float
            Size along the specified axis, units: meters.
        """
        if axis == "x":
            return self.x
        elif axis == "y":
            return self.y
        elif axis == "z":
            return self.z
        else:
            raise ValueError("Invalid axis. Must be 'x', 'y', or 'z'.")

    def assign_corner_coordinates(self, x, y, z):
        """
        Assigns corner coordinates for the element's position in a 3D model.

        Args:
        -----
        x : float
            X-coordinate of the corner, units: meters.
        y : float
            Y-coordinate of the corner, units: meters.
        z : float
            Z-coordinate of the corner, units: meters.
        
        Returns:
        --------
        None
        """
        self.x_corner = x
        self.y_corner = y
        self.z_corner = z

    def get_corner_coordinates(self):
        """
        Returns the corner coordinates of the element.

        Args: None
        Returns:
        --------
        tuple
            (x, y, z) coordinates of the element's corner, units: meters.
        """
        return (self.x_corner, self.y_corner, self.z_corner)
