from element import *

class Node:

    """
    A Node is a pair of element references.
    A list of such nodes describes all valid-only element interactions that must be computed every timestep.
    Nodes are used because it cuts out the extra computation time of re-observing invalid element interaction pairs
    """

    def __init__(self, el1:Element, el2:Element, orientation):
        self.el1 = el1
        self.el2 = el2
        self.rth = None

        if orientation == "x":
            self.rth = el1.rth_x_half + el2.rth_x_half
        elif orientation == "y":
            self.rth = el1.rth_y_half + el2.rth_y_half
        elif orientation == "z":
            self.rth = el1.rth_z_half + el2.rth_z_half
        else:
            raise ValueError("invalid orientation argument")

    def calculate_heattransfer(self, dt):
        """
        computes dE for one timestep for both elements linked by this node
        dE is storred in elements' buffers, awaiting more energy transfers, and eventually a dE->temperature change calculation

        Args: dt, timestep, s

        Returns: None
        """
        q_dot = (self.el2.T - self.el1.T) / self.rth
        q = q_dot * dt
        self.el1.accumulate_energy(q)
        self.el2.accumulate_energy(-q)