from simframe.utils.color import colorize

import numpy as np


class Boundary(object):
    '''Class for managing boundary conditions for gas and dust.'''

    def __init__(self, r, ri, S, condition=None, value=None):
        """Class that manages boundary conditions.

        Parameters
        ----------
        r : Field
            Radial grid
        ri : Field
            Radial grid cell interfaces
        S : Field
            Field on which the boundary condition shall be imposed
        condition: string or None, optional, default : None
            Type of boundary condition
            - None : No boundary condition imposed (default)
            - "const_grad" : Constant gradient
            - "const_pow" : constant power law
            - "const_val" : Constant value
            - "grad" : Custom gradient
            - "pow" : Custom power law with set exponent
            - "val" : Custom value

        Notes
        -----
        r and S have to be passed such that the boundary is at the lower end, i.e.,
        r[0] and S[0] are the boundary. If you want to modify the outer boundary
        you have to pass e.g. ``Boundary(sim.grid.r[::-1], sim.gas.Sigma[::-1].``"""
        self.setcondition(condition, value)
        self._r = r[:3]
        self._ri = ri[:3]
        self._S = S[:3]

    def __repr__(self):
        """Returns a meaningful description.

        Returns
        -------
        repr : string
            Description of boundary condition"""
        if self._condition == None:
            text = "No boundary condition set"
        elif self._condition == "const_val":
            text = "Constant value"
        elif self._condition == "const_pow":
            text = "Constant power law"
        elif self._condition == "const_grad":
            text = "Constant gradient"
        elif self._condition == "val":
            text = "Value"
        elif self._condition == "grad":
            text = "Gradient"
        elif self._condition == "pow":
            text = "Power law with set exponent"
        ret = "{}".format(text)
        return ret

    @property
    def condition(self):
        '''The type of boundary condition to be applied.'''
        return self._condition

    @condition.setter
    def condition(self, val):
        msg = "{} Use <Boundary>.setcondition() to set boundary condition.".format(
            colorize("Warning:", color="yellow"))
        print(msg)

    @property
    def value(self):
        '''The value of the boundary to be applied if applicable.'''
        return self._value

    @value.setter
    def value(self, value):
        self._value = value

    def setcondition(self, condition, value=None):
        """Function to set boundary condition.

        Parameters
        ----------
        condition : string
            Type of boundary conditon:
                - "const_grad" : constant gradient
                - "const_pow" : constant power law
                - "const_val" : constant value
                - "val" : custom value
                - "grad" : custom gradient
                - "pow" : custom power law with set exponent
                - None : Don't impose boundary condition (default)
        value : float or array, optional, default : None
            Value if needed for boundary condition"""
        if condition in [None, "const_val", "const_pow", "const_grad", "grad", "val", "pow"]:
            if condition in ["val", "grad", "pow"] and value is None:
                msg = "You have to give a value for condition '{:}'.".format(
                    condition)
                raise ValueError(msg)
            self._condition = condition
            self._value = value
        else:
            raise ValueError(
                "Unknown boundary condition '{}'.".format(condition))

    def _getboundary(self):
        """Function returns the calculated value at the boundary.

        Returns
        -------
        vb: float or array
            value at boundary"""
        # Helper
        r = self._r
        ri = self._ri
        S = self._S
        if self._condition is None:
            return self._S[0]
        elif self._condition == "const_grad":
            D = ri[1]/ri[2] * (r[1]-r[0]) / (r[2]-r[1])
            return (1. + D)*r[1]/r[0]*S[1] - r[2]/r[0]*D*S[2]
        elif self._condition == "const_pow":
            p = np.log(S[2]/S[1]) / np.log(r[2]/r[1])
            return S[1]*(r[0]/r[1])**p
        elif self._condition == "const_val":
            return self._S[1]
        elif self._condition == "pow":
            return S[1]*(r[0]/r[1])**self._value
        elif self._condition == "grad":
            return r[1]/r[0]*S[1] - self._value*ri[1]/r[0]*(r[1]-r[0])
        elif self._condition == "val":
            return self._value

    def setboundary(self):
        """Function sets the boundary value."""
        vb = self._getboundary()
        if vb is not None:
            self._S[0] = self._getboundary()
