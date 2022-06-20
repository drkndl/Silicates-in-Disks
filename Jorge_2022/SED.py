import numpy as np
import matplotlib.pyplot as plt


def Plancks(T, lamda):
    """
    Calculates Planck function for a particular temperature
    """
    I = 2*h*c**2 / (lamda**5 * np.exp( h*c/(lamda*k*T) ) - 1)
    return I

