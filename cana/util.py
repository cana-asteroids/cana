r"""Utility functions."""

from textwrap import dedent
from scipy.stats import truncnorm
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline


def verboseprint(text, underline=False, bold=False):
    r"""Adjust the style of the print."""
    text = text.splitlines()
    for line in text:
        print(line.strip())
        if underline:
            print('-'*len(line))
        # text = PrintStyle.UNDERLINE+text+PrintStyle.UNDERLINE
    # if bold:
        # text = PrintStyle.BOLD+text+PrintStyle.BOLD
    # print(text)


def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    r"""Return truncaded normal."""
    norm = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
    return norm


def find_nearest(array, value):
    r"""
    Find nearest value in array.

    Parameters
    ----------
    array : 1D float numpy array
        Array for the search
    value: float
        Value to be searched

    Returns
    -------
    idx , array[idx]: int, float
        Index and value of nearest point

    """
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]


def kwargupdate(defaults, options):
    r"""
    Modify or extend a kwarg dict.

    Parameters
    ----------
    defauts: dict
        The dictionary containing the default values

    options: dict
        The dictionary containing the modified or new values

    Returns
    -------
    A merged kwarg dictionary

    """
    if options is None:
        options = {}
    for key, value in defaults.items():
        if key not in options.keys():
            options[key] = value
    return options


def curvature(func, x, ftype='polynomial', order=4):
    r"""
    Measure the curvature of a function.

    Parameters
    ----------
    func: coefcients, spline or function
        The function for deriving the curvature. F(x)

    x:
        The x axis for the F(x) function.

    ftype: string
        'polynomial': func is the polynomial coefcients arrays
        'spline': func is a scipy.UnivariateSpline object
        'analytical': a self defined function

    order: int (optional)
        The order for fitting the analytical function.
        Only used if ftype='analytical'

    Returns
    -------
    The curvature array

    """
    if ftype == 'polynomial':
        # getting polynomial derivatives
        func_ = np.polyder(func)
        func_2 = np.polyder(func, m=2)
        # building arrays with derived polynomial
        y_ = np.polyval(func_, x)
        y_2 = np.polyval(func_2, x)
    if ftype == 'analytical':
        # fitting analytical function with a spline
        spl = UnivariateSpline(x, func(x), k=order)
        ftype = 'spline'
    if ftype == 'spline':
        # deriving function
        func_ = spl.derivative()
        func_2 = spl.derivative(n=2)
        y_ = func_(x)
        y_2 = func_2(x)
    # normalizing vectors according to the first derivative
    y_norm = np.median(y_)
    y_ = y_/y_norm
    y_2 = y_2/y_norm
    # building terms of the curvature
    aux = np.power((np.power(y_, 2) + 1), 1.5)
    aux2 = y_2
    # calculating curvature and radius of curvature
    r = aux/aux2
    k = 1/r
    return k


class Parameter(object):
    """

    """
    def __init__(self, value, uncertainty):
        r"""
        """
        self.DataFrame = pd.DataFrame()

    def __getitem__(self, item):
        assert isinstance(self.DataFrame, pd.core.frame.DataFrame)
        if item in self.DataFrame.columns:
            return self.DataFrame[item]
        elif item in self.DataFrame.index:
            return self.DataFrame.loc[item]
        else:
            raise ValueError('could not find %s' % (item))

    def to_latex(self):
        return self.DataFrame.to_latex()

    def __repr__(self):
        return self.DataFrame.__repr__()

    def to_csv(self, fname, sep=' '):
        self.DataFrame.to_csv(fname, sep=sep)
