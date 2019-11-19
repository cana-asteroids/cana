r"""Utility functions."""

import numpy as np
import pandas as pd


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
