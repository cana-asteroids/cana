r"""Handle Optical Constants."""

import os
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
from ..specdata import SpectralData


def read_constant(filename, label=None, **kwargs):
    r"""Read optical constant file.

    The data should have 3 columns in the format: wavelength, n, k

    Parameters
    ----------
    filename: string
        The path for the optical constant file

    label: string (optional)
        A label to used to identify the optical constant. If None will use the
        file basename as label. Default is None

    grainsize: float (optional)


    Returns
    -------
    Sample

    """
    if isinstance(filename, str):
        data = np.loadtxt(filename, dtype=[('w', np.float64),
                                           ('n', np.float64),
                                           ('k', np.float64)],
                          usecols=(0, 1, 2), **kwargs)
    elif isinstance(filename, np.ndarray):
        data = filename
    if label is None:
        if isinstance(filename, str):
            label = os.path.basename(filename)
        else:
            label = 'sdoc-constant'
    return OpticalConstant(data['w'], data['n'], data['k'], label=label)


def read_constant_batch(filelist, labels=None, **kwargs):
    r"""Read optical constant file.

    The data should have 3 columns in the format: wavelength, n, k

    Parameters
    ----------
    filename: string
        The path for the optical constant file

    label: string (optional)
        A label to used to identify the optical constant. If None will use the
        file basename as label. Default is None

    grainsize: float (optional)


    Returns
    -------
    Sample

    """
    if labels is None:
        out = [read_constant(f, label=None, **kwargs) for i, f in enumerate(filelist)]
    else:
        out = [read_constant(f, label=labels[i], **kwargs) for i, f in enumerate(filelist)]
    return out


@dataclass(repr=False)
class OpticalConstant(SpectralData):
    r"""Optical Constant Class."""

    w: np.ndarray
    n: np.ndarray
    k: np.ndarray
    label: str = field(default='asteroid')

    def __post_init__(self):
        r"""Inicialize class."""
        self.w = np.array(self.w, ndmin=1)
        self.n = np.array(self.n, ndmin=1)
        self.k = np.array(self.k, ndmin=1)
        assert self.w.size == self.n.size
        assert self.w.size == self.k.size
        dtype = ('w', 'n', 'k')
        SpectralData.__init__(self, self.w, dtype, self.label)

    def plot(self):
        r"""Visualization of Sample coeficients."""
        # Just a draft method for now
        _, ax1 = plt.subplots()
        ax1.set_xlabel('Wavelength (microns)', fontsize=14)

        ax1.plot(self.w, self.n, label='Real', color='teal', lw=2)
        ax1.tick_params(axis='y', labelcolor='teal')
        ax1.set_ylabel('Refractive Index (n)', color='teal', fontsize=14)

        ax2 = ax1.twinx()
        ax2.plot(self.w, self.k, label='Imaginary', color='firebrick', lw=2)
        ax2.tick_params(axis='y', labelcolor='firebrick')
        ax2.set_ylabel('Extinction coefficient (k)',
                       color='firebrick', fontsize=14)
        plt.show()
