r"""Core class to handle Spectral Data."""

import copy
import numpy as np


class SpectralData:
    r"""Base Class to handle Spectral Data."""

    def __init__(self, w, dtype, label, **kwargs):
        r"""Initialize the class.

        Parameters
        ----------
        w: np.ndarray
            The wavelentgh axis

        dtype: list or tuple
            The columns names of the Spectral data

        label: str
            The identifier of the data

        """
        self.w = w
        self.label = label
        self.dtype = dtype
        self.repr = self.make_repr()
        self.kwargs = kwargs

    def __getitem__(self, item):
        r"""Get item magic method."""
        if isinstance(item, (int, slice)):
            aux = self.__class__(*[self[d][item] for d in self.dtype],
                                 label=self.label, **self.kwargs)
        elif isinstance(item, (list, np.ndarray)):
            aux = self.__class__(*[self[d][item] for d in self.dtype],
                                 label=self.label, **self.kwargs)
        else:
            aux = self.__getattribute__(item)
        return aux

    def make_repr(self, maxlines=20, nlines=4):
        r"""Make representation of the Spectral Data."""
        basestr = self.__class__.__name__ + '(['
        whitespace = ' ' * len(basestr)
        if len(self) > maxlines:
            idss = list(range(nlines))
            idss.extend([-i-1 for i in idss[::-1]])
        else:
            idss = range(len(self))
        for n, i in enumerate(idss):
            line = ', '.join(map(str, [self[idx][i] for idx in self.dtype]))
            if i == 0:
                basestr += '('+line+'),\n'
            elif n == len(idss)-1:
                basestr += whitespace+'('+line+')'
            else:
                if (n == nlines) & (len(self) > maxlines):
                    basestr += whitespace + ' ...\n'
                basestr += whitespace + '(' + line + '),\n'
        basestr += '], \n'
        basestr += whitespace + 'columns=(' + ', '.join(self.dtype) + '),\n'
        basestr += whitespace + 'label="{0}"'.format(self.label)
        basestr += ')'
        return basestr

    def __repr__(self):
        r"""Representation of a Spectral Data."""
        return self.repr

    @property
    def shape(self):
        r"""Return the shape of the Spectral Data."""
        return (len(self.w), len(self.dtype))

    def __len__(self):
        r"""Return the length of the Spectral Data."""
        return len(self.w)

    def copy(self):
        r"""Return a copy of the object."""
        return copy.copy(self)

    def sort(self, order='w'):
        r"""Sort the Spectral data.

        Parameters
        ----------
        order: str
            The name of the
        """
        aux = np.argsort(self[order])
        return self.__class__(*[self[d][aux] for d in self.dtype],
                              label=self.label, **self.kwargs)

    def trim(self, w_low=1.0, w_up=2.3):
        r"""Trim to the desired wavelength range.

        Parameters
        ----------
        w_low: list [begin, end]
            Lower wavelength value to trim the data

        w_up: list [begin, end]
            Lower wavelength value to trim the data
        """
        aux = np.where((self.w > w_low) & (self.w < w_up))[0]
        args = [self[d][aux] for d in self.dtype]
        return self.__class__(*args, label=self.label, **self.kwargs)

    def rebase(self, baseaxis):
        r"""Interpolate the Spectral Data to a new wavelentgh axis.

        Parameters
        ----------
        basename: array
            New wavelength array

        Returns
        -------
        The interpolated object

        """
        args = [np.interp(baseaxis, self.w, self[d]) for d in self.dtype[1:]]
        return self.__class__(*args, label=self.label, **self.kwargs)
