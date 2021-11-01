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
        # self.repr = self.make_repr()
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
        return self.make_repr()

    @property
    def shape(self):
        r"""Return the shape of the Spectral Data."""
        return (len(self.w), len(self.dtype))

    def __len__(self):
        r"""Return the length of the Spectral Data."""
        return len(self.w)

    def copy(self, deep=True):
        r"""Return a copy of the object."""
        if deep:
            return copy.deepcopy(self)
        else:
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
        r"""Interpolate spectra to a new axis.

        Parameters
        ----------
        baseaxis: np.array
            The new wavelength axis for the Spectrum
        """
        args = [np.interp(baseaxis, self.w, self[d]) for d in self.dtype[1:]]
        new = self.__class__(baseaxis, *args, label=self.label,
                             **self.kwargs)
        return new

    def save(self, fname, fmt='%.5f', delimiter=' ', header='', footer='',
             comments='#', encoding=None):
        r"""
        Save the spectrum data into file.

        Parameters
        ----------
        fname : filename or file handle
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format.
            loadspec understands gzipped files transparently.

        fmt : str or sequence of strs, optional

        delimiter : str, optional
            String or character separating columns.

        newline : str, optional
            String or character separating lines.

        header : str, optional
            String that will be written at the beginning of the file.

        footer : str, optional
            String that will be written at the end of the file.

        comments : str, optional
            String that will be prepended to the header and footer strings,
            to mark them as comments.

        encoding : {None, str}, optional
            Encoding used to encode the outputfile. Does not apply to output
            streams.

        """
        arr = np.array([self[a] for a in self.dtype]).T
        np.savetxt(fname, arr, fmt=fmt, delimiter=delimiter, header=header,
                   footer=footer, comments=comments, encoding=encoding)
