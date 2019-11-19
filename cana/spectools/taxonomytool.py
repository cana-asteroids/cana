r"""Tool for spectral taxonomic classification."""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import pandas as pd
from .. import loadspec
from ..util import kwargupdate, find_nearest, Parameter

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))


def _load_demeo(fpath=PWD+'/../datasets/data/taxonomy/bus_demeo.tab'):
    r"""Load the demeo PDS templates file to an record array."""
    # Loading Demeo mean values for each tax classes
    demeo = np.loadtxt(fpath, dtype=[('wavelength', 'f'),
                                     ('A', 'f'), ('std_A', 'f'),
                                     ('B', 'f'), ('std_B', 'f'),
                                     ('C', 'f'), ('std_C', 'f'),
                                     ('Cb', 'f'), ('std_Cb', 'f'),
                                     ('Cg', 'f'), ('std_Cg', 'f'),
                                     ('Cgh', 'f'), ('std_Cgh', 'f'),
                                     ('Ch', 'f'), ('std_Ch', 'f'),
                                     ('D', 'f'), ('std_D', 'f'),
                                     ('K', 'f'), ('std_K', 'f'),
                                     ('L', 'f'), ('std_L', 'f'),
                                     ('O', 'f'), ('std_O', 'f'),
                                     ('Q', 'f'), ('std_Q', 'f'),
                                     ('R', 'f'), ('std_R', 'f'),
                                     ('S', 'f'), ('std_S', 'f'),
                                     ('Sa', 'f'), ('std_Sa', 'f'),
                                     ('Sq', 'f'), ('std_Sq', 'f'),
                                     ('Sr', 'f'), ('std_Sr', 'f'),
                                     ('Sv', 'f'), ('std_Sv', 'f'),
                                     ('T', 'f'), ('std_T', 'f'),
                                     ('V', 'f'), ('std_V', 'f'),
                                     ('X', 'f'), ('std_X', 'f'),
                                     ('Xc', 'f'), ('std_Xc', 'f'),
                                     ('Xe', 'f'), ('std_Xe', 'f'),
                                     ('Xk', 'f'), ('std_Xk', 'f')])
    # Classes std_
    tax_classes = ['A', 'B', 'C', 'Cb', 'Cg', 'Cgh', 'Ch', 'D', 'K', 'L', 'O',
                   'Q', 'R', 'S', 'Sa', 'Sq', 'Sr', 'Sv', 'T', 'V', 'X', 'Xc',
                   'Xe', 'Xk']
    return demeo, tax_classes


def _load_bus(fpath=PWD+'/../datasets/data/taxonomy/bus_demeo.tab'):
    r"""
    Load the demeo PDS templates and trim to bus wavelegth coverage.

    Also excluding classes that does not exist in Bus taxonomy.
    However, note ther it is not real bus templates.
    """
    demeo = _load_demeo()[0]
    # Classes std_
    tax_classes = ['A', 'B', 'C', 'Cb', 'Cg', 'Cgh', 'Ch', 'D', 'K', 'L', 'O',
                   'Q', 'R', 'S', 'Sa', 'Sq', 'Sr', 'T', 'V', 'X', 'Xc', 'Xe',
                   'Xk']
    tax_classes_aux = tax_classes[:]
    tax_classes_aux.extend(['std_'+i for i in tax_classes])
    tax_classes_aux.extend(['wavelength'])
    # aux = tax_classes.extend(tax_classes_aux)
    bus = demeo[tax_classes_aux]
    return bus, tax_classes


def taxonomy(spec, system='demeo', method='chisquared', return_n=3, norm=0.55,
             fitspec=True, speckwargs=None):
    r"""
    Perform taxonomic classification.

    Parameters
    ----------
        spec: Spectrum object
            The spec object in which the classification is being
            applied. It is set as None when the class is iniciated,
            and filled when the classify method is called.

        system: str
            The taxonomic system. Default is 'demeo', for
            a classification in the Bus-Demeo scheme.
            Options are: 'demeo'.

        method: str
            The classification method. Defalt is 'chi-squared'

        fitspec: boolean
            If required to fit the spectrum before interpolating the
            values to the comparison wavelengths. Default is True.

        return_n: int
            The number of classes to output. Default is 1, which will
            output the class with the lowest chi-squared value.

        norm: float
            The normalization point.

    Returns
    -------
    Taxonomic classification.
    """
    tax = Taxonomy(tax=system, norm=norm)
    speckwars_default = {'unit': 'micron'}
    speckwargs = kwargupdate(speckwars_default, speckwargs)
    if not isinstance(spec, list):
        if isinstance(spec, str):
            spec = loadspec(spec, **speckwargs)
        tclass = tax.classify(spec, cmethod=method, return_n=return_n,
                              fitspec=fitspec)
    elif isinstance(spec, list):
        tclass = pd.DataFrame(columns=['tax', 'chi'])
        for fsp in spec:
            sp = loadspec(fsp)
            tclass_aux = tax.classify(sp, cmethod=method, return_n=1,
                                      fitspec=fitspec)
            tclass.loc[tclass_aux.label] = tclass_aux.DataFrame.values[0]

    return tclass


class Taxonomy(object):
    r"""
    Class to handle spectral taxonomic classification.

    Parameters
    ----------
        tax: str
            The taxonomic system. Default is 'demeo', for
            a classification in the Bus-Demeo scheme.
            Options are: 'demeo'.

        dataset: numpy record array
            The mean values for the taxonomic classes.

        norm: float
            The normalization point.

    """

    def __init__(self, tax='demeo', norm=0.55):
        self.system = tax.lower()
        self.dataset, self.classes = self._load()
        self.norm = norm

    def _load(self):
        r"""Load the correspondent dataset of the Taxonomy system."""
        if self.system == 'demeo':
            return _load_demeo()
        if self.system == 'bus':  # -> not real bus templates
            return _load_bus()

    def classify(self, spec, cmethod='chi-squared', return_n=1, fitspec=True):
        r"""
        Classify a spectrum in the defined taxonomic system.

        Parameters
        ----------
        spec: Spectrum object
            The spec object in which the classification is being
            applied. It is set as None when the class is iniciated,
            and filled when the classify method is called.

        cmethod: str
            The classification method. Defalt is 'chi-squared'

        fitspec: boolean
            If required to fit the spectrum before interpolating the
            values to the comparison wavelengths. Default is True.

        return_n: int
            The number of classes to output. Default is 1, which will
            output the class with the lowest chi-squared value.

        Returns
        -------
        numpy record array
        dtype = ('tax', 'chi')
            tax: for the taxonomic classes
            chi: the chi-squared value for the class

        """
        # Selecting regions for comparison
        compspec, tax2comp = self._prep_spec(spec, fitspec)
        # separating single value classes
        chi = [(tcls, self.chisquared(compspec, tax2comp[tcls])) for tcls in
               self.classes]
        # sorting and outputing
        obj_class = np.array(chi, dtype=[('tax', 'U4'), ('chi', 'f')])
        obj_class = np.sort(obj_class, order=['chi'])
        tax = TaxClass(obj_class[:return_n], spec, compspec, tax2comp,
                       self.system, self.norm)
        return tax

    def chisquared(self, spec1, spec2):
        r"""
        Calculate the chi-squared between two spectra.

        The wavelenghts of the two spectra should be the same.

        Parameters
        ----------
        spec1: numpy array
            2D array corresponding to the wavelength and relative reflectance
            of an asteroid

        spec2: numpy array
            2D array corresponding to the wavelength and relative reflectance
            of the taxonomic class

        Returns
        -------
        The chi-squared value

        """
        # normalizing both specs on the secound point
        spec1 = spec1/spec1[1]
        spec2 = spec2/spec2[1]
        # calculating the chi-quared
        n_points = len(spec1)
        chi_aux = np.square(spec1-spec2)/spec1
        chi = np.sum(chi_aux)/n_points
        return chi

    def _prep_spec(self, spec, fitspec=True):
        r"""Prepare asteroid spectrum for comparison with classes templates."""
        # Searching the comparable region between the object and
        # taxonomic classes --> needs improvement
        comparable = np.argwhere((self.dataset['wavelength'] >= spec.w.min()-0.01) &
                                 (self.dataset['wavelength'] <= spec.w.max()+0.01))
        # Trimming the taxonomy dataset to the comparable wavelengths
        tax_comparable = self.dataset[comparable]
        norm_id = find_nearest(tax_comparable['wavelength'], self.norm)[0]
        for tclass in tax_comparable.dtype.names:
                if ('std' not in tclass) and (tclass != 'wavelength'):
                    tax_comparable[tclass] = tax_comparable[tclass] / \
                                             tax_comparable[tclass][norm_id]
        # Checks if desired to autpfit the spectra to resample the
        # spectrum in the comparable wavelengths
        if fitspec:
            _, fcoef = spec.autofit()
            spec_reflectances = np.polyval(fcoef, tax_comparable['wavelength'])
            spec_reflectances = spec_reflectances / spec_reflectances[norm_id]
        # If fispec is False, then will interpolate the values directly from
        # the spectrum ---> needs implementation

        return spec_reflectances, tax_comparable

    def plot_class(self, tclass, tax2comp=None, fax=None, axistitles=True,
                   show=True, legendkwargs=None, taxkwargs=None):
        r"""Plot taxonomic class templates."""
        taxsty_defaults = {'linestyle': '-', 'marker': 'o'}
        legendsty_defaults = {'loc': 'best'}
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        taxkwargs = kwargupdate(taxsty_defaults, taxkwargs)
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = fig.gca()
        if tax2comp is None:
            tax2comp = self.dataset
        # Checking if it is an auxiliary plot
        if isinstance(tclass, str) or isinstance(tclass, np.bytes_):
            fax.plot(tax2comp['wavelength'],
                     tax2comp[tclass],
                     label=tclass, **taxkwargs)
        else:
            for t in tclass:
                fax.plot(tax2comp['wavelength'],
                         tax2comp[t],
                         label=t, **taxkwargs)
        # else ---> implement! view all classes
        fax.legend(**legendkwargs)
        if axistitles:
            plt.xlabel('Wavelength')
            plt.ylabel('Normalized Reflectance')
        if show:
            plt.show()


class TaxClass(Taxonomy, Parameter):
    r"""A taxonomic class representation."""

    def __init__(self, obj_class, spec, compspec, tax2comp, system, norm,
                 label=None):
        r"""
        """
        super(TaxClass, self).__init__(system)
        self.tax = obj_class
        self.tax2comp = tax2comp
        self.spec = spec
        self.compspec = compspec
        if label is None:
            self.label = spec.label
        self.DataFrame = pd.DataFrame(obj_class, columns=['tax', 'chi'])
        self.norm = norm

    def is_primitive(self):
        r"""Check if the closest class is an Primitive class."""
        if self.system == 'demeo':
            if self.tax['tax'][0] in ('B', 'C', 'Cb', 'Ch', 'Cg', 'Cgh',
                                      'X', 'Xc', 'Xe', 'Xk', 'D'):
                return True
            else:
                return False

    def plot(self, fax=None, show=True, savefig=None, axistitles=True,
             speckwargs=None, legendkwargs=None, dotskwargs=None,
             taxkwargs=None):
        r"""
        Plot taxonomic classification.

        Parameters
        ----------
        fax (Optional): matplotlib.axes
            If desired to subplot image in a figure. Default is 'None', which
            will open a new plt.figure()

        show (Optional): boolean
            True if want to plt.show(). Default is True.

        savefig (Optional): str
            The path to save the figure. If set to None, wont save the figure.
            Default is None

        specparams: dict
            Arguments for matplotlib plot function, to be applied in the
            spectrum plot.

        fitparams: dict
            Arguments for matplotlib plot function, to be applied in the
            fitted spectrum plot (used for the classes comparison).

        tclassparams: dict
            Arguments for matplotlib plot function, to be applied in the
            plot for the classes.

        """
        specsty_defaults = {'c': '0.5', 'zorder': 0, 'label': 'Raw Spectrum'}
        specpointssty_defaults = {'c': 'darkred', 'label': 'Spectrum Fit',
                                  'linestyle': '-', 'marker': 'o'}
        taxsty_defaults = {'linestyle': '-', 'marker': 'o'}
        legendsty_defaults = {'loc': 'best'}

        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        taxkwargs = kwargupdate(taxsty_defaults, taxkwargs)
        dotskwargs = kwargupdate(specpointssty_defaults, dotskwargs)

        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = fig.gca()
        # Ploting the spec
        # if self.norm is None:
            # self.norm = find_nearest(self.dataset['wavelegth'], 0.55)[1]
        self.spec = self.spec.normalize(self.norm, window=0.03)
        self.spec.plot(fax=fax, axistitles=axistitles, show=False,
                       speckwargs=speckwargs)
        fax.scatter(self.tax2comp['wavelength'], self.compspec, **dotskwargs)
        # plotting the classes
        # generating the colors
        cmap = get_cmap('viridis')
        color_aux = np.linspace(0, 1, len(self.tax['tax']))
        for tindex, tcl in enumerate(self.tax):
            taxkwargs['c'] = cmap(color_aux[tindex])
            self.plot_class(tcl[0], tax2comp=self.tax2comp, fax=fax,
                            show=False, taxkwargs=taxkwargs)
        fax.legend(prop={'size': 8})
        # check if save the image
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
        # show in the matplotlib window?
        if show:
            plt.show()
