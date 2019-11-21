r"""Tool for measuring spectral slope."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from .uncertainties import SpecError
from .. import loadspec, Spectrum
from ..util import kwargupdate, Parameter


def slope(spec, wmin=0.4, wmax=0.9, norm=0.55, errormethod='rms',
          error_param=None, montecarlo=1000, speckwargs=None):
    r"""
    Calculate the spectral gradient.

    Parameters
    ----------
    spec: Spectrum, spectrum file, spectrum file list
        The input can be a Spectrum object, a spectrum file or
        a list of spectrum files

    wmin: float
        wavelength lower limit for the adjust.

    wmax: float
        wavelength upper limit for the adjust.

    norm: float
        The wavelength for normalizing the slope

    errormethod: 'rms', 'removal' or 'bin'
        The error methodology that will be applied for estimating the
        slope error. Default is 'rms'

    error_param: None or float
        The error methodoly parameter if needed. If errormethod='rms',
        then no value is necessary. If it is set for removal, the percentage
        of points to remove. For rebin, the param represents the binsize.

    montecarlo: integer
        Default is 1000

    Returns
    -------
    SlopeValue or Pandas.DataFrame
        For a single spectrum it will return a SlopeValue, for a list
        of spectra, returns a pandas.DataFrame with the results

    """
    slopemodel = Slope(wmin=wmin, wmax=wmax, norm=norm)
    speckwars_default = {'unit': 'microns'}
    speckwargs = kwargupdate(speckwars_default, speckwargs)
    error = SpecError(n=montecarlo, method=errormethod, param=error_param)
    if isinstance(spec, Spectrum):
        slp = slopemodel.measure(spec, error=error)
    elif isinstance(spec, str):
        spec = loadspec(spec, **speckwargs)
        slp = slopemodel.measure(spec, error=error)
    elif isinstance(spec, list):
        aux = []
        for spfile in spec:
            sp = loadspec(spfile)
            slp_aux = slopemodel.measure(sp, error=error)
            aux.append(slp_aux.DataFrame.T)
        slp = pd.concat(aux)
    return slp


class Slope(object):
    r"""
    Calculate  the spectral gradient of a spectrum.

    Methods
    -------
    measure

    """

    def __init__(self, wmin=0.4, wmax=0.9, norm=0.55):
        r"""
        Initialize the Slope model.

        Attributes
        ----------
            wmin: float
                wavelength lower limit for the adjust.

            wmax: float
                wavelength upper limit for the adjust.

            norm: float
                Normalization wavelength

        """
        self.wmin = wmin
        self.wmax = wmax
        self.norm = norm

    def measure(self, spec, error=SpecError(), label=None):
        r"""
        Calculate the slope in a spectral region.

        If the "error" is set as True, then uses a Monte-Carlo model
        to estimate the error in the slope.

        Parameters
        ----------
        spec: Spectrum
            The spectrum object

        error: SpecError
            The model to estimate the slope uncertainty.

        label: str (optional)
            The label of the spectrum. This is just to organize the output.

        Returns
        -------
        The slope value. If error==True, then also returns the
        uncertainty.

        """
        # Trimming the spectrum in the defined region
        slspec = spec.trim(self.wmin, self.wmax)
        # calculating the slope
        # slp, fspec = self._calc_slope(slspec)
        # filling class atributes
        if error is None:
            slp_value = self._calc_slope(slspec)
            slp_unc = None
        else:
            slp_value, slp_unc = error.estimate(slspec, self._calc_slope,
                                                error.n, error.param)
#            slp_value = np.mean(slp_aux)
#            slp_unc = np.std(slp_aux)
        slp = SlopeValue(self, slspec, slp_value, slp_unc=slp_unc, label=label)
        return slp

    def _calc_slope(self, spec):
        r"""Calculate the slope."""
        # fitting the spectra with first order polynomial
        fspec, fcoefs = spec.fit(order=1, ftype='polynomial')
        # Normalizing and calculating the slope
        norm_reflectance = np.polyval(fcoefs, self.norm)
        if spec.unit == 'micron':
            factor = 10
        if spec.unit == 'angstrom':
            factor = 100000
        slp = (fcoefs[0]/norm_reflectance) * factor
        return slp


class SlopeValue(Slope, Parameter):
    r"""
    Representation of a slope mesurement.

    Methods
    -------
    plot

    """

    def __init__(self, model, spec, slp, slp_unc=None, label=None):
        r"""
        Initialize the class.

        Parameters
        ----------
        model: Slope class
            The slope model

        spec: Spectrum
            The spectrum from which the slope was measured

        slp: float
            The measured slope

        slp_unc: float
            The slope uncertainty

        label: str
            The spectrum label, for outputing value

        """
        self.spec = spec
        self.slope = slp
        self.slope_unc = slp_unc
        self.DataFrame = pd.DataFrame(columns=['slope', 'slope_unc']).T
        self.label = label
        if label is None:
            self.label = spec.label
        self.DataFrame[self.label] = [slp, slp_unc]
        self.model = model

    def plot(self, fax=None, show=True, savefig=None,
             axistitles=True, speckwargs=None,
             slopekwargs=None, legendkwargs=None):
        r"""
        Plot the Spectrum and the measured slope.

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

        **kwargs: matplotlib plot kwargs

        """
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = fig.gca()
        # setting default values for image plot with matplotlib
        specsty_defaults = {'c': '0.3', 'lw': '1', 'zorder': 0}
        legendsty_defaults = {'loc': 'best'}
        label = '{0} $\pm$ {1}'.format(self.slope, self.slope_unc)
        slopesty_defaults = {'c': 'r', 'lw': '2', 'zorder': 1, 'label': label}
        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        slopekwargs = kwargupdate(slopesty_defaults, slopekwargs)
        # Ploting the spec
        self.spec = self.spec.normalize(self.model.norm, window=0.03)
        self.spec.plot(fax=fax, axistitles=axistitles, show=False,
                       speckwargs=speckwargs)
        # ploting the slope
        fspec, _ = self.spec.fit(order=1, ftype='polynomial')
        fspec = fspec.normalize(self.model.norm, window=0.03)
        fax.plot(fspec.w, fspec.r, **slopekwargs)
        # plot legend?
        if 'label' in slopekwargs:
            fax.legend(**legendkwargs)
        # check if save the image
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            matplotlib.use('TkAgg')
        # # show in the matplotlib window?
        if show:
            plt.show()
