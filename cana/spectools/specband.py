r"""Tool for measuring an absorptium band."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .uncertainties import SpecError
from .. import loadspec, Spectrum
from ..util import kwargupdate, find_nearest, Parameter


def depth(spec, wmin=0.54, wmax=0.88, cont_window=0.04, resolution='auto',
          errormethod='rms', error_param=None, montecarlo=1000,
          min_depth=1., theoric_min=0.7, max_dist=0.05, n_sigma=3,
          speckwargs=None):
    r"""
    Calculate the depth of an absorptium band.

    Parameters
    ----------
    spec: Spectrum, spectrum file, spectrum file list
        The imput can be

    wmin: float
        The wavelength for the beginning of the band.

    wmax: float
        The wavelength for the beginning of the band

    cont_window: float
        The wavelength window size for measuring the continuum

    resolution: 'auto' or int
        The spectrum wavelength resolution for measuring the center and depth
        of the band. If 'auto' will take the spec input resolution.

    errormethod: 'rms', 'removal' or 'bin'
        The error methodology that will be applied for estimating the
        depth error. Default is 'rms'

    error_param: None or float
        The error methodoly parameter if needed. If errormethod='rms',
        then no value is necessary. If it is set for removal, the percentage
        of points to remove. For rebin, the param represents the binsize.

    montecarlo: integer
        Default is 1000

    min_depth: float (optional)
        minimum depth for the band to be considered.
        Default is 1 percent (for hydration analysis).

    theoric_min: float (optional)
        theorical central wavelength of the band.
        Default is 0.7 (for hydration analysis)

    max_dist: (optional)
        Maximal distance from the therorical position for the
        band to be considered.
        Default is 0.05  (for hydration analysis).

    n_sigma: integer (optional)
        The minimum sigma level for the absorptium band detection.
        Default is 3.

    Returns
    -------
    DepthValue or Pandas.DataFrame
        For a single spectrum it will return a Depthvalue, for a list
        of spectra, returns a pandas.DataFrame with the results

    """
    error = SpecError(n=montecarlo, method=errormethod, param=error_param)
    cont = Continuum(lowerwindow=cont_window, upperwindow=cont_window)
    speckwars_default = {'unit': 'microns'}
    speckwargs = kwargupdate(speckwars_default, speckwargs)
    band = Depth(wmin, wmax, continuum=cont)
    if isinstance(spec, Spectrum):
        depth = band.measure(spec, resolution=resolution, error=error)
    elif isinstance(spec, str):
        spec = loadspec(spec, **speckwargs)
        depth = band.measure(spec, resolution=resolution, error=error)
    elif isinstance(spec, list):
        aux = []
        for spfile in spec:
            sp = loadspec(spfile)
            dth = band.measure(sp, resolution=resolution, error=error)
            if not dth.is_band():
                dth.DataFrame[dth.label] = ['-', '-', '-', '-']
            aux.append(dth.DataFrame.T)
        depth = pd.concat(aux)
    return depth


class Continuum(object):
    r"""
    Model for fitting and subtracting the continum.

    Methods
    -------
        fit
        remove
        plot_continuum
        plot_continuum_region


    Atributes:
    ----------
    cont_arr: numpy array
        2D numpy array of the fitted Continuum,
        in the same wavelenghs as the inputed spectrum.
        The cont_arr is set as None when the class is inicialized and
        will be rewriten when the 'fit' or 'remove' methods are called.

    """

    def __init__(self, lowerwindow=0.04, upperwindow=0.04):
        r"""
        Initialize Continuum fit model.

        Parameters
        ----------
            lowerwindow: float
                The wavelength windows for the region where the continuum
                will be fitted. E.G. if the AbsorptiumBand.wmin=0.55 and the
                lowerwindow=0.3, then the region 0.55-0.58 microns will be
                the lower window for the Continuum fitting.

            upperwindow: float
                The wavelength windows for the region where the continuum
                will be fitted. E.G. if the AbsorptiumBand.wmax=0.85 and the
                upperwindow=0.3, then the region 0.83-0.85 microns will the
                upper window for the Continuum fitting.

        """
        self.lowerwindow = lowerwindow
        self.upperwindow = upperwindow
        self.cont_arr = None
        self.bcoefs = None

    def fit(self, spec):
        r"""
        Fit the continuum.

        Parameters
        ----------
        spec_arr: numpy array
            The 2D (wavelength and Normalized reflectance) array

        Returns
        -------
        Continuum 2D numpy array, in the same wavelenghs as the input.

        """
        # Trimming the continum region
        lower_cont = spec.T[(spec.w < spec.w.min()+self.lowerwindow)].T
        upper_cont = spec.T[(spec.w > spec.w.max()-self.upperwindow)].T
        cont_region = np.hstack([lower_cont, upper_cont])
        self.bcoefs = np.polyfit(cont_region['w'], cont_region['r'], 1)
        cont_y = np.polyval(self.bcoefs, spec.w)
        self.cont_arr = np.array([spec.w, cont_y])
        return self.cont_arr, self.bcoefs

    def remove(self, spec):
        r"""
        Remove the continuum.

        Parameters
        ----------
        spec_arr: numpy array
            The 2D (wavelength and Normalized reflectance) array

        Returns
        -------
        2D numpy array of the spec_arr, with the continuum removed.

        """
        # Fitting continuum
        continuum = self.fit(spec)
        # Removing continuum
        spec.r = np.divide(spec.r, continuum[0][1])
        return spec

    def plot_continuum_region(self, spec, fax, cont_style=None):
        r"""
        Add the continuum line in a plot.

        Parameters
        ----------
        fax: matplotlib axes
            The current axes for ploting the continuum

        contstyle: None or dict
                Matplotlib arguments for plot function.
                If none, default values are:{'c':'k', 'linestyle':'--'}

        """
        if cont_style is None:
            cont_style = {'color': '0.7', 'linestyle': '--', 'zorder': 0}
        fax.axvspan(spec.w.min(), spec.w.min()+self.lowerwindow, **cont_style)
        fax.axvspan(spec.w.max(), spec.w.max()-self.upperwindow, **cont_style)

    def plot_continuum(self, fax, cont_style=None):
        r"""
        Add the continuum line in a plot.

        Parameters
        ----------
        fax: matplotlib axes
            The current axes for ploting the continuum

        contstyle: None or dict
                Matplotlib arguments for plot function.
                If none, default values are:{'c':'k', 'linestyle':'--'}

        """
        if cont_style is None:
            cont_style = {'c': 'k', 'linestyle': '--'}
        fax.plot(self.cont_arr[0], self.cont_arr[1], **cont_style)


class Depth(object):
    r"""
    Model for estimating the depth and center of a spectral absorptium band.

    Methods
    -------
    measure

    """

    def __init__(self, wmin=0.54, wmax=0.88, continuum=Continuum()):
        r"""
        Initialize Absorptium band model.

        Parameters
        ----------
        wmin: float
            The wavelength for the beginning region to estimate the band.

        wmax: float
            The wavelength for the end of the region to estimate of the band

        continuum: Continuum
            The object containing the information about the continuum fitting

        """
        self.wmin = wmin
        self.wmax = wmax
        self.cont = continuum

    def measure(self, spec, error=SpecError(), label=None, resolution='auto'):
        r"""
        Measure the absorptium band in a spectra.

        The error is estimated using an errormodel.

        Parameters
        ----------
        spec: Spectrum
            The spectrum object

        error: SpecError
            The model to estimate the slope uncertainty.

        label: None or str
            The spectrum label (for the output)

        resolution: 'auto' or int
           The spectrum wavelength resolution for measuring the center and
           depth of the band. If 'auto' will take the spec input resolution.

        Returns
        -------
            DepthValue

        """
        # Trimming spec to the band region
        bspec = spec.trim(self.wmin, self.wmax)
        # measuring the band using the error model
        band_aux = error.estimate(bspec, self._measure_band, error.n,
                                  resolution=resolution)
        # taking the mean and the std for both values
        minpos = (band_aux[0][0], band_aux[1][0])
        depth = (band_aux[0][1], band_aux[1][1])
        # formating
        minpos = (np.around(minpos[0], 4), np.around(minpos[1], 4))
        depth = (np.around(depth[0], 3), np.around(depth[1], 3))
        band = DepthValue(bspec, minpos, depth, cont=self.cont, label=label)
        return band

    def _measure_band(self, spec, resolution):
        r"""Auxialiary function for measuring the absorptium band."""
        # fitting the band
        fspec, fcoefs = spec.fit(order=4, ftype='spline')
        if isinstance(resolution, int):
            x = np.linspace(spec.w.min(), spec.w.max(), resolution)
            ref = fcoefs(x)
            fspec = Spectrum(x, ref)
        # removing continuum
        fspec_wo_cont = self.cont.remove(fspec)
        # finding the minimum
        band_min_index = fspec_wo_cont.r.argmin()
        # characterizing the band
        band_min = fspec_wo_cont.w[band_min_index]
        band_depth = float(1-fspec_wo_cont.r[band_min_index])*100
        return band_min, band_depth


class DepthValue(Depth, Parameter):
    r"""
    Representation of the a band depth and center value.

    Methods
    -------
    is_band
    plot

    """

    def __init__(self, spec, center, depth, cont, label=None):
        r"""
        Initialize class to vizualization of Depth value.

        Parameters
        ----------
        spec: Spectrum
            The spectrum information

        center: (float, float)
            The value and error of the band central position.

        depth: (float, float)
            The value and error of the band depth.

        cont: Continuum
            The continuum information.

        label: string
            The spectrum label. If None, it will get the spec.label.

        """
        self.spec = spec
        self.center = center
        self.depth = depth
        self.DataFrame = self._build_dataframe()
        self.label = label
        if label is None:
            self.label = spec.label
        self.DataFrame[self.label] = [depth[0], depth[1], center[0], center[1]]
        self.cont = cont

    @staticmethod
    def _build_dataframe():
        hcolumns = ['depth', 'depth_unc', 'center', 'center_unc']
        dataframe = pd.DataFrame(columns=hcolumns).T
        return dataframe

    def is_band(self, min_depth=1., theoric_min=0.7, max_dist=0.05,
                sigma=3):
        r"""
        Ask if the mesuared parameters can be considered an absorptium band.

        Parameters
        ----------
        min_depth: float (optional)
            minimum depth for the band to be considered.
            Default is 1 percent (for hydration analysis).

        theoric_min: float (optional)
            theorical central wavelength of the band.
            Default is 0.7 (for hydration analysis)

        max_dist: (optional)
            Maximal distance from the therorical position for the
            band to be considered.
            Default is 0.05  (for hydration analysis).

        n_sigma: integer (optional)
            The minimum sigma level for the absorptium band detection.
            Default is 3.

        Returns
        -------
        boolean

        """
        # check if minimum is in within distance
        if self.center[0] <= theoric_min:
            dist = (abs(self.center[0] - self.center[1] - theoric_min)
                    < max_dist)
        else:
            dist = (abs(self.center[0] + self.center[1] - theoric_min)
                    < max_dist)
        # check if band depth is higher than noise
        # measuring snr
        rms = self.spec.estimate_rms()
        # Restablishing the band depth unit to compare with snr
        band_snr = self.depth[0]-sigma*self.depth[1] > 100*rms
        # checking if band depth higher than a limit
        depth = (self.depth[0] > min_depth)
        # now checking the combined result
        ans_aux = [dist, band_snr, depth]
        ans = all(ans_aux)
        return ans

    def plot(self, fax=None, show=True, savefig=None,
             axistitles=True, speckwargs=None,
             bandkwargs=None, contkwargs=None,
             dotkwargs=None, legendkwargs=None):
        r"""
        Plot for the band.

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
        legendsty_defaults = {'loc': 'best', 'prop': {'size': 8}}
        dotsty_defaults = {'c': 'steelblue', 's': 70, 'zorder': 2}
        bandsty_defaults = {'c': 'steelblue', 'lw': '2', 'zorder': 1}
        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        bandkwargs = kwargupdate(bandsty_defaults,  bandkwargs)
        dotkwargs = kwargupdate(dotsty_defaults,  dotkwargs)
        # Ploting the spec
        # remove continuum
        spec = self.cont.remove(self.spec)
        spec.plot(fax=fax, axistitles=axistitles, show=False,
                  speckwargs=speckwargs)
        # ploting the band
        fspec = self.spec.fit(order=4, ftype='spline')[0]
        fax.plot(fspec['w'], fspec['r'], **bandkwargs)
        # ploting the continuum
        fax.axhline(y=1, c='k')
        self.cont.plot_continuum_region(self.spec, fax=fax)
        # plotting minimum
        min_original_aux = find_nearest(fspec.w, self.center[0])[0]
        min_original_ref = fspec.r[min_original_aux]
        label = """band minimum: {0}$\pm${1} \n band depth: {2}$\pm${3}
                """.format(self.center[0], self.center[1], self.depth[0],
                           self.depth[1])
        fax.scatter(self.center[0], min_original_ref, label=label, **dotkwargs)
        fax.legend(**legendkwargs)
        # check if save the image
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
        # # show in the matplotlib window?
        if show:
            plt.show()
