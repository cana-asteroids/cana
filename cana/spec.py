r"""Tools to handle asteroid spectra."""

from os.path import splitext, basename
from dataclasses import dataclass, field
from scipy.interpolate import UnivariateSpline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error
import numpy as np


import matplotlib.pyplot as plt
import matplotlib
from .util import find_nearest, kwargupdate
from .specdata import SpectralData

# Omiting the warning for bad polynomial fitting
# import warnings
# warnings.simplefilter('ignore', np.RankWarning)


def loadspec(filename, unit='micron', r_error_col=None,
             masknull=True, label=None, **kwargs):
    r"""
    Load a spectrum file in units of *unit*. Returns a new Spectrum object.

    Parameters
    ----------
    filename: str
        Path for the spectrum file

    unit: str
        The wavelength unit. Possible values are: 'micron', 'angstron',
        'nanometer'.

    r_error_col: None or integer (optional)
        The column for the errors in the reflectance.

    masknull: boolean (optional)
        If True removes points where the wavelength is zero.

    label: None or string (optional)
        The spectrum label. If None it will take the file basename.

    **kwargs: Other arguments. See numpy.loadtxt

    Returns
    -------
    spec: Spectrum object

    """
    # setting default values for loading the spectrum file
    default_kwargs = {'unpack': True}
    for k in default_kwargs:
        if k not in kwargs:
            kwargs[k] = default_kwargs[k]
    # Loading the file using numpy.loadtxt
    spec_data = np.loadtxt(filename, **kwargs)
    # masking zero values in the wavelength array
    if masknull:
        mask = np.argwhere(spec_data[0] == 0)
        spec_data = np.delete(spec_data, mask, axis=1)
    # inserting it in as a Spectrum object
    if r_error_col is not None:
        r_error_col = spec_data[r_error_col]
    # setting the label
    if label is None:
        label = basename(splitext(filename)[0])
    # creating Spectrum object with given data
    spec = Spectrum(spec_data[0], spec_data[1],
                    r_unc=r_error_col, unit=unit,
                    label=label)
    return spec


def stack_spec(tup):
    r"""
    Stack Spectra arrays.

    Parameters
    ----------
    tup : array of Spectrum objects
          A sequence of Spectrum that must be the same shape along all but
          the second axis. 1-D arrays can be any length.

    Returns
    -------
    stacked : Spectrum
        The Spectrum array formed by stacking the given spectra, sorted by the
        wavelength.

    """
    wave = np.hstack([i.w for i in tup])
    ref = np.hstack([i.r for i in tup])
    r_unc_aux = [type(i.r_unc) for i in tup]
    if type(None) not in r_unc_aux:
        r_unc = np.hstack([i.r_unc for i in tup])
    else:
        r_unc = None
    stacked = Spectrum(wave, ref, r_unc=r_unc, unit=tup[0].unit,
                       label='_'.join(t.label for t in tup))
    stacked = stacked.sort(order='w')
    return stacked


@dataclass(repr=False)
class Spectrum(SpectralData):
    r"""Create a spectrum object.

        A spectrum array is a subclass of a record or structured numpy array,
        with one axis representing the wavelength vector (w) and other the
        reflectance (r). The optional

    Attributes
    ----------
        w: numpy array
            array corresponding to the wavelength vector

        r: numpy array
            array corresponding to the relative reflectance of
            the asteroid

        unit: str
            The wavelength units. Default is 'microns'.

        r_unc: numpy array (optional)
            array corresponding to the relative reflectance
            uncertainty of the asteroid

        path: None or str (optional)
            The path of the spectrum file

        label: None or str
            The spectrum label

        res: float
            The spectrum resolution (number of points).

    Methods
    -------
    trim
    fit
    autofit
    estimate_rms
    clean_spec
    mad
    rebin
    normalize
    mask_region
    save
    plot

    """

    w: np.ndarray
    r: np.ndarray
    r_unc: np.ndarray = field(default=None)
    label: str = field(default='asteroid')
    unit: str = field(default='micron')

    def __post_init__(self):
        r"""Inicialize class."""
        self.w = np.array(self.w, ndmin=1)
        self.r = np.array(self.r, ndmin=1)
        assert self.w.size == self.r.size
        if self.r_unc is None:
            dtype = ('w', 'r')
        else:
            self.r_unc = np.array(self.r_unc, ndmin=1)
            assert self.r_unc.size == self.r.size
            dtype = ('w', 'r', 'r_unc')
        SpectralData.__init__(self, self.w, dtype, self.label, unit=self.unit)

    def angstrom2micron(self):
        r"""Convert wavenlength axis from angstrom to micron."""
        self.w = self.w / 10000.0
        self.unit = 'micron'
        return self

    def micron2angstrom(self):
        r"""Convert wavenlength axis from micron to angstrom."""
        self.w = self.w * 10000.
        self.unit = 'angstrom'
        return self

    def fit(self, order=4, ftype='spline'):
        r"""
        Fit the spectrum using a polynomial or a smooth spline.

        Parameters
        ----------
        order: int
            Order of the fitting.

        ftype: str
            Type of algorithm to use for the fitting.
            Options are: 'spline' or 'polynomial'.

        Returns
        -------
        fspec: Spectrum Object
            The fitted spectrum array

        fcoefs: array-like
            the fitting coefficients

        """
        # Performing the fit
        if ftype == 'spline':
            fcoefs = UnivariateSpline(self.w, self.r, k=order)
            fspec_y = fcoefs(self.w)
        elif ftype == 'polynomial':
            fcoefs = np.polyfit(self.w, self.r, order)
            fspec_y = np.polyval(fcoefs, self.w)
        y_err = np.abs(fspec_y - self.r)
        # building new array
        fspec = self.__class__(w=self.w, r=fspec_y, r_unc=y_err, unit=self.unit,
                               label=self.label + '_fit')
        return fspec, fcoefs

    def autofit(self, degree_min=1, degree_max=12):
        r"""
        Find the best order of the polynomial fitting of the spectrum.

        Parameters
        ----------
        degree_min: int
            The minimal order for a fit

        degree_max: int
            The maximal order for a fit

        Returns
        -------
        fspec: Spectrum Object
            The fitted spectrum array

        fcoefs: array-like
            the fitting coefficients
        """
        # degrees which we will test the fitting
        degrees = np.arange(degree_min, degree_max+1)
        # calculation cross-validation score and error for each degree
        cross, _ = np.array([self._autofit_aux(deg) for deg in degrees]).T
        # getting minimun cross validation score
        aux = np.argmin(cross)
        bestfit = degrees[aux]
        # fitting spectrum
        fspec, fcoefs = self.fit(order=bestfit, ftype='polynomial')
        return fspec, fcoefs

    def _autofit_aux(self, degree):
        r"""Auxiliary funtion for the autofit method."""
        # creating polynomial and reshaping array for model validation
        polynom = PolynomialFeatures(degree=degree, include_bias=False)
        x_train = self.w.reshape((-1, 1))
        x_train_trans = polynom.fit_transform(x_train)
        # Create the linear regression model and train
        model = LinearRegression()
        model.fit(x_train_trans, self.r)
        # Calculate the cross validation score
        cross_valid = cross_val_score(model, x_train_trans, self.r,
                                      scoring='neg_mean_squared_error', cv=10)
        # Training predictions and error
        y_train_predictions = model.predict(x_train_trans)
        training_error = mean_squared_error(self.r, y_train_predictions)
        return -np.mean(cross_valid), training_error

    def estimate_rms(self, ftype='auto', order=5):
        r"""
        Estimate the signal-to-noise ratio in a spectrum.

        Parameters
        ----------
        ftype: str
            Type of fitting for the snr estimation. Options are:
            'spline', 'polynomial', 'auto'. Default is 'auto'.

        order: int
            Order of the adjust. Ignored if ftype is 'auto'.

        Returns
        -------
        rms: float
            The estimated snr value

        """
        # fitting the spectrum
        if ftype == 'auto':
            fspec, _ = self.autofit()
        else:
            fspec, _ = self.fit(order, ftype)
        # Estimating the SNR
        std_arr = np.abs(self.r - fspec.r)
        rms = np.std(std_arr)
        return rms

    def clean_spec(self, method='sigmaclip', sigma=3, fit='auto'):
        r"""
        Remove outliers from the spectrum.

        Parameters
        ----------
        method: str
            Method for detecting outliers. Currently only 'sigmaclip' available
            Default is 'sigmaclip'.

        sigma: int
            Remove points higher than sigma.

        fit: 'auto' or integer
            The order of the polynomial fit. If auto it will try to find
            automaticaly. Default is 'auto'.

        Returns
        -------
        spec: Spectrum
            Spectrum object with outliers removed.

        """
        if fit == 'auto':
            fspec, _ = self.autofit()
        else:
            fspec, _ = self.fit(order=fit, ftype='polynomial')
        if method == 'sigmaclip':
            cspec = np.divide(self.r, fspec.r)
            cspec_index = [self._sigma_clip(val, sigma=sigma, cspec=cspec)
                           for val in cspec]
        aux = self[cspec_index]

        spec = self.__class__(aux.w, aux.r, r_unc=aux.r_unc, unit=self.unit,
                              label=self.label + '_cleaned')
        return spec

    def _sigma_clip(self, val, sigma, cspec):
        r"""Auxiliary method to perform sigma-clipping on array elements."""
        if (np.median(cspec) - self._mad(cspec) * sigma < val) and \
           (val < np.median(cspec) + self._mad(cspec)*sigma):
            return True
        return False

    def mad(self, axis=None):
        r"""
        Calculate the median absolute deviation.

        Parameters
        ----------
        axis: str
            'wave', 'ref' or None.
            It will return the mad in the defined axis.
            If None, than returns the mad in both axis

        Returns
        -------
        The median absolute deviation

        """
        if axis is not None:
            return self._mad(axis)
        return self._mad('w'), self._mad('r')

    @staticmethod
    def _mad(arr):
        r"""Auxiliary function for calculating the MAD."""
        return np.median(np.abs(arr - np.median(arr)))

    def rebin(self, binsize=11, method='median', std=True,
              rem_trend=False):
        r"""
        Rebin the spectrum.

        Parameters
        ----------
        binsize: int
            The number of points in the bin.

        method: str
            The method for the rebinning.
            Options are:'mean and 'median'.

        std: boolean
            If True, also returns the deviation.
            In the case of the median, returns the
            MAD (median absolute deviation).

        rem_trend: boolean

        Returns
        -------
        The rebined spectrum

        """
        spec_size = len(self.w)
        y_stop = (spec_size//binsize)*binsize
        wave_arr = self.w[:y_stop]
        ref_arr = self.r[:y_stop]
        if method == 'median':
            func = np.median
            std_func = np.std  # ->>>>>>change later
        if method == 'mean':
            func = np.mean
            std_func = np.std
        wave_reb = func(wave_arr.reshape(spec_size // binsize, binsize),
                        axis=-1).T
        ref_reb = func(ref_arr.reshape(spec_size // binsize, binsize),
                       axis=-1).T
        if rem_trend:
            fspec = self.autofit()[0]
            ref_arr = np.divide(ref_arr, fspec.r[:y_stop])
        if std:
            std = std_func(ref_arr.reshape(spec_size // binsize, binsize),
                           axis=-1).T
            std = np.array(std)
        else:
            std = None
        return self.__class__(w=wave_reb, r=ref_reb, r_unc=std, unit=self.unit,
                              label=self.label + '_binned')

    def ref_from_wavelength(self, w, interpolate=True):
        r"""
        Get the spectrum reflectance in a particular wavelength.

        Parameters
        ----------
        w: float
            Wavelength value
            If interpolate=False, The code will search the closest
            value.

        interpolate: boolean (optional)
            If interpolate=False, The code will search the closest
            value. If True it will interpolate the value of w.

        Returns
        -------
        The reflectance value

        """
        if not interpolate:
            aux = find_nearest(self.w, w)[0]
            ref = self.r[aux]
        else:
            ref = np.interp(w, self.w, self.r)
        return ref

    def normalize(self, wnorm=0.55, window=None, interpolate=True):
        r"""
        Normalize the spectrum in a particular wavelength.

        Parameters
        ----------
        wnorm: float
            Wavelength value to normalize the spectrum.
            If interpolate=False, The code will search the closest
            value.

        window: None or float (optional)
            The wavelenght window size for normalizing.
            If None it will normalize in the wnorm point only.

        interpolate: boolean (optional)
            If interpolate=False, The code will search the closest
            value. If True it will interpolate the value of wnorm.

        Returns
        -------
        The normalized Spectrum

        """
        if window is None:
            norm_factor = self.ref_from_wavelength(wnorm,
                                                   interpolate=interpolate)
        else:
            aux = np.argwhere((self.w > wnorm-window) &
                              (self.w < wnorm+window))
            norm_factor = np.mean(self.r[aux])
        self.r = self.r / norm_factor
        if self.r_unc is not None:
            self.r_unc = self.r_unc / norm_factor
        return self

    def mask_region(self, region=[(1.3, 1.45), (1.8, 1.95)]):
        r"""
        Exclude a region of the spectrum.

        Parameters
        ----------
        w_min: float
            Wavelength lower limit of the masked region

        w_max: float
            Wavelength upper limit of the masked region

        Returns
        -------
        masked_spec: Spectrum
            The Spectrum array without the masked region

        """
        if isinstance(region[0], (float, int)):
            masked_spec = self.mask_region_aux(self, wmin=region[0],
                                               wmax=region[1])
        else:
            masked_spec = self
            for rr in region:
                masked_spec = self.mask_region_aux(masked_spec,
                                                   wmin=rr[0],
                                                   wmax=rr[1])
        return masked_spec

    def mask_region_aux(self, spec, wmin, wmax):
        r"""
        Exclude a region of the spectrum.

        Parameters
        ----------
        w_min: float
            Wavelength lower limit of the masked region

        w_max: float
            Wavelength upper limit of the masked region

        Returns
        -------
        The Spectrum array without the masked region

        """
        aux = np.argwhere((spec.w > wmin) & (spec.w < wmax))
        mask = np.ones(len(spec.w), dtype=bool)
        mask[aux] = 0
        w = spec.w[mask]
        r = spec.r[mask]
        r_unc = spec.r_unc
        if r_unc is not None:
            r_unc = spec.r_unc[mask]
        return self.__class__(w=w, r=r, r_unc=r_unc, unit=spec.unit,
                              label=spec.label)


    def plot(self, fax=None, show=False, savefig=None,
             axistitles=True, speckwargs=None, legendkwargs=None):
        r"""
        Quick plot of the Spectrum.

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

        axistitles: boolean
            If True will label the axis. Default is True.

        speckwargs: dict
            Arguments for matplotlib plot function.
            default values: {'c':'0.9', 'lw':'1'}.

        legendkwargs: dict
            Arguments for matplotlib legend function.
            default values: {'loc':'best'}.

        Returns
        -------
        the matplotlib.axes of the figure

        """
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = plt.gca()
        # setting default values for image plot with matplotlib
        specsty_defaults = {'c': '0.1', 'lw': 1}
        legendsty_defaults = {'loc': 'best'}
        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        # Ploting the spec
        fax.errorbar(self.w, self.r, yerr=self.r_unc, **speckwargs)
        # Checking if desired to plot the axis labels
        if axistitles:
            if self.unit == 'micron':
                unit_label = '$\mu$m'
            fax.set_xlabel('Wavelength (%s)' % unit_label)
            fax.set_ylabel('Reflectance')
        # plot legend?
        if 'label' in speckwargs:
            fax.legend(**legendkwargs)
        # check if save the image
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            matplotlib.use('TkAgg')
        # show in the matplotlib window?
        if show:
            plt.show()
        # return fax
