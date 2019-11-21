r"""Tools to handle asteroid spectra."""

from os.path import splitext, basename
import numpy as np
from scipy.interpolate import UnivariateSpline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt
import matplotlib
from .util import find_nearest, kwargupdate
# Omiting the warning for bad polynomial fitting
# import warnings
# warnings.simplefilter('ignore', np.RankWarning)


def loadspec(filename, unit='micron', r_error_col=None,
             masknull=True, label=None, **kwargs):
    r"""
    Load a spectrum file and to an Spectrum object.

    Parameters
    ----------
    filename: str
        Path for the spectrum file

    unit: str
        The wavelength unit. Possible values are: 'micron', 'angstron',
        'nanometer'. Default is 'micron'.

    r_error_col: None or integer (optional)
        The collumn for the errors in the reflectance.
        Default is None

    masknull: boolean (optional)
        If True removes points where the wavelength is zero. Default
        is True

    label: None or string (optional)
        The spectrum label. If None it will take the file basename.
        Default is None.

    **kwargs: Other arguments for numpy.loadtxt

    Returns
    -------
    A Spectrum object

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
    spec = Spectrum(w=spec_data[0], r=spec_data[1],
                    r_unc=r_error_col, unit=unit,
                    path=filename, label=label)
    return spec


def stack_spec(tup):
    r"""
    Stack Spectra arrays.

    Parameters
    ----------
    tup : sequence of Spectrum
          The arrays must have the same shape along all but the second axis,
          except 1-D arrays which can be any length.

    Returns
    -------
    stacked : Spectrum
        The Spectrum array formed by stacking the given spectra, sorted by the
        wavelength.

    """
    wave = np.hstack((i.w for i in tup))
    ref = np.hstack((i.r for i in tup))
    spec = Spectrum(w=wave, r=ref, r_unc=None, unit=tup[0].unit,
                    path=None, label=tup[0].label)
    spec.sort(order='w')
    return spec


class Spectrum(np.recarray):
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

        r_unc: numpy array (optinal)
            array corresponding to the relative reflectance
            uncertainty of the asteroid

        path: None or str (optional)
            The path of the spectrum file

        label: None or str
            The spectrum label

    --Atributes
    ----------
        res: float
            The spectrum resolution (number of points).

    Methods
    -------
    trim
    fit
    autofit
    clean_spec
    estimate_rms
    rebin
    save
    plot

    """

    def __new__(cls, w, r, unit='micron', r_unc=None, path=None, label=None):
        r"""Generate Spectrum object."""
        # Asserting Spectrum is well structured
        assert len(w) == len(r)
        assert unit.lower() in ['micron', 'angstrom'], \
            'unit should be microns or angstroms'
        # Creating Spectrum
        dt_list = [('w', 'float'), ('r', 'float')]
        arr_list = [w, r]
        if (type(r_unc) is np.ndarray) and (r_unc.size == r.size):
            assert len(r) == len(r_unc)
            dt_list.append(('r_unc', 'float'))
            arr_list.append(r_unc)
        dt = np.dtype(dt_list)
        buffer = np.array(list(zip(*arr_list)), dtype=dt)
        obj = super(Spectrum, cls).__new__(cls, buffer.shape, dtype=dt)
        obj.w = buffer['w']
        obj.r = buffer['r']
        if r_unc is not None:
            obj.r_unc = buffer['r_unc']
        obj = obj.view(Spectrum)
        obj.unit = unit.lower()
        obj.res = len(obj.w)
        obj.path = path
        obj.label = label
        if r_unc is None:
            obj.r_unc = r_unc
        return obj

    def __array_finalize__(self, obj, context=None):
        if obj is None:
            return
        self.info = getattr(obj, 'info', None)

    def __array_prepare__(self, out_arr, context=None):
        return np.ndarray.__array_prepare__(self, out_arr, context)

    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def angstrom2micron(self):
        r"""Convert wavenlength axis from angstrom to micron."""
        self.w = self.w/10000.0
        self.unit = 'micron'
        return self

    def micron2angstrom(self):
        r"""Convert wavenlength axis from micron to angstrom."""
        self.w = self.w*10000.
        self.unit = 'angstrom'
        return self

    def trim(self, w_min=0.55, w_max=0.85):
        r"""
        Trim the spectrum in a defined region.

        Parameters
        ----------
        w_min: float
            Wavelength lower limit

        w_max: float
            Wavelength upper limit

        Returns
        -------
        self
            the trimmed spectrum.

        """
        aux = self[(self.w > w_min) & (self.w < w_max)]
        r_unc = self.r_unc
        if self.r_unc is not None:
            r_unc = aux['r_unc']
        return Spectrum(w=aux['w'], r=aux['r'], r_unc=r_unc, unit=self.unit,
                        path=self.path, label=self.label)

    def fit(self, order=4, ftype='spline'):
        r"""
        Fit the spectrum using a polynomial or a smooth spline.

        Parameters
        ----------
        order: int
            Order of the fitting.

        ftype: str
            Type of algorithm to use for the fitting. Options are: 'spline' or
            'polynomial'. Default is 'spline'.

        Returns
        -------
        The fitted spectrum array, the fitting coefcients

        """
        # Performing the fit
        if ftype == 'spline':
            fcoefs = UnivariateSpline(self.w, self.r, k=order)
            fspec_y = fcoefs(self.w)
        elif ftype == 'polynomial':
            fcoefs = np.polyfit(self.w, self.r, order)
            fspec_y = np.polyval(fcoefs, self.w)
        yerr = np.abs(fspec_y - self.r)
        # building new array
        fspec = Spectrum(w=self.w, r=fspec_y, r_unc=yerr, unit=self.unit,
                         path=self.path, label=self.label)
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
        The fitted spectrum array, the fitting coefcients

        """
        # degrees which we will test the fitting
        degrees = np.arange(degree_min, degree_max+1)
        # calculation cross-validation score and error for each degree
        cross, err = np.array([self._autofit_aux(deg) for deg in degrees]).T
        # getting minimun cross validation score
        aux = np.argmin(cross)
        bestfit = degrees[aux]
        # fitting spectrum
        fspec, fcoefs = self.fit(order=bestfit, ftype='polynomial')
        return fspec, fcoefs

    def _autofit_aux(self, degree):
        r"""Auxiliary funtion for the autofit method."""
        # creating polinomial and reshaping array for model valitation
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
        The Spectrum with removed "bad" points

        """
        if fit == 'auto':
            fspec, _ = self.autofit()
        else:
            fspec, _ = self.fit(order=fit, ftype='polynomial')
        cspec = np.divide(self.r, fspec.r)
        cspec_index = [self._sigma_clip(val, sigma=sigma, cspec=cspec)
                       for val in cspec]
        aux = self[cspec_index]
        return Spectrum(aux.w, aux.r, unit=self.unit,
                        path=self.path, label=self.label)

    def _sigma_clip(self, val, sigma, cspec):
        r"""Auxiliary method to perform sigma-clipping on array elements."""
        if (np.median(cspec) - self._mad(cspec)*sigma < val) and \
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
            The method for the rebinning. Options are:
            'mean and 'median'. Default is 'median'.

        std: boolean
            If True, also returns the deviation.
            In the case of the median, returns the
            MAD (median absolute deviation).

        rem_trend=False

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
        return Spectrum(w=wave_reb, r=ref_reb, r_unc=std, unit=self.unit,
                        path=self.path, label=self.label)

    def normalize(self, wnorm=0.55, window=None, interpolate=False):
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

        interpolate: boolean (optional) -> not implemented
            If interpolate=False, The code will search the closest
            value. If True it will interpolate the value of wnorm.

        Returns
        -------
        The normalized Spectrum

        """
        if window is None:
            aux = find_nearest(self.w, wnorm)[0]
            self.r = self.r / self.r[aux]
        else:
            aux = np.argwhere((self.w > wnorm-window) &
                              (self.w < wnorm+window))
            self.r = self.r / np.mean(self.r[aux])
        return self

    def mask_region(self, wmin, wmax):
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
        aux = np.argwhere((self.w > wmin) & (self.w < wmax))
        mask = np.ones(len(self.w), dtype=bool)
        mask[aux] = 0
        w = self.w[mask]
        r = self.r[mask]
        r_unc = self.r_unc
        if r_unc is not None:
            r_unc = self.r_unc[mask]
        return Spectrum(w=w, r=r, r_unc=r_unc, unit=self.unit,
                        path=self.path, label=self.label)

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
        np.savetxt(fname, self, fmt=fmt, delimiter=delimiter, header=header,
                   footer=footer, comments=comments, encoding=encoding)
        self.path = fname

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
            fax = fig.gca()
        # setting default values for image plot with matplotlib
        specsty_defaults = {'c': '0.1', 'lw': '1'}
        legendsty_defaults = {'loc': 'best'}
        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        # Ploting the spec
        fax.plot(self.w, self.r, **speckwargs)
        # Checking if desired to plot the axis labels
        if axistitles:
            fax.set_xlabel('Wavelength (%s)' % self.unit)
            fax.set_ylabel('Normalized Reflectance')
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
