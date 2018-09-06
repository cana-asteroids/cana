r'''
Tools to handle asteroid spectra
'''

from os.path import splitext, basename
import numpy as np
from scipy.interpolate import UnivariateSpline
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from cana import find_nearest, kwargupdate
# Omiting the warning for bad polynomial fitting
import warnings
warnings.simplefilter('ignore', np.RankWarning)


def loadspec(fname, unit='microns', r_error_col=None,
             masknull=True, label=None, **kwargs):
    r'''
    Loads a spectrum file and to an Spectrum object

    Parameters
    ----------
    fname: str
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
    A Spectrum 
    '''
    # setting default values for loading the spectrum file
    default_kwargs = {'unpack': True}
    for k in default_kwargs:
        if k not in kwargs:
            kwargs[k] = default_kwargs[k]
    # Loading the file using numpy.loadtxt
    spec_data = np.loadtxt(fname, **kwargs)
    # masking zero values in the wavelength array
    if masknull:
        mask = np.argwhere(spec_data[0] == 0)
        spec_data = np.delete(spec_data, mask, axis=1)
    # inserting it in as a Spectrum object
    if r_error_col is not None:
        r_error_col = spec_data[r_error_col]
    # setting the label 
    if label == None:
        label = basename(splitext(fname)[0])
    spec = Spectrum(w=spec_data[0], r=spec_data[1],
                    r_unc=r_error_col, unit=unit, 
                    path=fname, label=label)
    return spec

class Spectrum(np.recarray):
    r'''Create a spectrum object.

        A spectrum array is a subclass of a record or structured numpy array,
        with one axis representing the wavelength vector (w) and other the 
        reflectance (r). The optional        

    Parameters
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

    Atributes
    ----------
        res: float
            The spectrum resolution (number of points).

    Extended Methods
    ----------------
    trim
    fit
    autofit
    clean_spec
    estimate_rms
    rebin
    save
    plot
    '''
    def __new__(cls, w, r, unit='microns', r_unc=None, path=None, label=None):                  
        assert len(w) == len(r)
        assert unit.lower() in ['microns', 'angstrom'], 'unit should be microns or angstroms'
        dt_list = [('w', 'float'), ('r', 'float')]
        arr_list = [w, r] 
        if (type(r_unc) is np.ndarray) and (r_unc.size==r.size):
            assert len(r) == len(r_unc)
            dt_list.append(('r_unc', 'float'))
            arr_list.append(r_unc)
        dt = np.dtype(dt_list)
        buffer = np.array(zip(*arr_list),dtype=dt)
        obj = super(Spectrum, cls).__new__(cls, buffer.shape, dtype=dt)
        obj.w =  buffer['w']
        obj.r =  buffer['r']
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
   
    def __array_prepare__(self, out_arr, context=None):
        # print 'In __array_prepare__:'
        # print '   self is %s' % repr(self)
        # print '   arr is %s' % repr(out_arr)
        # then just call the parent
        return np.ndarray.__array_prepare__(self, out_arr, context)
    
    def angstron2micron(self):
        self.w = self.w/10000.
        self.unit = 'micron'
        return self


    def __array_wrap__(self, out_arr, context=None):
        # print 'In __array_wrap__:'
        # print '   self is %s' % repr(self)
        # print '   arr is %s' % repr(out_arr)
        # then just call the parent
        return np.ndarray.__array_wrap__(self, out_arr, context)
    
    def trim(self, w_min=0.55, w_max=0.85):
        r'''
        Trim the spectrum in a defined region

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
        '''
        aux = self[(self.w > w_min) & (self.w < w_max)]
        return Spectrum(w=aux['w'], r=aux['r'], r_unc=self.r_unc, unit=self.unit,
                        path=self.path, label=self.label)

    def fit(self, order=4, ftype='spline'):
        r'''
        Fits the spectrum using a polynomial or a smooth spline

        Parameters
        ----------
        order: int
            Order of the fitting.

        ftype: str
            Type of algorithm to use for the fitting. Options are: 'spline' or 'polynomial'.
            Default is 'spline'.

        Returns
        -------
        The fitted spectrum array, the fitting coefcients
        '''
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

    def autofit(self, degree_min=1, degree_max=12, test_set_fraction=0.40):
        r'''
        Method for finding the best order of the polynomial fitting of the
        spectrum.

        Parameters
        ----------
        degree_min: int
            The minimal order for a fit

        degree_max: int
            The maximal order for a fit

        test_set_fraction: float
            Percentage of points to randomly separate from the spectrum to set as
            the test data. The remaining is used as the training set.

        Returns
        -------
        The fitted spectrum array, the fitting coefcients

        '''
        # Dividing the spectrum between train and test data
        x_train, x_test, y_train, y_test = train_test_split(self.w, self.r,
                                                            test_size=test_set_fraction)
        # degrees which we will test the fitting
        degrees = np.arange(degree_min, degree_max+1)
        # calculating the rms
        rms = np.array([self._autofit_rms(deg, x_train, y_train, x_test, y_test)
                        for deg in degrees])
        # calculating the inflection point
        cont_coef = np.polyfit([degrees[0], degrees[-1]], [rms[0], rms[-1]], 1)
        cont_line = np.polyval(cont_coef, degrees)
        rms_aux = np.divide(rms, cont_line)
        rms_min = rms_aux.argmin()
        bestfit_deg = degrees[rms_min]
        fspec, fcoefs = self.fit(order=bestfit_deg, ftype='polynomial')
        return fspec, fcoefs

    @staticmethod
    def _autofit_rms(degree, x_train, y_train, x_test, y_test):
        r'''
        Auxiliary funtion for the autofit method.
        '''
        coef = np.polyfit(x_train, y_train, degree)
        # predicting model
        # y_predict_train = np.polyval(coef, x_train)
        y_predict_test = np.polyval(coef, x_test)
        # measuring the rms
        # rms_train = np.sqrt(np.sum(np.square(y_predict_train-y_train)))
        rms_test = np.sqrt(np.sum(np.square(y_predict_test-y_test)))
        # w_rms_train = rms_train/(len(x_train)-degree-1)
        w_rms_test = rms_test/(len(x_test)-degree-1)
        return w_rms_test #abs(w_rms_test -w_rms_train)


    def estimate_rms(self, ftype='auto', order=5):
        r'''
        Estimates the signal-to-noise ratio in a spectrum.

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
        '''
        # fitting the spectrum
        if ftype == 'auto':
            fspec, _ = self.autofit()
        else:
            fspec, _ = self.fit(order, ftype)
        # Estimating the SNR
        std_arr = np.abs(self.r - fspec.r)
        rms = np.std(std_arr)
        return rms

    def clean_spec(self, method='sigmaclip', sigma=3, fit=1):
        r'''
        Remove outliers from the spectrum

        Parameters
        ----------
        method: str
            Method for detecting outliers. Currently only 'sigmaclip' available
            Default is 'sigmaclip'

        sigma: int
            Remove points higher than sigma.
        
        fit: 'auto' or integer
            The order of the polynomial fit. If auto it will try to find automaticaly

        Returns
        -------
        The clean Spectrum
        '''
        if fit =='auto':
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
        r'''
        Auxiliary method to performs the sigma-clipping of array elements.
        '''
        if (np.median(cspec) - self._mad(cspec)*sigma < val
           ) and (val < np.median(cspec) + self._mad(cspec)*sigma):
            return True
        return False

    def mad(self, axis=None):
        r'''
        Calculates the median absolute deviation.

        Parameters
        ----------
        axis: str
            'wave', 'ref' or None.
            It will return the mad in the defined axis.
            If None, than returns the mad in both axis

        Returns
        -------
        The median absolute deviation
        '''
        if axis is not None:
            return self._mad(axis)
        return self._mad('w'), self._mad('r')

    @staticmethod
    def _mad(arr):
        r'''
        Auxiliary function for calculating the MAD
        '''
        return np.median(np.abs(arr- np.median(arr)))


    def rebin(self, binsize=11, method='median', std=True, 
              rem_trend=False):
        r'''
        Rebins the spectrum.

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
        '''
        spec_size = len(self.w)
        y_stop = (spec_size//binsize)*binsize
        wave_arr = self.w[:y_stop]
        ref_arr = self.r[:y_stop]
        if method == 'median':
            func = np.median
            std_func = np.std ###### ->>>>>>change later
        if method == 'mean':
            func = np.mean
            std_func = np.std
        wave_reb = func(wave_arr.reshape(spec_size// binsize, binsize),
                        axis=-1).T
        ref_reb = func(ref_arr.reshape(spec_size// binsize, binsize),
                       axis=-1).T
        if rem_trend:
            fspec = self.autofit()[0]
            ref_arr = np.divide(ref_arr, fspec.r[:y_stop])
        if std:
            std = std_func(ref_arr.reshape(spec_size// binsize, binsize),
                          axis=-1).T
        return Spectrum(wave_reb, ref_reb, r_unc=std, unit=self.unit,
                         path=self.path, label=self.label)

    def normalize(self, wnorm=0.55, window=None, interpolate=False):
        r''' Normalize the spectrum in a particular wavelength.

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
        '''
        if window == None:
            aux = find_nearest(self.w, wnorm)[0]
            self.r = self.r / self.r[aux]
        else:
            aux = np.argwhere((self.w>wnorm-window) & (self.w<wnorm+window))
            self.r = self.r / np.mean(self.r[aux])
        return self
    
    def save(self, fname, fmt='%.5f', delimiter=' ', header='', footer='', 
             comments='#', encoding=None):
        r'''
        
        fname : filename or file handle
            If the filename ends in .gz, the file is automatically saved in compressed gzip format.
            loadtxt understands gzipped files transparently.

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
            String that will be prepended to the header and footer strings, to mark them as comments.

        encoding : {None, str}, optional
            Encoding used to encode the outputfile. Does not apply to output streams. 

        '''
        np.savetxt(fname, self, fmt=fmt, delimiter=delimiter, header=header, footer=footer, 
             comments=comments, encoding=encoding)

    def plot(self, fax=None, show=True, savefig=None,
             axistitles=True, speckwargs=None, legendkwargs=None):
        r'''
        Method for the spectrum vizualization. This method is implemented for a 
        quick visualization of the spectrum. For more complexes viazualitions we
        we recommend 

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
        ------
        the matplotlib.axes of the figure

        '''
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = fig.gca()
        # setting default values for image plot with matplotlib
        specsty_defaults = {'c':'0.1', 'lw':'1'}
        legendsty_defaults = {'loc':'best'}
        # updating plot styles
        speckwargs = kwargupdate(specsty_defaults, speckwargs)
        legendkwargs = kwargupdate(legendsty_defaults, legendkwargs)
        # Ploting the spec
        fax.plot(self.w, self.r, **speckwargs)
        # Checking if desired to plot the axis labels
        if axistitles:
            plt.xlabel('Wavelength (%s)' % self.unit)
            plt.ylabel('Normalized Reflectance')
        # plot legend?
        if 'label' in speckwargs:
            fax.legend(**legendkwargs)
        # check if save the image
        if savefig != None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
        # show in the matplotlib window?
        if show:
            plt.show()
        return fax


