r'''
'''
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import pandas as pd
from .photoconfig import CONFIG
from .spec import Spectrum

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))


class PhotometryBase(object):
    r'''
    Basic information of photometric systems

    Parameters
    ----------
    system: str
        The Photometric system identifier: 'sdss', 'jplus', 'jwst.

    Attributes
    ----------
    name: str
        The Photometric system identifier

    system: pandas.DataFrame
        A DataFrame with the photometric filters info

    center: pandas.DataFrame
        The central wavelength for each filter

    ids: list
        A list with the filters identifiers

    ids_err: list
        A list for the filters errors

    '''
    def __init__(self, system):
        self.name = system
        self.system = self._build_info()
        self.center = self.system['center']
        self.fwhm = self.system['fwhm']
        self.ids = self.system.index.tolist()
        self.ids_err = [iid+'_err' for iid in self.ids]

    def _build_info(self):
        r'''
        '''
        con = CONFIG[self.name]
        info = pd.DataFrame.from_dict(con, orient='columns')
        info = info.sort_values(by='center')
        return info

    def __repr__(self):
        return self.system.T.__repr__()

    def __getitem__(self, item):
        if item in self.system.columns:
            return self.system[item]
        elif item in self.system.index:
            return self.system.loc[item]
        else:
            raise ValueError('could not find %s' % (item))

    def plot_system(self):
        r'''
        '''


class Photometry(PhotometryBase):
    r'''
    Representation of a Photometric System

    Parameters
    ----------
    system: str
        The Photometric system identifier: 'sdss', 'jplus', 'jwst.

    tcurve: boolean
        True if there is

    Attributes
    ----------
    transmission_curves: numpy.recarray
        A numpy array containing the transmission curves for each
        filter.

    name: str
        The Photometric system identifier

    system: pandas.DataFrame
        A DataFrame with the photometric filters info

    center: pandas.DataFrame
        The central wavelength for each filter

    ids: list
        A list with the filters identifiers

    ids_err: list
        A list for the filters errors
    '''
    def __init__(self, system='sdss', tcurve=True):
        self.name = system.lower()
        super(Photometry, self).__init__(system=system)
        if (system == 'dawn') or (system == 'standard'):
            tcurve = False
        self.convol_type = tcurve
        if self.convol_type:
            self.transmission_curves = self._load()


    def _load(self):
        fsystems_dir = PWD+'/datasets/data/photometry/'
        if self.name == 'sdss':
            photofile = 'sdss/sdss.tab'
            ids = ['w', 'u', 'g', 'r', 'i', 'z']
        if self.name == 'ecas':
            pass
        if self.name == 'osiris':
            pass
        if self.name == 'jplus':
            photofile = 'jplus/jplus.tab'
            ids = ['w', 'F378', 'F861', 'F430', 'r', 'F410', 'F395', 'i', 'g',
                   'z', 'F660', 'F348', 'F515']
        if self.name == 'dawm': #-> need implementation
            pass
        if self.name == 'jwst':
            photofile = 'jwst/jwst.tab'
            ids = ['w', 'F070W', 'F090W', 'F115W', 'F140M', 'F150W2',
                   'F150W', 'F162M', 'F164N', 'F182M', 'F187N', 'F200W',
                   'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2',
                   'F323N', 'F335M', 'F356W', 'F360M', 'F405N', 'F410M',
                   'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
        if self.name == '2mass':
            photofile = '2mass/2mass.tab'
            ids = ['w', 'j', 'h', 'k']
        if self.name == 'spitzer':
            photofile = 'spitzer/spitzer.tab'
            ids = ['w', 'ch1', 'ch2', 'ch3', 'ch4']
        dtype = [(ii, np.float64) for ii in ids]
        transmission_curves = np.loadtxt(fsystems_dir+photofile,
                                         dtype=dtype)
        return transmission_curves

    @staticmethod
    def _fields_view(arr, fields):
        dtype2 = np.dtype({name: arr.dtype.fields[name] for name in fields})
        return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

    def bandpass_array(self, fid):
        r"""Transimission curve of a photomotric filter.

        Parameters
        ----------
        fid: str
            The filter identifier

        Returns
        -------
        numpy.ndarray
            The filter transimission curve in the relevant wavelengths
        """
        aux = np.argwhere(self.transmission_curves[fid] > 0)
        bpass = self._fields_view(self.transmission_curves, ['w', fid])[aux]
        return bpass

    def convol(self, spec, filters=None, label=None, interp=False):
        r'''Apply a convolution of a spectrum with the photometric system

        Parameters
        ----------
        spec: Spectrum
            A Spectrum object

        filters: None, str or list
            The filters to apply the convolution. None will try the convolution
            in all filters

        label: string
            The label for the output

        Returns
        -------
        PhotoDataframe with colvoluted spectrophotometry
        '''
        ## -> check here the bandpass
        # checking which bands to convolve
        if (filters is not None) and (not isinstance(filters, list)):
            cfilters = filters.split()
        else:
            cfilters = self.ids
        # making the collumns for the errors
        if self.convol_type:
            convoluted = [self._transmission_convol_aux(spec, fil, interp) for fil in cfilters]
        else:
            convoluted = [self._box_convol_aux(spec, fil, interp) for fil in cfilters]
        if label == None:
            label = spec.label
        # making output
        photodata = PhotoDataframe([convoluted], system=self.name, columns=cfilters,
                                   index=[label], autofill=False)
        return photodata

    def _transmission_convol_aux(self, spec, fband, interp):
        r'''
        '''
        if not interp:
            # interpolating the the transimission curve
            ref_f = np.interp(spec.w, self.transmission_curves['w'],
                              self.transmission_curves[fband])
            # performing the convolution
            cval = np.trapz(ref_f*spec.r*spec.w, spec.w) /     \
                   np.trapz(ref_f*spec.w, spec.w)
        else:
            wave = np.linspace(self.center[fband] - self.fwhm[fband]/2.,
                               self.center[fband] + self.fwhm[fband]/2.,
                               interp)
            ref = np.interp(wave, spec.w, spec.r)
            ref_f = np.interp(wave, self.transmission_curves['w'],
                              self.transmission_curves[fband])
            # performing the convolution
            cval = np.trapz(ref_f*ref*wave, wave) /     \
                            np.trapz(ref_f*wave, wave)
        return cval

    def _box_convol_aux(self, spec, fband, interp):
        r'''
        '''
        if not interp:
            interp = 200 #-> work on how to pass this param
        wave = np.linspace(self.center[fband] - self.fwhm[fband]/2.,
                           self.center[fband] + self.fwhm[fband]/2.,
                           interp)
        ref = np.interp(wave, spec.w, spec.r)
        ref_f = np.ones(interp)
        # performing the convolution
        cval = np.trapz(ref_f*ref*wave, wave) /     \
               np.trapz(ref_f*wave, wave)
        return cval

    def _check_bandpass(self):
        r'''
        Checks if bandpass is covered though the spectrum
        when performing a convolution.
        '''
        return True


    def plot_curves(self, fax=None, show=True, legend=True,
                    linestyle=None, colormap='viridis'):
        r'''
        View the transmission curves

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

        linestyle: dict
            Arguments for matplotlib plot function, to be applied in the
            spectrum plot.
        '''
        # creates the defautl values for plotting
        if linestyle is None:
            linestyle = {'lw':'1'}
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure()
            fax = fig.gca()
        # making colors for plot
        cmap = get_cmap(colormap)
        if self.convol_type:
            color_aux = np.linspace(0, 1, len(self.transmission_curves.dtype.names[1:]))
            for nid, nam in enumerate(self.transmission_curves.dtype.names[1:]):
                linestyle['c'] = cmap(color_aux[nid])
                fax.plot(self.transmission_curves['w'], self.transmission_curves[nam],
                         label=nam, **linestyle)
        else:
            pass #-> make gaussian based on available info
        if legend:
            fax.legend(prop={'size': 8})
        if show:
            plt.show()


class PhotoDataframe(pd.DataFrame, PhotometryBase):
    r'''
    A subclass of a pandas.DataFrame to handle spectrophotometric tables

    Parameters
    ----------
    data: list or pandas.DataFrame

    system: None or string
        The photometric system beeing used. Options are: sdss, jplus, jwst, ecas and osiris

    missing_columns: list

    autofill: boolean
        If True will fill collumns based on the photometric system

    Extended Attributes
    -------------------
    system: PhotometryBase
        A simplified dataframe with photometric system information

    Exetended Methods
    -----------------
    insert_errors
    plot
    flux2mag
    mag2flux
    color
    to_spec
    '''

    @property
    def _constructor(self):
        return PhotoDataframe._internal_ctor

    _metadata = ['system', 'missing_columns', 'autofill']

    @classmethod
    def _internal_ctor(cls, *args, **kwargs):
        kwargs['missing_columns'] = None
        kwargs['autofill'] = False
        return cls(*args, **kwargs)

    def __init__(self, data=None, system=None, missing_columns=None,
                 autofill=True, columns=None, index=None, copy=True):
        if isinstance(data, pd.core.internals.BlockManager):
            system = self._recognize_system(data.items)
        # getting photometric system information
        if isinstance(system, str):
            self.system = PhotometryBase(system)
        # building columns names
        if system is not None:
            if autofill:
                cols = self._columns_ids()
                # certifying data has proper dimensions
                if data is not None:
                    if not isinstance(data, pd.core.internals.BlockManager):
                        data = np.array(data, ndmin=2)
                        # checking missing columns
                        if columns is not None:
                            missing_columns = [cnam for cnam in cols if cnam not in columns]
                        if (data.shape[1] == len(cols)/2) and missing_columns is None:
                            missing_columns = 'err'
                        if missing_columns == 'err':
                            missing_columns = self.system.ids_err
                        if missing_columns is not None:
                            mid_aux = [cols.index(mid) for mid in missing_columns]
                            for mid in mid_aux:
                                data = np.insert(data, mid, 0, axis=1)
                    super(PhotoDataframe, self).__init__(data=data,
                                                     index=index,
                                                    columns=cols,
                                                    copy=copy)
                else:
                    super(PhotoDataframe, self).__init__(index=index,
                                                        columns=cols,
                                                        copy=copy)
            else:
                super(PhotoDataframe, self).__init__(data=data,
                                    index=index,
                                    columns=columns,
                                    copy=copy)
        else:
            super(PhotoDataframe, self).__init__(data=data,
                                    index=index,
                                    columns=columns,
                                    copy=copy)
    @staticmethod
    def _recognize_system(block):
        aux_ids =  [fid for fid in list(block) if '_err' not in fid]
        for k,v in CONFIG.iteritems():
            if sorted(aux_ids) == sorted(v['center'].keys()):
                return k

    def _columns_ids(self):
        cols = []
        for nfil, fil in enumerate(self.system.ids):
            cols.extend([fil, self.system.ids_err[nfil]])
        return cols


    def insert_error(self, err, filters=None):
        r''' Inserts error on filters errors collumns

        Parameters
        ----------
        err: float or list
            A error value or a list of errors respectively to the filters

        filters: None, string or list
            The filters identifiers. None will apply on all filters

        Returns
        -------
        PhotoDataframe
        '''
        if isinstance(filters, str):
            filters = filters.split()
        elif filters is None:
            filters = self.system.ids
        if isinstance(err, float):
            err = [err for i in filters]
        ids = [fil+'_err' for fil in filters]
        for i, ii in enumerate(ids):
            self[ii] = [err[i] for _ in xrange(self.shape[0])]
        return self


    def plot(self):
        pass

    def flux2mag(self):
        pass

    def mag2flux(self):
        pass

    def color(self):
        pass

    def tospec(self, idx):
        aux = [i for i in self.columns if '_err' not in i]
        aux_err = [i for i in self.columns if '_err' not in i]

        w = self.system.system.loc[aux]['center']
        r = self.loc[idx][aux]
        if aux_err:
            r_unc = self.loc[idx][aux_err]
        r_unc = None
        return Spectrum(w, r, r_unc, label=idx)
