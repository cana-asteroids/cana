r'''
'''

import numpy as np
from cana.spec import Spectrum


class SpecError(object):
    r'''
    Monte-Carlo method for estimating the Spectrum parameters errors.
    All resampling methods assume a gaussian error for the reflectances. 

    Methods:
    --------
    rms: Estimates the rms of the spectrum and resample the reflectances
        considering the rms as the errors in the reflectances values
    
    removal: Resample the spectrum by randomly removing a percentage of
        the points.
    
    rebin: Rebin the Spectrum and takes the bin deviation as the error
        for resampling the reflectance values. 


    Parameters
    ----------
    n: integer
        Number of iterations for Monte-Carlo model.
        Default 1000.
    
    method: 'rms', 'removal', 'rebin'
        The resampling method.
    
    param: float or integer
        The error methodoly parameter if needed. If errormethod='rms', 
        then no value is necessary. If it is set for removal, the percentage
        of points to remove. For rebin, the param represents the binsize.
    '''
    def __init__(self, n=1000, method='rms', param=5):
        r'''
        '''
        self.n = n
        self.method = method
        self.param = param
        # sppeding up some methods
        self.resample_vec = np.vectorize(self._ref_resample_aux)


    def binmethod(self, spec, binsize):
        rspec = spec.rebin(binsize=binsize, std=True, rem_trend=True)
        spec_resamp = self.resample_vec(rspec.r, rspec.r_unc)
        return Spectrum(rspec.w, spec_resamp, r_unc=rspec.r_unc, unit=rspec.unit)

    
    def removalmethod(self, spec, rem_per):
        r'''
        '''
        rem_size = int(spec.res * rem_per)
        rid = np.random.random_integers(0, spec.res, size=rem_size)
        wave_copy = np.delete(spec.w, rid)
        ref_copy = np.delete(spec.r, rid)
        return Spectrum(wave_copy, ref_copy, r_unc=spec.r_unc, unit=spec.unit)

    def rmsmethod(self, spec):
        r'''
        '''
        rms = spec.estimate_rms()
        resampled_spec = self._rms_resample(spec, rms)
        return resampled_spec

    def _rms_resample(self, spec, rms):
        r'''
        '''
        spec_resamp = self.resample_vec(spec.r, rms)
        return Spectrum(spec.w, spec_resamp,r_unc=spec.r_unc, unit=spec.unit)

    @staticmethod
    def _ref_resample_aux(point, rms):
        aux = np.random.normal(point, rms, 1)
        return aux[0]

    def resample(self, spec):
        r'''
        '''
        if self.method == 'rms':
            return self.rmsmethod(spec)
        if self.method == 'removal':
            return self.removalmethod(spec, self.param)
        if self.method == 'rebin':
            return self.binmethod(spec, self.param)


    def distribution(self, spec, func):
        out = []
        for _ in xrange(self.n):
            sp = self.resample(spec)
            out_aux = func(sp)
            out.append(out_aux)
        return np.array(out)
