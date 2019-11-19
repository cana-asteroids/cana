r"""Tools for estimating parameters uncertainty."""

import numpy as np
import sys
from .. import Spectrum


class SpecError(object):
    r"""
    Monte-Carlo method for estimating the Spectrum parameters errors.

    All resampling methods assume a gaussian error for the reflectances.

    Methods
    -------
    rms: Estimates the rms of the spectrum and resample the reflectances
        considering the rms as the errors in the reflectances values

    removal: Resample the spectrum by randomly removing a percentage of
        the points.

    rebin: Rebin the Spectrum and takes the bin deviation as the error
        for resampling the reflectance values.

    """

    def __init__(self, n=100, method='rms', param=None):
        r"""
        Initialize the error model.

        Parameters
        ----------
        n: integer
            Number of iterations for Monte-Carlo model.
            Default 1000.

        method: 'rms', 'removal', 'rebin'
            The resampling method.

        param: float or integer
            The error methodoly parameter if needed. If errormethod='rms',
            then no value is necessary. If it is set for removal,
            the percentage of points to remove. For rebin, the param represents
             the binsize.

        """
        self.n = n
        self.method_name = method
        self.param = param
        if self.method_name == 'rms':
            self.resample = rms_resampling
            self.estimate = rms_error_estimate
        if self.method_name == 'rebin':
            self.resample = rebin_resampling
            self.estimate = rebin_error_estimate


def rebin_resampling(spec, binsize=5):
    r"""
    Resample the spectra considering the error from the rebinned spectrum.

    Parameters
    ----------
    spec: spectrum
        The spectrum object

    bisize: int
        The size of the bin

    Returns
    -------
    The resampled binned Spectrum

    """
    rspec = spec.rebin(binsize=binsize, std=True, rem_trend=True)
    spec_resamp = ref_resample(rspec.r, rspec.r_unc)
    return Spectrum(rspec.w, spec_resamp, r_unc=rspec.r_unc, unit=rspec.unit)


def rebin_error_estimate(spec, func, n, param=5, **kwargs):
    r"""
    Estimate error by a monte-carlo with spectra resampling from the rebining.

    Parameters
    ----------
    spec: Spectrum
        The Spectrum object.

    func: method
        The method for estimating a parameter from the spectrum.

    n: integer
        Number of iterations for Monte-Carlo model.
        Default 1000.

    param: integer
        Binsize for the rebining

    """
    rspec = spec.rebin(binsize=param, std=True, rem_trend=True)
    bin_err = rspec.r_unc
    out = []
    out_append = out.append
    i = 0
    while i < n:
        sp = spec_resample(rspec, bin_err)
        vl = func(sp, **kwargs)
        out_append(vl)
        i += 1
    out = np.array(out)
    return out.mean(axis=0), out.std(axis=0)


def rms_resampling(spec):
    r"""
    Resample the spectra considering the reflectance rms.

    Parameters
    ----------
    spec: Spectrum
        The Spectrum object

    Returns
    -------
    The resampled Spectrum

    """
    rms = spec.estimate_rms()
    rspec = spec_resample(spec, rms)
    return rspec


def rms_error_estimate(spec, func, n, param=None, **kwargs):
    r"""
    Estimate error by a monte-carlo with spectra resampling from the rms.

    Parameters
    ----------
    spec: Spectrum
        The Spectrum object.

    func: method
        The method for estimating a parameter from the spectrum.

    n: integer
        Number of iterations for Monte-Carlo model.
        Default 1000.

    param: None
        Not necessary

    """
    rms = spec.estimate_rms()
    out = []
    out_append = out.append
    i = 0
    while i < n:
        sp = spec_resample(spec, rms)
        vl = func(sp, **kwargs)
        out_append(vl)
        i += 1
    out = np.array(out)
    return out.mean(axis=0), out.std(axis=0)


def spec_resample(spec, rms):
    r"""Auxiliary method for resampling the spectrum."""
    spec_resamp = ref_resample(spec.r, rms)
    return Spectrum(spec.w, spec_resamp, r_unc=spec.r_unc, unit=spec.unit)


@np.vectorize
def ref_resample(point, rms):
    r"""Auxiliary method for resampling the reflectance."""
    aux = np.random.normal(point, rms, 1)
    return aux[0]
