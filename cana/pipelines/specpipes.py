
import pandas as pd
import matplotlib.pyplot as plt

from .. import loadspec, Spectrum, Depth, Slope, Taxonomy
from ..spectools.uncertainties import SpecError
from ..util import kwargupdate



def primitive_visible(spec, params='default', error=SpecError(),
                      isband=None, speckwargs=None, outplot=None,
                      verbose=False):
    r"""
    Pipeline for the analysis of primitive asteroid visible spectrum.

    Parameters
    ----------
    spec: Spectrum, spectrum file, spectrum file list
        The input can be a Spectrum object, a spectrum file or
        a list of spectrum files

    params: list of methodologies classes
        A list containing all methodologies to be applyied.
        If is set as 'default', it will run the default pipeline
        for pimitive asteroids with characterization of slope,
        taxonomy and hydration band.

    error: SpecError
        The model to estimate the slope uncertainty.

    isband: None or dictionary (optional)
        The parameters to consider a real absorption band.
        The default values are{'min_depth': 1, 'theoric_min': 0.7,
                               'max_dist': 0.5, 'sigma': 3}

    speckwargs: dictionary (optional)
        Kwargs for cana.loadspec method.
        It is only used if spec is str or list of str.
        Default is {'unit': 'microns'}

    outplot: None, 'show' or path (str)
        None will not make plots. 'show' will call plt.show().
        Any other string will be considered a directory path for
        plt.savefig to save the plots. In this case filenames will
        be read from the spectra file
    """
    # Checking which parametrization we will apply
    if params == 'default':
        slope = Slope()
        depth = Depth()
        tax = Taxonomy()
        params = [slope, depth, tax]
    # Stablishing default values for the analysis
    speckwars_default = {'unit': 'micron'}
    speckwargs = kwargupdate(speckwars_default, speckwargs)
    isband_default = {'min_depth': 1., 'theoric_min': 0.7, 'max_dist': 0.05,
                      'sigma': 3}
    isband = kwargupdate(isband_default, isband)
    # applying for Spectrum object
    if isinstance(spec, Spectrum):
        return primitive_aux(spec, params, isband=isband, error=error,
                             outplot=outplot)
    # for single spectrum path
    elif isinstance(spec, str):
        spec = loadspec(spec, **speckwargs)
        return primitive_aux(spec, params, isband=isband, error=error,
                             outplot=outplot)
    # for list of spectra path
    elif isinstance(spec, list):
        aux = []
        for fsp in spec:
            if verbose is not None:
                print(fsp)
            sp = loadspec(fsp, **speckwargs)
            aux.append(primitive_aux(sp, params, isband=isband, error=error,
                                     outplot=outplot))
        out = pd.concat(aux)
        return out


def primitive_aux(spec, params, isband, error, outplot=None):
    r"""Auxiliary function to run primitive asteroid visibile pipeline."""
    for par in params:
        if isinstance(par, Slope):
            # calculating the slope
            specslope = par.measure(spec)
        if isinstance(par, Taxonomy):
            spectax = par.classify(spec, return_n=1)
            taxaux = pd.DataFrame(spectax.DataFrame.values,
                                  columns=spectax.DataFrame.columns,
                                  index=[spec.label])
        if isinstance(par, Depth):
            specband_ = par.measure(spec)
            if specband_.is_band(**isband):
                specband = specband_.DataFrame.T
            else:
                specband = pd.DataFrame(['-', '-', '-', '-'],
                                        index=['depth', 'depth_unc',
                                               'center', 'center_unc'],
                                        columns=[spec.label]).T
        else:
            specband = pd.DataFrame(['-', '-', '-', '-'],
                                    index=['depth', 'depth_unc',
                                           'center', 'center_unc'],
                                    columns=[spec.label]).T
    # Creating output image
    if outplot is not None:
        _, ax = plt.subplots(2, 2, figsize=(10, 10))
        spec.plot(fax=ax[0][0], show=False)
        if any(isinstance(par, Slope) for par in params):
            specslope.plot(fax=ax[0][1], show=False)
        if any(isinstance(par, Taxonomy) for par in params):
            spectax.plot(fax=ax[1][0], show=False)
        if any(isinstance(par, Depth) for par in params):
            if not specband_.is_band(**isband):
                specband_.plot(fax=ax[1][1], show=False,
                               dotkwargs={'c': 'r', 's': 70, 'zorder': 2},
                               bandkwargs={'c': 'r', 'lw': '2', 'zorder': 1})
            else:
                specband_.plot(fax=ax[1][1], show=False)
        plt.tight_layout()
        # Saving figure
        if outplot != 'show':
            plt.savefig('{0}/{1}.png'.format(outplot, spec.label))
        else:
            plt.show()
    # Returning the results
    out = pd.concat([specslope.DataFrame.T, taxaux, specband],
                    axis=1, join='inner')
    return out
