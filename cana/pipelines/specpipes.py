
import cana
import pandas as pd
import matplotlib.pyplot as plt

def primitive_visible(spec, tax=None, slope=None, band=None,
              isband=None, error=cana.SpecError(),
              speckwargs=None, outplot=None, verbose=False):
    r'''
    Pipeline method to run standard methods for the analysis of primitive asteroid
    visible spectrum.

    Parameters
    ----------


    '''
    speckwars_default = {'unit': 'microns'}
    speckwargs = cana.kwargupdate(speckwars_default, speckwargs)
    isband_default = {'min_depth': 1., 'theoric_min': 0.7, 'max_dist': 0.05,
                      'sigma':3}
    isband = cana.kwargupdate(isband_default, isband)
    # applying for Spectrum object
    if isinstance(spec, cana.Spectrum):
        return _primitive(spec, tax=tax, slope=slope, band=band,
              isband=isband, error=cana.SpecError(), outplot=outplot)
    # for single spectrum path
    elif isinstance(spec, basestring):
        spec = cana.loadspec(spec, **speckwargs)
        return _primitive(spec, tax=tax, slope=slope, band=band,
              isband=isband, error=cana.SpecError(), outplot=outplot)
    #for list of spectra path
    elif isinstance(spec, list):
        aux = []
        for fsp in spec:
            if verbose is not None:
                print(fsp)
            sp = cana.loadspec(fsp, **speckwargs)
            aux.append(_primitive(sp, tax=tax, slope=slope, band=band,
                       isband=isband, error=cana.SpecError(), outplot=outplot))
        out = pd.concat(aux)
        return out

def _primitive(spec, tax=None, slope=None, band=None,
               isband=None, error=cana.SpecError(), outplot=None):
    r'''
    Auxiliary function to run primitive asteroid visibile analysis pipeline
    '''
    # calculating the slope
    if slope is not None:
        specslope = slope.measure(spec)
    # Applying taxonomical classification
    if tax is not None:
        spectax = tax.classify(spec, return_n=1)
        taxaux = pd.DataFrame(spectax.DataFrame.values,
                                columns=spectax.DataFrame.columns,
                                index=[spec.label])

    if band is not None:
        if tax is not None or spectax.is_primitive():
            specband1 = band.measure(spec)
            if specband1.is_band(**isband):
                specband = specband1.DataFrame.T
            else:
                specband = pd.DataFrame(['-','-','-','-'],
                                        index=['depth', 'depth_unc', 'center', 'center_unc'],
                                        columns=[spec.label]).T
        else:
            specband = pd.DataFrame(['-','-','-','-'],
                                    index=['depth', 'depth_unc', 'center', 'center_unc'],
                                    columns=[spec.label]).T
    if outplot is not None:
        fig, ax = plt.subplots(2,2, figsize=(10,10))
        spec.plot(fax=ax[0][0], show=False)

        specslope.plot(fax=ax[0][1], show=False)

        spectax.plot(fax=ax[1][0], show=False)

        if not specband1.is_band(**isband):
            specband1.plot(fax=ax[1][1], show=False, dotkwargs = {'c':'r', 's':70, 'zorder':2},
                           bandkwargs = {'c':'r', 'lw':'2', 'zorder':1})
        else:
            specband1.plot(fax=ax[1][1], show=False)

        plt.savefig('{0}/{1}.png'.format(outplot, spec.label))

    out = pd.concat([specslope.DataFrame.T,taxaux, specband],axis=1, join='inner')
    return out
