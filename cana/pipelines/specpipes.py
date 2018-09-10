
import cana
import pandas as pd


def primitive(spec, tax=None, slope=None, band1=None,
              isband=None, error=cana.SpecError(), 
              speckwargs=None):
    r'''
    '''
    speckwars_default = {'unit':'microns'}
    speckwargs = cana.kwargupdate(speckwars_default, speckwargs)
    isband_default = {'min_depth':1., 'theoric_min':0.7, 'max_dist':0.05,
                      'sigma':3}
    isband = cana.kwargupdate(isband_default, isband)
    # applying for Spectrum object
    if isinstance(spec, cana.Spectrum):
        return _primitive(spec, tax=tax, slope=slope, band1=band1,
              isband=isband, error=cana.SpecError())
    # for single spectrum path
    elif isinstance(spec, basestring):
        spec = cana.loadspec(spec, **speckwargs)
        return _primitive(spec, tax=tax, slope=slope, band1=band1,
              isband=isband, error=cana.SpecError())
    #for list of spectra path
    elif isinstance(spec, list):
        aux = []
        for fsp in spec:
            sp = cana.loadspec(fsp, **speckwargs)
            aux.append(_primitive(sp, tax=tax, slope=slope, band1=band1,
                       isband=isband, error=cana.SpecError()))
        out = pd.concat(aux)
        return out

def _primitive(spec, tax=None, slope=None, band1=None,
              isband=None, error=cana.SpecError()):
    r'''
    '''
    if slope is not None:
        specslope = slope.measure(spec)
    if tax is not None:
        spectax = tax.classify(spec, return_n=1)
        taxaux = pd.DataFrame(spectax.DataFrame.values,
                                columns=spectax.DataFrame.columns,
                                index=[spec.label])
    if band1 is not None:
        if tax is not None and spectax.is_primitive():
            specband1 = band1.measure(spec)
            if specband1.is_band(**isband):
                specband1 = specband1.DataFrame.T
            else:
                specband1 = pd.DataFrame(['-','-','-','-'],
                                        index=['depth', 'depth_unc', 'center', 'center_unc'],
                                        columns=[spec.label]).T
        else:
            specband1 = pd.DataFrame(['-','-','-','-'],
                                    index=['depth', 'depth_unc', 'center', 'center_unc'],
                                    columns=[spec.label]).T

    out = pd.concat([specslope.DataFrame.T,taxaux, specband1],axis=1, join='inner')
    return out