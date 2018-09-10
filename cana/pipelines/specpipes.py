
import cana
import pandas as pd

def primitive(spec, tax=None, slope=None, band1=None, band2=None,
              error=cana.SpecError()):
    r'''
    '''
    if not isinstance(spec, list):
        if isinstance(spec, basestring):
            spec = cana.loadspec(spec)
        if slope is not None:
            specslope = slope.measure(spec)
        if tax is not None:
            spectax = tax.classify(spec)
        if band1 is not None:
            specband1 = band1.measure(spec)
        if band2 is not None:
            specband2 = band2.measure(spec)
        
