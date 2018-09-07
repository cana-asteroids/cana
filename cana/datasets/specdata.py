
import os
from glob import glob
from cana.spec import loadspec
import pandas as pd

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))

def getspectrum(aid, ref='smass'):
    r'''
    Gets a spectrum from the available datasets

    Parameters
    ----------
    sid: integer
        The asteroid number identification
    
    ref: str
        The reference for the spectrum: 'smass' or 's3os2'

    Returns
    --------
    Spectrum
    '''
    ref = ref.lower()
    if not isinstance(aid, basestring):
        sid = str(aid)
    specpath = PWD + '/data/test/spectra/{0}/{1}.tab'.format(ref, aid)
    spec = loadspec(specpath, unit='micron')
    return spec


def listdata():
    r'''
    '''