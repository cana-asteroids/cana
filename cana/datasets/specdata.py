r"""Functions for data access."""

import os
from .. import loadspec

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))


def getspectrum(aid, ref='primass'):
    r"""
    Get a spectrum from the available datasets.

    Parameters
    ----------
    sid: integer
        The asteroid number identification

    ref: str
        The reference for the spectrum: 'primass'

    Returns
    -------
    Spectrum

    """
    ref = ref.lower()
    if not isinstance(aid, str):
        aid = str(aid)
    specpath = PWD + '/data/testdata/spectra/{0}/{1}.tab'.format(ref, aid)
    spec = loadspec(specpath, unit='micron')
    return spec


def listdata():
    r"""List available data for testing cana methods."""
    print('Not yet implemented')
