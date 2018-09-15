
import pandas as pd
import numpy as np
from cana import PhotoError

def splitandcheck(ref, filters):
    r'''
    '''
    if isinstance(filters, basestring):
        filters = filters.split()
    # checking if selected filters are in the database
    aux = ref.columns.tolist()
    assert set(filters).issubset(aux), \
    """selected filters are not a subset of the dataframe. {0} not in {1}""".format(filters, aux)
    return filters


def photoslope(photodata, norm='g', filters='g r i z', error=PhotoError()):
    r'''
    '''
    out = pd.DataFrame(columns=['slope', 'slope_err'])
    filters = splitandcheck(photodata, filters)
    for index, ast in photodata.iterrows():
        print index, len(photodata)
        slope_arr = error.distribution(photodata.system['center'], ast, filters, _slope_aux, norm=norm)
        out.loc[index] = [slope_arr.mean(), slope_arr.std()]
    return out

def _slope_aux(ref, wave, filters, norm):
    r'''
    '''
    ref[filters] = ref[filters]/ref[norm]
    coefs = np.polyfit(wave[filters].values, ref[filters].values, 1)
    slp = (coefs[0]/ref[norm]) * 10
    return slp