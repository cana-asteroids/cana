
import os
from astropy.io import ascii
import numpy as np
from cana.photo import PhotoDataframe

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))

def read_sdss(fname='https://sbnarchive.psi.edu/pds3/non_mission/EAR_A_I0035_5_SDSSTAX_V1_1/data/sdsstax_obs_table.tab'):
    r'''
    Reads SDSS-based taxonomic classificationfor all asteroids observed with SDSS.

    The reference is: 
        Hasselmann, P. H., Carvano, J. M., and Lazzaro, D., SDSS-based Asteroid Taxonomy V1.1.
        EAR-A-I0035-5-SDSSTAX-V1.1. NASA Planetary Data System, 2012.

    This function is oriented for the sdsstax_obs_table.tab

    Parameters
    ----------
    fname: string
        Path of the sdsstax_obs_table.tab file. 
        Default gets it directly from PDS Dust Archive
    
    Returns
    -------
    PhotoDataFrame
    '''
    sdss = ascii.read(fname, format='fixed_width_no_header',
                      names=('aid', 'ast_name', 'prov_desig',
                             'tax', 'score', 'moid', 'bad',
                             'u', 'u_err', 'g', 'g_err', 'r', 'r_err',
                             'i', 'i_err', 'z', 'z_err'),
                      col_starts=(0, 7, 24, 36, 39, 44, 52, 55, 61, 68, 74,
                                  81, 87, 94, 100, 107, 113))
    not_num = np.where(sdss['aid'] == 0)[0]
    sdss['aid'] = sdss['aid'].astype('str')
    sdss['aid'][not_num] = sdss['prov_desig'][not_num]
    sdss = sdss.to_pandas()
    psdss = PhotoDataframe(sdss, 'sdss', columns=sdss.columns, autofill=False)
    return psdss
