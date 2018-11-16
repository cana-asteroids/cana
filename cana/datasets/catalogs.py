
import os
from astropy.io import ascii
import numpy as np
from cana.photo import PhotoDataframe
import pandas as pd

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))

def read_sdss(fname='https://sbnarchive.psi.edu/pds3/non_mission/EAR_A_I0035_5_SDSSTAX_V1_1/data/sdsstax_obs_table.tab'):
    r'''
    Reads SDSS-based taxonomic classificationfor all asteroids observed with SDSS.

    The reference is: 
        Hasselmann, P. H., Carvano, J. M., and Lazzaro, D., SDSS-based Asteroid Taxonomy V1.1.
        EAR-A-I0035-5-SDSSTAX-V1.1. NASA Planetary Data System, 2012.

    This function is oriented for the sdsstax_obs_table.tab

    Note: We are removing the logarithm scale in the flux.

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
    g_min = 10**(sdss['g'] - sdss['g_err']) 
    g_max = 10**(sdss['g'] + sdss['g_err'])
    g = np.mean([g_min, g_max])
    
    for i in ['u', 'g', 'r', 'i', 'z']:
        aux_min = 10**(sdss[i]) /g
        aux_max = 10**(sdss[i]) /g
        sdss[i+'_err'] = 10**(sdss[i+'_err']) -1
    psdss = PhotoDataframe(sdss, 'sdss', columns=sdss.columns, autofill=False)
    return psdss


def read_astdys(fname='http://hamilton.dm.unipi.it/~astdys2/propsynth/all.syn'):
    r'''
    Reads AstDys Proper Elements Database
     
    References
    ----------

    ANALYTIC PROPER AND MEAN ELEMENTS
    A. Milani and Z. Knezevic: Asteroid mean elements: higher order and iterative theories,
                             Celestial Mechanics, 71, 55-78, 1999.

    PROPER ELEMENTS MAIN BELT
    Z. Knezevic and A. Milani: Synthetic proper elements for outer main belt asteroids,
                             CMDA, 78, 17-46, 2000.
    Z. Knezevic, A. Lemaitre and A. Milani: 2002. The determination of asteroid proper elements.
                             In: Asteroids III (W. Bottke, A. Cellino, P. Paolicchi and R.P. Binzel, Eds.),
                             Univ. Arizona Press and LPI, 603-612.
    Z. Knezevic and A. Milani: 2003. Proper element catalogs and asteroid families.
                             Astron. Astrophys. 403, 1165-1173.

    PROPER ELEMENTS TROJANS

    Milani, A. The Trojan asteroid belt: proper elements, stability, chaos and families,
                             CMDA 57, 59-94, 1993. abstract 
    Beauge, C. and Roig, F. A Semianalytical Model for the Motion of the Trojan Asteroids:
                             Proper Elements and Families. 
                             Icarus, Volume 153, Issue 2, pp. 391-415, 2001 (abstract) 
    
    Returns
    -------
    pandas.DataFrame
    '''
    milani = ascii.read(fname, format='fixed_width_no_header', fast_reader=True,
                        names=('aid', 'mag_prop', 'a_prop',
                               'e_prop', 'sini_prop', 'n', 'g', 's', 'lce',
                               'my'),
                        col_starts=(0, 10, 18, 29, 38, 50, 64, 76, 90, 97),
                        comment='%')
    milani = milani.to_pandas()
    return milani

def read_wise(fname='https://sbnarchive.psi.edu/pds3/non_mission/EAR_A_COMPIL_5_NEOWISEDIAM_V1_0/data/neowise_mainbelt.tab'):
    r'''
    '''
    wise = pd.read_csv(fname, names=('aid', 'prov_desig', 'mpc_packed', 'H', 'G',
                               'mean_jd', 'N_W1', 'N_W2', 'N_W3', 'N_W4', 'fitcode',
                               'diameter', 'diameter_err',  'v_albedo',  'v_albedo_err',
                               'ir_albedo', 'ir_albedo_err', 'beaming', 'beaming_err',
                               'stacked_flag', 'reference'),
)
    return wise
