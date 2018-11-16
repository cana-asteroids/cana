
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad
from cana import PhotoError, PhotoDataframe, Photometry, loadspec
from cana.photo import PhotometryBase

# Absolute directory path for this file
PWD = os.path.dirname(os.path.abspath(__file__))

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

def measureband(func, out, photodata, method='joint', filters='r i z', method_class=0., error=PhotoError(), label=None):
    # operating for each asteroid
    for index, ast in photodata.iterrows():
        print index, len(photodata)
        ## applying methodology
        if method == 'threshold':
            # Calculating the depth distribution
            depth_arr = error.distribution(photodata.system.center, ast, filters, func)
            # if threshold not parsed, calculating threshold
            if not isinstance(method_class, float):
                err = [ast[col] for col in ast.index if '_err' in col]
                # If no threshold and no class is beeing passed,
                # calculate the threshold using default values
                if method_class == 'auto':
                    method_class = BandThreshold(photodata.system.name)
                threshold = method_class.measure(err, func, filters)
            else: 
                threshold = method_class
            band_prob = len(np.argwhere(depth_arr>threshold))/float(len(depth_arr))
            depth = depth_arr.mean()
        # checking if method is joint 
        elif method == 'joint':
            joint = method_class
            # checking if joint class is beeing passed
            if not isinstance(joint, JointProbability):
                joint = JointProbability(photodata.system.name)
            depth = func(ast, photodata.system.center, filters)
            band_prob = joint.measure(depth, ast, func, filters )
        # returning probability
        if label is not None:
            index = ast[label]
        out.loc[index] = [band_prob*100, depth]
        out = out.round(3)

    return out

def threepointsband(photodata, method='joint', filters='r i z', method_class=0., error=PhotoError(), label=None):
    r'''
    '''
    prob = pd.DataFrame(columns=['probability', 'depth'])
    filters = splitandcheck(photodata, filters)
    prob = measureband(_threepointsband, prob, photodata, method, filters, method_class, error, label)
    return prob

def _threepointsband(ref, wave, filters):
    r'''
    '''
    # the filters for fitting the continuum are the first and last filters
    contfilters = [filters[0], filters[-1]]
    # print contfilters, wave[filters[1]],filters[1]
    contcoef = np.polyfit(wave[contfilters], ref[contfilters], 1)
    contref = np.polyval(contcoef, wave[filters[1]])
    # print ref[filters[1]],  contref
    depth = 1- ref[filters[1]] / contref
    return depth 


def photoarea(photodata, method='joint', filters='r F660 i F861', method_class=0., error=PhotoError()):
    prob = pd.DataFrame(columns=['probability', 'depth'])
    filters = splitandcheck(photodata, filters)
    prob = measureband(_area_aux, prob, photodata, method, filters, method_class, error)
    return prob


def _area_aux(ref, wave, filters):
    # measuring continumn
    contfilters = [filters[0], filters[-1]]
    contcoef = np.polyfit(wave[contfilters], ref[contfilters], 1)
    contref = np.polyval(contcoef, wave[filters])
    #subtracting continumn
    ref = np.divide(ref[filters], contref)
    above = np.argwhere(ref[filters[1:-1]] > 1)
    if len(above) == 0:
        area = PolyArea(wave[filters], ref[filters])
    elif len(above) == len(ref[filters[1:-1]]):
        area = - PolyArea(wave[filters], ref[filters])
    else:
        areaabv = area_above(ref, wave, filters, above[0][0])
        bellow = np.argwhere(ref[filters[1:-1]] < 1)
        areablw = area_above(ref, wave, filters, bellow[0][0])
        # area = PolyArea(wave[filters], ref[filters])

        area = areablw - areaabv
        pass
    return area

def area_above(ref, wave, filters, idabv):
    r'''
    '''
    left_filters = filters[idabv: idabv+2]
    laux = np.polynomial.Polynomial.fit(wave[left_filters], ref[left_filters], 1)
    left_corner = (laux -1).roots()[0]
    #right side
    right_filters = filters[idabv+1: idabv+3]
    raux = np.polynomial.Polynomial.fit(wave[right_filters], ref[right_filters], 1)
    right_corner = (raux -1).roots()[0]
    #corners
    x = [left_corner, wave[filters[idabv+1]], right_corner]
    y = [1, ref[filters[idabv+1]], 1]
    area = PolyArea(x,y)
    return area

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))



def fitband(photodata, method='joint', filters='r F660 i F861', method_class=0., error=PhotoError()):
    prob = pd.DataFrame(columns=['probability', 'depth'])
    filters = splitandcheck(photodata, filters)
    prob = measureband(_fitband, prob, photodata, method, filters, method_class, error)
    return prob

def _fitband(ref, wave, filters):
    # the filters for fitting the continuum are the first and last filters
    contfilters = [filters[0], filters[-1]]
    contcoef = np.polyfit(wave[contfilters], ref[contfilters], 1)
    contref = np.polyval(contcoef, wave[filters])
    # print ref[filters[1]],  contref
    depth = 1- ref[filters[1]] / contref
    return depth 

def calculate_threshold(photodata, system, filters,func,
                       flat_ref=1,
                       band_ref='/home/mario/projetos/astertools/test/testdata/spectra/91.tab',
                       errormodel=PhotoError()):
    r'''
    '''
    err = [photodata[col] for col in photodata.index if '_err' in col]
    band_threshold = BandThreshold(system,
                                   flat_ref=flat_ref,
                                   band_ref=band_ref,
                                   errormodel=errormodel)
    threshold = band_threshold.measure(err, func, filters)
    return threshold


class BandProbability(object):

    def __init__(self, system, flat_ref=1,
                 band_ref='/home/mario/projetos/astertools/test/testdata/spectra/91.tab',
                 errormodel=PhotoError()):
        r'''
        '''

        self.system = PhotometryBase(system)


        self.data = PhotoDataframe(np.full(len(self.system.ids), flat_ref),
                                   system = system,
                                   index=['flat'])
        self.data = self.data.append(self.make_band_ref(band_ref))
        self.errormodel = errormodel


    def gen_distributions(self, err, func, filters):
        r'''
        '''
        filters = filters
        if isinstance(filters, basestring):
            filters = filters.split()
        
        self.data.insert_error(err, filters)
        flat_dist = self.errormodel.distribution(self.system['center'],
                                                 self.data.loc['flat'],
                                                 filters, func)
        band_dist = self.errormodel.distribution(self.system['center'],
                                                  self.data.loc['band'],
                                                  filters, func)
        return band_dist, flat_dist

    def make_band_ref(self, band_ref):
        r'''
        '''
        photometry = Photometry(self.system.name)
        spec = loadspec(band_ref)
        photospec = photometry.convol(spec, label='band')
        return photospec

    @staticmethod
    def gaussian(x, mean, std, a=1):
        val = a * np.exp(-(x - mean)**2 / (2*std**2))
        return val


class BandThreshold(BandProbability):

    def __init__(self, system, flat_ref=1,
                 band_ref='/home/mario/projetos/astertools/test/testdata/spectra/91.tab',
                 errormodel=PhotoError()):
        r'''
        '''
        super(BandThreshold, self).__init__(system, flat_ref, band_ref, errormodel)

    def measure(self, err, func, filters):
        r'''
        '''
        # Calculating the distributions for each 
        band_dist, flat_dist = self.gen_distributions(err, func, filters)
        # treshold should be at the intersection between the distrbutions
        thresh = (band_dist.mean() - flat_dist.mean()) / 2 
        return thresh



class JointProbability(BandProbability):

    def __init__(self, system, flat_ref=1,
                 band_ref='/home/mario/projetos/astertools/test/testdata/spectra/91.tab',
                 errormodel=PhotoError()):
        r'''
        '''
        super(JointProbability, self).__init__(system, flat_ref, band_ref, errormodel)
        self.joint = np.vectorize(self.joint_aux)
        self.flat = None
        self.band = None

    def joint_aux(self, x):
        r'''
        '''
        thresh = (self.band.mean() - self.flat.mean()) / 2 
        if x > thresh:
            return self.gaussian(x, self.flat.mean(), self.flat.std())
        if x <= thresh:
            return self.gaussian(x, self.band.mean(), self.band.std())


    def integrate(self, depth):
        r'''
        '''
        return quad(self.joint_aux, -np.inf, depth)[0] / quad(self.joint_aux, -np.inf, np.inf)[0]

    def measure(self, depth, ast, func, filters):
        r'''
        '''
        # err_fil = [i+'_err' for i in filters]
        err = [ast[i+'_err'] for i in filters]
        self.band, self.flat = self.gen_distributions(err, func, filters)

        # plt.figure()
        # plt.hist(self.flat, color='r',alpha=0.5, bins=40, density=True)
        # plt.axvline(depth)
        # plt.hist(self.band, color='b',alpha=0.5, bins=40, density=True)
        # plt.xlim(-0.06, 0.04)
        # plt.show()

        # print self.band.mean(), self.band.std()
        # print self.flat.mean(), self.flat.std()
        # print

        prob = self.integrate(depth)
        return prob
