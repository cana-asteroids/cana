"""
Calculating the center and depth of an absortium band
-----------------------------------------------------
"""

from cana.datasets import getspectrum
import cana

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.


spec = getspectrum('000752', ref='primass')

# Define region to measure the band
# The wavelength where the band starts
wmin = 0.55
# The wavelength where the band ends
wmax = 0.85
# The wavelength windows to calculate the continuun.
# A value of 0.03 means that the continuum will be measure
#  at 0.55-0.58 and 0.83-0.85 range
c_window = 0.03

# Calculating the band parameters
# Defaults: wmin=0.54, wmax=0.88, cont_window=0.04, resolution='auto',
#          errormethod='rms', error_param=None, montecarlo=1000,
#          min_depth=1., theoric_min=0.7, max_dist=0.05, n_sigma=3,
#          speckwargs=None
band_depth = cana.depth(spec, wmin=wmin, wmax=wmax, cont_window=c_window)

# print the calculated band parameters
print(band_depth)
# visualizing the band
band_depth.plot()

# check if the values can be considered an absortium band.
# is_band Defaults: min_depth=1., theoric_min=0.7, max_dist=0.05, sigma=3
print('Is an absortion band detected?')
print(band_depth.is_band())


# Note: cana.depth can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.

###############################################
# You can also try a object-oriented approach
###############################################

# Defining methodoly for the continumn
continumnmodel = cana.Continuum(lowerwindow=c_window, upperwindow=c_window)

# creating the depth model
depthmodel = cana.Depth(wmin=wmin, wmax=wmax, continuum=continumnmodel)

# creating the depth model for estimating the parameters uncertainties

errormodel = cana.spectools.SpecError(method='rms', n=1000)

# Calculating the band parameters
band_depth = depthmodel.measure(spec, error=errormodel)
