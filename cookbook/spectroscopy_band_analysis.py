"""
Spectroscopy: Calculating the center and depth of an absortium band
"""

from cana.datasets import getspectrum
import cana.bandutils as canaband
from cana.errormodels import SpecError

# First load an spectrum
# you can do: spec = cana.spec.loadspec('path to your spectrum file')
#             or simple pass the file path or a list of files.

# For the example, we just gonna use one from the available datasets.
# To see what else is available at the datasets check the datasets documention.

spec = getspectrum('38', ref='smass')

# For a simple approach
band_starts = 0.55       #-> the wavelength where the band starts
band_ends = 0.85         #-> the wavelength where the band ends
continumn_window = 0.03  #-> The wavelength windows to calculate the continumn
errormethod = 'rms'      #-> the method for estimating the parameters uncertainties
errorrepetitions = 1000  #-> The number of repetitions in the Monte-Carlo method to estimate the parameters uncertainties

# Calculating the band parameters
band_depth = canaband.depth(spec, wmin=band_starts, wmax=band_ends, cont_window=continumn_window,
                            errormethod=errormethod, montecarlo=errorrepetitions)

# print the calculated band parameters
print band_depth
# visualizing the band
band_depth.plot()

# check if the values can be considered an absortium band.
minimum_depth = 1.             # The minimum depth value for it to be considered a band
theoric_minimum_position = 0.7 # The central wavelengh of the theorical band
maximum_distance = 0.0         # The maximum distance from the calculated center to the theorical center
sigma_level = 3                # The sigma level from the depth and the spetrum noise

print band_depth.is_band(min_depth=minimum_depth, theoric_min=theoric_minimum_position,
                         max_dist=maximum_distance, sigma=sigma_level)


###############################################
# You can also try a object-oriented approach
###############################################

# Defining methodoly for the continumn
continumnmodel = canaband.Continuum(lowerwindow=continumn_window, upperwindow=continumn_window)

# creating the depth model
depthmodel = canaband.Depth(wmin=band_starts, wmax=band_ends, continuum=continumnmodel)

# creating the depth model for estimating the parameters uncertainties

errormodel = SpecError(method='rms', n=1000)

# Calculating the band parameters
band_depth = depthmodel.measure(spec, error=errormodel)

# print the calculated band parameters
print band_depth
# visualizing the band
band_depth.plot()
# check if the values can be considered an absortium band.
print band_depth.is_band(min_depth=minimum_depth, theoric_min=theoric_minimum_position,
                         max_dist=maximum_distance, sigma=sigma_level)

###############################################
# TIP:
# Check out the *cana.spec.depth_list* method
# for a direct application in a list of spectra files


