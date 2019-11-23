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

## Define region to measure the band
# The wavelength where the band starts
wmin = 0.55
# The wavelength where the band ends
wmax = 0.85
# The wavelength windows to calculate the continuun.
# A value of 0.03 means that the continuum will be measure at 0.55-0.58 and 0.83-0.85 range
continumn_window = 0.03


# Methodology to estimate uncertainty
errormethod = 'rms'
iterations = 1000

# Calculating the band parameters
band_depth = cana.depth(spec, wmin=wmin, wmax=wmax, cont_window=continumn_window,
                        errormethod=errormethod, montecarlo=iterations)

# print the calculated band parameters
print(band_depth)
# visualizing the band
band_depth.plot()

# check if the values can be considered an absortium band.
minimum_depth = 1.             # The minimum depth value for it to be considered a band
theoric_minimum_position = 0.7 # The central wavelengh of the theorical band
maximum_distance = 0.05        # The maximum distance from the calculated center to the theorical center
sigma_level = 3                # The sigma level from the depth and the spetrum noise

print('Is an absortion band detected?')
print(band_depth.is_band(min_depth=minimum_depth, theoric_min=theoric_minimum_position,
                         max_dist=maximum_distance, sigma=sigma_level))


# Note: cana.depth can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.

###############################################
# You can also try a object-oriented approach
###############################################

# Defining methodoly for the continumn
continumnmodel = cana.Continuum(lowerwindow=continumn_window, upperwindow=continumn_window)

# creating the depth model
depthmodel = cana.Depth(wmin=wmin, wmax=wmax, continuum=continumnmodel)

# creating the depth model for estimating the parameters uncertainties

errormodel = cana.spectools.SpecError(method='rms', n=1000)

# Calculating the band parameters
band_depth = depthmodel.measure(spec, error=errormodel)
