"""
Spectroscopy: Calculating the spectral gradient
"""

from cana.datasets import getspectrum
from cana.errormodels import SpecError
import cana

# First load an spectrum
# you can do: spec = cana.spec.loadspec('path to your spectrum file')
#             or simple pass the file path or a list of files.

# For the example, we just gonna use one from the available datasets.
# To see what else is available at the datasets check the datasets documention.

spec = getspectrum('38', ref='smass')

# Defining parameters for calculating the slope
wavemin = 0.4
wavemax = 0.9
wavenorm = 0.55
errormethod = 'rms'      #-> the method for estimating the parameters uncertainties
errorrepetitions = 1000  #-> The number of repetitions in the Monte-Carlo method to estimate the parameters uncertainties


# For a simple aproach
astslope = cana.slope(spec, wmin=wavemin, wmax=wavemax, norm=wavenorm, errormethod=errormethod, montecarlo=errorrepetitions)

# print the results
print astslope
# plot the results
astslope.plot(show=True) 
