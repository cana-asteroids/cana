"""
Spectroscopy: Calculating the spectral gradient
"""

from cana.datasets import getspectrum
import cana.spec.slopeutils as canaslp
from cana.spec.errormodels import MonteCarlo

# First load an spectrum
# you can do: spec = cana.spec.loadspec('path to your spectrum file')

# For the example, we just gonna use one from the available datasets.
# To see what else is available at the datasets check the datasets documention.

spec = getspectrum('38', ref='smass')

# Defining parameters for calculating the slope
wavemin = 0.4
wavemax = 0.9
wavenorm = 0.55
errormodel = MonteCarlo()

# For a simple aproach

slp = canaslp.slope(spec, wmin=wavemin, wmax=wavemax, norm=wavenorm, error=errormodel)

# print the results
print slp
# plot the results
slp.plot(show=True)


###############################################
# TIP:
# Check out the *cana.spec.slope_list* method
# for a direct application in a list of spectra files
