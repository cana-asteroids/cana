"""
Calculating the spectral gradient
---------------------------------
"""

from cana.datasets import getspectrum
import cana

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.


spec = getspectrum('000334', ref='primass')

# Defining parameters for calculating the slope
wavemin = 0.4
wavemax = 0.9
wavenorm = 0.55


# Methodology to estimate uncertainty
errormethod = 'rms'
iterations = 1000

# Measuring the slope
slope = cana.slope(spec, wmin=wavemin, wmax=wavemax, norm=wavenorm,
                   errormethod=errormethod, montecarlo=iterations)

# print the results
print(slope)

# plot the results
slope.plot(show=True)

# Note: cana.slope can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.
