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

# Measuring the slope
# Defaults: wmin=0.4, wmax=0.9, norm=0.55, errormethod='rms',
#           error_param=None, montecarlo=1000, speckwargs=None
slope = cana.slope(spec)

# print the results
print(slope)

# plot the results
slope.plot(show=True)

# Note: cana.slope can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.
