"""
Autofit a Spectra
------------------
"""

import cana
import matplotlib.pyplot as plt

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.

spec = cana.datasets.getspectrum('000334', ref='primass')

# fitting the spectrum, you can give the minimal and maximal order that could
# be fitted. Default is 1 and 12.
specfit, coefs = spec.autofit(degree_min=1, degree_max=12)

# Ploting
plt.plot(spec.w, spec.r, c='0.3', lw=1, label='Spectrum')
plt.plot(specfit.w, specfit.r, c='r', label='Fit')
plt.legend()
plt.show
