"""
Autofit a Spectra
------------------
"""

import cana
import matplotlib.pyplot as plt

# First load an spectrum, we will just gonna use one from the available datasets.
# you can do: spec = cana.loadspec('path to your spectrum file')

# See spec.py Spectrum class for spec attributes
spec = cana.datasets.getspectrum('000334', ref='primass')

# fitting the spectrum, you can give the minimal and maximal order that could
# be fitted. Default: degree_min=1, degree_max=12.
specfit, coefs = spec.autofit()

# Plotting
plt.plot(spec.w, spec.r, c='0.3', lw=1, label='Spectrum')
plt.plot(specfit.w, specfit.r, c='r', label='Fit')
plt.legend()
plt.show