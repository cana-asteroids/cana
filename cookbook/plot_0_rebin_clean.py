"""
Rebining and cleaning the spectrum from bad points
--------------------------------------------------
"""

import matplotlib.pyplot as plt
from cana.datasets import getspectrum
import cana

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.

fig, ax = plt.subplots(1,3, figsize=(15,5), sharey=True)

spec = getspectrum('000752b', ref='primass')

spec.plot(fax=ax[0], speckwargs={'label':'Raw'})

# Rebining spectrum

spec_rebin = spec.rebin(binsize=11)

spec_rebin.plot(fax=ax[1], speckwargs={'label':'Rebinned'})

# Using sigma clip to remove bad points from spectrum

spec_clean = spec.clean_spec(method='sigmaclip', sigma=3, fit='auto')

spec_clean.plot(fax=ax[2], speckwargs={'label':'Cleaned'})

plt.show()




