"""
Spectroscopy: Calculating the slope from photometric data
"""

from cana.datasets import getspectrum
import cana
from cana.errormodels import SpecError

# First load an spectrum
# you can do: spec = cana.spec.loadspec('path to your spectrum file')
#             or simple pass the file path or a list of files.

# For the example, we just gonna use one from the available datasets.
# To see what else is available at the datasets check the datasets documention.

spec = getspectrum('38', ref='smass')



slope = cana.slope(spec)

