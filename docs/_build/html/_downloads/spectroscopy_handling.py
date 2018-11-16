"""
Spectroscopy: Loading and Visualizing a Spectrum
"""

from cana.spec import loadspec
from cana.datasets import getspec

#Loading an spectra from the dataset
spec = getspec('38', ref='smass')

spec = loadspec()
