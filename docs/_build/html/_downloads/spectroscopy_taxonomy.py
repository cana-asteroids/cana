"""
Spectroscopy: Taxonomic classification
"""

from cana.datasets import getspectrum
import cana.spec.taxutils as canatax

# First load an spectrum
# you can do: spec = cana.spec.loadspec('path to your spectrum file')

# For the example, we just gonna use one from the available datasets.
# To see what else is available at the datasets check the datasets documention.

spec = getspectrum('38', ref='smass')

# Defining parameters for taxonomic classitication
taxonomy_system = 'demeo'
classification_method = 'chisquared'
return_n = 3
fitspec = True

# For a simple aproach

tax = canatax.taxonomy(spec, system=taxonomy_system, method=classification_method,
                       return_n=return_n, fitspec=fitspec)

# print the results
print tax
# plot the results
tax.plot()


###############################################
# TIP:
# Check out the *cana.spec.taxonomy_list* method
# for a direct application in a list of spectra files


