"""
Taxonomic classification
------------------------
"""

import cana

# First load an spectrum, we will just gonna use one from the available datasets.
# you can do: spec = cana.loadspec('path to your spectrum file')

spec = cana.datasets.getspectrum('000334', ref='primass')

# Perform taxonomic classification
# Defults: system='demeo', method='chisquared', return_n=3, norm=0.55,
#             fitspec=True, speckwargs=None
tax = cana.taxonomy(spec)

# print the results
print(tax)

# plot the results
tax.plot()

# Note: cana.slope can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.
