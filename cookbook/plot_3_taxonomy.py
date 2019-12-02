"""
Taxonomic classification
------------------------
"""

import cana

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.

spec = cana.datasets.getspectrum('000334', ref='primass')

# Defining parameters for taxonomic classitication
taxonomy_system = 'bus'
classification_method = 'chisquared'
return_n = 3
fitspec = True


# Perform taxonomic classification
tax = cana.taxonomy(spec, system=taxonomy_system, method=classification_method,
                    return_n=return_n, fitspec=fitspec)

# print the results
print(tax)

# plot the results
tax.plot()

# Note: cana.slope can also get the direct spectrum file path,
# or a list of spectra file. In this case it would return the a
# pandas table with results.
