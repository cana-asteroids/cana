"""
Pipeline for Visible Spectra
----------------------------
"""

import cana

# First load an spectrum
# you can do: spec = cana.loadspec('path to your spectrum file')

# For the example, we will just gonna use one from the available datasets.

spec = cana.datasets.getspectrum('000752', ref='primass')

# Defining parameters that will be used in the pipeline
# We are explicitating here, but all values above are actually set as default
# and can simply be setted as params='default'

# The Slope method
slope = cana.Slope(wmin=0.5, wmax=0.9, norm=0.55)
# The Taxonomy method
tax = cana.Taxonomy(system='bus')
# The band method
continuum = cana.Continuum(lowerwindow=0.03, upperwindow=0.03)
band = cana.Depth(wmin=0.55, wmax=0.85, continuum=continuum)
isband = {'min_depth': 1, 'theoric_min': 0.7,
          'max_dist': 0.5, 'sigma': 3}

# Running the pipeline
# if spec is a list of spectra path, outplot can be setted as
# a directory to save the plots

pipe = cana.primitive_visible(spec, params=[slope, tax, band],
                              isband=isband, outplot='show')

# Printing results
print(pipe)
