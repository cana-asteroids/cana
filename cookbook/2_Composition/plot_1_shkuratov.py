"""
Create a spectrum using Shkuratov model
---------------------------------------
"""

import sdoc
from cana.composition import read_constant, Shkuratov
import numpy as np

# First load an spectrum, we will just gonna use one from the SDOC database.
# See more about SDOC at https://https://github.com/depra/sdoc
# Alternatively, you can read your own optical constant using:
# oc = cana.read_constant('path to your optical constant file')

# First initialize SDOC
sdb = sdoc.SDOC(mode='r')

# To see what is inside the database: sdb.contents
# Here we will just show an example for a titan tholin
oc_label, oc_data = sdb.get_constant('T_0')

oc = read_constant(oc_data, label=oc_label)

# Initializing the model
shkuratov = Shkuratov(sample=oc, grainsize=30, porosity=0.5)

# Build the wavelentgh axis
w = np.linspace(0.4, 2.3, 200)

spec, albedo = shkuratov.build_spec(wavelengths=w)
# Plot the optical constant
spec.plot()

# Print the spectrum and albedo
print(spec)
print('Albedo: ', albedo)
