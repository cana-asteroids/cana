"""
Create an Intimate Mixture using Shkuratov model
------------------------------------------------
"""

import sdoc
from cana.composition import read_constant_batch, IntimateMixture
import numpy as np

# First load an spectrum, we will just gonna use one from the SDOC database.
# See more about SDOC at https://https://github.com/depra/sdoc
# Alternatively, you can read your own optical constant using:
# oc = cana.read_constant('path to your optical constant file')

# First initialize SDOC
sdb = sdoc.SDOC(mode='r')

# To see what is inside the database: sdb.contents
# Here we will just show an example for a titan tholin and ice tholin
ocs_label, ocs_data = sdb.get_constant_batch(['T_0', 'I_0'])

ocs = read_constant_batch(ocs_data, label=ocs_label)

# Initializing the model
mixture = IntimateMixture(samples=ocs,
                          grainsizes=[30, 100],
                          porosity=0.5)

# Build the wavelentgh axis
w = np.linspace(0.4, 2.3, 200)

spec, albedo = mixture.make(wavelengths=w)
# Plot the optical constant
spec.plot()

# Print the spectrum and albedo
print(spec)
print('Albedo: ', albedo)
