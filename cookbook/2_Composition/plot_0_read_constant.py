"""
Read an optical constant
------------------------
"""

import sdoc
from cana.composition import read_constant


# First load an spectrum, we will just gonna use one from the SDOC database.
# See more about SDOC at https://https://github.com/depra/sdoc
# Alternatively, you can read your own optical constant using:
# oc = cana.read_constant('path to your optical constant file')

# First initialize SDOC
sdb = sdoc.SDOC(mode='r')

# To see what is inside the database: sdb.contents
# Here we will just show an example for a ice tholin
oc_label, oc_data = sdb.get_constant('T_0')

oc = read_constant(oc_data, label=oc_label)

print(oc)

# Plot the optical constant
oc.plot()
