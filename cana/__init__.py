from .spec import loadspec, Spectrum, stack_spec
from . import spectools
from .spectools import slope, Slope, taxonomy, Taxonomy, \
                       depth, Depth, Continuum
from .pipelines import primitive_visible
from .util import curvature

from .photo import PhotoDataframe
from .convolution import convolution, interp_spec
from . import datasets
from .composition import OpticalConstant, read_constant, Shkuratov, \
                         IntimateMixture, read_constant_batch \

from .thermal import NEATM
