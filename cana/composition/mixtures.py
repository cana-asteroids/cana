r"""Mixtures used in the Shkuratov model."""

import copy
import numpy as np

from .shkuratov import Shkuratov
from .. import interp_spec


def create_inclusion(constant1, constant2, concentration):
    r"""
    """
    #interpolating second constant
    constant2_ = constant2.rebase(constant1.w)
    ime = (constant1.k * constant2_.n) - (constant1.n * constant2_.k) / \
          (constant2_.n**2 + constant2_.k**2)
    # ime =
    ka = 1.5 * concentration * constant1.n * ime
    constant1_ = constant1.copy()
    constant1_.k += ka
    return constant1_


class IntimateMixture(object):
    r"""Mixture of samples."""

    def __init__(self, samples, grainsizes=[30,50], proportions=[0.5, 0.5],
                 porosity=0.5, model='shkuratov'):
        r"""Initialize the Initimate Mixture class.

        Parameters
        ----------
        samples: list
            List containing the optical constants Sample

        proportions: list
            List of proportions in the same order of the samples List

        porosity: float
            The value for porosity

        model: string
            The name of the reflectance model to be used

        """
        self.modelname = model.capitalize()
        self.samples = samples
        self.proportions = proportions
        self.porosity = porosity
        self.samples_ids = ','.join([sam.label for sam in self.samples])
        self.grainsizes = grainsizes
        self.model = self.select_model()

    def __str__(self):
        r"""Representation of the IntimateMixture class."""
        grains = ','.join([str(g) for g in self.grainsizes])
        prop = ','.join([str(p) for p in self.proportions])
        return '''IntimateMixture(samples:({0}), grainsizes:({1}),
                proportions=({2}), porosity:{3})'''.format(self.samples_ids,
                                                           grains, prop,
                                                           self.porosity)

    def __repr__(self):
        r"""Representation of the IntimateMixture class."""
        return self.__str__()

    def copy(self):
        r"""Return a copy of the object."""
        return copy.copy(self)

    def select_model(self):
        r"""Select the reflectance model that will be used."""
        if self.modelname == 'Shkuratov':
            model = Shkuratov
        return model

    def make(self, albedo_w=0.55, wavelengths=None,
             datatype=None):
        r"""Make the reflectance spectrum of the mixture.

        Parameters
        ----------
        rebase: boolean
            True if you want to rebase the wavelength axis of the samples.
            If False, the provided samples should have equal wavelentghs.

        interpmethod:
            Method that will be used to rebase the samples.
            Default is 'interpolate'.

        baseaxis: None or array
            The wavelength array that will be used if rebase=True.
            If baseaxis=None and rebase=True, it will select the wavelength
            of the first sample as a reference. Not used if rebase=False.

        """
        if wavelengths is None:
            wavelengths = self.samples[0].w
        # building coef structure
        coef = np.zeros(len(wavelengths),
                        dtype=[('b', np.float64), ('f', np.float64)])
        # calculating coeficients
        for s, sam in enumerate(self.samples):
            # no need to calculate coeficients if proportion is zero
            if self.proportions[s] != 0.:
                # print(sam,self.proportions[s],self.grainsizes[s])
                sam = sam.rebase(baseaxis=wavelengths)
                model = self.model(sample=sam, grainsize=self.grainsizes[s],
                                   porosity=self.porosity)
                coef_aux = model.scattering_coef()
                # print(coef_aux)
                # summing coeficients
                coef['b'] += coef_aux['b']*self.proportions[s]
                coef['f'] += coef_aux['f']*self.proportions[s]
        spec, albedo = model.build_spec(coef=coef,
                                        albedo_w=albedo_w)
        spec = interp_spec(spec, baseaxis=wavelengths, datatype=datatype)
        return spec, albedo
