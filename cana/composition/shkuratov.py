r"""Shkuratov Model."""

import numpy as np
from .. import Spectrum
from numba import jit


class Shkuratov(object):
    r"""Shkuratov Reflectance Model.

    Adapted from Shkuratov (1999).

    """

    def __init__(self, sample=None, grainsize=30, porosity=0.5):
        r"""Init Shkuratov model class.

        Parameters
        ----------
        sample: Sample object
            The sample containing the information of the optical constant
            and grainsize

        porosity: float
            Positiy of the sample

        """
        self.sample = sample
        self.porosity = porosity
        self.grain = grainsize

    @property
    def sample(self):
        """Contain Optical constant data array."""
        return self._sample

    @sample.setter
    def sample(self, val):
        self._sample = val

    def optical_density(self):
        r"""Calculate Optical Density."""
        tau = (4 * np.pi * self.sample.k * self.grain) / self.sample.w
        return tau

    def fresnel_coef(self):
        r"""Calculate Fresnel coeficient."""
        r_0 = (self.sample.n - 1) ** 2 / (self.sample.n + 1) ** 2
        return r_0

    def scattering_coef(self):
        r"""Scatering coeficients based on Shkuratov (1999) approximations."""
        # calculating internal_particle_reflection
        r_b, r_f = self.scattering_coef_aux(
            self.sample.n, self.sample.k, self.grain, self.sample.w
        )
        # creating output
        scat = np.zeros(len(r_b), dtype=[("b", np.float64), ("f", np.float64)])
        scat["b"] = r_b
        scat["f"] = r_f
        return scat

    @staticmethod
    @jit(nopython=True)
    def scattering_coef_aux(n, k, grain, w):
        r"""Calculate the Scatering coeficient."""
        tau = (4 * np.pi * k * grain) / w
        r0 = (n - 1) ** 2 / (n + 1) ** 2
        Ri = 1.04 - (1 / n ** 2)
        # Using approximations for reflection
        Re = r0 + 0.05
        Rb = (0.28 * n - 0.2) * Re
        Rf = r0 + 0.05 - (0.28 * n - 0.2) * Re
        # calculating transmision coeficients
        Te = 1 - Re
        Ti = 1 - Ri
        # average light scattering indicatrix of a particle
        rb = Rb + 0.5 * Te * Ti * Ri * np.exp(-2 * tau) / (1 - Ri * np.exp(-1 * tau))
        rf = (
            Rf
            + Ti * Te * (np.exp(-1 * tau))
            + (0.5 * Te * Ti * Ri * np.exp(-2 * tau)) / (1 - Ri * np.exp(-1 * tau))
        )
        # print(rb,rf)
        return rb, rf

    def build_spec(self, coef=None, albedo_w=0.55, wavelengths=None):
        r"""Build the modeled spectra from the optical constants.

        Returns
        -------
        Spectrum, Geometric Albedo
        """
        if wavelengths is not None:
            self.sample = self.sample.rebase(wavelengths)
        if coef is None:
            coef = self.scattering_coef()
        # albedo indicatrix
        rho_b = self.porosity * coef["b"]
        rho_f = (self.porosity * coef["f"]) + 1 - self.porosity
        # building albedo
        aux = (1 + rho_b ** 2 - rho_f ** 2) / (2 * rho_b)
        albedo = aux - np.sqrt(aux ** 2 - 1)
        # building Spectrum
        spec = Spectrum(self.sample.w, albedo, label="modeled_spec", unit="micron")
        # return albedo
        if albedo_w is not None:
            geom_alb = np.interp(albedo_w, self.sample.w, albedo)
        else:
            geom_alb = None
        return spec, geom_alb
