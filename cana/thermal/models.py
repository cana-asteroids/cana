r"""Thermal Models."""

import numpy as np
from scipy.integrate import dblquad
from .. import Spectrum
sig = 5.670373 * (10**-8)          # Stefan-Boltzmann Constant [W/m^2/K^4]
SUNc = 1367.567          # solar constant at 1 AU in [W/m^2]
AU = 1.49597870 * (10**8)        # AU in [km]
k_bol = 1.380622 * (10**-23)     # Boltzmann Konstante k = 1.4e-16 erg/grad
c_light = 2.99792458 * (10**8)  # speed of light [m/s]
h_planck = 6.6260755 * (10**-34)	 # Planck constante Js


class ThermalCore(object):
    r"""
    """

    def __init__(self, albedo, eta=0.9, epsilon=0.9, G=0.15):
        r"""
        """
        self.eta = eta
        self.epsilon = epsilon
        self.albedo = albedo
        self.G = G

    def phase_integral(self):
        r"""Phase Integral.

        See Bowell et al. (1989) for a reference.

        Parameters
        ----------
        gslope: float
            The phase slope

        Return
        ------
        phase integral
        """
        return 0.29 + 0.684 * self.G

    def bolometric_albedo(self):
        r"""Bolometric Bond Albedo.

        Uses the geometric albedo and the phase integral to calculate The
        bolometric bond albedo.
        """
        return self.albedo * self.phase_integral()

    def subsolar_temperature(self, r):
        r"""Calculate the sub-solar point temperature."""
        A = self.bolometric_albedo()
        Tss = (SUNc * (1-A) / ((r**2) * self.epsilon * sig * self.eta))**0.25
        return Tss


class NEATM(ThermalCore):
    r"""The Near Earth Asteroid Thermal Model.

    Reference: Harris (1998, Icarus, 131, 291-301).
    """

    def __init__(self, r, delta, albedo, diameter, phase_angle=0,
                 eta=0.9, epsilon=0.9, G=0.15):
        r"""Initialize the model.

        Parameters
        ----------
        r: float
            Heliocentric  distance in unit of AU.

        delta: float
            Geocentric distance in unit of AU.

        albedo: float
            Geometric albedo

        diameter: float
            Diameter of the asteroid in unit of km

        phase_angle: float
            Phase angle in unit of degrees

        eta: float
            the NEATM eta parameter
        """
        ThermalCore.__init__(self, albedo, eta=eta,
                             epsilon=epsilon, G=G)
        self.r = r
        self.D = diameter
        self.delta = delta
        self.phase_angle = phase_angle
        self.alpha = np.radians(phase_angle)
        self.T0 = self.subsolar_temperature(r)
        self._integral = np.vectorize(self._integral_aux)

    @staticmethod
    def planckian(wavelength, temp):
        r"""Calculate the planck function.

        Parameters
        ----------
        wavelength: float
            The lambda in unit of microns.

        temp: float
            The temperature in unit of Kelvin.

        Returns
        -------
        The black-body spectral radiance
        """
        wavelength = wavelength * 1e-6
        aux = np.exp(h_planck * c_light / k_bol / wavelength / temp)
        B = 2 * h_planck * c_light / (wavelength**3 * (aux - 1.0))
        return B

    def _integral_aux(self, wavelength):
        r"""Auxiliary method to calculate the integral per wavelength."""
        aux = dblquad(self.flux_aux, 0.0, np.pi/2.0,
                      lambda phi: -np.pi/2.0 + phi, lambda phi: np.pi/2.0,
                      args=(wavelength, self.alpha))[0]
        return aux

    def flux_aux(self, phi, theta, wavelength, alpha):
        r"""Auxiliary method to calculate the flux."""
        T = self.T0 * np.cos(phi)**0.25 * np.cos(theta)**0.25
        B = self.planckian(wavelength, T)
        a = B * np.pi * np.cos(phi)**2 * np.cos(theta - self.alpha)
        return a

    def flux(self, wavelengths=[0.5, 2.5, 4.5]):
        r"""Calculate the thermal flux.

        Parameters
        ----------
        wavelengths: list or numpy.ndarray
            The wavelength vector to measure the thermal flux

        Returns
        -------
        ThermalSpec
            The thermal spectrum array
        """
        flux = self._integral_aux(wavelengths)
        # fixing units
        norm = (self.epsilon * (self.D / self.delta)**2
                / np.pi / 2.0) * 4.468370499519714 * 1e9

        flux *= norm

        return Spectrum(wavelengths, flux, label='thermal-flux')


class ThermalSpec(Spectrum):
    r""" Thermal Spectrum."""
