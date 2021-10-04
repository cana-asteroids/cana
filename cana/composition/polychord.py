import os
from functools import partial
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
import corner

import pypolychord
from pypolychord.settings import PolyChordSettings
from .mixtures import IntimateMixture



class PolyChord:
    r"""
    """

    def __init__(self, constants, grains=[12, 20], outdir='model',
                 nlive=1000, do_clustering=True, read_resume=False):
        r"""
        """
        self.constants = constants
        self.grains = grains
        self.proportions = [0, 1]
        self.ndims = 2 * len(self.constants)
        self.nderived = 1
        self.settings = PolyChordSettings(self.ndims, self.nderived)
        self.settings.file_root = outdir
        self.settings.nlive = nlive
        self.settings.do_clustering = do_clustering
        self.settings.read_resume = read_resume
        self.metric = None

    @property
    def nconstants(self):
        r"""Return the number of optical constants."""
        return len(self.constants)

    @staticmethod
    def chi_like(spec, mspec):
        r"""Calculate the likelihood from the chi-square."""
        res_aux = (spec.r - mspec.r)**2
        err = 1. / spec.r_unc**2
        return -0.5*(np.sum((res_aux*err) - np.log(err)))

    @staticmethod
    def weighted_chi(spec, mspec, wranges=None, weigths=None):
        r"""Calculate a weighted chi square.

        Parameters
        ----------
        spec: Spectrum
            The spectrum to be modeled

        mspec: Spectrum
            The spectrum of a mixture

        wranges: list
            A list of wavelength ranges that will be used to apply a weigth

        weigths:list
            The list of weights corresponding to the wavelength
            ranges (wranges)
        """
        weight_arr = np.ones(len(mspec.w))
        weight_sum = np.sum(weigths)
        for i, wave in enumerate(wranges):
            aux = np.ravel(np.argwhere((mspec.w >= wave[0]) & (mspec.w <= wave[1])))
            weight_arr[aux] = np.array([weigths[i]/weight_sum for _ in range(len(aux))])
        weight_arr = (weight_arr / np.sum(weight_arr))*len(mspec.w)
        res_aux = ((spec.r - mspec.r)**2) * weight_arr
        res_aux = (spec.r - mspec.r)**2
        err = 1. / spec.r_unc**2 * weight_arr
        return -0.5*(np.sum((res_aux*err) - np.log(err)))

    def likelihood(self, theta, spec, norm, albedo_lim, **kwargs):
        """Calculate the Likelihood."""
        r_sqr = sum(theta**2)
        proportions = theta[:self.nconstants]
        grains_ = theta[self.nconstants:]
        grains = (grains_ * (self.grains[1] - self.grains[0])) + self.grains[0]
        mix = IntimateMixture(self.constants,
                              proportions=proportions,
                              grainsizes=grains)
        mspec, alb_mix = mix.make(wavelengths=spec.w)
        mspec = mspec.normalize(norm)
        if (alb_mix >= albedo_lim[0]) & (alb_mix <= albedo_lim[1]):
            log_like = self.metric(spec, mspec, **kwargs)
        else:
            log_like = -np.inf
        return log_like, [r_sqr]

    def prior(self, hypercube):
        """Return prior."""
        # print(hypercube, self.prior_transform(hypercube))
        return self.prior_transform(hypercube)

    @staticmethod
    def dumper(live, dead, logweights, logZ, logZerr):
        print("Last dead point:", dead[-1])

    def prior_transform(self, values):
        r"""Transform the prior."""
        props = values[:self.nconstants]

        props = -np.log(props)
        theta_p = props / np.sum(props)
        grains = np.array(values[self.nconstants:])
        theta = np.concatenate([theta_p, grains])
        return theta

    def _likelihood(self, spec, norm, albedo_lim, **kwargs):
        r""" """
        return partial(self.likelihood, spec=spec, norm=norm,
                       albedo_lim= albedo_lim)

    def run(self, spec, norm=0.5, albedo_lim=None, metric='chi',
           percentiles=[16, 50, 84], sigma=1, **kwargs):
        r"""
        """
        if metric == 'chi':
            self.metric = self.chi_like
        elif metric == 'wchi':
            self.metric = partial(self.weighted_chi, **kwargs)
        likelihood = self._likelihood(spec, norm, albedo_lim, **kwargs)

        output = pypolychord.run_polychord(likelihood,
                                           self.ndims,
                                           self.nderived,
                                           self.settings,
                                           self.prior,
                                           self.dumper)
        return PolyChordOut(output, self.constants, self.grains,
                            percentiles, sigma)
        # return output


class PolyChordOut:
    r"""
    """

    def __init__(self, output, constants, grains, percentiles, sigma):
        r"""
        """
        self.output = output
        self.grains = grains
        self.constants = constants
        self.samples, self.labels = self.rename_cols()
        self.constants_labels = self.labels[:self.nconstants]
        self.grains_labels = self.labels[self.nconstants:]
        self.convert_grains()

        self.parameters = self.set_result(percentiles)

    def __repr__(self):
        return self.parameters.__repr__()

    def _repr_html_(self):
        return self.parameters._repr_html_()

    def convert_grains(self):
        r"""
        """
        self.samples[self.grains_labels] = (self.samples[self.grains_labels] *
                                            (self.grains[1] - self.grains[0])) + self.grains[0]



    @property
    def nconstants(self):
        return len(self.constants)

    def rename_cols(self):
        r"""
        """
        params = [f'p{i}' for i in range(len(self.constants)*2)]
        labels = [i.label for i in self.constants]
        labels_grains = [f'grain_{i.label}' for i in self.constants]
        labels.extend(labels_grains)
        cols_aux = zip(params, labels)
        new_cols = dict(cols_aux)
        samples = self.output.samples.rename(columns=new_cols)
        return samples, labels


    def set_result(self, percentiles=[16, 50, 84]):
        r"""
        """
        values_aux = self.samples[self.labels].to_numpy()
        values_aux = values_aux.reshape(-1, 2 * self.nconstants)
        values = list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                          zip(*np.percentile(values_aux,
                                             percentiles, axis=0))))

        best_values = pd.DataFrame(columns=['best_fit', 'sigma_low', 'sigma_upp'])
        for i, val in enumerate(values):
            best_values.loc[self.labels[i]] = val
        return best_values

    def spec_from_index(self, index, wavelengths, **kwargs):
        r"""
        """
        proportions = self.samples[self.labels[:self.nconstants]].loc[index]
        grains = self.samples[self.labels[self.nconstants:]].loc[index]

        mix = IntimateMixture(self.constants,
                              proportions=proportions,
                              grainsizes=grains)
        mspec, alb = mix.make(wavelengths=wavelengths, **kwargs)
        return mspec, alb

    def select_sigma(self, sigma=1):
        r"""
        """
        out = self.samples.copy()
        for label in self.labels:
            out = out[(out[label] > self.parameters.loc[label, 'best_fit'] -
                       sigma * self.parameters.loc[label, 'sigma_low']) &
                      (out[label] < self.parameters.loc[label, 'best_fit'] +
                       sigma * self.parameters.loc[label, 'sigma_upp'])]
        return out

    def plot_sigma_spectrum(self, sigma=1, wavelengths=None,
                            norm=0.5, **kwargs):
        r"""
        """
        ax = plt.gca()
        sigma_mixes = self.select_sigma(sigma)
        for index in sigma_mixes.index:
            spec, alb = self.spec_from_index(index, wavelengths, **kwargs)
            if norm is not None:
                spec = spec.normalize(norm)
            spec.plot(fax=ax, show=False, speckwargs={'c':'0.1', 'alpha':0.5})

    def corner_proportions(self, **kwargs):
        r"""
        """
        proportions = self.samples[self.labels[:self.nconstants]].to_numpy()
        fig = corner.corner(proportions,
                            labels=self.labels[:self.nconstants],
                            **kwargs)
        return fig

    def corner_grains(self, **kwargs):
        r"""
        """
        grains = self.samples[self.grains_labels].to_numpy()
        fig = corner.corner(grains,
                            labels=self.grains_labels,
                            **kwargs)
        return fig

    def corner_all(self, **kwargs):
        r"""
        """
        params = self.samples[self.labels[self.nconstants:]].to_numpy()
        fig = corner.corner(params,
                            labels=self.labels[:self.nconstants],
                            **kwargs)
        return fig
