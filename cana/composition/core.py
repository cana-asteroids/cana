r"""Compositional model core classes."""

from scipy.stats import norm, chi2
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from .mixtures import IntimateMixture


class BaseCompositionalModel(object):
    r"""Compositional Model core class."""

    def __init__(self, samples, grainsizes=[30, 50],
                 proportions=[0, 1], porosity=0.5,
                 mixtype='intimate',
                 model='shkuratov', metric='chi2-error'):
        r"""Initialize Compositional model.

        Parameters
        ----------
        samples: list
            List containing the optical constants Sample

        grainsizes: list
            Grainsize Prior. List the range of which the grainsize can vary.
            Default is [30, 50].

        proportions: list
            Proportions prior. List range of which proportions can vary.
            Default is [0, 1], which means that pure materials can tested.

        porosity: float
            The value for porosity. #-> Need to modify to prior later

        mixtype: string
            Type of mixture that will be used to model the spectrum.

        model: string
            The name of the reflectance model to be used

        method: string
            The name of the metric that will be used to evaluate the fitness
            of the model. Default is chisquare (scipy.stats.chisquare)

        """
        # defining mixtures parameters
        self.refmodel = model
        self.mixtype = mixtype
        self.porosity = porosity
        self.samples = [sam.copy() for sam in samples]
        self.samples_ids = [sam.label for sam in samples]
        self.grains = grainsizes
        self.grains = np.asarray(self.grains)

        self.proportions = proportions
        self.proportions = np.asarray(self.proportions)

        # model selection criteria
        self.metric = self.fitness_method(metric)
        # defining mixture type
        self.mixtype = self.setmixturetype(mixtype)

    def eval(self, spec, mspec, **kwargs):
        r"""Evaluate spectrum fitness according to the defined metric."""
        val = self.metric(spec, mspec, **kwargs)
        return val

    def fitness_method(self, method_name):
        r"""Set the method to evaluate the fitness."""
        if method_name == 'chi2-error':
            method = self.chi2_error
        if method_name == 'weighted-chi2-error':
            method = self.weighted_chi2
        if method_name == 'chi2':
            method = self.chisquare
        return method

    @staticmethod
    def setmixturetype(mixtype):
        r"""Set the class for making the mixtures."""
        if mixtype == 'intimate':
            mix = IntimateMixture
        return mix

    def chi2_error(self, spec, mspec):
        r"""Calculate the Reduced chi-square.

        In this case the spec have errorbars."""
        res_aux = (spec.r - mspec.r)**2
        res = np.sum(res_aux/spec.r_unc**2) / (2*len(self.samples))
        # print(res_aux)
        # calculating p-value
        pte = 1 - chi2.cdf(res, len(spec.w)-2*len(self.samples))
        return res, pte

    def chisquare(self, spec, mspec):
        r"""Calculate the Reduced chi-square."""
        res_aux = (spec.r - mspec.r)**2
        res = np.sum(res_aux) / (2*len(self.samples))
        # calculating p-value
        pte = 1 - chi2.cdf(res, len(spec.w)-2*len(self.samples))
        return res, pte

    def weighted_chi2(self, spec, mspec, wranges, weigths):
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
        for i, w in enumerate(wranges):
            aux = np.ravel(np.argwhere((mspec.w >= w[0]) & (mspec.w <= w[1])))
            weight_arr[aux] = np.array([weigths[i]/weight_sum for _ in range(len(aux))])
        weight_arr = (weight_arr / np.sum(weight_arr))*len(mspec.w)
        res_aux = ((spec.r - mspec.r)**2) * weight_arr
        res = np.sum(res_aux/(spec.r_unc**2)) / (2*len(self.samples))
        pte = chi2.sf(res, len(spec.w)-2*len(self.samples)-1)
        return res, pte

    @staticmethod
    def residual(spec, mspec):
        r"""Calculate the residual."""
        res_aux = np.abs(spec.r - mspec.r)
        res = np.sum(res_aux)
        return res


class ModelOutput(object):
    r"""Class to visualize the output of the compositional model."""

    def __init__(self, spec, albedo_lim, model,
                 normalize_spec, wnorm):
        r"""Initialize output for models.

        Parameters
        ----------
        spec: cana.Spectrum
            The modeled Spectrum.

        samples:
            The samples used to model the spectrum.

        model:
            The reflectance model that was used to model.

        """
        self.spec = spec
        self.albedo_lim = albedo_lim
        self.wnorm = wnorm
        self.normalize_spec = normalize_spec
        self.model = model
        # building DataFrame
        self.allmixes = self.build_allmixes_dataframe(model.samples_ids)
        # placeholders
        self.bestmix = -1
        self.bestalbedo = -1
        self.bestscore = -1
        self.bestspec = -1
        self.valid_mixes = -1
        self.result = -1

    @staticmethod
    def build_allmixes_dataframe(samples, albedo=True):
        r"""Build empty dataframe to vizualize evaluated mixtures."""
        columns = samples.copy()
        grain_aux = [sam+'_grain' for sam in samples]
        columns.extend(grain_aux)
        columns.append('porosity')
        if albedo:
            columns.append('albedo')
        columns.extend(['score', 'p-value', 'iter'])
        out = pd.DataFrame(columns=columns)
        return out

    def build_result(self):
        r"""
        """
        result = self.valid_mixes.describe()
        result.loc['Bestfit'] = self.allmixes.iloc[0]
        result = result.drop(['count', '25%', '50%', '75%'])
        result = result.drop(['score', 'iter', 'p-value'], axis=1)
        self.result = result
        return result

    def __repr__(self):
        if not isinstance(self.result, pd.DataFrame):
            repr_aux = 'Model have not yet produced results'
        else:
            repr_aux = self.result.T.__repr__()
        return repr_aux

    def make_spec(self, wavelengths=None):
        r"""Make the spectrum of the bestfit model."""
        if wavelengths is None:
            wavelengths = self.spec.w
        spec, albedo = self.bestmix.make(wavelengths=wavelengths)
        if self.normalize_spec:
            spec = spec.normalize(self.wnorm)
        return spec, albedo

    def set_bestfit(self, index=None):
        r"""Set which is the bestfit."""
        self.sort_mixes()
        if index is None:
            best = self.allmixes.iloc[0]
            index = best.name
        self.bestmix = self.mixfromindex(index=index)
        self.bestspec, self.bestalbedo = self.make_spec()
        self.bestscore = self.allmixes.loc[index, 'score']
        return self

    def mixfromindex(self, index=0):
        r"""Make the mixture from the results table."""
        aux = self.allmixes.loc[index]
        props = list(aux[[sam for sam in self.model.samples_ids]])
        grains = list(aux[[sam+'_grain' for sam in self.model.samples_ids]])
        samples = self.model.samples
        mix = self.model.mixtype(samples, proportions=props,
                                 grainsizes=grains,
                                 porosity=self.model.porosity,
                                 model=self.model.refmodel)
        return mix

    def sort_mixes(self, col='score'):
        r"""Sort the rows of results table.

        Parameters
        ----------
        col: string
            The name of the column used to sort the table.
            Default is 'score'.

        """
        self.allmixes = self.allmixes.sort_values(by=[col])
        if self.albedo_lim is not None:
            aux = self.allmixes[(self.allmixes['albedo'] >= self.albedo_lim[0]) &
                                (self.allmixes['albedo'] <= self.albedo_lim[1])]
            aux2 = list(aux.index)
            aux3 = [i for i, v in self.allmixes.iterrows() if i not in aux2]
            aux2.extend(aux3)
            self.allmixes = self.allmixes.loc[aux2]
        return self.allmixes

    def select_valid(self, sigma=3):
        r"""Select valid models.

        Parameters
        ----------
        sigma: float
            The sigma value used to select valid models.

        Returns
        -------
        pandas.DataFrame
            Valid models
        """
        if self.albedo_lim is not None:
            valid = self.allmixes[(self.allmixes['albedo'] >= self.albedo_lim[0]) &
                                  (self.allmixes['albedo'] <= self.albedo_lim[1])]
        conf_value = 1 - 2*norm.sf(sigma)
        self.valid_mixes = valid[valid['p-value'] > conf_value]
        return self.valid_mixes

    def append_mix(self, val, index=None):
        r"""Append the mixture evaluation result to the table.

        Parameters
        ----------
        val: list
            List of mixture parameter.

        index: none or integer
            The index that will be used in the table.
            If None, will autoincrement index.
            Default is None.

        """
        if index is not None:
            self.allmixes.loc[index] = val
        else:
            aux = len(self.allmixes)
            self.allmixes.loc[aux] = val
        return self.allmixes

    def plot_bestfit_residual(self, fax=None, show=False, savefig=None,
                              axistitles=True, speckwargs=None, mixwargs=None, reskwargs=None,
                              legendkwargs=None):
        r"""Vizualize the bestfit."""
        # checking if plot in another frame
        if fax is None:
            fig = plt.figure(1)
            plt.tight_layout()
            fax1 = fig.add_axes((.1, .3, .8, .6))
        # Ploting the spec
        self.spec.plot(fax=fax1, show=False, speckwargs={'c': 'k', 'lw': 2,
                                                         'label': 'spectrum'})
        # Ploting mix
        self.bestspec.plot(fax=fax1, show=False,
                           speckwargs={'c': 'firebrick', 'lw': 2,
                                       'label': 'bestfit'})
        fax1.legend()
        # Ploting residuals
        fax2 = fig.add_axes((.1, .1, .8, .2))
        fax2.plot(self.spec.w, self.spec.r - self.bestspec.r)
        fax2.set_ylabel('Residuals')
        # check if save the image
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()


    def plot_valid_proportions(self, savefig=None, show=True):
        r"""
        """
        aux = self.valid_mixes[self.model.samples_ids]
        aux.boxplot(rot=45)
        # plt.savefig('varuna-spitizer-values.png')
        plt.ylim(self.model.proportions[0], self.model.proportions[1])
        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            # matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()

    def plot_valid_models(self, savefig=None, show=True):
        r"""
        """
        f, ax = plt.subplots()

        for i, v in self.valid_mixes.iterrows():
            b = self.mixfromindex(index=i)
            sp, alb = self.make_spec()
            sp.normalize(0.55)
            ax.plot(sp.w, sp.r, lw=0.5, c='0.7')
        ax.errorbar(self.spec.w, self.spec.r, yerr=self.spec.r_unc, zorder=1000000)

        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            # matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()

    def plot_proportions(self, shape=[3,2],figsize=(10,15), step=1,
                         fax=None, sharex=True, savefig=None, show=True):
        r"""Vizualize the bestfit."""
        # checking if plot in another frame
        if fax is None:
            _, ax = plt.subplots(*shape, figsize=figsize, sharex=sharex)
            plt.tight_layout()

        ax = np.ravel(ax)
        for n, i in enumerate(self.model.samples_ids):
            ax[n].scatter(self.allmixes[i][::step],
                          np.log10(self.allmixes['score'][::step]),
                          c='0.7', s=1, label=i)
            ax[n].scatter(self.valid_mixes[i], np.log10(self.valid_mixes['score']),
                          c='b', s=1)
            ax[n].set_ylabel('score')
            ax[n].set_xlabel('proportion')
            ax[n].legend()
            # ax[n].axhline(np.log10(aux['score'].max()))

        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            # matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()

    def plot_grains(self, shape=[3,2], figsize=(10,15), step=1,
                    sharex=True, savefig=None, show=True):
        r"""Vizualize the bestfit."""
        # checking if plot in another frame
        _, ax = plt.subplots(*shape, figsize=figsize, sharex=sharex)
        plt.tight_layout()

        ax = np.ravel(ax)
        for n, i in enumerate(self.model.samples_ids):
            ax[n].scatter(self.allmixes[i + '_grain'][::step],
                          np.log10(self.allmixes['score'][::step]),
                          c='0.7', s=1, label=i)
            ax[n].scatter(self.valid_mixes[i + '_grain'],
                          np.log10(self.valid_mixes['score']),
                          c='b', s=1)
            ax[n].set_ylabel('score')
            ax[n].set_xlabel('grainsize')
            ax[n].legend()
            # ax[n].axhline(np.log10(aux['score'].max()))

        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            # matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()
