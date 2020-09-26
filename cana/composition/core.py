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
                 model='shkuratov', metric='chisquare'):
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
        if method_name == 'chisquare':
            method = self.chi2
        if method_name == 'weighted-chisquare':
            method = self.weighted_chi2
        return method

    @staticmethod
    def setmixturetype(mixtype):
        r"""Set the class for making the mixtures."""
        if mixtype == 'intimate':
            mix = IntimateMixture
        return mix

    def chi2(self, spec, mspec):
        r"""Calculate the Reduced chi-square."""
        res_aux = (spec.r - mspec.r)**2
        res = np.sum(res_aux/spec.r_unc**2) / (2*len(self.samples))
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
        res = np.sum(res_aux/spec.r_unc**2) / (2*len(self.samples))
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
        self.table = self.build_dataframe(model.samples_ids)
        # placeholders
        self.bestmix = -1
        self.bestalbedo = -1
        self.bestscore = -1
        self.bestspec = -1
        self.valid = -1

    @staticmethod
    def build_dataframe(samples, albedo=True):
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

    def make_spec(self, baseaxis=None):
        r"""Make the spectrum of the bestfit model."""
        if baseaxis == None:
            baseaxis = self.spec.w
        spec, albedo = self.bestmix.make(baseaxis=baseaxis)
        if self.normalize_spec:
            spec = spec.normalize(self.wnorm)
        return spec, albedo

    def set_bestfit(self, index=None):
        r"""Set which is the bestfit."""
        if not self.sorted:
            self.sort()
        if index is None:
            best = self.table.iloc[0]
            index = best.name
        self.bestmix = self.mixfromindex(index=index)
        self.bestspec, self.bestalbedo = self.make_spec()
        self.bestscore = self.table.loc[index, 'score']
        return self

    def mixfromindex(self, index=0):
        r"""Make the mixture from the results table."""
        aux = self.table.loc[index]
        props = list(aux[[sam for sam in self.model.samples_ids]])
        grains = list(aux[[sam+'_grain' for sam in self.model.samples_ids]])
        samples = self.model.samples
        mix = self.model.mixtype(samples, proportions=props,
                                 grainsizes=grains,
                                 porosity=self.model.porosity,
                                 model=self.model.refmodel)
        return mix

    def sort(self, col='score'):
        r"""Sort the rows of results table.

        Parameters
        ----------
        col: string
            The name of the column used to sort the table.
            Default is 'score'.

        """
        self.table = self.table.sort_values(by=[col])
        if self.albedo_lim is not None:
            aux = self.table[(self.table['albedo'] >= self.albedo_lim[0]) &
                             (self.table['albedo'] <= self.albedo_lim[1])]
            aux2 = list(aux.index)
            aux3 = [i for i, v in self.table.iterrows() if i not in aux2]
            aux2.extend(aux3)
            self.table = self.table.loc[aux2]
        return self.table

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
            valid = self.table[(self.table['albedo'] >= self.albedo_lim[0]) &
                               (self.table['albedo'] <= self.albedo_lim[1])]
        conf_value = 1-norm.pdf(sigma)
        self.valid = valid[valid['p-value'] > conf_value]
        return self.valid

    def __str__(self):
        r"""
        """
        pass

    def append(self, val, index=None):
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
            self.table.loc[index] = val
        else:
            aux = len(self.table)
            self.table.loc[aux] = val
        return self.table

    def __getattr__(self, name):
        r"""
        """

    def plotbestfit(self, residual=True, fax=None, show=False, savefig=None,
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
        if residual:
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

    def plot_proportions(self, shape=[3,2],figsize=(10,15), step=1,
                         fax=None, savefig=None, show=True):
        r"""Vizualize the bestfit."""
        # checking if plot in another frame
        if fax is None:
            _, ax = plt.subplots(*shape, figsize=figsize)
            plt.tight_layout()

        ax = np.ravel(ax)
        for n, i in enumerate(self.model.samples_ids):
            ax[n].scatter(self.table[i][::step],
                          np.log10(self.table['score'][i][::step]),
                          c='0.7', s=1, label=i)
            ax[n].scatter(self.valid[i], np.log10(self.valid['score']),
                          c='b', s=1)
            ax[n].set_ylabel('score')
            ax[n].set_xlabel('proportion')
            # ax[n].axhline(np.log10(aux['score'].max()))

        if savefig is not None:
            plt.savefig(savefig)
            if not show:
                plt.clf()
            matplotlib.use('TkAgg')
            # show in the matplotlib window?
        if show:
            plt.show()
