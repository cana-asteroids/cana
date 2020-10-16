r"""Random search to find best fit model."""

import itertools
import numpy as np
from .core import BaseCompositionalModel, ModelOutput
from ..util import verboseprint, get_truncated_normal
from .. import interp_spec
# import time

class RandomSearch(BaseCompositionalModel):
    r"""Grid Model for compositional analysis."""

    def __init__(self, samples, grainsizes=[5, 200],
                 proportions=[0, 1], porosity=0.5,
                 population=1000, alpha_scale=20,
                 niterations=5, nwalkers=4,
                 verbose=2, initpopulation=None,
                 mixtype='intimate',
                 model='shkuratov', metric='weighted-chisquare'):
        r"""Initialize GridModel.

        Parameters
        ----------
        samples: list
            List of optical constants

        grainsizes: list or list of lists
            List of limits for the grainsizes

        proportions: list
            List of limits for the proportions

        porosity: float of list
            Value of the porosity or list of limits for the porosity

        population: int
            Number of mixtures that will be tested in each iteration

        niterations: int
            Number of iterations

        nwalkers: int
            Number of walkers that will be used in each iteration.

        verbose: 0,1 or 2
            Verbose level: 0 to avoid any print while running modelling
                           1 to print minimal info
                           2 to print higher level info

        mixtype: string
            Type of mixture that will be used to model the spectrum.

        model: string
            The name of the reflectance model to be used

        method: string
            The name of the metric that will be used to evaluate the fitness
            of the model. Default is reduced chisquare

        """
        BaseCompositionalModel.__init__(self, samples, grainsizes,
                                        proportions, porosity,
                                        mixtype=mixtype, model=model,
                                        metric=metric)
        self.n = niterations
        self.iter = -1
        self.population = population
        self.nwalkers = nwalkers
        self.verbose = verbose
        self.alpha_scale = alpha_scale
        if initpopulation is None:
            self.initpopulation = population
        else:
            self.initpopulation = initpopulation

    def gen_initial_population(self):
        r"""Generate initial population."""
        # Creating parameters for initial population
        prop = np.random.dirichlet(np.ones(len(self.samples)) /
                                   len(self.samples),
                                   size=self.initpopulation)
        if self.grains.shape == (2,):
            grains = np.random.randint(self.grains[0], self.grains[1]+1,
                                       size=(self.initpopulation,
                                       len(self.samples)))
        else:
            grains = np.zeros((self.initpopulation, len(self.samples)))
            for n, g in enumerate(self.grains):
                grains[:, n] = np.random.randint(g[0], g[1]+1,
                                                 size=self.initpopulation)
        # creating mixtures
        for i in range(self.initpopulation):
            yield grains[i], prop[i]

    def gen_population_iter(self, base):
        r"""Generate population for the iteration."""
        prop_base = np.array(base.proportions)
        population = int(self.population/self.nwalkers)
        aux = np.argwhere(prop_base < 0.03)
        if len(aux) > 0:
            prop_base[aux] = 0.03
        prop_base = prop_base * self.alpha_scale
        # print(prop_base)
        prop = np.random.dirichlet(prop_base, size=population)
        grains = np.zeros((population, len(self.samples)))
        for n, g in enumerate(self.grains):
            gaux = get_truncated_normal(base.grainsizes[n], np.std([g[0], g[1]])/ 10, #2**(self.iter-1),
                                        low=g[0], upp=g[1]+1)
            grains[:, n] = gaux.rvs(population)
        for i in range(population):
            yield grains[i], prop[i]

    def mixturegrid(self, base=None):
        r"""Make mixture grid."""
        # makind mixture grid
        if self.iter == -1:
            grid = self.gen_initial_population()
        else:
            # running iterations and walkers
            aux = base.allmixes.copy()
            aux2 = aux[self.samples_ids]
            aux2 = aux2.drop_duplicates()
            aux2 = aux2.iloc[:self.nwalkers]
            grid = []
            if self.verbose > 1:
                verboseprint("Number of Walkers {0}".format(self.nwalkers))
            for _, i in enumerate(aux2.index):
                # base = base.set_bestfit(index=i)
                mix = base.mixfromindex(index=i)
                grid_aux = self.gen_population_iter(mix)
                grid.append(grid_aux)
            grid = itertools.chain(*grid)
        for mix in grid:
            grains, prop = mix
            m_aux = self.mixtype(self.samples, proportions=prop,
                                 grainsizes=grains, porosity=self.porosity,
                                 model=self.refmodel)
            yield m_aux

    def model(self, spec, albedo_lim=None, albedo_w=0.55,
              normalize_spec=True, wnorm=1.2, normvalue=0.127, wavelengths=None,
              sigma=1, datatype=None, filterfunc=None, **kwargs):
        r"""Find the best fit to the spectrum.

        Parameters
        ----------
        baseaxis: None or array
            The wavelength array that will be used if rebase=True.
            If baseaxis=None it will select the wavelength
            of the first sample as a reference. Not used if rebase=False.

        """
        if self.verbose >= 1:
            verboseprint("Running grid search model for spectrum {0}".format(spec.label),
                         underline=True)
        if self.verbose == 2:
            verboseprint("""\
                         Using optical constans: {0}
                         with grainsize ranges: {1}
                         and proportion ranges: {2}
                         """.format(','.join(self.samples_ids),
                                    ','.join(map(str, self.grains)),
                                    ','.join(map(str, self.proportions))))

        result = ModelOutput(spec, albedo_lim,
                             normalize_spec=normalize_spec, wnorm=wnorm,
                             model=self)
        if wavelengths is None:
            wavelengths = spec.w
        # timeinit = time.time()
        self.iter = -1
        for i in range(1, self.n+1):
            # time1 = time.time()
            # printing info
            if self.verbose == 2:
                verboseprint('Running Iteration: {0}'.format(i), underline=True)
            if self.verbose == 1:
                verboseprint('Running Iteration: {0} ...'.format(i))
            # generating mixtures
            if self.iter == -1:
                mixtures = self.mixturegrid()
            else:
                mixtures = self.mixturegrid(result)
            # running model
            for _, mix in enumerate(mixtures):
                if albedo_lim is not None:
                    mspec, malbedo = mix.make(albedo_w=albedo_w,
                                              wavelengths=None)
                    mspec = interp_spec(mspec, wavelengths, datatype)
                # normalizing
                if normalize_spec:
                    mspec = mspec.normalize(wnorm)
                    # if normvalue is not None:
                    #     mspec.r = mspec.r*normvalue
                # evaluate mixture
                score, p_value = self.eval(spec, mspec, **kwargs)
                # making output
                out = list(np.round(mix.proportions, 4))
                out.extend(mix.grainsizes)
                out.append(self.porosity)
                if albedo_lim is not None:
                    out.append(malbedo)
                out.extend([score, p_value, i])
                result.append_mix(out, index=None)
            # preparing for next iteration
            result.sort_mixes()
            if filterfunc is not None:
                result = filterfunc(result, spec)
            self.iter = i
            # timeend =  time.time()
            # print("--- %s seconds ---" % (timeend - timeinit))
            # print("--- %s seconds ---" % (timeend - time1))
        result = result.set_bestfit()
        result.select_valid(sigma)
        result.build_result()

        if self.verbose >= 1:
            verboseprint("""\
                         Best mixture with score {0}:
                         {1}
                         """.format(result.bestscore,
                                    result.bestmix.__str__()))

        return result
