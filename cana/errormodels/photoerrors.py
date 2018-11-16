import numpy as np
# import numba as nb

class PhotoError(object):
    
    def __init__(self, n=1000):
        r'''
        '''
        self.n = n
        # sppeding up some methods
        self.resample_vec = np.vectorize(self._ref_resample_aux)
        # self.distribution = distribution

    @staticmethod
    def _ref_resample_aux(point, rms):
        aux = np.random.normal(point, rms, 1)
        return aux[0]
    
    def distribution(self, wave, ref, filters, func, **funckwargs):
        r'''
        '''
        out = []
        filters_err = [fil+'_err' for fil in  filters]
        ref_aux = ref.copy()
        for _ in xrange(self.n):
            ref_aux[filters] = self.resample_vec(ref[filters], ref[filters_err])
            out_aux = func(ref_aux, wave, filters, **funckwargs)
            out.append(out_aux)
        return np.array(out)

    # def distribution2(self, wave, ref, filters, func, **funckwargs):
    #     r'''
    #     '''
    #     out = []
    #     filters_err = [fil+'_err' for fil in  filters]
    #     # ref_aux = ref.copy()
    #     for _ in xrange(self.n):
    #         out.append(self.dist(ref, wave, filters, filters_err, func, **funckwargs))
    #     return np.array(out)

    # def dist(self, ref, wave, filters, filters_err, func):
    #     ref_aux = ref.copy()
    #     # print ref
    #     ref_aux[filters] = self.resample_vec(ref[filters], ref[filters_err])
    #     # print ref
    #     func = nb.njit(func)
    #     out_aux = func(ref_aux, wave, filters)
    #     return out_aux

# @np.vectorize
# def resample_vec(point, rms):
#     aux = np.random.normal(point, rms, 1)
#     return aux[0]

# @nb.jit
# def distribution(n, wave, ref, filters, func):
#         r'''
#         '''
#         out = []
#         filters_err = [fil+'_err' for fil in  filters]
#         ref_aux = ref.copy()
#         for m in range(n):
#             # print ref
#             ref_aux[filters] = resample_vec(ref[filters], ref[filters_err])
#             # print ref
#             out_aux = func(ref_aux, wave, filters)
#             out.append(out_aux)
#         return np.array(out)