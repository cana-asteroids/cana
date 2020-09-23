import numpy as np


class SpectralData:
    r"""
    """

    def __init__(self):
        self.repr = self.make_repr()

    def __getitem__(self, item):
        if isinstance(item, (int, slice)):
            return self.__class__(*[self[d][item] for d in self.dtype],
                                  label = self.label)
        elif isinstance(item, (list, np.ndarray)):
            return self.__class__(*[self[d][item] for d in self.dtype],
                                  label = self.label)
        else:
            return self.__getattribute__(item)

    def make_repr(self,maxlines=20, nlines=4):
        r"""Make representation of the Spectral Data"""
        basestr = self.__class__.__name__ + '(['
        whitespace = ' ' * len(basestr)
        if len(self)> maxlines:
            idss = list(range(nlines))
            idss.extend([-i-1 for i in idss[::-1]])
        else:
            idss = range(len(self))
        for n,i in enumerate(idss):
            line = ', '.join(map(str, [self[idx][i] for idx in self.dtype]))
            if i == 0:
                basestr += '('+line+'),\n'
            elif n ==len(idss)-1:
                basestr += whitespace+'('+line+')'
            else:
                if (n == nlines) & (len(self)> maxlines):
                    basestr += whitespace + ' ...\n'
                basestr += whitespace + '(' + line + '),\n'
        basestr +='], \n'
        basestr += whitespace + 'columns=(' + ', '.join(self.dtype) + '),\n'
        basestr += whitespace + 'label="{0}"'.format(self.label)
        basestr += ')'
        return basestr

    def __repr__(self):
        return self.repr

    @property
    def shape(self):
        return (len(self.w), len(self.dtype))

    def __len__(self):
        return len(self.w)


    def sort(self, order='w'):
        r"""
        """
        aux = np.argsort(self[order])
        return self.__class__(*[self[d][aux] for d in self.dtype],
                              label = self.label)
