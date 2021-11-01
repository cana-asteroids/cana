r"""Methods for convolution between a spectrum and a photometric system."""

import matplotlib.pyplot as plt
from .spec import loadspec, stack_spec
from .util import kwargupdate
from .photo import Photometry, PhotoDataframe


def interp_spec(spec, baseaxis, datatype=None):
    r"""Interpolate the Spectral Data to a new wavelentgh axis.

    Parameters
    ----------
    basename: array
        New wavelength array

    Returns
    -------
    The interpolated object

    """
    if baseaxis is None:
        new_spec = spec
    if datatype is None:
        new_spec = spec.rebase(baseaxis)
    else:
        aux = []
        for key, val in datatype.items():
            if key == 'spec':
                baxis = baseaxis[(baseaxis >= val[0]) & (baseaxis <= val[1])]
                spc = spec.rebase(baxis)

            else:
                aux2 = convolution(spec, system=key, filters=val)
                spc = aux2.tospec(spec.label)
            aux.append(spc)
        new_spec = stack_spec(aux)
    return new_spec


def convolution(spec, system='sdss', filters=None, label=None, speckwargs=None):
    r"""
    Perform the convolution between an Spectrum and a the filter system.

    Parameters
    ----------
    spec: Spectrum object

    system: str
        The photometric system name. Options are: sdss, jplus, osiris

    filters:
        The filters to apply the convolution. None will try the convolution
        in all filters

    label: string
        The label for the output. Used only for a single file

    speckwargs: dict
        Options for reading files with loadspec. Default is {'unit':'microns'}.
        For more options check cana.loadspec documentation.

    Returns
    -------
    PhotoDataframe with colvoluted spectrophotometry
    """
    if isinstance(system, str):
        photo = Photometry(system)
    else:
        photo = system
    speckwars_default = {'unit': 'micron'}
    speckwargs = kwargupdate(speckwars_default, speckwargs)
    if not isinstance(spec, list):
        if isinstance(spec, str):
            spec = loadspec(spec, label=label, **speckwargs)
        pspec = photo.convol(spec, filters, spec.label)
    # appling on a list of spectra
    elif isinstance(spec, list):
        pspec = PhotoDataframe(system=system)
        for fsp in spec:
            spc = loadspec(fsp, **speckwargs)
            ps_aux = photo.convol(spc, filters)
            pspec = pspec.append(ps_aux)
    return pspec



def plot_convolution(spec, photo, fax=None, savefig=None, axistitles=True,
                     show=True, speckwargs=None, dotkwargs=None):
    r'''
    '''
    # setting default values for image plot with matplotlib
    specsty_defaults = {'c':'0.3', 'lw':'1', 'zorder':0}
    specsty = kwargupdate(specsty_defaults, speckwargs)
    # setting default values for image plot with matplotlib
    dotsty_defaults = {'color':'steelblue', 's':70, 'zorder':10}
    dotkwargs = kwargupdate(dotsty_defaults, dotkwargs)
    # checking if plot in another frame
    if fax is None:
        fig = plt.figure()
        fax = fig.gca()
    # Ploting the spec
    spec.plot(fax=fax, axistitles=axistitles, show=False, speckwargs=speckwargs)
    # photo = photo.dropna(axis='columns')
    fils = [i for i in photo.columns if '_err' not in i]
    for i in fils:
        fax.scatter(photo.system['center'][i], photo[i], **dotkwargs)
    # check if save the image
    if savefig is not None:
        plt.savefig(savefig)
        if not show:
            plt.clf()
    # show in the matplotlib window?
    if show:
        plt.show()
