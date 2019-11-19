r"""Tests for Spectrum core module."""

import numpy as np
import sys
sys.path.insert(0, "..")
import cana


def testloadspec():
    r"""Test if able to load Spectrum."""
    spec = cana.datasets.getspectrum('000334', ref='primass')
    assert len(spec.w) == 1788
    assert np.mean(spec.r) > 0.


def testunitconversion():
    r"""Test if unit transformation is acting properly."""
    spec = cana.datasets.getspectrum('000334', ref='primass')
    assert spec.unit == 'micron'
    assert np.mean(spec.w) < 10
    spec = spec.micron2angstrom()
    assert spec.unit == 'angstrom'
    assert np.mean(spec.w) > 10
    spec = spec.angstrom2micron()
    assert spec.unit == 'micron'
    assert np.mean(spec.w) < 10


def testspectrim():
    r"""Test if is trimming spectrum properly."""
    spec = cana.datasets.getspectrum('000334', ref='primass')
    assert spec.w.min() < 0.5
    assert spec.w.max() > 0.8
    spec = spec.trim(0.55, 0.75)
    assert spec.w.min() > 0.5
    assert spec.w.max() < 0.8


def testrebin():
    r"""Test if is rebing spectrum properly."""
    spec = cana.datasets.getspectrum('000334', ref='primass')
    assert len(spec.w) == 1788
    spec = spec.rebin(binsize=6)
    assert len(spec.w) == 298


def testautofit():
    r"""Test if is spectrum automatic fitting is working properly."""
    spec = cana.datasets.getspectrum('000334', ref='primass')
    aux = spec.autofit(degree_min=1, degree_max=12)
    order = len(aux[1])
    assert order == 5
