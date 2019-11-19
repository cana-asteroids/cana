r"""Tests for Spectrum tools module."""

import sys
sys.path.insert(0, "..")
import cana


def testslope():
        r"""Test computing the slope."""
        spec = cana.datasets.getspectrum('000334', ref='primass')
        slp = cana.slope(spec)
        assert 4.50 <= slp.slope <= 5.5


def testtaxonomy():
        r"""Test taxonomic classification."""
        spec = cana.datasets.getspectrum('000334', ref='primass')
        tax = cana.taxonomy(spec, system='bus')
        assert tax['tax'][0] == 'Xc'


def testband():
    r"""Test band depth."""
#         spec = cana.datasets.getspectrum('000334', ref='primass')
