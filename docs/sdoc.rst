SDOC: Small Database of Optical Constants
==========================================

SDOC is a database of selected optical constants to be used in compositional modeling of solar system minor bodies,
maintained separately from cana. Note that SDOC is not required to run cana.


All data of the optical constants are stored in a single HDF file, for easy sharing and access of the data.

More information about SDOC and its code is accessible at `github.com/depra/sdoc <https://github.com/depra/sdoc>`_.

**This database is currently under construction, the contents of the database are incomplete and used for testing purposes only.
A version of database with carefully selected optical constants and their metadata will be available soon.**


HDF Database Schema
--------------------

.. image:: https://raw.githubusercontent.com/depra/sdoc/master/docs/images/sdoc.png?token=ABGQYMX3AAZIGIHN7JL7FL27GP7W2
   :align: center


Dependencies
------------

- h5py
- pandas
- numpy


How to install
--------------

If you have `Anaconda <https://www.anaconda.com/distribution/>`_ or `pip <https://pypi.org/project/pip/>`_ installed:

::

   pip install cana-sdoc

Usage
-----
For a example of how to access the database and search for a optical constant, take a look at `sdoc notebooks <https://github.com/depra/sdoc/blob/master/notebooks/accesing_the_database.ipynb>`_.,
or at the `cookbook recipes related to compositional modeling <gallery/index.html#compositional-models>`_.

Cite
----

A paper with the description of the database and the compositional modeling techniques is in preparation. However, if you use the optical constants inserted in the db, **make sure to cite the references from where the optical constants were extracted.**
