.. CANA documentation master file, created by
   sphinx-quickstart on Mon Nov 18 23:11:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://raw.githubusercontent.com/depra/cana/readme/docs/_static/Cana_logo_small.png
   :align: center
   :scale: 50
   :alt: Codes for ANalisis of Asteroids
    


What is CANA?
-------------

The CANA package stands for "*Codes for ANalysis of Asteroids*". The tool was designed to perform scientific analysis for asteroids spectroscopic data.

As part of the PRIMitive Spectroscopic Survey (PRIMASS), we have collected and anaylized hundreds of spectra of primitive asteroids among the last years.
In this context we collected the routines that we typically use for the analysis and parametrization of these data to make a package.
As we are preparing to make the data library (PRIMASS-L) public to the comunity through 
`Small Bodie Node of the the Planetary Data Science (PDS) <https://pds-smallbodies.astro.umd.edu/>`_, it is our aim to make our science reproductible and of easy access. 


Contents
--------
The package is currently under development. At the moment CANA counts with:

Spectroscopic analysis tools:

* Handling tools, such as load, trim, fit (and autofit), estimate SNR, clean spectrum
* Slope calculation
* Taxonomic Classification
* Hydration band analysis

We will soon make available new parametrization methods and tools for handling Photometric data and Compositional modeling!

Dependencies
------------
As Python 2.7 is coming to an end, we have updated all code for Python >=3.6, and do not provide support for Python 2.

We recommend using the `Anaconda Python distribution <https://www.anaconda.com/distribution/>`_ to ensure you have all dependencies installed.

If you have any problems with the instalation, check if all dependencies are installed.
The dependencies can be seen below:

- `Numpy <http://www.numpy.org/>`__
- `Pandas <https://pandas.pydata.org/>`_
- `Scipy <https://www.scipy.org/>`_
- `Matplotlib <https://matplotlib.org/>`_
- `Sklearn <http://scikit-learn.org/stable/>`_

Installing
-----------

You can get the latest stable version
`here <https://github.com/depra/cana/releases>`_, or clone the repository using git:

::

      git clone https://github.com/depra/cana.git

Unpack it if necessary, and go into the directory "*cana-master*", then run the below commands on a terminal shell:

:: 

   python setup.py install


Get Started
-----------

See the documentation page to get started!


Cite
----
CANA is a open tool. You can use, adapt and modify it as wanted. But if you any of these, please cite us!

   De Pra, M., Carvano, J., Morate, D., Licandro, J. Pinilla-Alonso, N. (2018). CANA: An open-source python tool to study hydration in the Solar System. 

See `cite <cite.html>`_ to get the bibtex entry for the citation.


