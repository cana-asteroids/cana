.. Cana documentation master file, created by
   sphinx-quickstart on Thu May 17 20:15:55 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: CANA: Codes for Analysis of Asteroids

.. raw:: html

    <div class="container-fluid banner">
        <div class="container">

            <div class="row site-title">
                <div class="col-lg-3 col-sm-2">
                </div>
                <div class="col-lg-6 col-sm-8">
                    <img src="_static/CanaColor1.png" width="80%">
                </div>
                <div class="col-lg-3 col-sm-2">
                </div>
            </div>
        </div>
    </div>
    <div class="row slogan">
            <p class="package-name">Codes for Analysis of Asteroids
            </p>
            <p class="slogan">An open-source Python library for handling asteroids 
            spectroscopic and spectrophotometric data	
            </p>
    </div>


.. note::

    The Cana package is a work in progress. We need **your help** to make it
    better! You can contribute by providing an feedback of your experience and sharing
    code of analysis that can be incorporated in the package.
    Write to the `Mário De Prá  <mariondepra@gmail.com>`__
    or send us a
    `pull request on Github <https://github.com/depra/cana/pulls>`__.

Introduction
------------

Cana is an open-source python package designed for handling the analysis of asteroids spectroscopic and spectrophotometrc data. We aim to make our science more **open** and **reproductible**. 

Don't reinvent the wheel! 
We just released an beta version! The code is currently under development, there are still many
modules we want to include.



Contents
---------

.. raw:: html

    <div class="home-row">

        <div class="row">
            <div class="col-md-6 home-overview">
                <h3><a href="api.html">Spectroscopy</a></h3>
                <p>
                Handling tools, Taxonomic Classification, Slope, Absorptium bands analysis.
                </p>
                <em>
                 The package is mainly designed for the indentification
                of features commonly present in the spectra of primitive asteroids,
                such as the hydration band near 0.7 microns and the turn-off point
                near the 0.5 micron region.
                </em>
            </div>

            <div class="col-md-6 home-overview">
                <h3><a href="api.html">Spectrophotometry</a></h3>
                <p>
                Spectral convolution, Slope and Absorptium bands probability.
                </p>
                <em>
                The Spectral convolution to ECAS, SDSS and J-PLUS. 
                The package is mainly designed for the indentification
                of features commonly present in the spectra of primitive asteroids,
                such as the hydration band near 0.7 microns and the turn-off point
                near the 0.5 micron region.
                </em>
            </div>
        </div>

        <div class="row">
            <div class="col-md-6 home-overview">
                <h3><a href="api.html">Compositional Analysis</a></h3>
                <p>
                To be implemented
                </p>
                <em>
                We will soon implement methods for performing compositional analysis
                for different types of Solar System small body objects.
                </em>
        </div>
    </div>


Get Started
------------
See the `install instructions <install.html>`_ to set up your computer and install Cana.

Once you have everything installed, we recommend you take a look at our recipes at the
`Cookbook <cookbook.html>`_, to get an idea of what cana can do. For more detailed 
documentation go to the `Documentation page <api.html>`_ to get the complete API description. 

Cana is still a work in progress, if you are have trouble getting it to work, 
or you want to give us a feedback or colaborate, you can:

* Write to `Mário De Prá  <mariondepra@gmail.com>`_
* Report bugs, enhancements or else through `Github <https://github.com/depra/cana/issues>`_.
 

Cite
----
CANA is **made by scientists**.  Citing us help us justify the effort that
goes into building, maintaining and making it an open tool.

   De Pra, M., Carvano, J., Morate, D., Licandro, J. Pinilla-Alonso, N. (2018). CANA: An open-source python tool to study hydration in the Solar System. 

See `cite <cite.html>`_ to get the bibtex entry for the citation.

|

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   docs.rst
   install.rst
   cookbook.rst
   contribute.rst
   cite.rst
