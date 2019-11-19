.. CANA documentation master file, created by
   sphinx-quickstart on Mon Nov 18 23:11:00 2019.
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
                    <img src="https://raw.githubusercontent.com/depra/cana/spec/docs/_static/Cana_logo.png" width="80%">
                </div>
                <div class="col-lg-3 col-sm-2">
                </div>
            </div>
        </div>
    </div>
    <div class="row slogan">
            <p class="package-name">Codes for ANalysis of Asteroids
            </p>
            <p class="slogan">An open-source python tool for asteroid science	
            </p>
    </div>

What is CANA?
-------------

The CANA package stands for "*Codes for ANalysis of Asteroids*". The tool was designed to perform scientific analysis for asteroids spectroscopic data.

As part of the PRIMitive Spectroscopic Survey (PRIMASS), we have collected and anaylized hundreds of spectra of primitive asteroids among the last years.
In this context we collected the routines that we typically use for the analysis and parametrization of these data to make a package.
As we are preparing to make the data library (PRIMASS-L) public to the comunity through 
`Small Bodie Node of the the Planetary Data Science (PDS) <https://pds-smallbodies.astro.umd.edu/>`_, it is our aim to make our science reproductible and of easy access. 

Click on the links bellow, if you would like to know more about PRIMASS or PRIMASS-L:

.. raw:: html


        <div class="primass">

            <div class="col-md-6 primass-home">
               <a href=about.html#primitive-asteroid-spectroscopic-survey><img src="https://raw.githubusercontent.com/depra/cana/spec/docs/_static/primass.png" width="60%"></a>
            </div>

            <div class="col-md-6 primass-home">
                <a href=about.html#primass-l><img src="https://raw.githubusercontent.com/depra/cana/spec/docs/_static/primass-l.png" width="40%"></a>
            </div>
        </div>




Contents
--------
The package is currently under development. At the moment CANA counts with:

Spectroscopic analysis tools:

* Handling tools, such as load, trim, fit (and autofit), estimate SNR, clean spectrum
* Slope calculation
* Taxonomic Classification
* Hydration band analysis

We will soon make available new parametrization methods and tools for handling Photometric data and Compositional modeling!

Get Started
------------

See the `install instructions <install.html>`_ to set up your computer and install Cana.

Once you have everything installed, we recommend you take a look at our recipes at the
`Cookbook <cookbook.html>`_, to get an idea of what cana can do. For more detailed 
documentation go to the `Documentation page <api.html>`_ to get the complete API description. 

Cana is still a work in progress, if you are have trouble getting it to work, 
or you want to give us a feedback or colaborate, you can:

* Report bugs, enhancements or else through `Github <https://github.com/depra/cana/issues>`_.

Cite
----
CANA is a open tool. You can use, adapt and modify it as wanted. But if you any of these, please cite us!

   De Pra, M., Carvano, J., Morate, D., Licandro, J. Pinilla-Alonso, N. (2018). CANA: An open-source python tool to study hydration in the Solar System. 

See `cite <cite.html>`_ to get the bibtex entry for the citation.

