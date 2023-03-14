Overview
========

`NEMO <https://www.nemo-ocean.eu/>`_ (Nucleus for European Modelling of the Ocean) is a state-of-the-art modelling framework for
research activities and forecasting services in ocean and climate sciences.

`FABM <https://github.com/fabm-model/fabm>`_ (Framework for Aquatic Biogeochemical Models) is a programming framework for biogeochemical 
models of marine and freshwater systems.

This repository adds the coupling within the NEMO-TOP module to run biogeochemical models through
FABM for NEMO4.0. The starting NEMO branch is `version 4.0.4, revision 14301 <http://forge.ipsl.jussieu.fr/nemo/browser/NEMO/releases/r4.0/r4.0.4?rev=14301>`_. 
Details of the changes made, along with an example compilation can be seen on the repository `wiki <https://github.com/pmlmodelling/NEMO4.0-FABM/wiki>`_. 

NEMO relies on `XIOS <https://forge.ipsl.jussieu.fr/ioserver>`_ to provide input/output. A brief overview of 
the format of XIOS xml files is available `here <https://github.com/pmlmodelling/NEMO4.0-FABM/wiki/Using-XIOS>`_.

Optimisation
========

As part of the ARCHER2 eCSE programme an extensive report has been produced on optimisation of the NEMO4.0-FABM codebase. This report examines the scalability of the code, both through a change in core count and a change in the number of tracers, and assesses the functions with the largest overhead when making these changes. There are a number of recommendations for node placement and I/O environment settings when running on ARCHER2 that may be relevant elsewhere. 

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7732984.svg
   :target: https://doi.org/10.5281/zenodo.7732984


