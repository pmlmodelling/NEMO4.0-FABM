#################
Quick Start Guide
#################

.. For headings markup, this convention is recommended from Python’s Style Guide
   # with overline, for parts
   * with overline, for chapters
   =, for sections
   -, for subsections
   ^, for subsubsections
   ", for paragraphs

   'global.rst' contains a list of roles, substitutions and URL links for the guide.
   It is loaded for each file with the 'rst_epilog' setting in 'conf.py'

.. toctree::
   :hidden:

   todos

.. Only displayed with 'make drafthtml'

.. toctree::
   :hidden:
   :caption: Unpublished READMEs
   :glob:
   :titlesonly:

   unpub/*

.. Only displayed with 'make drafthtml'

.. toctree::
   :hidden:
   :caption: Getting started

   install
   cfgs
   tests

.. toctree::
   :hidden:
   :caption: Setup your configuration

   setup
   diags
   tools

.. toctree::
   :hidden:
   :caption: Advanced use

   zooms
   cplg
   da
   tracers

.. toctree::
   :hidden:
   :caption: Miscellaneous

   cite
   contrib
   changes
   acro

:Release:  |release|

.. only:: draft

   .. Only on draft

   :*Date*:     |today|
   :*SVN rev*:  |revision|

.. include:: readme.rst
