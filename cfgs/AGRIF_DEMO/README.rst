**************
Embedded zooms
**************

.. todo::



.. contents::
   :local:

Overview
========

AGRIF (Adaptive Grid Refinement In Fortran) is a library that
allows the seamless space and time refinement over rectangular regions in NEMO.
Refinement factors can be odd or even (usually lower than 5 to maintain stability).
Interaction between grid is "two-ways" in the sense that
the parent grid feeds the child grid open boundaries and
the child grid provides volume averages of prognostic variables once
a given number of time step is completed.
These pages provide guidelines how to use AGRIF in NEMO.
For a more technical description of the library itself, please refer to AGRIF_.

Compilation
===========

Activating AGRIF requires to append the cpp key ``key_agrif`` at compilation time:

.. code-block:: sh

   ./makenemo [...] add_key 'key_agrif'

Although this is transparent to users,
the way the code is processed during compilation is different from the standard case:
a preprocessing stage (the so called ``conv`` program) translates the actual code so that
saved arrays may be switched in memory space from one domain to an other.

Definition of the grid hierarchy
================================

An additional text file :file:`AGRIF_FixedGrids.in` is required at run time.
This is where the grid hierarchy is defined.
An example of such a file, here taken from the ``ICEDYN`` test case, is given below

.. literalinclude:: ../../../tests/ICE_AGRIF/EXPREF/AGRIF_FixedGrids.in

The first line indicates the number of zooms (1).
The second line contains the starting and ending indices in both directions on the root grid
(``imin=34 imax=63 jmin=34 jmax=63``) followed by the space and time refinement factors (3 3 3).
The last line is the number of child grid nested in the refined region (0).
A more complex example with telescoping grids can be found below and
in the :file:`AGRIF_DEMO` reference configuration directory.

.. todo::

   Add some plots here with grid staggering and positioning?

When creating the nested domain, one must keep in mind that
the child domain is shifted toward north-east and
depends on the number of ghost cells as illustrated by
the *attempted* drawing below for ``nbghostcells=1`` and ``nbghostcells=3``.
The grid refinement is 3 and ``nxfin`` is the number of child grid points in i-direction.

.. image:: _static/agrif_grid_position.jpg

Note that rectangular regions must be defined so that they are connected to a single parent grid.
Hence, defining overlapping grids with the same refinement ratio will not work properly,
boundary data exchange and update being only performed between root and child grids.
Use of east-west periodic or north-fold boundary conditions is not allowed in child grids either.
Defining for instance a circumpolar zoom in a global model is therefore not possible.

Preprocessing
=============

Knowing the refinement factors and area,
a ``NESTING`` pre-processing tool may help to create needed input files
(mesh file, restart, climatological and forcing files).
The key is to ensure volume matching near the child grid interface,
a step done by invoking the :file:`Agrif_create_bathy.exe` program.
You may use the namelists provided in the :file:`NESTING` directory as a guide.
These correspond to the namelists used to create ``AGRIF_DEMO`` inputs.

Namelist options
================

Each child grid expects to read its own namelist so that different numerical choices can be made
(these should be stored in the form :file:`1_namelist_cfg`, :file:`2_namelist_cfg`, etc...
according to their rank in the grid hierarchy).
Consistent time steps and number of steps with the chosen time refinement have to be provided.
Specific to AGRIF is the following block:

.. literalinclude:: ../../namelists/namagrif
   :language: fortran

where sponge layer coefficients have to be chosen according to the child grid mesh size.
The sponge area is hard coded in NEMO and applies on the following grid points:
2 x refinement factor (from ``i=1+nbghostcells+1`` to ``i=1+nbghostcells+sponge_area``)

.. rubric:: References

.. bibliography:: zooms.bib
   :all:
   :style: unsrt
   :labelprefix: A
   :keyprefix: a-
