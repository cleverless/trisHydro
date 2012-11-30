trisHydro
=========

Three-dimensional Relativistic Israel-Stewart Hydrodynamics

Hydrodynamics code designed to study the effects of viscosity and non-trivial longitudinal dynamics in relativistic heavy ion collisions.  Built to be light-weight enough to run on a single-machine with multiple threads but is not stable in the presence of shocks at this time. 

Paramter input is done by ParameterMap which is part of the utility package with CoRaL. Freely available for download at http://www.pa.msu.edu/~pratts/freecodes/home.html

Makefile included allows for linking against the boost threads library for multithread support and for high density output through xdmf.  Both features are removable by disabling preprocessor flags in src/hydrodef.h and removing the associated linking from the loader command.

The makefile has a bunch of linking crap in it that needs to be fixed, as it assumes that everything is getting compiled as part of a larger repository that is not yet publicly available.  For now this code is provided for reference and storage purposes. I can be available to help with installation and can be reached at: cleverless(at)gmail(dot)com.