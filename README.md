# MATLAB Maxwell-Solver (2D)

## Description
A 2D-Maxwell-solver using the Finite Difference Time Domain Method combined with an explicit Heat-Diffusion solver.
Its main purpose is providing a first impression of how an ultrashort laser pulse might interact with a material of some specific geometry.
The demo illustrates the partial shielding of an aluminum target by previously ablated liquid droplets.
Note: All the optical, thermophysical and transport-properties are modeled as constants here.
In general they do depend at least on the electron- and lattice-temperatures, the density and the degree of ionization.

Dispersion is taken into account by means of the Drude-Lorentz model having a single Lorentz-pole using the Auxiliary-Differential-Equation (ADE) method.
Reflected waves are absorbed at the boundaries of the simulation domain using Berenger's split field method.
Both, the TMZ and the TEZ modes are solved simultaneously.
A soft source distributes the excitation energy among all modes evenly.
The model material investigated in this script is aluminum.

![demo](<fdtd.gif>)

## Features
- user defined, arbitrary 2D geometry
- TMZ and TEZ modes calculated simultaneously
- Drude-Lorentz Dispersion with single Lorentz-pole
- computes the absorbed power density
- coupled to an explicit two-temperature heat-diffusion solver
- easily adjustable spatial and temporal laser source profles (currently both gaussian)
- contour plots for the absorbed power density, the electron- and the lattice temperature and all the E- and H-field components

## Usage

The main input-parameters are the following:
- *xdim,ydim*: simulation box domain in units of FD-cells
- *delta*: the width of a single FD-cell in meters (same for x- and y-direction)
- *lambda*: wavelength of the incident radiation in units of meters
- *I0*: laser peak intensity in units of W/m^2
- *tmax*: maximum simulation time in seconds
- *t_FWHM*: full width at half maximum duration of the gaussian pulse
- *srcx*: x-coordinate of the laser source in units of FD-cells
- *srcw*: beam radius of the gaussian pulse in units of FD-cells
- *t0*: time of peak intensity of the gaussian pulse
- *material*: a 2D-matrix describing the geometry of the sample, where 1=material and 0=vacuum
- *bound_width*: the width of the perfeclty matched layers in units of FD-cells
