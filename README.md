# MATLAB Maxwell-Solver (2D)

## Description
A 2D-Maxwell-solver using the Finite Difference Time Domain Method combined with an explicit Heat-Diffusion solver.
Its main purpose is providing a first impression of how an ultrashort laser pulse might interact with a material of some specific geometry.
The demo illustrates the partial shielding of an aluminum target by previously ablated liquid droplets.
Note: All the optical, thermophysical and transport-properties are modeled as constants here.
In general they do depend at least on the electron- and lattice-temperatures, the density and the degree of ionization.

![example](<fdtd.gif>)

## Features
- user defined, arbitrary density profile
- you can adjust the angle of incidence and the polarizatio of the incoming light (s- or p-polarization)
- plots the profile for the calculated absorbed power density
- outputs the integral absorption, reflection and transmission
- the material's permittivity can be modeled as a function of wavelenght, density and temperature

## Usage

In order to compute the absolute power density in units of W/m^3, simply multiply the relative power density by the incident laser intensity.
The main input-parameters are the following:
- *m_polar*: either 1 (for s-polarization) or 2 (for p-polarization)
- *lambda*: wavelength of the incident radiation in units of meters
- *theta*: angle of incidence
- *nelements*: nr of piecewise constant material elements
- *delta*: the widht of each element in units of nano-meters
- *dprof*-vector: an adjustable vector containing the density profile of the material. The same can be done for the temperature.
- *getEpsilon(lambda,Te,Ti,rho)*: a function, giving the relative permittivity as a function of electron temperature, ion-temperature, density and wavelength. 
