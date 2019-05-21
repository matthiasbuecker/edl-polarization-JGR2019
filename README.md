# Numerical and Analytical Modelling of Stern- and Diffuse-Layer Polarization Mechanisms in Porous Media

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3066177.svg)](https://doi.org/10.5281/zenodo.3066177)

This repository contains Matlab scripts, Comsol Multiphysics models, and numerical simulation data used to generate the plots in the manuscript

Bücker, M., Flores Orozco, A., Undorf, S., and Kemna, A., 2019, *On the Role of Stern- and Diffuse-Layer Polarization
Mechanisms in Porous Media*, submitted to JGR: Solid Earth.

If you find this data useful in your own research, please mention the this manuscript.

## Summary

> Water-saturated porous media exhibit a low-frequency (<1 MHz) dispersion of the electrical conductivity caused by the polarization of the electrical double layer (EDL) coating the charged solid-liquid interface. We develop a mathematical framework describing the polarization caused by field-induced perturbations of the ion densities in the Stern and the diffuse layer of the EDL for two different geometrical configurations of solid and liquid phase. For spherical grains immersed in an electrolyte we derive an improved analytical description by combining suitable models for diffuse- and Stern-layer polarization. The selected models differ from those usually used in geophysical literature and improve the agreement with the corresponding finite-element (FE) solution significantly. We then employ the validated FE model to examine the EDL in a pore-constriction geometry, which is often used to study membrane polarization. Here, a suitable analytical model can only be set up for a pure diffuse-layer polarization. The results for the coupled Stern- and diffuse-layer polarization in both geometries indicate that (1) the polarization of the Stern layer is much stronger than the polarization of the diffuse layer as long as the EDL is not connected at the system scale; (2) this dominance of the Stern-layer polarization can be observed in both geometries, but (3) the contribution of the diffuse layer increases with increasing compaction as represented by the pore-constriction geometry; and (4) the contributions of both parts of the EDL reach similar levels, when the EDLs on different surfaces are interconnected at the system scale.

## Instructions

* Plots for grain and membrane (pore-constriction) models are generated separately
* Access directory (*grain* or *membrane*)
* Open and run *parameters.m* first, then run *spectra.m*

Comsol models are those with the *.mph* extension. Yet, all numerical data used in the manuscript is available without running the Comsol simulation.
