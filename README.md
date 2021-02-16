# Optometrika

![optometrika](https://user-images.githubusercontent.com/46988982/51661552-0c009200-1f66-11e9-8d38-79f35f6ac8d8.png)

## Overview

Optometrika MATLAB library implements analytical and iterative ray tracing approximation to optical image formation using Snell’s and Fresnel’s laws of refraction and reflection.

Currently, the library implements refractive and reflective general surfaces, aspheric (conic) surfaces with astigmatism, Fresnel surfaces, cones and cylinders (elliptic too), planes, circular and ring-shaped apertures, rectangular flat screens, spheroidal screens, and a realistic model of the human eye with accommodating lens and spheroidal retina. See example*.m files for examples of ray tracing in general (user-defined shape) lenses, aspheric lenses, Fresnel lenses, prisms, mirrors, and human eye. 

The library traces refracted rays, including intensity loss at the refractive surface. Reflected rays are currently traced for mirrors and also for a single total internal reflexion or double refraction, if it happens. Note that the Bench class object is not a real physical bench, it is only an ordered array of optical elements, and it is your responsibility to arrange optical objects in the right order. In particular, if you need to trace rays passing through the same object multiple times, you have to add the object multiple times to the bench array in the order the object is encountered by the rays. For example, double refraction/reflection for cylindrical and conical surfaces can be calculated by adding the surface twice to the bench (see example10). 

The library is compact and fast. It was written using Matlab classes and is fully vectorized. It takes about 2 seconds to trace 100,000 rays through an external lens and the human eye (8 optical surfaces) on a 3 GHz Intel Core i7 desktop. Fresnel lens tracing is somewhat slower due to looping through the Fresnel cones describing the lens surface. Tracing through user-defined (general) surfaces is significantly slower due to iterative search of ray intersections with the surface. 

Thank you for downloading Optometrika, enjoy it!


## List of examples:

example1.m: tests the basic functionality of the Optometrika library

example2.m: demonstrates the Optometrika's optical model of the human eye

example3.m: demonstrates accommodation of the human eye by minimizing the retinal image

example4.m: tests a ring lens with the cosine surface profile defined in coslens.m

example5.m: tests planar mirrors

example6.m: tests planar and parabolic mirrors (a Newtonian refractor telescope)

example7.m: tests a Fresnel lens

example8.m: tests a lens with polynomial aspheric terms

example9.m: tests cone mirrors

example10.m: tests cylinder and cone surfaces with double refraction

example11.m: demonstrates ray tracing for rays originating inside the human eye

example12.m: draws a lens and determines its front surface, back surface, and total height. Makes an animated gif of the lens and an engineering drawing of the lens.

example13.m: tests refraction through the lens edge and backward rays refraction (sub-aperture Maksutov-Cassegrain telescope)

example14.m: tests refraction through a lens with astigmatism (different vertical and horizontal radii of curvature)

example15.m: simulates a hexagonal array of spherical micro lenses

example16.m: demonstrates STL export of various lenses

example17.m: demonstrates tracing of optical pathlength (opl) for three different lenses

example18.m: demonstrates optical path length through wedged optical windows

example19.m: tests precision of the numerical solver in default and precision mode


## To acknowledge use of this software, please cite the following publication:

A. Schultze and Y. Petrov. OptoMetrika -Open source ray tracer for MATLAB. 2021.DOI:10.13140/RG.2.2.10027.98081/1.

We also welcome opportunities to collaborate as co-author(s) on project applications that make use of OptoMetrika, and we aim to continue supporting and developing LoadDef based on community input. If you have any questions about the software or would like to discuss a particular project idea further, please get in touch! Thanks!
