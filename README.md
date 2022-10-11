Codes in Fortran developed in 2012 during my first year as a postdoctoral fellow at Escola Politécnica, Universidade de São Paulo, in a [project](https://bv.fapesp.br/en/bolsas/132418/computer-simulations-of-the-microstructural-evolution-of-fe-ni-c-alloys/) funded by [São Paulo Research Foundation](https://fapesp.br/en).

The first program, *grain_builder*, allows you to create a three-dimensional simulation box for molecular dynamics containing at least eight grains using Voronoi tessellation. Euler angles defining each grain orientation can be determined by the user or at random. It is also possible to generate a two-dimensional, columnar microstructure, also using Voronoi tessellation, with *columnar_builder*, in which case only tilt grain boundaries are present. Both codes can be compiled by a standart Fortran 90-compliant compiler, such as *gfortran* (recommended).

Usage is quite straightforward, just provide the input needed by any of the programs while running. Feel free to modify the codes according to your needs and translate them to other computer languages.

If you find them useful and if you don't mind, please cite my papers below in your papers in which you use any of the codes:

* [Stability of nanocrystalline Ni-based alloys: coupling Monte Carlo and molecular dynamics simulations](https://doi.org/10.1088/1361-651X/aa83ef)
* [Thermal conductivity of nanocrystalline SiGe alloys using molecular dynamics simulations](http://dx.doi.org/10.1063/1.4826526)

Questions, comments or suggestions? Feel free to send an [e-mail](mailto:roberto.veiga@ufabc.edu.br).
