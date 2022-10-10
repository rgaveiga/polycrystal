Codes in Fortran developed in 2012 during my first year as a postdoctoral fellow at Escola Politécnica, Universidade de São Paulo, in a [project](https://bv.fapesp.br/en/bolsas/132418/computer-simulations-of-the-microstructural-evolution-of-fe-ni-c-alloys/) funded by [São Paulo Research Foundation](https://fapesp.br/en).

The first program, *grain_builder*, allows you to create a three-dimensional simulation box for molecular dynamics containing at least eight grains using Voronoi tessellation. Euler angles defining each grain orientation can be determined by the user or at random. It is also possible to generate a two-dimensional, columnar microstructure, also using Voronoi tessellation, with *columnar_builder*, in which case only tilt grain boundaries are present. Both codes can be compiled by a standart Fortran 90-compliant compiler, such as *gfortran* (recommended).

Usage is quite straightforward, just provide the input needed by any of the programs while running.

Questions, comments or suggestions? Feel free to send an [e-mail](mailto:roberto.veiga@ufabc.edu.br).
