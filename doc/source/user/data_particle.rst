****************
Enzo-E Particles
****************

*[ This page is under development ]*

This section describes all particle types in Enzo-E, including
particle attributes, where each particle type is used and modified,
and groups of related particles.

* dark
* trace
* star
* sink

IMPORTANT: If running a simulation with particles that have mass, and
you want to analyse the output with yt, you must include a "mass_is_mass"
parameter in the "Particle" group in the input parameter file. It can be
set equal to anything, so long as it can be read by yt. It is recommended that
the value is set to be "true".

In addition, you will need to install the latest development version of yt from source.
This can be done with:

git clone https://github.com/yt-project/yt.git

cd yt

pip install -e .

