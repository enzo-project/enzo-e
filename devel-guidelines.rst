.. include:: roles.incl

************************
Enzo-E Coding Guidelines
************************

**[ This page is under development:** *last updated 2019-03-26* **]**


This page describes current coding guidelines for Enzo-E development.  It is
a working document and always open to suggestions for changes.  

-----------------------

============
Naming files
============

All Enzo-E-related files should be in the ``src/Enzo`` subdirectory of
the Enzo-E source code directory.  File names have the form
*component*\ ``_``\ *class*\ ``.[hc]pp``, where *component* is the
high-level code component (for Enzo-E code it is always "``enzo``"),
*class* is the class-name (see below for Enzo-E class naming
guidelines), and extension ``.hpp`` for header (declaration) files and
``.cpp`` for source (definition) files.

+---------------------+---------------------------------------------+
| File type           |  Enzo-E file name                           |
+=====================+=============================================+
| Physics method      |  ``enzo_EnzoMethod``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Linear solver       |  ``enzo_EnzoSolver``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Initial conditions  | ``enzo_EnzoInitial``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Boundary conditions | ``enzo_EnzoBoundary``\ *Name*\ ``.[hc]pp``  |
+---------------------+---------------------------------------------+
| Interpolation       | ``enzo_EnzoProlong``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Restriction         | ``enzo_EnzoRestrict``\ *Name*\ ``.[hc]pp``  |
+---------------------+---------------------------------------------+

==============
Naming classes
==============

Enzo-E classes all begin with ``Enzo``, followed by the (capitalized)
general type indicating the type of Cello class ("``Method``",
"``Solver``", "``Initial``" for initial conditions, etc.), followed by
the specific name for the class type.  For example, the Conjugate
Gradient (CG) linear solver class is named ``EnzoSolverCg``.  Note
that all code implementing an Enzo-E class should (generally) live in
the files ``enzo_``\ *class-name*\ ``.hpp`` and ``enzo_``\ *class-name*\ ``.cpp``.

+---------------------+---------------------------+------------------+
| Class type          |  Enzo-E class name        | Cello base class |
+=====================+===========================+==================+
| Physics method      |  ``EnzoMethod``\ *Name*   | ``Method``       |
+---------------------+---------------------------+------------------+
| Linear solver       |  ``EnzoSolver``\ *Name*   | ``solver``       |
+---------------------+---------------------------+------------------+
| Initial conditions  | ``EnzoInitial``\ *Name*   | ``Initial``      |
+---------------------+---------------------------+------------------+
| Boundary conditions | ``EnzoBoundary``\ *Name*  | ``Boundary``     |
+---------------------+---------------------------+------------------+
| Interpolation       | ``EnzoProlong``\ *Name*   | ``Prolong``      | 
+---------------------+---------------------------+------------------+
| Restriction         | ``EnzoRestrict``\ *Name*  | ``Restrict``     |
+---------------------+---------------------------+------------------+

====================
Naming class methods
====================

public methods thing_1()
private methods thing_2_()
entry methods p_blah()
reduction entry methods r_reduce()
Naming variables
Array dimensions mx,my,mz
Active region size nx,ny,nz
Ghost zone depth gx,gy,gz
Loop variables ix,iy,iz

====================
Accessing Field data
====================

Field field = block->data()->field();
id = field.field_id("density");
field.dimensions(id, &mx, &my, &mz);
field.size (&nx, &ny, &nz);
field.ghost_depth(id &gx, &gy, &gz);
double * d = (double *) field.values(id);
for (int iz=gz; iz<gz+nz; iz++) {
for (int iy=gy; iy<gy+ny; iy++) {
for (int ix=gx; ix<gx+nx; ix++) {
int i = ix + mx*(iy + my*iz);
d[i] = tiny ;
}
}
}

===================
Accessing Particles
===================
