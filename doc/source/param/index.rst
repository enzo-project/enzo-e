.. include:: ../roles.incl
	     
*************************
Enzo-E / Cello Parameters
*************************

This page documents all current parameters implemented in Enzo-E /
Cello.  Each parameter is summarized, its type or types are listed,
and the default value (if any) is provided.  The scope of the
parameter is also listed, which is either "Cello" or "Enzo", depending
on whether the parameter is associated with Cello framework or Enzo-E
application.  Any assumptions associated with a parameter are also
listed; for example, a parameter may only be valid if some other
parameter is set to a certain value.

If you find any errors in the documentation, or have any specific
suggestions, please contact the documentation maintainer at
jobordner@ucsd.edu.

-----
Adapt
-----

Adapt parameters define how the mesh hierarchy dynamically adapts to
the solution.  It is closely related to the Mesh parameters, which
defines the root grid size, number of blocks in the root grid, and
size of blocks.

.. include:: adapt.incl


-------
Balance
-------

Parameters for controlling dynamic load balancing are enclosed within
the :p:`Balance` group.  Currently only one :p:`Balance` parameter is
available, which is used to control how frequently load balancing
is performed.

.. include:: balance.incl
	     
--------
Boundary
--------

:p:`Boundary` group parameters define boundary conditions.  For simple
(non-mixed) boundary conditions, only the :p:`type` parameter is
required, e.g. :p:`Boundary` { :p:`type` = :t:`"periodic"`; }.  For more
complicated boundary conditions, the :p:`list` parameter is used to
define :p:`Boundary` subgroups, where each subgroup
specifies boundary conditions for some subset of the domain.  The
:p:`axis` and :p:`face` parameters are available to restrict boundary
conditions to a subset of faces, whereas the :p:`mask` parameter is
available for even finer control of mixed boundary conditions, which
may be time-dependent.  Inflow boundary conditions use the :p:`value`
parameter to specify field values at the boundary.

.. include:: boundary.incl

------
Domain
------

Domain parameters specify the lower and upper extents of the computational domain, using the :p:`lower` and :p:`upper` parameters.

.. include:: domain.incl

-----
Field
-----

Fields and their properties are defined using the :p:`Field` parameter
group.  All fields must be explicitly defined using the :p:`list` Field
parameter, and must match the names expected by the respective
Methods.  Properties include the number of ghost zones, precision, and whether a
field is centered or lies on some face, edge, or corner.  Some
performance-related parameters are available as well, including
alignment in memory, and memory padding between fields.

.. include:: field.incl

-----
Group
-----

.. _groups:

The :p:`Group` parameter group is used to differentiate between
different types of Field's and Particles.  For example, field groups
may include "color" and "temporary", and particle groups may include
"dark_matter" and "star".

:: 

    Group {

       list = ["color", "temporary"];

       color {
          field_list = ["species_HI", "species_HII" ]; 
       } 

       temporary {
          field_list = ["pressure", "temperature"]; 
       }

    }

Field and Particle groups can analogously be defined in the respective
Field and Particle parameter groups:

:: 

    Field {

       list = ["density", "velocity_x", "species_HI"];

       species_HI {

          group_list = ["temporary"]; 

       }

    }

Groups allow Cello applications to differentiate between these
different types of fields and particles using the ``Grouping`` class
(see ``src/Cello/data_Grouping.?pp``).

.. include:: group.incl
	     
-------
Initial
-------

The :p:`Initial` group is used to specify initial conditions.  :p:`cycle` specifies the initial cycle number (usually 0), :p:`list` specifies a list of initial conditions, which may include ``"value"`` for initializing fields directly, or other problem-specific initial condition generators.

.. include:: initial.incl
	     
------
Memory
------

Parameters in the :p:`Memory` group are used to define the behavior
of Cello's dynamic memory allocation and deallocation.

.. include:: memory.incl

----
Mesh
----

.. include:: mesh.incl


------
Method
------

.. include:: method.incl

-------
Monitor
-------

.. include:: monitor.incl

------
Output
------

Output parameters are used to specify what types of disk output to
perform and on what schedule.

.. include:: output.incl

--------
Particle
--------

Cello supports any number of particle types--e.g. `"dark"` for dark
matter particles, or `"trace"` for tracer particles.  Each particle
type in turn may have any number of attributes--e.g. `"x"` or
`"position_x"` for position, `"vx"` or `"velocity_x"` for velocity,
`"mass"`, `"id"`, etc.  Attributes can have any basic floating-point
or integer type.

All particle types must have at least attributes for position, defined
using the `position` parameter.  This allows Cello to know whether
particles have moved off of a Block, and if so to relocate them to
the correct new block.

Particle positions may be defined as integer types instead of
floating-point.  When a particle position attribute is defined as an
integer, then the coordinate value is defined relative to the enclosed
Block instead of a global coordinate system.  This can be useful both
to reduce memory usage, and to simultaneously improve accuracy--it
avoids possible catastrophic cancellation errors that are especially
large in "deep" Blocks in an AMR hierarchy whose position is far
from 0.  When positions are defined as integers, 0 is defined to be
the center of the block, and [ -*min-int* / 2 , *max-int* / 2) are the
bounds of the Block, where *min-int* is the minimum value of the
signed integer of the corresponding size.  Integer types allowed
include `"int8"`, `"int16"`, `"int32"`, and `"int64"`.  Two byte
integers `"int16"` should be sufficient for most simulations: it
has a range of [ -16384, 16384 ) within the particle's containing
Block, and ranges [-32768, -16384) and [16384, 32768) on either side
of the associated Block.

Particles are allocated and operated on in "batches".  The
`batch_size` parameter defines how many particles are in a batch.  By
operating on particles in batches, the frequency of memory operations
is greatly reduced, and functions operating on particle attributes can
be more efficient due to reduced overhead.  It should also simplify
writing particle methods to be executed on accelerators, such as
NVIDIA or AMD GPU's.

Just as with fields, particle types can be assigned to groups_.

.. include:: particle.incl

-------   
Physics
-------

.. include:: physics.incl

-------
Restart
-------

.. include:: restart.incl
   
--------
schedule
--------

"schedule" is a parameter *subgroup* that defines when to do something, such as
perform output, apply a method, or to apply the dynamic load balancer.
Schedules can be specified as a :p:`list` of values, or as an interval of
values specified using some subset of :p:`start`, :p:`stop`, and
:p:`step`.  The associated variable, set using :p:`var`, can be "cycle",
"time", or "seconds".  Here "time" refers to simulation time, and
"seconds" to wall-clock time.  At each cycle, all schedules are
checked to see if the cycle number, simulation time or wall-clock
seconds match the list or interval of values.  If there is a match,
the associated output or is performed; otherwise, it is skipped.

Note that when simulation "time" is specified, then the simulation's
time step may be reduced so that the corresponding output occurs
exactly at the specified time.

:: 

    Output {


       list = ["check", "dump", "image"];
    
       check {

          # **** write a checkpoint every 100.0 seconds ****

          schedule {
             var = "seconds";
             start = 100.0;
             step =  100.0;
          }
           ...
       }

       dump {

          # **** perform a data dump every 50 cycles until cycle 1000 ****

          schedule {
             var = "cycle";
             step =   50;
             stop = 1000;
           }
            ...
       }

       image {

          # **** write an image at times t = 1.0,  2.0, and 5.0 ****

          schedule {
             var = "time";
             list = [1.0, 2.0, 5.0];
           }
            ...
       }
    }
            
.. include:: schedule.incl

------
Solver
------

.. include:: solver.incl

--------
Stopping
--------

.. include:: stopping.incl

-------
Testing
-------

.. include:: testing.incl

-----   
Units
-----   

.. include:: units.incl
