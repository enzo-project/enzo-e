Parameter files
---------------

Whereas Enzo uses a flat parameter list with relatively simple data
types, Cello uses a more "hierarchical" parameter list, and allows
more complicated data types. The rationale is that the power of the
application is ultimately limited by the expressiveness of its input;
conversely, a rich and powerful input language can allow the
application itself to be more powerful.

This page describes how Cello parameter files are structured, which is
similar to configuration files readable by `libconfig
<http://www.hyperrealm.com/libconfig/libconfig.html>`_.  For
documentation on specific parameters in Enzo-E / Cello, please see the
`Enzo-E / Cello parameter reference
<http://cello-project.org/doc/parameters-list.html>`_ page.

Groups
******

A **group** is a named collection of related parameters or subgroups enclosed in braces:

  ::

   <group> ::= <group-name> "{" <parameter-assignment-list> "}"

All group names are capitalized, and all parameters must be contained in exactly
one group. Below, ``Domain`` and ``Mesh`` are groups, and ``lower``, ``upper``,
and ``root_size`` are parameters contained in the respective groups.

  ::

     Domain { 
       lower = [0.0, 0.0, 0.0];
       upper = [1.0, 1.0, 1.0];
     } 

     Mesh { 
       root_size = [128,128,128];
     }

Parameters within a group can be specified in different parts of a file,
or even in separate files.  Ordering of parameters within a group only matters
when a given parameter is listed more than once, in which case only its
last value is used.  Below, boundary conditions will be initialized to
be "``reflecting``".

  ::

     Boundary { 
       type = "periodic";
     } 

     Mesh { 
       root_size = [128,128,128];
     }

     Boundary { 
       type = "reflecting";
     } 

Subgroups
*********

A **subgroup** is a related collection of parameters within a group.
Subgroups typically begin with a lower case letter, and may be nested.

Below, ``Initial`` is a group, ``value`` is a subgroup, and
``cycle`` and ``density`` are parameters.

 ::

  Initial {

     value {
        density = [ x + y, x - y < 0.5, 1.0 ];
     }

     cycle = 10;

  }

The shorthand used for naming parameters in the documentation is
<group> : <subgroup> : <parameter>, e.g. "``Mesh : root_size``" or
"``Initial : value : density``".  Note that This shorthand is used
only for documentation--it is not valid syntax within the parameter file.

Parameters
**********

A **parameter** is a variable associated with a group and (optionally)
a subgroup. A parameter assignment is the assignment of a value to a
parameter. The scope of a parameter is its enclosing group and
subgroup.  Parameters may be either required or optional, and all
optional parameters have default values that are specified in the API
function call used to access the parameter value in the code.  Some
parameters may have more than one type, e.g. a scalar, or a list of
scalars.  Parameters may be assigned to multiple times, in which case
the last value is used.  All parameter assignments end with a
semi-colon ``";"``

  ::

    <parameter-assignment> ::= <parameter-name> '=' <value> ';'

Data types
**********

Parameters can be assigned different *values*, depending on the
parameter name and the group the parameter belongs to.

  ::

    <value> ::= <value-integer>
              | <value-scalar>
              | <value-string>
              | <value-variable>
              | <scalar-expression>
              | <logical-expression>
              | <list>

A parameter value is one of several different basic data types:

  ::

    ==================	=============================
    Type         	Example
    ==================	=============================
    integer 	        42
    scalar 	        6.673e-8
    string         	"velocity_x"
    variable 	        x, y, z, or t
    scalar expression 	2.0*sin(pi*x) + cos(pi*y)
    logical expression 	(x < y) || (y < z)
    list 	        [ "velocity_x", 3.14, x < y ]
    ==================	=============================

**Integer types** are integers, and must be representable using a
32-bit integer.    

    
**Scalar types** are any floating point or integral numerical values.  
The constant 'pi' is also recognized as a scalar.

   *Note that floating-point and integers are not interchangeable: if a
   floating point type is expected, one cannot use an integer.*

**String types** are enclosed in double-quotes. 

**Variables** represent the position coordinates in space (x, y, and z) and time
(t).

**Scalar expressions** are any "C-like" expression evaluating to a
Scalar, and involving Scalar's, Variable's, operations '+' '-' '*'
'/', '^' (for power), parenthesis, and (almost all) standard functions
in math.h. Scalar expressions are typically used for specifying initial or
boundary conditions, etc.

   *Note that "-" when used for subtraction must have blank space
   after it:* ``x-1.0`` *will not be parsed correctly, but* ``x -
   1.0`` *will.  Similarly, "-" when used for negation must not have a
   blank space after it.*

**Logical expressions** are any "C-like" expressions that evaluate to
"true" or "false", and involve Scalars, Variabless, and at least one
relational operator ``==`` ``!=`` ``>`` ``<`` ``<=`` ``>=``. Logical
operators ``&&`` and ``||`` are also recognized.  Logical expressions
are typically used for defining subregions of the domain for
initial  or boundary conditions.

**Lists** represent an ordered sequence of values of mixed types,
separated by commas.  Lists can be assigned a value, e.g. ``list = ["dark","star"];``,
or can be appended to, e.g. ``list += ["trace"];``  Appending to a parameter
that has not been assigned to yet is allowed, and equivalent to assignment.

Comments
********

Comments begin with # and extend to the end of the line.

Include files
*************

The ``include`` directive is used to include other parameter
declarations from other files. For example, one can have a file of
parameters for AMR that is maintained separately from problem specific
declarations:

::

   include "amr_defaults.incl"
   include "hydro_defaults.incl"

The advantage of using ``include`` is that repetition between
different parameter files can be reduced.  However, a disadvantage is
that parameters for a given run can be spread out among different
files, making it difficult to understand what parameters are defined
and their values.  Because of this, Cello writes out its parameters to
a single file ``"parameters.out"``.  Since it is a valid parameter file itself,
it can even be used to rerun the simulation, though it should be renamed first
to avoid it being overwritten.


Examples
********

Below is a list of sample input files used for developing Enzo-E
parameters. Individual parameters are expected to evolve, though the
underlying grammar and syntax are relatively fixed.

  ::

     Boundary {
         type = "reflecting";
     }

     Domain {
         lower = [ 0.0, 0.0 ];
         upper = [ 0.3, 0.3 ];
     }

     Field {

         list = [ "density", "velocity_x", "velocity_y",
                  "total_energy", "internal_energy", "pressure" ];

         courant = 0.8;
         gamma = 1.4;
         ghost_depth = 4;
     }

     Initial {
        list = ["value"];
        value {
           type = "value";
           density = [ 0.125, ( x  +  y ) <  0.1517 , 1.0 ];
           total_energy = [ 2.8, ( x  +  y ) <  0.1517 , 2.5 ];
           velocity_x = 0.0;
           velocity_y = 0.0;
        }     
     }

     Adapt {
         list = [ "SLOPE" ];

         SLOPE {
             field_list = [ "density" ];
             min_refine = 2.0;
             max_coarsen = 0.5;
             type = "slope";
         }

     }

     Mesh {

         max_level = 3;
         root_blocks = [ 2, 2 ];
         root_rank = 2;
         root_size = [ 40, 40 ];
     }

     Method {

         list = [ "ppm" ];

     }

     Output {

         list = [ "DENSITY", "MESH" ];

         DENSITY {
             name = [ "implosion-d-%03d.png", "count" ];
             type = "image";
             image_type = "data";
             field_list = [ "density" ];
             colormap = [ 0.0, 0.0, 0.0,
                          1.0, 0.0, 0.0,
                          1.0, 1.0, 0.0,
                          1.0, 1.0, 1.0 ];
             schedule {
                 step = 10;
                 type = "interval";
                 var = "cycle";
             }
         }

         MESH {
             name = [ "implosion-mesh-%03d.png", "count" ];
             type = "image";
             image_type = "mesh";
             image_reduce_type = "max";
             image_size = [ 513, 513 ];
             colormap = [ 0.0, 0.0, 0.0,
                          0.0, 0.0, 1.0,
                          0.0, 1.0, 1.0,
                          0.0, 1.0, 0.0,
                          1.0, 1.0, 0.0,
                          1.0, 0.0, 0.0 ];
             schedule {
                 step = 10;
                 type = "interval";
                 var = "cycle";
             }
         }
     }

     Stopping {
         cycle = 100;
         time = 2.50;
     }

----

2020-04-10: Updated with corrections from Joshua Smith.
     
