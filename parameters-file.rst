Parameter File Syntax
---------------------

Whereas Enzo uses a flat parameter list with relatively simple data
types, Cello uses a more "hierarchical" parameter list and allows more
complicated data types. The rationale is that the power of the
application is limited by the expressiveness of its input; conversely,
a rich and powerful input language can allow the application itself to
be more powerful.

Groups
******

A **group** is a named collection of related parameters or subgroups enclosed in braces:

  ::

   <group> ::= <group-name> "{" <parameter-assignment-list> "}"

All group names are capitalized, and all parameters belong to exactly
one group. Below, ``Domain`` and ``Mesh`` are groups, and ``extent``
and ``root_size`` are parameters contained in the respective groups.

  ::

     Domain { 
       extent = [0.0, 0.3, 0.0, 0.3] 
     } 

     Mesh { 
       root_size = [400,400];
     }
      
Subgroups
*********

A **subgroup** is a related collection of parameters within a group.
Subgroups begin with a lower case letter, and may be nested.

   *Whereas groups do not requires semi-colons* ``";"`` *after their
   definition, subgroups do*.

Below, ``Initial`` is a group, ``density`` is a subgroup, and
``cycle`` and ``value`` are parameters.

 ::

  Initial {

     density {
        value = [ x + y, x - y < 0.5, 1.0 ];
     };

     cycle = 10;

  }

The shorthand used for naming parameters in the documentation is
<group> : <subgroup> : <parameter>, e.g. ``Mesh : root_size`` or
``Initial : density : value``.  

   *This shorthand is used only for documentation--it is not valid
   syntax in the parameter file.*

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

    <value> ::= <value-scalar>
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
    scalar 	        6.673e-8
    string         	"velocity_x"
    variable 	        x, y, z, or t
    scalar expression 	2.0*sin(x) + cos(y)
    logical expression 	(x < y) || (y < z)
    list 	        [ "velocity_x", 3.14, x < y ]
    ==================	=============================

**Scalar types** are any floating point or integral numerical values.  
The constant 'pi' is also recognized.

   *Note that floating-point and integers are not interchangeable: if a
   floating point type is expected, one cannot use an integer.*

**String types** are enclosed in double-quotes. 

**Variables** represent the position in space (x, y, and z) and time
(t).

**Scalar expressions** are any "C-like" expression evaluating to a
Scalar, and involving Scalar's, Variable's, operations '+' '-' '*'
'/', '^' (for power), parenthesis, and (almost all) standard functions
in math.h. Scalar expressions are used for specifying initial or
boundary conditions, etc.

   *Note that "-" when used for subtraction must have blank space
   after it:* ``x-1.0`` *will not be parsed correctly, but* ``x -
   1.0`` *will.  Similarly, "-" when used for negation must not have a
   blank space after it.*

**Logical expressions** are any "C-like" expressions that evaluate to
"true" or "false", and involve Scalar's, Variables's, and at least one
relational operator == != > < <= >=. Logical expressions are used for
defining subregions of the domain for initializing Field's.

**Lists** represent an ordered sequence of values of mixed types, including other lists, separated by commas.

    **NEW!** *In addition to assigning a list of values to a parameter,     as in* ``list = ["dark","star"];``, *one may also append to an existing     list.  For example,* ``list = ["dark","star"];`` *followed by* ``list += ["trace"];`` *would be equivalent to* ``list = ["dark","star","trace"]``.   *Appending to a parameter that has not been assigned to yet is acceptable, and  equivalent to assignment.*

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
different parameter files can be reduced; however, a disadvantage is
that parameters for a given run can be spread out among different
files.  Because of this, Cello writes out its parameters to the file
``"parameters.out"``, which can be used to compare parameters used
with those expected.  Since it is a valid parameter file itself, it
can even be used to rerun the simulation.


Examples
********

Below is a list of sample input files used for developing Enzo-P
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
                   "total_energy", "internal_energy" ];

          courant = 0.8;
          gamma = 1.4;
          ghosts = 4;
      }

      Initial {
          density {       value = [ 0.125, ( x  +  y ) <  0.1517 , 1.0 ]; };
          total_energy {  value = [ 2.8, ( x  +  y ) <  0.1517 , 2.5 ]; };
          velocity_x {    value = 0.0; };
          velocity_y {    value = 0.0; };
      }

      Mesh {

          list = [ "SLOPE" ];

          SLOPE {
              field_list = [ "density" ];
              max_refine = 10.0;
              min_coarsen = 4.0;
              type = "slope";
          };

          max_level = 4;
          root_blocks = [ 2, 2 ];
          root_rank = 2;
          root_size = [ 48, 48 ];
      }

      Method {

          list = [ "ppm" ];

          ppm {
              diffusion = true;
              dual_energy = false;
              flattening = 3;
              steepening = true;
          };
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
              };
          };

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
              };
          };

      }

      Stopping {
          cycle = 20;
          time = 2.50;
      }
