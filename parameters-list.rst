  .. role:: p
  .. role:: g
  .. role:: t
  .. role:: s
  .. role:: d
  .. role:: e
  .. role:: o
  .. role:: c
  .. role:: z

  .. |H+| replace:: H\ :sup:`+`
  .. |D+| replace:: D\ :sup:`+`
  .. |H++| replace:: H\ :sup:`++`
  .. |++| replace:: \ :sup:`++`
  .. |H-| replace:: H\ :sup:`-`
  .. |He+| replace:: He\ :sup:`+`
  .. |H2+| replace:: H\ :sub:`2`:sup:`+`
  .. |H| replace:: H
  .. |H2| replace:: H\ :sup:`2`
  .. |He| replace:: He
  .. |e-| replace:: e\ :sup:`-`
  .. |aij| replace:: a\ :sub:`i,j`		    
  .. |a00| replace:: a\ :sub:`0,0`		    
  .. |a10| replace:: a\ :sub:`1,0`		    
  .. |am0| replace:: a\ :sub:`m,0`		    
  .. |a01| replace:: a\ :sub:`0,1`		    
  .. |a0n| replace:: a\ :sub:`0,n`		    
  .. |a0j| replace:: a\ :sub:`0,j`		    
  .. |a1j| replace:: a\ :sub:`1,j`		    
  .. |amj| replace:: a\ :sub:`m,j`		    
  .. |ai0| replace:: a\ :sub:`i,0`		    
  .. |ai1| replace:: a\ :sub:`i,1`		    
  .. |ain| replace:: a\ :sub:`i,n`		    
  .. |amn| replace:: a\ :sub:`m,n`		    
  .. |am1| replace:: a\ :sub:`m,1`
  .. |a1n| replace:: a\ :sub:`1,n`

***************
Parameters List
***************

This page documents all current parameters implemented in Enzo-P /
Cello.  Each parameter is summarized, its type or types are listed,
and the default value (if any) is provided.  The scope of the
parameter is also listed, which is either "Cello" or "Enzo", depending
on whether the parameter is associated with Cello code or Enzo-P code.
Some parameters also may have specified assumptions associated with
them; for example, a parameter may only be valid if some
other parameter is set to a certain value.

The last parameters review was performed on 2015-09-10 with code
revision 3836.  If you find any errors in the documentation, or have
any specific suggestions, please contact the documentation maintainer
at jobordner@ucsd.edu.

-----
Adapt
-----

Adapt parameters define how the mesh hierarchy dynamically adapts to
the solution.  It is closely related to the Mesh parameters, which
defines the root grid size, number of blocks in the root grid, and
size of blocks.

----

:Parameter: :p:`Adapt` : :p:`interval`
:Summary:   :s:`Number of cycles between adapt steps`
:Type:      :t:`integer`
:Default:   :d:`1`
:Scope:     :c:`Cello`

:e:`The interval parameter is used to set the number of root-level cycles between mesh adaptation.  The default is 1.`

----

:Parameter: :p:`Adapt` : :p:`max_level`
:Summary:   :s:`Maximum level in the adaptive mesh hierarchy`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`This parameter specifies the level of the most highly refined Block in the mesh hierarchy.  The default is 0, meaning there is no refinement past the initial root-level grid.`

----

:Parameter: :p:`Adapt` : :p:`min_level`
:Summary:   :s:`Minimum level in the adaptive mesh hierarchy`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`This parameter specifies the coarsest level of "sub-root" Blocks, and must non-positive.  This is used primarily for multigrid methods, such as in the` :t:`"gravity_mg0"` :e:`method.  The default is 0, meaning no sub-root Blocks are created.  If multigrid is used, then both` :p:`Adapt` : :p:`min_level` :e:`and` :p:`Method` : :p:`gravity_mg` : :p:`min_level` :e:`must be set.`

----

:Parameter: :p:`Adapt` : :p:`list`
:Summary:   :s:`List of refinement criteria`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`List of mesh refinement criteria, each of which has its own associated` ``Adapt : <criteria> :`` :e:`parameters.  When multiple criteria are used, if all refinement criteria evaluate to "coarsen", then the block will coarsen; if any refinement criteria evaluate as "refine", then the block will refine.`

:e:`The items in the list need not be the same as the (required)` :p:`Adapt` : :g:`<criterion>` : :p:`type` :e:`parameter; they are solely used to identify and distinguish between different criteria in the simulation.  This allows the user to use multiple criteria of the same type but with different parameters, e.g. "mask" with different masks:`

   ::

       Adapt {
          list = ["criterion_1", "criterion_2"];
          criterion_1 {
             type = "shock";
          };
          criterion_2 {
             type = "shear";
          }
       }

----

:Parameter:  :p:`Adapt` : :p:`min_face_rank`
:Summary:    :s:`Minimum rank of Block faces to check for 2:1 refinement restriction`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`Many numerical methods require a 2:1 refinement restriction on adaptive meshes, such that no Block in level i is adjacent to another Block in a level j with |i - j|>1.  This assumption may be required across corners and edges as well as 2D faces.  This parameter specifies the minimum rank (dimensionality) of Block faces across which to enforce the 2:1 refinement restriction.`

----

:Parameter: :p:`Adapt` : :g:`<criterion>` : :p:`field_list`
:Summary:   :s:`List of field the refinement criterion is applied to`
:Type:        [ :t:`string` | :t:`list` ( :t:`string` ) ]
:Default:     :d:`[]` ( all fields )
:Scope:     :c:`Cello`

:e:`This parameter specifies the fields that the refinement criteria is applied to.  For example, if type = "slope" and field_list = ["density"], then the "refine by slope" refinement criterion is applied to the density field.  This is only used for refinement by slope.`


----

:Parameter: :p:`Adapt` : :g:`<criterion>` : :p:`level_exponent`
:Summary:   :s:`Level exponent parameter`
:Type:        :t:`float`
:Default:     :d:`0.0`
:Scope:     :c:`Cello`
:Assumes:   :g:`<criterion>` is of :p:`type` :t:`"mass"`

:e:`The level exponent parameter is used in the "mass" refinement criterion type only.  It is used as a scaling factor for the refinement criteria for different mesh levels.`


----

:Parameter: :p:`Adapt` : :g:`<criterion>` : :p:`max_coarsen`
:Summary:   :s:`Cutoff value for coarsening a block`
:Type:        [ :t:`float` | :t:`list` ( :t:`float` ) ]
:Default:     :d:`0.15`
:Scope:     :c:`Cello`

:e:`A block may coarsen if the refinement criterion applied to the block is smaller than this value everywhere in the block.   A list is used for the` :t:`"shock"` :e:`refinement criterion type, in which case the first value is for pressure and the second is for the energy ratio.`

----

:Parameter: :p:`Adapt` : :g:`<criterion>` : :p:`include_ghosts`
:Summary:   :s:`Whether to include ghost zones when applying the refinement criterion`
:Type:      :t:`logical`
:Default:   :d:`false`
:Scope:     :c:`Cello`

:e:`When applying a mesh refinement criterion, this parameter specifies whether to apply it to ghost zones in the block as well as non-ghost zones.`

----

:Parameter: :p:`Adapt` : :g:`<criterion>` : :p:`min_refine`
:Summary:   :s:`Cutoff value for refining a block`
:Type:        [ :t:`float` | :t:`list` ( :t:`float` ) ]
:Default:     :d:`0.3`
:Scope:     :c:`Cello`

:e:`A block must refine if the refinement criterion applied to the block is larger than this value anywhere in the block.  A list is used for the` :t:`"shock"` :e:`refinement criterion type, in which case the first value is for pressure and the second is for the energy ratio.`

----

:Parameter:  :p:`Adapt` : :g:`<criterion>` : :p:`output`
:Summary:    :s:`Name of a field in which to store the result of the refinement criterion`
:Type:    :t:`string`
:Default: :d:`""`
:Scope:     :c:`Cello`

:e:`In addition to evolved field values, one may also output the refinement criteria.  This may be  useful for example for debugging or for finding appropriate values for max_coarsen and min_refine.  A value of -1 specifies coarsening, +1 for refining, and 0 for staying the same.`

----

:Parameter:  :p:`Adapt` : :g:`<criterion>` : :p:`type`
:Summary:    :s:`Type of mesh refinement criteria`
:Type:    :t:`string`
:Default: :d:`"unknown"`
:Scope:     :c:`Cello`

:e:`Type of mesh refinement criteria.  This is a required parameter, and must be one of "slope", "shear", "mask", "mass", or "shock".`
 

-------
Balance
-------

Parameters for controlling dynamic load balancing are enclosed within
the :p:`Balance` group.  Currently only one :p:`Balance` parameter is
available, which is used to control how frequently load balancing
is performed.

----

:Parameter:  :p:`Balance` : :p:`schedule` : :g:`<schedule_parameter>`
:Summary:    :s:`Scheduling parameters for dynamic load balancing`
:Type:       :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`Dynamic load balancing is scheduled according to` :p:`schedule` :e:`parameters.  Scheduling parameters---including` :p:`var`, :p:`list`, :p:`start`, :p:`stop`, and :p:`step` :e:`---are documented in the` `schedule`_ :e:`section.`

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

----

:Parameter:  :p:`Boundary` : :p:`list`
:Summary:    :s:`List of boundary condition subgroups`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`For mixed boundary conditions, the` :p:`list` :e:`parameter specifies the list of names of subgroups that define boundary conditions on each portion of the domain boundary.  Boundary conditions in each subgroup are applied in the order listed.  In the example below, two subgroups` :t:`"one"` :e:`and` :t:`"two"` :e:`are defined, which specify reflecting boundary conditions along the x-axis and outflow boundary conditions along the y-axis:`

  ::

       Boundary {
          list = ["one", "two"];
          one {
             type = "reflecting";
             axis = "x";
          };
          two {
             type = "outflow";
             axis = "y";
          }
       }

----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`type`
:Summary:    :s:`Type of boundary condition`
:Type:    :t:`string`
:Default: :d:`"undefined"`
:Scope:     :c:`Cello`

:e:`Boundary conditions in Enzo-P include` :t:`"reflecting"` :e:`,` :t:`"outflow"` :e:`,` :t:`"inflow"` :e:`, and` :t:`"periodic"`.  :e:`Other boundary condition types can be implemented by either a) modifying the existing` :p:`EnzoBoundary` :e:`class or b) creating a new class inherited from the` :p:`Boundary` :e:`base class.`  :t:`"inflow"` :e:`boundary conditions additionally require` :p:`value` :e:`and` :p:`field_list` :e:`parameters.`

----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`axis`
:Summary:    :s:`Axis along which boundary conditions are to be enforced`
:Type:    :t:`string`
:Default: :d:`"all"`
:Scope:     :c:`Cello`

:e:`The`  :p:`axis` :e:`parameter restricts the boundary conditions to the face orthogonal to the specified axis.`  :p:`axis` :e:`must be` :t:`"x"` , :t:`"y"` , :t:`"z"` :e:`or` :t:`"all"`.  :e:`The` :p:`axis` :e:`parameter may be used in conjunction with the` :p:`face` :e:`parameter, or by itself.`


----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`face`
:Summary:    :s:`Face along which boundary conditions are to be enforced`
:Type:    :t:`string`
:Default: :d:`"all"`
:Scope:     :c:`Cello`

:e:`The` :p:`face` :e:`parameter can restrict the boundary conditions to be applied only to the` :p:`upper` :e:`or` :p:`lower` :e:`faces.  face orthogonal to the given face.`  :p:`face` :e:`must be` :t:`"x"` , :t:`"y"` , :t:`"z"` :e:`or` :t:`"all"`.  :e:`The` :p:`face` :e:`parameter may be used in conjunction with the` :p:`axis` :e:`parameter, or by itself.`

----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`mask`
:Summary:    :s:`Subregion in which boundary conditions are to be enforced`
:Type:    :t:`logical-expr`
:Default: :d:`none`
:Scope:     :c:`Cello`

:e:`The`  :p:`mask` :e:`parameter specifies the subregion of the boundary on which to apply the boundary conditions.  The logical expression  may be a function of x, y, z, and t, and boundary conditions are restricted to where (and when) it evaluates to true`::

       Boundary {
          ...
          OUT {
             type = "outflow";
             mask = (x >= 4.0) || 
                    (y >= 1.0 && (x >= 0.744017 + 11.547* t));
          };
       }


----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`value`
:Summary:    :s:`Value for inflow boundary conditions`
:Type:    :t:`float`
:Type:    :t:`float-expr`
:Type:    :t:`list` ( :t:`float-expr` [, :t:`logical-expr`, :t:`float-expr` [, ... ] ] )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`The` :p:`value` :e:`parameter is used to specify field values for` :p:`inflow` :e:`type boundary conditions.  The` :p:`value` :e:`parameter is used in conjunction with the` :p:`field_list` :e:`parameter.` :p:`value` :e:`may be of type` :t:`float`, :t:`float-expr`, :e:`or a list of alternating` :t:`float-expr` :e:`and` :t:`logical-expr` :e:`types`.  :t:`float-expr` :e:`may be a function of x, y, z, and t.  When a list is specified, the` :t:`logical-expr` :e:`is treated as a mask, similar to an 'if-then-else' clause`

   ::

       Boundary {
          ...
          VELOCITY_Y {
             type = "inflow";
             field_list = "velocity_y";
             value = [ -8.25*0.5,
                        ((x <= 0.166667) && (y <= 0.0) ) ||
                         (x <= 0.0) ||
                         ((x < 0.744017 + 11.547*t) && (y >= 1.0))
                      ];
           };
       }


----

:Parameter:  :p:`Boundary` : :g:`<condition>` : :p:`field_list`
:Summary: :s:`List of fields to apply boundary conditions to`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`The` :p:`field_list` :e:`parameter is used to restrict boundary conditions to the specified fields.  An empty list, which is the default, is used to specify all fields.`

------
Domain
------

Domain parameters specify the lower and upper extents of the computational domain, using the :p:`lower` and :p:`upper` parameters.

----

:Parameter:  :p:`Domain` : :p:`lower`
:Summary: :s:`Lower domain extent`
:Type:    :t:`list` ( :t:`float` )
:Default: :d:`[0.0, 0.0, 0.0]`
:Scope:     :c:`Cello`

:e:`Lower extent of the computational domain,` [x\ :sub:`min`], [ x\ :sub:`min`\, y\ :sub:`min`], :e:`or` [ x\ :sub:`min`\, y\ :sub:`min`\, z\ :sub:`min`].

----

:Parameter:  :p:`Domain` : :p:`upper`
:Summary: :s:`Upper domain extent`
:Type:    :t:`list` ( :t:`float` )
:Default: :d:`[1.0, 1.0, 1.0]`
:Scope:     :c:`Cello`

:e:`Upper extent of the computational domain,` [x\ :sub:`max`], [ x\ :sub:`max`\, y\ :sub:`max`], :e:`or` [ x\ :sub:`max`\, y\ :sub:`max`\, z\ :sub:`max`].



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


----

:Parameter:  :p:`Field` : :p:`list`
:Summary: :s:`List of fields`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`All fields must be explicitly listed in the` :p:`list` :e:`parameter.  Field names depend on the Method(s) used; e.g., PPM uses` :t:`"density"`,  :t:`"velocity_x"`, :t:`"velocity_y"`, :t:`"total_energy"`, :e:`and`  :t:`"internal_energy"`.
  

----

:Parameter:  :p:`Field` : :p:`gamma`
:Summary: :s:`Adiabatic exponent`
:Type:    :t:`float`
:Default: :d:`5.0 / 3.0`
:Scope:     :z:`Enzo`
:Todo:  :o:`perhaps move this to a different group, e.g. Physics`

:p:`gamma` :e:`specifies the ratio of specific heats for the ideal gas used by the PPM hydrodynamics solver.`


----

:Parameter:  :p:`Field` : :p:`alignment`
:Summary: :s:`Force field data on each block to start on alignment bytes`
:Type:    :t:`integer`
:Default: :d:`8`
:Scope:     :c:`Cello`

:e:`Depending on the computer architecture, variables can be accessed from memory faster if they have at least 4-byte or 8-byte alignment.  This parameter forces each field block array to have an address evenly divisible by the specified number of bytes.`

----

:Parameter:  :p:`Field` : :g:`<field>` : :p:`centering`
:Summary: :s:`Specify the position of the given field variable within the computational cell.`
:Type:    :t:`list` ( :t:`logical` )
:Default: :d:`[ true, true, true ]`
:Scope:     :c:`Cello`

:e:`By default, variables are centered within a computational cell.  Some methods expect some variable, e.g. velocity components, to be positioned on a cell face.  The effect of this parameter is to increase the dimension of the field block by one along each axis with a value of "false".  Numerical method implementations like PPML that assume (NX,NY,NZ) sized blocks even for offset variables, as opposed to e.g. (NX+1,NY,NZ), should still define the variable as centered.`

----

:Parameter:  :p:`Field` : :g:`<field>` : :p:`group_list`
:Summary: :s:`Specify a list of groups that the Field belongs to`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[ ]`
:Scope:     :c:`Cello`

:e:`Different Fields may belong to any number of different "groups".  For example, Enzo uses "color fields", which Enzo-P implements as defining color fields to belong to the group "color".`

----

:Parameter:  :p:`Field` : :p:`courant`
:Summary: :s:`Courant safety factor for fields`
:Type:    :t:`float`
:Default: :d:`0.6`
:Scope:     :c:`Cello`
:Todo:    :o:`Rename?`

:e:`Courant safety factor for all fields.  This is a multiplication factor for the time step as determined by the respective Method(s) used.  This parameter can be updated on restart using the` `Restart : file` :e:`restart parameter file.`

----

:Parameter:  :p:`Field` : :p:`ghost_depth`
:Summary: :s:`Field ghost zone depths`
:Type:    [ :t:`integer` | :t:`list` ( :t:`integer` ) ]
:Default: :d:`[ 0, 0, 0 ]`
:Scope:     :c:`Cello`

:e:`The default storage patch / block ghost zone depths [gx, gy, gz] along each axis for fields.  If an integer, then the same ghost zone depth is used for each axis.  Currently this value needs to be $4$ for PPM when AMR is used.`

----

:Parameter:  :p:`Field` : :p:`padding`
:Summary: :s:`Add padding of the specified number of bytes between fields on each block.`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`If block sizes are large and a power of two, and if the computer's cache has low associativity, performance can suffer due to cache thrashing.  This can be avoided by introducing padding between fields.  A value of twice the cache line width is recommended.  Since field blocks are usually small, this should not usually be an issue.`

----

:Parameter:  :p:`Field` : :p:`precision`
:Summary: :s:`Default field precision`
:Type:    :t:`string`
:Default: :d:`"default"`
:Scope:     :c:`Cello`

:e:`Default precision for all fields.  Supported precisions include "single" (32-bit) and "double" (64-bit).  "quadruple" is accepted, but not implemented by most numerical methods (e.g. PPM).  "default" is for compatibility with Enzo, and corresponds to either "single" or "double" depending on the CELLO_PREC configuration flag setting.  This precision parameter must not conflict with the CELLO_PREC setting.`

----

:Parameter:  :p:`Field` : :p:`prolong`
:Summary: :s:`Type of prolongation (interpolation)`
:Type:    :t:`string`
:Default: :d:`"linear"`
:Scope:     :c:`Cello`

:e:`For adaptive mesh refinement, field values may need to be transferred from coarser to finer blocks, either from coarse neighbor blocks in the refresh phase, or to fine child blocks during refinement in the adapt phase.  Valid values include` :t:`"linear"` :e:`; other values accepted but not implemented include` :t:`"enzo"` :e:`and` :t:`"MC1"` :e:` ; which are unfinished implementations of Enzo's` :t:`"InterpolationMethod"` :e:`functionality.`

----

:Parameter:  :p:`Field` : :p:`restrict`
:Summary: :s:`Type of restriction (coarsening)`
:Type:    :t:`string`
:Default: :d:`"linear"`
:Scope:     :c:`Cello`

:e:`For adaptive mesh refinement, field values may need to be transferred from finer to coarser blocks, either from fine neighbor blocks in the refresh phase, or to the parent block during coarsening in the adapt phase.  Valid values include` :t:`"linear"` :e:`\; ;other values accepted but not implemented include` :t:`"enzo"`.

----

:Parameter:  :p:`Field` : :p:`interpolation_method`
:Summary: :s:`Type of "enzo" interpolation and coarsening`
:Type:    :t:`string`
:Default: :d:`"SecondOrderA"`
:Scope:     :c:`Cello`
:Status:  **Not accessed**

:e:`For the "enzo"` :p:`prolong` :e:`or` :p:`restrict` :e:`Field parameters, this parameter defines the specific interpolation method used.  It is analogous to the` ``InterpolationMethod`` :e:`parameter in Enzo.  Valid values include` ``"ThirdOrderA"`` ,   ``"SecondOrderA"`` ,    ``"SecondOrderB"``, ``"SecondOrderC"`` , :e:`and` ``"FirstOrderA"``.

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
       }; 

       temporary {
          field_list = ["pressure", "temperature"]; 
       };

    }

Field and Particle groups can analogously be defined in the respective
Field and Particle parameter groups:

::

    Field {

       list = ["density", "velocity_x", "species_HI"];

       species_HI {

          group_list = ["temporary"]; 

       };

    }

Groups allow Cello applications to differentiate between these
different types of fields and particles using the ``Grouping`` class
(see ``src/Cello/data_Grouping.?pp``).

----

:Parameter:  :p:`Group` : :p:`list`
:Summary: :s:`List of groups`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`This parameter defines all groups.`

----

:Parameter:  :p:`Group` : :g:`<group>` : :p:`field_list`

:Summary: :s:`List of fields belonging to the group`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`This parameter is used to assign fields to a given group.`

----

:Parameter:  :p:`Group` : :g:`<group>` : :p:`particle_list`

:Summary: :s:`List of particle types belonging to the group`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`This parameter is used to assign particle groups to a given group.`

-------
Initial
-------

The :p:`Initial` group is used to specify initial conditions.  :p:`cycle` specifies the initial cycle number (usually 0), :p:`type` specifies the type of initial conditions, either ``"value"`` for initializing fields directly, or other problem-specific initial condition generators.

----

:Parameter:  :p:`Initial` : :p:`cycle`
:Summary: :s:`Initial cycle number`
:Type:    :t:`list` ( :t:`integer` )
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`Initial value for the cycle number.`


----

:Parameter:  :p:`Initial` : :p:`type`
:Summary: :s:`Identifier specifying the type of initial conditions`
:Type:    :t:`string`
:Default: :d:`"value"`
:Scope:     :c:`Cello`

:e:`This parameter specifies how the field variables in the simulation will be initialized.  The default is` ``"value"`` :e:`, in which case field variables are initialized directly in the input file, e.g.`

    ::

       Initial {
          type = "value";
	  density { value = [ 0.125, x + y < 0.5, 1.0 ]; };
          ...
       }

:e:`Here, the` :p:`density` :e:`field is initialized to be 0.125 in the subregion of the domain where x + y < 0.5, and is initialized to 1.0 elsewhere.  Other valid values for` :p:`type` :e:`include` ``"implosion_2d"``, ``"sedov_array_2d"``, ``"sedov_array_3d"``, :e:`and` ``grackle_test``, :e:`which are problem-specific initializers analogous to those in the original Enzo application.`

----

:Parameter:  :p:`Initial` : :p:`time`
:Summary: :s:`Initial time`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     :c:`Cello`

:e:`Initial time in code units.`

----

:Parameter:  :p:`Initial` : :g:`<field>` : :p:`value`
:Summary: :s:`Initialize field values`
:Type:    :t:`list` ( :t:`float-expr`, [ :t:`logical-expr`, :t:`float-expr`, [ ... ] ] )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`This parameter is used to initialize fields when the` :p:`type` :e:`parameter is` ``"value".``  :e:`The first element of the list must be a` :t:`float` :e:`expression, and may include arithmetic operators, variables "x", "y", "z", and most functions in the POSIX math library /include/math.h.  The second optional list element is a logical expression, and  serves as a "mask" of the domain.  The third` :t:`float` :e:`expression parameter is required if a mask is supplied, and serves as the "else" case.  Multiple such mask-value pairs may be used.  Example: [ sin ( x + y ), x - y < 0, 1.0 ] is read as "sin ( x + y ) where x - y < 0, 1.0 elsewhere".`

sedov
-----

:Parameter:  :p:`Initial` : :p:`sedov` : :p:`array`
:Summary: :s:`Size of array of Sedov blasts`
:Type:    :t:`list` ( :t:`integer` )
:Default: :d:`[ 1, 1, 1 ]`
:Scope:   :z:`Enzo`

:e:`This parameter defines the size of the array of Sedov blast waves.  The default is a single blast.`

----

:Parameter:  :p:`Initial` : :p:`sedov` : :p:`radius_relative`
:Summary: :s:`Initial radius of the Sedov blast`
:Type:    :t:`float`
:Default: :d:`0.1`
:Scope:   Enzo  
:Todo:    :o:`write`

----

:Parameter:  :p:`Initial` : :p:`sedov` : :p:`pressure_in`
:Summary: :s:`Pressure inside the Sedov blast`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     Enzo  
:Todo:    :o:`write`

----

:Parameter:  :p:`Initial` : :p:`sedov` : :p:`pressure_out`
:Summary: :s:`Pressure outside the Sedov blast`
:Type:    :t:`float`
:Default: :d:`1.0e-5`
:Scope:     Enzo  
:Todo:    :o:`write`

----

:Parameter:  :p:`Initial` : :p:`sedov` : :p:`density`
:Summary: :s:`Density for the Sedov blast array problem`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     Enzo  
:Todo:    :o:`write`

turbulence
----------

:Parameter:  :p:`Initial` : :p:`turbulence` : :p:`density`
:Summary: :s:`Initial density for turbulence initialization and method`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     Enzo  

:e:`Initial density for initializing the turbulence problem.`

----

:Parameter:  :p:`Initial` : :p:`turbulence` : :p:`pressure`

:Summary: :s:`Initial pressure for turbulence initialization and method`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     Enzo  

:e:`Initial pressure for initializing the turbulence problem.  Default is 0.0, meaning it is not used.  Either` `pressure` :e:`or` `temperature` :e:`should be defined, but not both.`

----

:Parameter:  :p:`Initial` : :p:`turbulence` : :p:`temperature`
:Summary: :s:`Initial temperature for turbulence initialization and method`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     Enzo  

:e:`Initial temperature for initializing the turbulence problem.  Default is 0.0, meaning it is not used.  Either` `pressure` :e:`or` `temperature` :e:`should be defined, but not both.`

------
Memory
------

Parameters in the :p:`Memory` group are used to define the behavior
of Cello's dynamic memory allocation and deallocation.

----

:Parameter:  :p:`Memory` : :p:`active`
:Summary: :s:`Whether to track memory usage`
:Type:    :t:`logical`
:Default: :d:`true`
:Scope:     :c:`Cello`

:e:`This parameter is used to turn on or off Cello's build-in memory tracking.  By default it is on, meaning it tracks the number and size of memory allocations, including the current number of bytes allocated, the maximum over the simulation, and the maximum over the current cycle.  Cello implements this by overloading C's new, new[], delete, and delete[] operators.  This can be problematic on some systems, e.g. if an external library also redefines these operators, in which case this parameter should be set to false.`

----
Mesh
----

:Parameter:  :p:`Mesh` : :p:`root_blocks`
:Summary: :s:`Number of Blocks used to tile the coarsest Patch`
:Type:    :t:`list` ( :t:`integer` )
:Default: :d:`[ 1, 1, 1 ]`
:Scope:     :c:`Cello`

:e:`This parameter specifies the number of Blocks along each axis in the mesh "forest".  The product must not be smaller than the number of processors used.`

----

:Parameter:  :p:`Mesh` : :p:`root_rank`
:Summary: :s:`Physical dimensionality of the problem`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`Number of physical dimensions in the problem, 1, 2, or 3.`

----

:Parameter:  :p:`Mesh` : :p:`root_size`
:Summary: :s:`Coarsest Patch size`
:Type:    :t:`list` ( :t:`integer` )
:Default: :d:`[ 1, 1, 1 ]`
:Scope:     :c:`Cello`

:e:`This parameter specifies the total size of the root-level mesh.  For example, [400, 400] specifies a two dimensional root-level discretization of 400 x 400 zones, excluding ghost zones.`

------
Method
------

:Parameter:  :p:`Method` : :p:`list`
:Summary: :s:`Sequence of numerical methods to apply.`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`none`
:Scope:     :c:`Cello`

:e:`This parameter specifies the list of numerical methods to use.  Each method in the list is applied in the order specified.  Possible values include:`

  *  :t:`"ppm"` :e:`for Enzo-P's PPM hydrodynamics method`
  *  :t:`"ppml"` :e:`for the PPML ideal MHD solver`
  *  :t:`"heat"` :e:`for the forward-Euler heat-equation solver, which is used primarily for demonstrating how new Methods are implemented in Enzo-P`
  *  :t:`"grackle"` (not fully implemented yet)  :e:`for heating and cooling methods in the Enzo Grackle library`
  * :t:`"gravity_[cg|bicgstab|mg0]"` :e:`for solving for the gravitational potential using the conjugate gradient method (CG), BiCG-STAB, or multigrid on the root-grid.`
  * :t:`"null"` :e:`for "no solver", which is used to specify time step size for testing the AMR meshing infrastructure without undue computation.`

  :e:`Parameters specific to individual methods are specified in subgroups, e.g.`::

     Method {
        list = ["ppm"];
        ppm {
           diffusion   = true;
           flattening  = 3;
           steepening  = true;
           dual_energy = false;
        }
     }


----

:Parameter:  :p:`Method` : :p:`courant`
:Summary: :s:`Global Courant safety factor`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     :c:`Cello`

:e:`The global Courant safety factor is a multiplication factor for the time step applied on top of any Field or Particle specific Courant safety factors.`

cosmology
---------

Currently cosmology parameters are not accessed.

----

:Parameter:  :p:`Method` : :p:`cosmology`
:Summary: :s:`Turn on or off cosmology machinery`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Turn on or off cosmology machinery.`

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`comoving_box_size`
:Summary: :s:`Enzo's CosmologyComovingBoxSize parameter`
:Type:    :t:`float`
:Default: :d:`64.0`
:Scope:     :z:`Enzo`

:e:`Enzo's` CosmologyComovingBoxSize :e:`parameter.`

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`hubble_constant_now`
:Summary: :s:`Hubble constant for Z=0`
:Type:    :t:`float`
:Default: :d:`0.701`
:Scope:     :z:`Enzo`

:e:`Hubble constant for Z=0.`  

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`initial_redshift`
:Summary: :s:`Enzo's CosmologyInitialRedshift parameter.`
:Type:    :t:`float`
:Default: :d:`20.0`
:Scope:     :z:`Enzo`

:e:`Enzo's` CosmologyInitialRedshift :e:`parameter.`

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`max_expansion_rate`
:Summary: :s:`Maximum expansion rate`
:Type:    :t:`float`
:Default: :d:`0.01`
:Scope:     :z:`Enzo`

:e:`Maximum expansion rate.`

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`omega_lamda_now`
:Summary: :s:`Omega lambda for Z=0`
:Type:   :t:`float`
:Default: :d:`0.721`
:Scope:     :z:`Enzo`

:e:`Omega lamda for Z=0.`

----

:Parameter:  :p:`Method` : :p:`cosmology` : :p:`omega_matter_now`
:Summary: :s:`Omega matter for Z=0`
:Type:    :t:`float`
:Default: :d:`0.279`
:Scope:     :z:`Enzo`

:e:`Omega matter for Z=0.`

gravity_bicgstab
----------------

:Parameter:  :p:`Method` : :p:`gravity_bicgstab` : :p:`iter_max`
:Summary: :s:`Iteration limit for the BiCGStab solver`
:Type:    :t:`int`
:Default: :d:`100`
:Scope:     :z:`Enzo`

:e:`Maximum number of BiCGStab iterations to take.`

----

:Parameter:  :p:`Method` : :p:`gravity_bicgstab` : :p:`res_tol`
:Summary: :s:`Residual norm reduction tolerance for the BiCGStab solver`
:Type:    :t:`float`
:Default: :d:`1e-6`
:Scope:     :z:`Enzo`

:e:`Stopping tolerance on the 2-norm of the residual relative to the initial residual, i.e. BiCGStab is defined to have converged when ||R_i ||`:sub:`2` `/ ||R_0 ||`:sub:`2` `< res_tol.`

----

:Parameter:  :p:`Method` : :p:`gravity_bicgstab` : :p:`grav_const`
:Summary: :s:`Gravitational constant`
:Type:    :t:`float`
:Default: :d:`6.67384e-8`
:Scope:     :z:`Enzo`

:e:`Gravitational constant used in place of G.  The default is G in cgs units.`


----

:Parameter:  :p:`Method` : :p:`gravity_bicgstab` : :p:`diag_precon`
:Summary: :s:`Whether to apply diagonal preconditioning`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Whether to diagonally precondition the linear system A*X = B in BiCGStab by 1.0 / (h^2).`


----

:Parameter:  :p:`Method` : :p:`gravity_bicgstab` : :p:`monitor_iter`
:Summary: :s:`How often to display progress`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`The current iteration, and minimum, current, and maximum relative residuals, are displayed every monitor_iter iterations.  If monitor_iter is 0, then only the first and last iteration are displayed.`

gravity_cg
----------

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`iter_max`
:Summary: :s:`Iteration limit for the CG solver`
:Type:    :t:`int`
:Default: :d:`100`
:Scope:     :z:`Enzo`

:e:`Maximum number of CG iterations to take.`

----

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`res_tol`
:Summary: :s:`Residual norm reduction tolerance for the CG solver`
:Type:    :t:`float`
:Default: :d:`1e-6`
:Scope:     :z:`Enzo`

:e:`Stopping tolerance on the 2-norm of the residual relative to the initial residual, i.e. CG is defined to have converged when ||R_i ||`:sub:`2` `/ ||R_0 ||`:sub:`2` `< res_tol.`

----

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`grav_const`
:Summary: :s:`Gravitational constant`
:Type:    :t:`float`
:Default: :d:`6.67384e-8`
:Scope:     :z:`Enzo`

:e:`Gravitational constant used in place of G.  The default is G in cgs units.`


----

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`diag_precon`
:Summary: :s:`Whether to apply diagonal preconditioning`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Whether to diagonally precondition the linear system A*X = B in EnzoMethodGravityCg by 1.0 / (h^2).`


----

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`monitor_iter`
:Summary: :s:`How often to display progress`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`The current iteration, and minimum, current, and maximum relative residuals, are displayed every monitor_iter iterations.  If monitor_iter is 0, then only the first and last iteration are displayed.`

gravity_mg
----------

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`iter_max`

:Summary: :s:`Maximum number of multigrid cycles.`
:Type:    :t:`int`
:Default: :d:`10`
:Scope:     :z:`Enzo`

:e:`Maximum number of cycles of the multigrid solver.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`res_tol`

:Summary: :s:`Residual norm reduction limit for the multigrid solver`
:Type:    :t:`float`
:Default: :d:`1e-6`
:Scope:     :z:`Enzo`

:e:`Stopping tolerance on the 2-norm of the residual relative to the initial residual, i.e. multigrid is defined to have converged when ||R_i ||`:sub:`2` `/ ||R_0 ||`:sub:`2` `< res_tol.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`grav_const`
:Summary: :s:`Gravitational constant`
:Type:    :t:`float`
:Default: :d:`6.67384e-8`
:Scope:     :z:`Enzo`

:e:`Gravitational constant used in place of G.  The default is G in cgs units.`

----

:Parameter:  :p:`Method` : :p:`gravity_cg` : :p:`monitor_iter`
:Summary: :s:`How often to display progress`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`The current iteration, and minimum, current, and maximum relative residuals, are displayed every monitor_iter iterations.  If monitor_iter is 0, then only the first and last iteration are displayed.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`smooth`

:Summary: :s:`Multigrid smoother`
:Type:    :t:`string`
:Default: :d:`"jacobi"`
:Scope:     :z:`Enzo`

:e:`The Compute object to use for smoothing the residual in the multigrid method.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`smooth_weight`

:Summary: :s:`Multigrid smoother weighting`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     :z:`Enzo`

:e:`The weighting for the multigrid smoother.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`smooth_pre`

:Summary: :s:`Number of multigrid pre-smoothings`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`Number of applications of the smoother for pre-smoothings`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`smooth_post`

:Summary: :s:`Number of multigrid post-smoothings`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`Number of applications of the smoother for post-smoothings`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`smooth_coarse`

:Summary: :s:`Number of multigrid smoothings for coarse solver`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :z:`Enzo`

:e:`Number of applications of the smoother for approximating the solution on the coarsest level.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`restrict`

:Summary: :s:`Multigrid restrict operation`
:Type:    :t:`string`
:Default: :d:`"linear"`
:Scope:     :z:`Enzo`

:e:`The Restrict type to use for transferring residuals to parent (coarser) blocks.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`prolong`

:Summary: :s:`Multigrid prolong operation`
:Type:    :t:`string`
:Default: :d:`"linear"`
:Scope:     :z:`Enzo`

:e:`The Prolong type to use for transferring corrections to child (finer) blocks.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`min_level`

:Summary: :s:`The coarsest multigrid level`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :z:`Enzo`

:e:`The coarsest level in the multigrid algorithm, which is the level at which the coarse grid solver is applied.  This number should be negative.`

----

:Parameter:  :p:`Method` : :p:`gravity_mg` : :p:`max_level`

:Summary: :s:`the finest multigrid level`
:Type:    :t:`integer`
:Default: :d:`Adapt:max_level`
:Scope:     :z:`Enzo`

:e:`The finest level in the multigrid algorithm.  Could be less than the finest level in the mesh hierarchy to improve computational speed at the cost of reduced accuracy.`


grackle
-------

"Grackle is a chemistry and radiative cooling library for astrophysical
simulations. It is a generalized and trimmed down version of the
chemistry network of the Enzo simulation code."

Most of the descriptions of the parameters come from the `Grackle documentation <http://grackle.readthedocs.org/en/grackle-1.0/index.html>`_; for the
most up-to-date description of Grackle parameters, see the `Grackle parameters <http://grackle.readthedocs.org/en/grackle-1.0/Parameters.html#id1>`_ section of the website.

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`density_units`

:Summary: :s:`Units for the density field`
:Type:    :t:`float`
:Default: :d:`1.67e-24  (1 m_H/cc)`
:Scope:     :z:`Enzo`

:e:`Units of density for the Grackle chemistry and cooling solver library.`


----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`length_units`

:Summary: :s:`Units for distance`
:Type:    :t:`float`
:Default: :d:`3.086e21 (1 kpc)`
:Scope:     :z:`Enzo`

:e:`Units of length for the Grackle chemistry and cooling solver library.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`time_units`

:Summary: :s:`Units for time`
:Type:    :t:`float`
:Default: :d:`3.15569e13 (1 Myr)`
:Scope:     :z:`Enzo`

:e:`Units of time for the Grackle chemistry and cooling solver library.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`a_units`
:Summary: :s:`Units for the cosmological expansion factor`
:Type:    :t:`float`
:Default: :d:`1.0`
:Scope:     :z:`Enzo`

:e:`Units of the cosmological expansion factor for the Grackle chemistry and cooling solver library.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`gamma`
:Summary: :s:`The ratio of specific heats for an ideal gas`
:Type:    :t:`float`
:Default: :d:`5/3`
:Scope:     :z:`Enzo`

:e:`The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if` :p:`primordial_chemistry` :e:`> 1.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`with_radiative_cooling`
:Summary:    :s:`Include radiative cooling`
:Type:       :t:`logical`
:Default:    :d:`true`
:Scope:     :z:`Enzo`

:e:`Flag to include radiative cooling and actually update the thermal energy during the chemistry solver. If off, the chemistry species will still be updated. The most common reason to set this to off is to iterate the chemistry network to an equilibrium state.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`primordial_chemistry`
:Summary: :s:`Flag to control which primordial chemistry network is used`
:Type:    :t:`logical`
:Default:  :d:`false`
:Scope:     :z:`Enzo`

:e:`Flag to control which primordial chemistry network is used.`

  **0:** :e:`no chemistry network. Radiative cooling for primordial species is solved by interpolating from lookup tables calculated with Cloudy. A simplified set of functions are available (though not required) for use in this mode. For more information, see` `Pure Tabulated Mode <http://grackle.readthedocs.org/en/grackle-1.0/Integration.html#tabulated-mode>`_.

  **1:** :e:`6-species atomic H and He. Active species:` |H|, |H+|, |He|, |He+|, |++|, |e-|.

  **2:** :e:`9-species network including atomic species above and species for molecular hydrogen formation. This network includes formation from the` |H-| :e:`and` |H2+| :e:`channels, three-body formation` ( |H| + |H| + |H|  :e:`and`  |H| + |H| + |H2|), |H2| :e:`rotational transitions, chemical heating, and collision-induced emission (optional). Active species: above +` |H-|, |H2|, |H2+|.

  **3:** :e:`12-species network include all above plus HD rotation cooling. Active species: above plus D,` |D+|, :e:`HD.`

  **Note:** :e:`In order to make use of the non-equilibrium chemistry network (primordial_chemistry options 1-3), you must add and advect baryon fields for each of the species used by that particular option.`


----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`metal_cooling`
:Summary:  :s:`Flag to enable metal cooling using the Cloudy tables`
:Type:     :t:`logical`
:Default:  :d:`false`
:Scope:     :z:`Enzo`

:e:`Flag to enable metal cooling using the Cloudy tables. If enabled, the cooling table to be used must be specified with the Grackle` :p:`data_file` :e:`parameter.`

**Note:** :e:`In order to use the metal cooling, you must add and advect a metal density field.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`h2_on_dust`
:Summary:     :s:`Flag to enable H2 formation`
:Type:        :t:`logical`
:Default:     :d:`false`
:Scope:     :z:`Enzo`

:e:`Flag to enable H2 formation on dust grains, dust cooling, and dust-gas heat transfer follow Omukai (2000). This assumes that the dust to gas ratio scales with the metallicity.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`cmb_temperature_floor`
:Summary:    :s:`Flag to enable an effective CMB temperature floor.`
:Type:       :t:`logical`
:Default:    :d:`true`
:Scope:     :z:`Enzo`

:e:`Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`data_file`
:Summary:     :s:`Path to the data file containing the metal cooling and UV background tables.`
:Type:        :t:`string`
:Default:     :d:`""`
:Scope:     :z:`Enzo`

:e:`Path to the data file containing the metal cooling and UV background tables.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`three_body_rate`
:Summary:      :s:`Flag to control which three-body H2 formation rate is used.`
:Type:        :t:`integer`
:Default:     :d:`0`
:Scope:     :z:`Enzo`
:Status:  **Not accessed**


:e:`Flag to control which three-body H2 formation rate is used.`

   **0:** `Abel, Bryan & Norman (2002) <http://adsabs.harvard.edu/abs/2002Sci...295...93A>`_

   **1:** `Palla, Salpeter & Stahler (1983) <http://adsabs.harvard.edu/abs/1983ApJ...271..632P>`_

   **2:** `Cohen & Westberg (1983) <http://adsabs.harvard.edu/abs/1983JPCRD..12..531C>`_

   **3:** `Flower & Harris (2007) <http://adsabs.harvard.edu/abs/2007MNRAS.377..705F>`_

   **4:** `Glover (2008) <http://adsabs.harvard.edu/abs/2008AIPC..990...25G.>`_

   :e:`These are discussed in` `Turk et. al. (2011) <http://adsabs.harvard.edu/abs/2011ApJ...726...55T>`_

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`cie_cooling`
:Summary:    :s:`Flag to enable |H2| collision-induced emission cooling`
:Type:        :t:`logical`
:Default:     :d:`false`
:Scope:     :z:`Enzo`

:e:`Flag to enable` |H2| :e:`collision-induced emission cooling from` `Ripamonti & Abel (2004) <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`_.

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`h2_optical_depth_approximation`
:Summary: :s:`Flag to enable |H2| cooling attenuation`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`
   
:e:`Flag to enable H2 cooling attenuation from` `Ripamonti & Abel (2004) <http://adsabs.harvard.edu/abs/2004MNRAS.348.1019R>`_.

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`photoelectric_heating`
:Summary: 
:Type:    
:Default:
:Scope:     :z:`Enzo`

:e:`Flag to enable a spatially uniform heating term approximating photo-electric heating from dust from Tasker & Bryan (2008)http://adsabs.harvard.edu/abs/2008ApJ...673..810T.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`photoelectric_heating_rate`
:Summary: 
:Type:    
:Default:  :d:`8.5e-26`
:Scope:     :z:`Enzo`

:e:`If` :p:`photoelectric_heating` :e:`is enabled, the heating rate in units of erg cm-3 s-1.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`UVbackground`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo:  :o:`write`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`UVbackground_redshift_on`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo:       :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`UVbackground_redshift_off`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo:       :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`UVbackground_redshift_fullon`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

   

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`UVbackground_redshift_drop`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

   


----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`Compton_xray_heating`
:Summary: 
:Type:    
:Default:   :d:`0`
:Scope:     :z:`Enzo`


:e:`Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999)http://adsabs.harvard.edu/abs/1999ApJ...517L...9M.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`LWbackground_intensity`
:Summary: 
:Type:    
:Default: :d:`0`
:Scope:     :z:`Enzo`

:e:`Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field in units of 10-21 erg s-1 cm-2 Hz-1 sr-1.`

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`LWbackground_sawtooth_suppression`
:Summary: 
:Type:    
:Default: :d:`0`
:Scope:     :z:`Enzo`


:e:`Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000)http://adsabs.harvard.edu/abs/2000ApJ...534...11H.`


----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`HydrogenFractionByMass`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`DeuteriumToHydrogenRatio`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`SolarMetalFractionByMass`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`NumberOfTemperatureBins`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`ih2co`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`ipiht`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`TemperatureStart`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`TemperatureEnd`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`comp_xray`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`temp_xray`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`CaseBRecombination`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`NumberOfDustTemperatureBins`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`DustTemperatureStart`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`DustTemperatureEnd`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

----

:Parameter:  :p:`Method` : :p:`grackle` : :p:`cloudy_electron_fraction_factor`
:Summary: 
:Type:    
:Default: 
:Scope:     :z:`Enzo`
:Todo: :o:`write`
:Status:  **Not accessed**

heat
----

:Parameter:  :p:`Method` : :p:`heat` : :p:`alpha`
:Summary:    :s:`Parameter for the forward euler heat equation solver`
:Type:       :t:`float`
:Default:    :d:`1.0`
:Scope:     :z:`Enzo`

:e:`Thermal diffusivity parameter for the heat equation.`

null
----

:Parameter:  :p:`Method` : :p:`null` : :p:`dt`
:Summary:    :s:`Set the time step for the "null" Method`
:Type:       :t:`float`
:Default:    :d:`max (float)`
:Scope:     :z:`Enzo`

:e:`Sets the time step for the` :p:`null` :e:`Method.  This is typically used for testing the AMR meshing infrastructure without having to use any specific method.  It can also be used to add an additional maximal time step value for other methods.`

ppm
---

:p:`Method : ppm` parameters are used to initialize parameters for
Enzo-P's PPM hydrodynamics method.

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`density_floor`
:Summary: :s:`Lower limit on density`
:Type:   :t:`float`
:Default: :d:`1.0e-6`
:Scope:     :z:`Enzo`

:e:`Density floor, which replaces Enzo's "tiny_number".`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`diffusion`
:Summary: :s:`PPM diffusion parameter`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`PPM diffusion parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`dual_energy`

:Summary: :s:`Whether to use dual-energy formalism`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Whether to use the dual-energy formalism.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`dual_energy_eta_1`
:Summary: :s:`Dual energy parameter eta 1`
:Type:   :t:`float`
:Default: :d:`0.001`
:Scope:     :z:`Enzo`

:e:`First dual-energy formalism parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`dual_energy_eta_2`
:Summary: :s:`Dual energy parameter eta 2`
:Type:   :t:`float`
:Default: :d:`0.1`
:Scope:     :z:`Enzo`

:e:`Second dual-energy formalism parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`flattening`
:Summary: :s:`PPM flattening parameter`
:Type:   :t:`integer`
:Default: :d:`3`
:Scope:     :z:`Enzo`

:e:`PPM flattening parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`minimum_pressure_support_parameter`
:Summary: :s:`Enzo's MinimumPressureSupportParameter`
:Type:   :t:`integer`
:Default: :d:`100`
:Scope:     :z:`Enzo`

:e:`Enzo's`  MinimumPressureSupportParameter :e:`parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`number_density_floor`
:Summary: :s:`Lower limit on number density`
:Type:   :t:`float`
:Default: :d:`1.0e-6`
:Scope:     :z:`Enzo`

:e:`Number density floor, which replaces Enzo's "tiny_number".`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`pressure_floor`
:Summary: :s:`Lower limit on pressure`
:Type:   :t:`float`
:Default: :d:`1.0e-6`
:Scope:     :z:`Enzo`

:e:`Pressure floor, which replaces Enzo's "tiny_number".`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`pressure_free`
:Summary: :s:`Pressure-free flag`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Pressure-free flag.` 

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`steepening`
:Summary: :s:`PPM steepening parameter`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`PPM steepening parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`temperature_floor`
:Summary: :s:`Lower limit on temperature`
:Type:   :t:`float`
:Default: :d:`1.0e-6`
:Scope:     :z:`Enzo`

:e:`Temperature floor, which replaces Enzo's "tiny_number".`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`use_minimum_pressure_support`
:Summary: :s:`Minimum pressure support`
:Type:   :t:`logical`
:Default: :d:`false`
:Scope:     :z:`Enzo`

:e:`Enzo's` UseMinimumPressureSupport :e:`parameter.`

----

:Parameter:  :p:`Method` : :p:`ppm` : :p:`mol_weight`
:Summary: :s:`Mean molecular mass`
:Type:   :t:`float`
:Default: :d:`0.6`
:Scope:     :z:`Enzo`

:e:`Mean molecular mass used in computing temperature.`


turbulence
----------

----

:Parameter:  :p:`Method` : :p:`turbulence` : :p:`edot`
:Summary: :s:`Initial value for edot for turbulence Method`
:Type:    :t:`float`
:Default: :d:`-1.0`
:Scope:     :z:`Enzo`
:Todo: :o:`write`

----

:Parameter:  :p:`Method` : :p:`turbulence` : :p:`mach_number`
:Summary: :s:`Value for Mach number in turbulence problem`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     :z:`Enzo`
:Todo: :o:`write`

-------
Monitor
-------

:Parameter:  :p:`Monitor` : :p:`debug`
:Summary: :s:`Whether to display debugging output`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :c:`Cello`

:e:`If true, then process DEBUG() statements, writing the output to both stderr and appending to files out.debug.<proc>, where <proc> is the (physical) process rank.  Note that out.debug.<proc> files are not erased at the start of a run. This parameter is not scalable and is inefficient since output files are continually opened and closed by each process.`

----

:Parameter:  :p:`Monitor` : :p:`verbose`
:Summary: :s:`Whether to display "verbose" output`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :c:`Cello`

:e:`If true, then output requests with Monitor::verbose() will be called.  This will generally produce more detailed output, such as which specific Blocks are refining and coarsening, etc.`

------
Output
------

Output parameters are used to specify what types of disk output to
perform and on what schedule.

----

:Parameter:  :p:`Output` : :p:`list`
:Summary: :s:`List of output file sets`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`List of active file sets, each of which has its own associated Output : <file_set> : parameters.  Any file set parameters associated with a file set not in the` `list` :e:`parameter are ignored.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`axis`
:Summary: :s:`Axis of projections for image output`
:Type:    :t:`string`
:Default: :d:`none`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`For the "image" output type, the axis along which to project the data for 3D problems.  Values are` `"x", "y", :e:`or` "z".  :e:`See the associated type parameter.` 

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`colormap`
:Summary: :s:`Color map for image output`
:Type:    :t:`list` ( :t:`float` )
:Default: :d:`[]`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`For the "image" output type, a list of the form` [r\ :sub:`0`\, g\ :sub:`0`\, b\ :sub:`0`\, r\ :sub:`1`\, g\ :sub:`1`\, b\ :sub:`1`\, ...], :e:`where` 0.0 ≤ r\ :sub:`i`\,g\ :sub:`i`\,b\ :sub:`i`\ ≤ :e:`1.0 are RGB values.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`field_list`
:Summary: :s:`List of fields to output`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`List of fields for this output file set.  For "image" field types, the field list must contain exactly one field.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`particle_list`
:Summary: :s:`List of particle types to output`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`List of particles types for this output file set..`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`name`
:Summary: :s:`File names`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`""`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is *not* of :p:`type` :t:`"restart"`

:e:`This parameter specifies the names of files in the corresponding file_group.  The first element is the file name, which may contain printf-style formatting fields.  Subsequent values correspond to variables for the formatting fields, which may include "cycle", "time", "count" (a counter incremented each time output is performed), and "proc" (the physical processor rank).  The file name should include an appropriate extension, e.g. ".png" for "image" output, and ".h5" or ".h5" for "data" output.  Example: ["projection-%04d.png", "cycle"].`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`dir`
:Summary: :s:`Name of the directory for restart dumps`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`""`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"restart"`

:e:`This parameter specifies the names of output restart parameter files.  The first element is the file name, which may contain printf-style formatting fields.  Subsequent values correspond to variables for the formatting fields, which may include "cycle", "time", "count" (a counter incremented each time output is performed), and "proc" (the physical processor rank).  Example: ["Restart-%02d", "count"].`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`stride`
:Summary: :s:`Subset of processors to perform write`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"data"`
:Status:    **Broken: see bug** # 13_

.. _13: http://client64-249.sdsc.edu/cello/bug/show_bug.cgi?id=13

:e:`This parameter allows for a strict subset  of physical processors to output data, which is especially helpful for large process counts  to reduce the load on parallel file systems.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`type`
:Summary: :s:`Type of output files`
:Type:    :t:`string`
:Default: :d:`"unknown"`
:Scope:     :c:`Cello`

:e:`The type of files to output in this output file set.  Supported types include "image" (PNG file of 2D fields, or projection of 3D fields) and "data".  For "image" files, see the associated colormap and axis parameters.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_min`
:Summary: :s:`Data value associated with the first color in the colormap`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`This parameter specifies the Field value associated with the first color in the file set's colormap.` **This value is only used if the** :p:`image_specify_bounds` **parameter is** :p:`true`.  :e:`If` :p:`image_specify_bounds` :e:`is` :p:`false`, :e:`then the minimum global value of the field is used instead.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_max`
:Summary: :s:`Data value associated with the last color in the colormap`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`This parameter specifies the Field value associated with the last color in the file set's colormap.` **This value is only used if the** :p:`image_specify_bounds` **parameter is** :p:`true`.  :e:`If` :p:`image_specify_bounds` :e:`is` :p:`false`, :e:`then the maximum global value of the field is used instead.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_specify_bounds`
:Summary: :s:`Whether to use` :p:`image_min` :s:`and` :p:`image_max`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`This parameter determines whether to use the` :p:`image_min` :e:`and` :p:`image_max` :e:`parameters for mapping the Field data to the color map, or to use the Field data's minimum and maximum.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_ghost`
:Summary: :s:`Whether to include ghost zones in the image`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`Setting the` :p:`image_ghost` :e:`to true will include ghost zone values in the image output.  This is typically used only when debugging.  The default is false.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_reduce_type`
:Summary: :s:`How to handle 3D field data orthogonal to the image`
:Type:    :t:`string`
:Default: :d:`"sum"`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`When images are generated for 3D problems, multiple data values will be associated with each pixel in the image.  This parameter defines how to handle these multiple values, including` :t:`"sum"`, :t:`"min"`, :t:`"max"`, :e:`and`, :t:`"avg"`.  :e:`For field data the default of` :t:`"sum"` :e:`is appropriate, though for images of meshes` :t:`"max"` :e:`should be used`.

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_face_rank`
:Summary: :s:`Whether to include neighbor markers in the mesh image output`
:Type:    :t:`integer`
:Default: :d:`3`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`This parameter is primarily used for debugging.  Internally, each node in the mesh keeps track of the mesh level of its neighbors.  This parameter includes a marker on each face colored according to the neighbor's level.  The value of this parameter specifies the lower limit on the face "rank" (0 for corners, 1 for edges, 2 for faces).  The default of 3 means no markers are displayed.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_size`
:Summary: :s:`Set the size of the image`
:Type:    :t:`list` ( :t:`integer` )
:Default: :d:`[0,0]`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`Specify the size of the output image.  By default it is sized to be one pixel per field value at the finest mesh level.  This is useful to keep images from being to big for large problems, or too small for small problems (e.g. for mesh images which could otherwise be too small).`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_log`
:Summary: :s:`Whether to output the log of the data`
:Type:    :t:`logical`
:Default: :d:`<false>`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`If true, then the natural logarithm of the field value is used for mapping values to the colormap, otherwise use the original field value.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_type`
:Summary: :s:`Type of image to write`
:Type:    :t:`string`
:Default: :d:`"data"`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`This parameter is used to control whether field values are used to generate the image, whether it's an image of the mesh structure, or a combination of both.  Valid values are` :t:`"data"`, :t:`"mesh"`, :e:`or` :t:`"data+mesh"`.

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_block_size`
:Summary: :s:`Number of pixels for fine-level blocks in a mesh image`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`For images of meshes, this parameter defines how many pixels wide each finest-level block is in the image.  This parameter and the image_size parameter should not both be set.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`image_mesh_color`
:Summary: :s:`How to color blocks in a mesh image`
:Type:    :t:`string`
:Default: :d:`"level"`
:Scope:     :c:`Cello`
:Assumes:   :g:`<file_set>` is of :p:`type` :t:`"image"`

:e:`By default, blocks in mesh images are colored according to the level of the block.  In addition to` :t:`"level"`, :e:`other possible ways to assign colors to blocks include` :t:`"process"` :e:`and` :t:`"age"`.

   
schedule
--------

The schedule parameter subgroup defines when to do something, such as
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

       check {

          # **** write a checkpoint every 100.0 seconds ****

          schedule {
             var = "seconds";
             start = 100.0;
             step =  100.0;
          }
           ...
       };

       dump {

          # **** perform a data dump every 50 cycles until cycle 1000 ****

          schedule {
             var = "cycle";
             step =   50;
             stop = 1000;
           }
            ...
       };

       image {

          # **** write an image at times t = 1.0,  2.0, and 5.0 ****

          schedule {
             var = "time";
             list = [1.0, 2.0, 5.0];
           }
            ...
       };
    }
            
----

:Parameter:    :p:`Output` : :g:`<file_set>` : :p:`schedule` : :p:`var`
:Summary: :s:`Variable associated with scheduling for the given file set`
:Type:    :t:`string`
:Default: :d:`"none"`
:Scope:     :c:`Cello`

:e:`The` :p:`var` :e:`parameter specifies what value is checked at each cycle, which may be` :t:`"cycle"`, :t:`"time"`, :e:`or` :t:`"seconds"` :e:`Here "time" refers to simulation time, and "seconds" to wall-clock time.  Note that when simulation "time" is specified, the simulation's time step may be reduced such that the corresponding output occurs exactly at the specified time.`

---- 

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`schedule` : :p:`list`
:Summary: :s:`List of scheduled values for the specified variable`
:Type:    [ :t:`list` ( :t:`integer` ) | :t:`list` ( :t:`float` ) ]
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`This parameter specifies a list of values to check against for output with respect to cycle, time, or seconds.  If the "var" parameter associated with the schedule is "cycle", then` :p:`value` :e:`must be a list of integers; otherwise,` :p:`value` :e:`must be a list of` :t:`float`:e:`'s  The default is an empty list.`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`schedule` : :p:`start`
:Summary: :s:`Starting value for scheduled interval`
:Type:    [ :t:`integer` | :t:`float` ]
:Default: :d:`0 | 0.0`
:Scope:     :c:`Cello`
:Todo:    :o:`write`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`schedule` : :p:`stop`
:Summary: :s:`Last value for scheduled interval`
:Type:    [ :t:`integer` | :t:`float` ]
:Default: :d:`max (integer) | max (double)`
:Scope:     :c:`Cello`
:Todo:    :o:`write`

----

:Parameter:  :p:`Output` : :g:`<file_set>` : :p:`schedule` : :p:`step`
:Summary: :s:`Stepping increment for interval`
:Type:    [ :t:`integer` | :t:`float` ]
:Default: :d:`1 | 1.0`
:Scope:     :c:`Cello`
:Todo:    :o:`write`

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
Block, and ranges [-32768, 16384) and [16384, 32768) on either side
of the associated Block.

Particles are allocated and operated on in "batches".  The
`batch_size` parameter defines how many particles are in a batch.  By
operating on particles in batches, the frequency of memory operations
is greatly reduced, and functions operating on particle attributes can
be more efficient due to reduced overhead.  It should also simplify
writing particle methods to be executed on accelerators, such as
NVIDIA or AMD GPU's.

Just as with fields, particle types can be assigned to groups_.

----

:Parameter:  :p:`Particle` : :p:`list`
:Summary: :s:`List of particle types`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`Cello allows arbitrary parameter types (dark matter particles, tracer particles, star particles, etc.), each with arbitrary attributes (position, velocity, etc.).  The` :p:`list` :e:`parameter defines which types of particles to use.`

  ::

    Particle {

        list = ["dark", "trace"];

    }

----

:Parameter:  :p:`Particle` : :p:`batch_size`
:Summary: :s:`Number of particles in a "batch" of particles`
:Type:    :t:`integer`
:Default: :d:`1024`
:Scope:     :c:`Cello`

:e:`Particles are allocated and operated on in` *batches*.  :e:`The number of particles in a batch is set using the` :p:`batch_size` :e:`parameter.  The default batch size is 1024.`

----

:Parameter:  :p:`Particle` : :g:`particle_type` : :p:`attributes`
:Summary: :s:`List of attribute names and data types`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`none`
:Scope:     :c:`Cello`

:e:`Each particle type can have multiple attributes of varying types, which are defined by the` :p:`attributes` :e:`parameter.  The` :p:`attributes` :e:`parameter is a list of strings, alternating between the name of the parameter, and its type.  Names may include` :t:`"position_x"`, :t:`"velocity_z"`, :t:`"mass"`,
:t:`"id"`, :e:`etc.  Types may include` :t:`"single"`, :t:`"double"`, :t:`"quadruple"`, :t:`"int8"`, :t:`"int16"`, :t:`"int32"`, or :t:`"int64"`.  :e:`Ordering of attributes in memory is as in the` :p:`attributes` :e:`parameter.`

   ::

    Parameter {

        list = ["trace", "dark"];

            trace {

                attributes = ["id", "int64",
                              "x",  "single",
                              "y",  "single",
                              "z",  "single"];
             }

             dark {

                 attributes = ["id",         "int64",
                               "mass",       "double",
                               "velocity_x", "single",
                               "velocity_y", "single",
                               "velocity_z", "single",
                               "position_x", "int16",
                               "position_y", "int16",
                               "position_z", "int16"];
           }
      }

:e:`Note that when attributes of multiple sizes are included in the same parameter type, it can be helpful to order the attributes so that larger-sized attributes are listed first, followed by smaller-sized attributes.  This can help prevent allocating more memory than necessary, since attributes may be padded with unused bytes for correct memory alignment.`

----

:Parameter:  :p:`Particle` : :g:`particle_type` : :p:`interleaved`
:Summary: :s:`Format of output files`
:Type:    :t:`logical`
:Default: :d:`false`
:Scope:     :c:`Cello`

:e:`Particle attributes within a batch of particles may be stored in memory either particle-by-particle, or "interleaved" (attribute-by-attribute).  If` |aij| :e:`represents the jth attribute of particle i, then with` :p:`interleaved = false`, :e:`attributes would be stored as` |a00| ... |am0|, |a01| ... |am1| ... |a0n| ... |amn|. :e:`If, however,` :p:`interleaved = true`, :e:`then attributes would be stored as`   |a00| ... |a0n|, |a10| ... |a1n| ... |am0| ... |amn|. :e:`Non-interleaved particle attributes have array accesses of stride 1 and minimal storage overhead, but may not utilize cache well.  Interleaved particle attributes` *may* :e:`have improved cache utilization, but will have stride > 1, and may require memory padding for correct alignment of attributes in memory.  The default is` :t:`false.`

  

----

:Summary: :s:`Specify a list of groups that the Particle type belongs to`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[ ]`
:Scope:     :c:`Cello`


:e:`Different Particle types may belong to any number of different "groups", which allows simulation code to loop over multiple related particle types.`

  ::

    Particle {
        list = ["trace","dark","star"];

        dark { group_list = ["has_mass"]; }
        star { group_list = ["has_mass"]; }
    }


:e:`This example can be rewritten as follows, which is completely equivalent:`

  ::


    Particle 
        list = ["trace","dark","star"];
    }

    Group {
        list = ["has_mass"];
        has_mass {
           particle_list = ["dark","star"];
        }
    }

----

:Parameter:  :p:`Particle` : :g:`particle_type` : :p:`position`
:Summary: :s:`Format of output files`
:Type:    :t:`string`
:Default: :d:`""`
:Scope:     :c:`Cello`

:e:`Cello needs to know which particle attributes represent position, so that it can determine when particles migrate out of a Block and need to be moved to a neighboring Block.  This is done using the` :p:`position` :e:`parameter:`

  ::

    Particle {

       list = ["trace"];

       trace {

          attributes = ["id",
                        "x","single",
	                "y","single",
	                "z","single"];

          position = ["x","y","z"];
       }
    }

----

:Parameter:  :p:`Particle` : :g:`particle_type` : :p:`velocity`
:Summary: :s:`Format of output files`
:Type:    :t:`string`
:Default: :d:`""`
:Scope:     :c:`Cello`

:e:`Enzo may need to know which particle attributes represent velocity, for example for kick() or drift() operations.  This is done using the` :p:`velocity` :e:`parameter, whose usage is analogous to the` :p:`position` :e:`parameter.  While specifying position is required, specifying velocity is optional.`

  ::

    Particle {

       list = ["dark"];

       trace {

          attributes = [ "x","single",   "y","single",   "z","single",
                        "vx","single",  "vy","single",  "vz","single",
			"mass","single"];

          velocity = ["vx","vy","vz"];
       }
    }

-----------
Performance
-----------

:Parameter:  :p:`Performance` : :p:`name`
:Summary: :s:`Format of output files`
:Type:    :t:`string`
:Default: :d:`""`
:Scope:     :c:`Cello`

:e:`This parameter specifies the format of output files, e.g. "cello-perf-%06d.data".  An` :t:`integer` :e:`specifier must be included for processor number.  Default is "", which means no performance data is written to disk.`

----

:Parameter:  :p:`Performance` : :p:`stride`
:Summary: :s:`Output file processor stride`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :c:`Cello`

:e:`Writing performance data to files can overload parallel file systems if one file is written per process on runs with high processor counts.  This parameter allows for a subset of processors to write to files.  E.g., if Performance:stride=3, then only processors 0, 3, 6, etc. write performance data files.`

----

:Parameter:  :p:`Performance` : :p:`warnings`
:Summary: :s:`Whether to output performance-related warnings`
:Type:    :t:`logical`
:Default: :d:`true`
:Scope:     :c:`Cello`

:e:`If calls to the Performance API are incorrect, e.g. if stop_region() is called on a region that has not been started, then this parameter specifies whether or not to display warning messages`

----

:Parameter:  :p:`Performance` : :p:`papi` : :p:`counters`
:Summary: :s:`List of PAPI counters`
:Type:    :t:`list` ( :t:`string` )
:Default: :d:`[]`
:Scope:     :c:`Cello`

:e:`List of PAPI hardware performance counters to trace, e.g. 'counters = ["PAPI_FP_OPS", "PAPI_L3_TCA"];'.`



-------
Restart
-------

:Parameter:  :p:`Restart` : :p:`file`
:Summary: :s:`Parameter file to read on restart`
:Type:    :t:`string`
:Default: :d:`""`
:Scope:     :c:`Cello`

:e:`This optional parameter is used to specify the name of a parameter file to read on restart.  Its purpose is to allow a simulation to be restarted with slightly different parameter values, e.g. with a smaller courant number.  Currently, very few parameters are supported in the restart parameter file, just` :p:`Field` : :p:`courant` :e:`and` :p:`Testing` : :p:`time_final`.

--------
Stopping
--------

:Parameter:  :p:`Stopping` : :p:`cycle`
:Summary: :s:`Stopping cycle`
:Type:    :t:`integer`
:Default: :d:`max (` :t:`integer` :d:`)`
:Scope:     :c:`Cello`

:e:`Stopping cycle.`

----

:Parameter:  :p:`Stopping` : :p:`time`
:Summary: :s:`Stopping time`
:Type:    :t:`float`
:Default: :d:`max (` :t:`double` :d:`)`
:Scope:     :c:`Cello`

:e:`Stopping time.`

----

:Parameter:  :p:`Stopping` : :p:`seconds`
:Summary: :s:`Stop after this number of seconds (wall-clock time)`
:Type:    :t:`float`
:Default: :d:`max (` :t:`double` :d:`)`
:Scope:     :c:`Cello`

:e:`End the calculation after this many seconds of wall-clock time.`

----

:Parameter:  :p:`Stopping` : :p:`interval`
:Summary: :s:`Stopping interval`
:Type:    :t:`integer`
:Default: :d:`1`
:Scope:     :c:`Cello`

:e:`Number of cycles between applying the stopping criteria.`


-------
Testing
-------

:Parameter:  :p:`Testing` : :p:`cycle_final`
:Summary: :s:`Enzo-P unit test parameter for expected final cycle number`
:Type:    :t:`integer`
:Default: :d:`0`
:Scope:     :c:`Cello`

:e:`Enzo-P unit test parameter for expected final cycle number.`

----

:Parameter:  :p:`Testing` : :p:`time_final`
:Summary: :s:`Enzo-P unit test parameter for expected final time`
:Type:    :t:`float`
:Default: :d:`0.0`
:Scope:     :c:`Cello`

:e:`Enzo-P unit test parameter for expected final time.`

----

:Parameter:  :p:`Testing` : :p:`time_tolerance`
:Summary: :s:`Tolerance on the absolute error between actual final time and time_final`
:Type:    :t:`float`
:Default: :d:`1.0e-6`
:Scope:     :c:`Cello`

:e:`Enzo-P unit test parameter for tolerance on the expected final time.`
