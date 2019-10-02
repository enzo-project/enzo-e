.. include:: ../roles.incl

**********************
Flux Correction Design
**********************

============
Requirements
============

Flux correction involves updating field values along faces between
neighboring blocks to ensure that conserved quantities are indeed
conserved.  This correction step isn't required at all block faces: if
both mesh refinement levels and time-steps are the same between two
adjacent blocks, then fluxes along the interfaces are expected to
already match to machine precision.  However, if neighboring blocks
have different spacial or temporal resolution, the computed solution
will not in general be conservative, so some type of correction due to
the inconsistent fluxes will be required.

.. note:: We will henceforth use the term "resolution jump" to include
          a difference in time steps as well as in mesh resolution,
          "coarse" and "fine" to include relatively long and short
          time steps as well as relativly large and small block cell
          widths.

The basic update involves adding correction factors to all coarse
values that lie along coarse-fine interfaces, where the correction
factor is computed using the coarse and fine fluxes across the
interface.  While the basic update is straightforward, care must be
taken to ensure that the correction factors are computed correctly,
especially for non-centered field variables or when adaptive
time-stepping is used.  Additional corrections will be required for
MHD, and expansion terms in cosmological problems may also need to be
considered.

Basic operations involved include the following:
   
#. allocating and deallocating flux data
#. setting and accessing flux values
#. communicating fluxes between neighboring blocks (fine-to-coarse)
#. computing correction factors given coarse and fine fluxes
#. correcting coarse-level field values given correction factors

.. note:: Relative to ENZO's structured AMR grids, flux correction for
          Enzo-E / Cello is simpler due to Cello's array-of-octree
          refinement: blocks only share fluxes at block faces (Enzo-E
          blocks do not contain sub-blocks), and the topology of fine-
          and coarse-level block intersections is simpler (Enzo-E
          coarse-level block faces are adjacent to exactly four
          fine-level block faces).  The first simplification removes
          the need to loop over sub-grids or store fluxes internal to
          a block (i.e. the "projection step" in ENZO), and the second
          simplifies loop indexing.

As a precondition we assume fluxes for all conserved fields are
provided by the hydrodynamics solver on all required block faces at
each time-step.

Below we discuss the software requirements specific to the Enzo-E and
Cello layers.

.. note:: As this is a working document, specific requirements are
          expected to be added, modified, or shifted between layers as
          needed.
   
     
Enzo-E requirements
===================

Enzo-E is responsible for storing, accessing, communicating, and
operating with fluxes, to implement flux correction, with support
provided by Cello.  Specific requirements are listed below:

- **RE-1.** Initialize and store fluxes, provided by the hydro solver,
  at each time-step for all conserved fields across all Block faces
  that lie along a resolution jump
- **RE-2.** Request required fluxes from neighboring Blocks
- **RE-3.** Compute correction factors for all required conserved field
  face values given block face and neighbor block fluxes
- **RE-4.** Correct conserved field face values given computed
  correction factors
- **RE-5.** Identify which fields require flux correction
- **RE-6.** Support adaptive time-stepping
- **RE-7.** Support MHD
- **RE-8.** Ensure conservation is not lost through any other operation, e.g. interpolation

Cello requirements
==================

Cello's primary responsibilities are to provide Enzo-E with sufficient
support for storing fluxes, computing flux correction factors, and
applying the flux correction to conserved field values.  Technical
details, especially those involving parallel computing and data
communication, should be isolated from the Enzo-E software layer as
much as possible.

Specific requirements include the following:
   
- **RC-1.** Store fluxes of conserved fields that lie along block faces
- **RC-2.** Store associated fluxes computed on adjacent blocks
- **RC-3.** Communicate fluxes between adjacent blocks when (and ideally only when) needed
- **RC-4.** Store the time interval along with each collection of fluxes
- **RC-5.** Allow multiple fluxes to be stored for the same field but
  different time-steps
- **RC-6.** Provide support for a block to compute correction factors
  along a block face given the block's fluxes and the corresponding
  neighboring block fluxes
- **RC-7.** Provide support for correcting field values along a block
  face given the computed correction factors
- **RC-8.** Allow for dynamic allocation and deallocation of fluxes (optional)
- **RC-9.** Ensure conservative inter-grid interpolation and coarsening

======
Design
======

Our design is developed top-down, starting with a Cello Method for
implementing flux correction in Enzo-E, *EnzoMethodFluxCorrect*.
Using a Cello Method to implement flux correction is a natural
approach, since Methods in Enzo-E are analagous to steps in the
high-level :p:`EvolveLevel()` function in ENZO, and flux correction is
one such step.  (Specifically, flux correction is performed in
ENZO---along with communication and a projection substep---in
:p:`UpdateFromFinerGrids()`, which is called by :p:`EvolveLevel()`).

To support coding :p:`EnzoMethodFluxCorrect`, our design introduces a
Cello class *FluxData*, whose responsibility is to collect together
all flux data required for flux correction on a block into a single
object.  This :p:`FluxData` class will be analagous to the existing
:p:`FieldData` and :p:`ParticleData` classes, which contain all field
data and particle data associated with a Block.  Like :p:`FieldData`
and :p:`ParticleData` objects, :p:`FluxData` will be contained in the
Blocks :p:`Data` object.  Communication will be handled by the refresh
mechanism by augmenting the existing communication of field and
particle data between Blocks to include communicating fluxes.

While flux data could be implemented directly as arrays
(e.g. std::vector, :p:`EnzoArray`, etc.) other attributes need to be
associated with each collection of flux data, specifically the Block
the fluxes were computed on, which Block face the fluxes are
associated with, the time interval for the fluxes, etc.  Thus, we
introduce a lower-level class *FaceFluxes* to store the array of flux
data for an individual face, along with its defining attributes.

We develop the interfaces for these classes below, including prototype
implementations of the higher-level classes to aid in developing the
design of the lower-level classes.


EnzoMethodFluxCorrect class design
==================================

The :p:`EnzoMethodFluxCorrect` class is a Cello :p:`Method`, whose
main virtual method is :p:`compute(Block)`.  This operates on some
subset of Fields and Particle types on a Block, and will include new
operations on fluxes as well.  The :p:`EnzoMethodFluxCorrect` method
communicates fluxes between neighboring blocks to ensure consistency
between processes, and computes and applies appropriate
flux-correction operations to required Fields values along block
interfaces.

-------------------------------
EnzoMethodFluxCorrect interface
-------------------------------

Since the :p:`EnzoMethodFluxCorrect` class is inherited from the Cello
:p:`Method` class, the public interface for this class is already
prescribed.  The two methods in the interface are the constructor
:p:`EnzoMethodFluxCorrect()` used to initialize the required
communication, and :p:`EnzoMethodFluxCorrect::compute(Block)` which
implements flux correction on a given Block.  (The other virtual
function in the :p:`Method` interface is :p:`timestep()`, which is not
required for flux-correction.)

----

.. glossary::

   ``EnzoMethodFluxCorrect()``
   
      *Create a new EnzoMethodFluxCorrect object, and define its refresh communication requirements*
  
----

.. glossary::

   ``virtual void compute (Block * block)``
   
      *Request Cello to refresh its flux data, then apply flux correction*

      * **block**: *Block that flux correction is being applied to*
   
----

------------------------------------
EnzoMethodFluxCorrect implementation
------------------------------------

Below is prototype code for various steps of flux-correction needed,
including initializing for the communication step to refresh fluxes
between adjacent blocks, performing the refresh communication, and
performing the flux correction given the refreshed fluxes.

register fluxes refresh
-----------------------

While most Methods require Field and/or Particle data to be consistent
between neighboring blocks, the flux correction step requires fluxes
to be consistent.  Refresh phases in Methods are typically registered
in the Method's constructor, so we register a refresh phase with all
:p:`FluxData` block data in :p:`EnzoMethodFluxCorrect()`.  Although
code will need to be updated in the Cello layer to support refreshing
flux data between blocks, the prototype code below should be
relatively complete.

.. code-block:: C++

   EnzoMethodFluxCorrect::EnzoMethodFluxCorrect()
   {
      // Register a refresh phase to update fluxes
      
      ir_flux_ = add_new_refresh_();
      Refresh & refresh = new_refresh(ir_flux_);
      refresh.add_all_fluxes ();

      refresh.set_callback(CkIndex_EnzoBlock::p_method_flux_correct());
   }
   

perform fluxes refresh
----------------------

Performing the communication step to refresh fluxes between adjacent
blocks is analagous to similar code in other Enzo-E Methods used to
refresh field and particle data along faces.  Again, while code will
need to be updated in the Cello layer to support refreshing flux data
between blocks, the prototype code below should be relatively
complete.

.. code-block:: C++

   void EnzoMethodFluxCorrect::compute(Block * block)
   {
   
      // Start refresh phase to update fluxes
      
      block->new_refresh(ir_flux_).set_active(is_finest_(block));

      block->new_refresh_start
           (ir_flux_, CkIndex_EnzoBlock::p_method_flux_correct());
   }                           

   void EnzoBlock::p_method_flux_correct()
   {  static_cast<EnzoMethodFluxCorrect*> (solver())->update_field_faces(this);   }  

   void EnzoMethodFluxCorrect::update_field_faces(Block* block) throw()
   {
       // continue with updated fluxes

       // ...
   }


perform flux correction
-----------------------

.. code-block:: C++

   void EnzoMethodFluxCorrect::update_field_faces(EnzoBlock* block) throw()
   {
       // Get this block's Field and FluxData objects

       Index index        = block->index();
       FluxData flux_data = block->data()->flux_data();
       
       // Loop over all fields in the "conservative" field group

       for (int index_field = ...) {

          // Loop over all faces of this block
          
          for (int axis=0; axis<rank; axis++) {
             for (int face=0; face<2; face++) {

                // get "coarse" fluxes
                   
                FaceFluxes fluxes_coarse = flux_data.face_fluxes
                    (index_field,axis,face,cycle_start,cycle_stop,
                    flux_type_mine);

                if (face_requires_correction_(block,axis,face)) {

                   // Loop over all adjacent fine-level faces
                   
                   for (std::pair child = ...) {

                      // Loop over fine-grid time steps

                      FaceFluxes fluxes_fine;
                      for (int cycle=cycle_begin; cycle<cycle_end; cycle++) {
                      
                         fluxes_fine += flux_data.face_fluxes
                            (index_field,axis,face,cycle_start,cycle_stop,
                             flux_type_neighbor,child);
                      }

                      // Compute and apply flux correction on a face

                      update_field_face_
                         (Block, index_field, axis,face,child,
                          fluxes_coarse,fluxes_fine);
                }
             }
          }
       }
   }

.. note:: When adaptive time-stepping is used, multiple time intervals
          of "fine-grid" fluxes will be required for a given face.  To
          avoid issues with comparing floating-point values, we assume
          an integer-valued "finest-grid" cycle number is available.
          However, until adaptive time-stepping is implemented, this
          reverts to the block's current cycle for non-adaptive
          (constant over all blocks per cycle) time-stepping.

.. code-block:: C++

   EnzoMethodFluxCorrect::update_field_face_
      (Block, index_field, axis,face,child,
       fluxes_coarse,fluxes_fine_sum);
   {
   
      // "attach" face_values to the Field's face (like subview of multiarray)
      
      FieldFace field_face = FieldFace(field,axis,face,child);

      int i_start[2],i_stop[2];
      fluxes_coarse.get_loop_limits(i_start,i_stop);
      
      for (int i1=i_start[1]; i1<i_stop[1]; i1++) {
         for (int i0=i_start[0]; i0<i_stop[0]; i0++) {

            // Apply flux correction at (i0,i1)

            field_face[i0][i1] += (fluxes_fine[i0][i1] - fluxes_coarse[i0][i1]);
         }
      }
   }

.. note:: The update is over-simplified, e.g. the field will have to
          be multiplied by density to get conservative form, fine-grid
          fluxes need to be weighted by relative volume, etc.
   
FluxData class design
=====================

The :p:`FluxData` class is used to store all flux data required by
a Block at any given time.

------------------
FluxData interface
------------------

Below we develop the interface for the :p:`FluxData` object for storing
and manipulating fluxes associated with a Block.

----

.. glossary::      

   ``FaceFluxes(index_block, index_field, axis, face, child, cycle, dt)``

      * **index_block**: *Index of the Block that fluxes are associated with*
      * **index_field**: *Index of the Field that fluxes are associated with*
      * **axis**: *Axis (0,1,2) of the face the fluxes are associated with*
      * **face**: *Lower (0) or upper (1) face the fluxes are associated with*
      * **child**: *Tuple identifying which quadrant of the face neighboring fine-block fluxes are associated with*
      * **cycle**: *Fine-level computational cycle of the fluxes*
      * **dt**: *Time step associated with the fluxes*
   
     *Return an alias to the fluxes for the given field and face.  This    should be usable as an lvalue so that updates to the fluxes will    be reflected in the values stored in the :p:`FluxData` object.*

----

.. glossary::

   ``~FluxData ()``

      *Delete the FluxData object*

----

.. glossary::      

   ``FaceFluxes & fluxes (index_block, index_field, axis, face, child, cycle)``

     *Return an alias to the fluxes for the given field and face.  This    should be usable as an lvalue so that updates to the fluxes will    be reflected in the values stored in the :p:`FluxData` object.*

     arguments:
     
      * **index_block**: *Index of the Block associated with the fluxes*
      * **index_field**: *Index of the Field associated with the fluxes*
      * **axis,face**: *Face of the Block associated with the fluxes, identified by axis (0,1,2) and face (0,1)*
      * **child**: *Tuple identifying which quadrant of the face neighboring fine-block fluxes are associated with*
      * **cycle**: *Fine-level computational cycle of the fluxes*
   
----

.. glossary::

  ``operator += (face_fluxes)``

     *Add the FaceFluxes to the existing FaceFluxes object.  Copies if uninitialized.  Cycle numbers of appended* face_fluxes *objects must be consecutive.  The resulting FaceFluxes object will internally store start and stop cycles, and flux values will be summed.*

     * **face_fluxes**: *FaceFluxes object being appended*
     
----

-----------------------
FluxData implementation
-----------------------

FaceFluxes class design
=======================

Create a new uninitialized FaceFluxes object

--------------------
FaceFluxes interface
--------------------

.. glossary::

   ``void FaceFluxes(index_block, index_field, axis, face, child, cycle, dt)``

      *Create a FaceFluxes object for the given Block, Field, axis, and face.*
   
      * **index_block**: *Index of the Block associated with the fluxes*
      * **index_field**: *Index of the Field associated with the fluxes*
      * **axis,face**: *Face of the Block associated with the fluxes, identified by axis (0,1,2) and face (0,1)*
      * **child**: *Tuple identifying which quadrant of the face neighboring fine-block fluxes are associated with*
      * **cycle**: *Fine-level computational cycle of the fluxes*
      * **dt**: time step for the fluxes
  
-------------------------
FaceFluxes implementation
-------------------------

   * **index_block**: index of the Block associated with the fluxes
   * **index_field**: index of the Field associated with the fluxes
   * **axis,face**: face of the Block, identified by axis (0,1,2) and face (0,1)
   * **child**: *Tuple identifying which quadrant of the face neighboring fine-block fluxes are associated with*
   * **cycle_begin**: starting cycle associated with fluxes
   * **cycle_end**: ending cycle associated with fluxes
   * **dt**: time step for the fluxes

Communication
=============

Communicating fluxes between adjacent blocks is similar to communicating
field face data and migrating particles, so it makes sense to
augment the existing refresh mechanism for communicating :p:`FieldData`
and :p:`ParticleData` to also refresh :p:`FluxData`.  We briefly
outline a couple possible differences below.

First, flux data is generally only communicated from coarse blocks to
neighboring fine blocks, and from blocks with a smaller time step to
neighboring blocks with larger time steps.  Thus an additional
communication pattern and associated synchronization could be
introduced.  Examples of existing synchronization patterns are
:p:`sync_neighbor` type for communicating between all pairs of
adjacent leaf blocks, and :p:`sync_level` between all adjacent blocks in the
same level.  For fluxes, we could introduce a :p:`sync_fluxes` synchronization
type, which by definition includes

   1. leaf blocks only
   2. *usually* only d-1 faces (may need block edge- or
      corner-adjacent flux data for MHD)
   3. only from finer resolution to coarser resolution (temporal as
      well as spacial)

.. figure:: flux-refresh-1.png
    :width: 400px
    :align: center
    :alt: figure illustrating that fine block fluxes are coarsened before being communicated
    :figclass: align-center

    Communicating fluxes assuming constant time steps.               
      
.. figure:: flux-refresh-2.png
    :width: 400px
    :align: center
    :alt:  figure illustrating that fine block fluxes with shorter time steps are accumulated in the receiving coarse grid neighbor fluxes
    :figclass: align-center

    Communicating fluxes assuming adaptive time steps.               
      
Testing
=======

In Enzo-E / Cello we strive to use test-driven development, which
helps keep development cycles short and helps keep the code base from
accumulating untested code.

We will include unit tests for the lower level :p:`FaceFluxes` and
:p:`FluxData` classes, and use application testing for
:p:`EnzoMethodFluxCorrect`.

Hydrodynamics application tests will include increasingly difficult
test problems, with the easiest being a 2D two-level AMR (L=2)
single-processor (P=1) single-cycle (C=1) problem.  This should be
simple enough to allow a debuggable problem, where the count of
flux-correction steps is finite.  A similar progression of
increasingly difficult problems will be used for MHD.  The "hardest"
problem should be many levels, many processes (multi-node at the
least) and long-duration.

Similar series of tests including self-gravity and cosmolgy will be
developed.

- **DT-1.** :p:`FaceFluxes` unit tests
- **DT-2.** :p:`FluxData` unit tests
- **DT-H.** Hydrodynamics application tests
- **DT-H1.** 2D L=2 P=1 C=1 hydro (easiest hydro)
- **DT-H2.** 3D L=2 P=1 C=1 hydro
- **DT-H3.** 3D L>2 P=1 C=1 hydro
- **DT-H4.** 3D L>2 P>1 C=1 hydro
- **DT-H5.** 3D L>2 P>1 C>1 hydro (hardest hydro)
- **DT-G.** Self-gravity application tests
- **DT-C.** Cosmology application tests

Documentation
=============

- **DD-1.** *Design*: add flux correction design to :p:`design/design-flux.rst` 
- **DD-2.** *Method*: add :p:`EnzoMethodFluxCorrect` method documentation to :p:`user/problem_method.rst`
- **DD-3.** *Testing*: add :p:`testing/testing_flux.rst` test documentation
   
==========
Milestones
==========

- **M-1.** :p:`EnzoMethodFluxCorrect` demonstrated  working for hydrodynamics
- **M-2.** :p:`EnzoMethodFluxCorrect` demonstrated  working for MHD

=====
Tasks
=====

- **T-1.** :p:`FaceFluxes` class design and implementation
- **T-2.** :p:`FluxData` class design and implementation
- **T-3.** :p:`EnzoMethodFluxCorrect` class design and implementation

=====
Notes
=====
