.. include:: ../roles.incl
  
**********************
Flux Correction Design
**********************
.. toctree::
     
============
Requirements
============

Flux correction involves updating field values along faces between
neighboring blocks to ensure that conserved quantities are 
conserved across jumps in spacial and temporal resolution.
The update involves adding correction factors to all coarse
values that lie along coarse-fine interfaces, where the correction
factors are computed using the coarse and fine fluxes across the
interface.  While the basic update is straightforward, care must be
taken to ensure that the correction factors are computed correctly,
especially for non-centered field variables or when adaptive
time stepping is used.  Additional corrections will be required for
MHD, and expansion terms in cosmological problems may also need to be
considered.

Basic operations involved include the following:
   
#. allocating and deallocating flux data
#. setting and accessing flux values
#. communicating fluxes between neighboring blocks (fine-to-coarse)
#. computing correction factors given coarse and fine fluxes
#. correcting coarse-level field values given correction factors

Relative to ENZO's structured AMR grids, flux correction for Enzo-E /
Cello is simpler due to Cello's array-of-octree refinement: blocks
only share fluxes at block faces (Enzo-E blocks do not contain
sub-blocks), and the topology of fine- and coarse-level block
intersections is simpler (Enzo-E coarse-level block faces are adjacent
to exactly four fine-level block faces).  The first simplification
removes the need to loop over sub-grids or store fluxes internal to a
block (i.e. the "projection step" in ENZO), and the second removes the
need to explicitly store loop indices for flux arrays.


We assume fluxes for all conserved fields are provided by the
hydrodynamics solver along all required block faces at each time step.

Method requirements
===================

The Method component is responsible for storing, accessing,
communicating, and operating with fluxes, to implement flux correction
using support provided by Cello flux classes.  Specific requirements
are listed below:

- **RE-1.** Initialize and store fluxes, provided by the hydro solver,
  at each time step for all conserved fields across all coarse-fine
  interfaces
- **RE-2.** Request required fluxes from neighboring Blocks
- **RE-3.** Compute correction factors for all required conserved
  field face values given a block face and fluxes from both the block
  and its face-sharing neighbor
- **RE-4.** Correct conserved field face values given computed
  correction factors
- **RE-5.** Identify which fields require flux correction
- **RE-6.** Support adaptive time stepping
- **RE-7.** Support MHD
- **RE-8.** Ensure conservation is not lost through any other operation, e.g. interpolation

Data requirements
=================

Cello's flux data classes provide the flux-correction method with
sufficient support for storing fluxes, computing flux correction
factors, and applying the flux correction to conserved field values.

Specific requirements include the following:
   
- **RC-1.** Store fluxes of conserved fields that lie along block faces
- **RC-2.** Store associated fluxes computed on adjacent blocks
- **RC-3.** Communicate fluxes between adjacent blocks when (and ideally only when) needed
- **RC-4.** Store the time interval along with each collection of fluxes
- **RC-5.** Allow multiple fluxes to be stored for the same field but
  different time steps
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
implementing flux correction, *MethodFluxCorrect*.

To support coding :p:`MethodFluxCorrect`, our design uses a *FluxData*
class, whose responsibility is to manage all flux data on a block.
This :p:`FluxData` class is analagous to the existing
:p:`FieldData` and :p:`ParticleData` classes.  Like :p:`FieldData`
and :p:`ParticleData` objects, :p:`FluxData` is contained in the
Blocks :p:`Data` object.  Communication is handled by the refresh
mechanism by augmenting the existing communication of field and
particle data between Blocks to include communicating fluxes.

While flux data could be implemented directly as arrays
(e.g. std::vector, :p:`EnzoArray`, etc.) other attributes need to be
associated with each collection of flux data, specifically the Block
the fluxes were computed on, which Block face the fluxes are
associated with, the time interval for the fluxes, etc.  For this we
introduce a lower-level class *FaceFluxes* to store the array of flux
data for an individual face, along with its defining attributes.
The :p:`FaceFluxes` class is in turn implemented usingt a simple *Face*
class that defines the block face on which the fluxes are defined.

Interfaces
==========

We develop the interfaces for the flux-related classes below, starting
with the top-level :p:`MethodFluxCorrect`, then the progressively
higher-level :p:`Face`, :p:`FaceFluxes`, and :p:`FluxData` classes.


---------------------------
MethodFluxCorrect class
---------------------------

The :p:`MethodFluxCorrect` class is a Cello :p:`Method`, whose
main virtual method is :p:`compute(Block)`.  This operates on some
subset of data types on a Block.  The :p:`MethodFluxCorrect` method
initiates communication 
between processes, and computes and applies appropriate
flux-correction operations on required Field values along block
interfaces.

Since the :p:`MethodFluxCorrect` class is inherited from the Cello
:p:`Method` class, the public interface for this class is already
prescribed.  The two methods in the interface are the constructor
:p:`MethodFluxCorrect()` used to initialize the required
communication, and :p:`MethodFluxCorrect::compute(Block)` which
implements flux correction on a given Block.  (The other virtual
function in the :p:`Method` interface is :p:`timestep()`, which is not
required for flux-correction.)

----

.. glossary::

   ``MethodFluxCorrect::MethodFluxCorrect()``
   
      *Create a new MethodFluxCorrect object, and define its refresh communication requirements*
  
----

.. glossary::

   ``virtual void MethodFluxCorrect::compute (Block * block)``
   
      *Request Cello to refresh its flux data, then apply flux correction*

      * **block**: *Block that flux correction is being applied to*
   
----------
Face class
----------

A block :p:`Face` is any facet, edge, or corner of a block.  For flux
correction in hydrodynamics, one typically only deals with facets, or
(d-1)-dimensional faces; however, for MHD, edge faces may be used as
well.

A face is determined by its "center" (ix,iy,iz), -1 <= ix,iy,iz <= +1,
assuming the block corners are at (+/-1,+/-1,+/-1).  As
examples, the positive Y-axis facet is (0,+1,0), the edge along X=+1
and Z=-1 is (+1,0,-1), and the entire block as a 3-D "face" is
(0,0,0).

Each Face is associated with a normal vector defining
the direction of the fluxes.  This direction is given by the axis (x=0,
y=1, z=2), and face (lower=0, upper=1).

----

.. glossary::

   ``Face::Face (int ix, int iy, int iz,  int axis, int face)``
    *Create a Face object for a Block associated with face (ix,iy,iz),  -1 <= ix,iy,iz <= +1, with fluxes in the direction (axis, face), 0 <= axis < rank, 0 <= face <= 1.*

----

.. glossary::

   ``void Face::get_face (int *ix, int *iy, int *iz)``
    *Return the tuple (ix,iy,iz), -1 <= ix,iy,iz <= 1, identifying the block's face, which may be a corner, edge, facet, or the entire block.*

----

.. glossary::

   ``int Face::axis()``
    *Return the axis associated with the normal direction: x=0, y=1, or z=2.*

----

.. glossary::

   ``int Face::face()``
    *Return whether the normal direction is towards the lower (0) or upper (1) face direction.*

----------------     
FaceFluxes class
----------------     

Face fluxes represent an array of fluxes of a given conserved Field
through a Block's face or subset of a face.  Operations available for
fluxes including *coarsening*, for summing fluxes in a finer block to
match the resolution of a neighboring coarser block, and *summing*,
for accumulating a sequence of fluxes associated with a block with a
finer time step to match a neighboring block's coarser time step.

----

..
   comment # [*]

.. glossary::
   
   ``FaceFluxes::FaceFluxes (Face face, int index_field, int nx, int ny, int nz, int cx, int cy, int cz)``
    *Create a FaceFluxes object for the given face, field, and block size. Optionally include centering adjustment (0 <= cx,cy,cz <= 1) for facet-, edge-, or corner-located field values*
   
----

.. glossary::

   ``void FaceFluxes::allocate()``
    *Allocate the flux array and initialize values to 0.0.*

----

.. glossary::

  ``void FaceFluxes::deallocate()``
   *Deallocate the flux array.*

----

.. glossary::

  ``void FaceFluxes::clear()``
   *Set flux array values to 0.0.*

----

.. glossary::

  ``Face FaceFluxes::face()``
   *Return the face associated with the FaceFluxes.*

----

.. glossary::

  ``int FaceFluxes::index_field()``
   *Return the associated field index.*


----

.. glossary::

  ``void FaceFluxes::get_size (int * mx, int * my, int * mz)``
   *Return the array dimensions of the flux array, including any adjustments for centering.  Indexing is ix + mx\*(iy +my\*iz).*
   
----

.. glossary::

  ``void FaceFluxes::set_flux_array ( std::vector<double> array, int dx, int dy, int dz)``
   *Copy flux values from an array to the FluxFaces flux array. Array element array[ix\*dx + iy\*dy + iz\*dz] should correspond to flux value (ix,iy,iz), where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).*
  
----

.. glossary::

  ``std::vector<double> & FaceFluxes::flux_array (int * dx=0, int * dy=0, int * dz=0)``
   *Return the array of fluxes and associated strides (dx,dy,dz) such that the (ix,iy,iz) flux value is fluxes[ix\*dx + iy\*dy + iz\*dz], where (0,0,0) <= (ix,iy,iz) < (mx,my,mz).*
  
----

.. glossary::

  ``void FaceFluxes::coarsen(int cx, int cy, int cz, int rank)``
   *Used for coarsening fine-level fluxes to match coarse level fluxes. Arguments (cx,cy,cz) specify the child indices of the block within its parent (not to be confused with centering (cx_,cy_,cz_); flux array size is kept the same, with offset determined by child indices.*
  
----

.. glossary::

  ``void FaceFluxes::accumulate (FaceFluxes & ff, int cx, int cy, int cz, rank)``
   *Add the FaceFluxes object to this one. Used for accumulating fluxes with finer time steps until they match the coarser time step. Assumes spacially-conforming FaceFluxes objects.*
  
----

.. glossary::

  ``FaceFluxes & FaceFluxes::operator *= (double weight)``
   *Scale the fluxes array by a scalar constant.*
  
--------------
FluxData class
--------------

The :p:`FluxData` class defines a collection of all :p:`FluxFaces`
required by a Block to perform flux corrections.  This includes all
flux arrays on :p:`Faces` whose neighboring Block differs in either
mesh refinement level or time step.  :p:`FluxFaces` for faces that
require flux-correction come in conforming pairs, one set of fluxes
corresponding to the Block, and one corresponding to the block's
neighbors. Flux arrays for the neighboring block are received in the
flux refresh operation.  Support for coarsening, adding, and
differencing fluxes is the responsibility of the :p:`FaceFluxes`
class; :p:`FluxData` is primarily a container.


.. glossary::      

   ``FluxData::FluxData()``
     *Create an empty FluxData() object*

----

.. glossary::

   ``void FluxData::allocate(int nx, int ny, int nz,  std::vector<int> field_list, std::vector<int> * cx_list=nullptr, std::vector<int> * cy_list=nullptr, std::vector<int> * cz_list=nullptr)``
    *Allocate all flux arrays for each field in the list of field indices.  Optional arrays to indicate the centering of fields may also be provided.*

----

.. glossary::

  ``void FluxData::deallocate()``
   *Deallocate all face fluxes for all faces and all fields.*

----

.. glossary::

  ``int FluxData::num_fields()``
   *Return the number of field indices.*

----

.. glossary::

  ``int FluxData::index_field(int i_f)``
   *Return the i'th field index.*

----

.. glossary::      

   ``void FluxData::block_fluxes(int axis, int face, int i_f)``
    *Return the face fluxes object associated with the given facet and field. Note 0 <= i_f  < num_fields() is an index into the field_list vector, not the field index itself.*

----

.. glossary::      

   ``void FluxData::neighbor_fluxes(int axis, int face, int i_f)``
    *Return the neighboring block's face fluxes associated with the given facet and field. Note 0 <= i_f  < num_fields() is an index into the field_list vector, not the field index itself.*

----

.. glossary::      

   ``void FluxData::set_block_fluxes(FaceFluxes * ff, int axis, int face, int i_f)``
    *Set the block's face fluxes associated with the given facet and field. Note 0 <= i_f  < num_fields() is an index into the field_list vector, not the field index itself.*

----

.. glossary::      

   ``void FluxData::set_neighbor_fluxes(FaceFluxes * ff, int axis, int face, int i_f)``
    *Set the neighboring block's face fluxes associated with the given facet and field. Note 0 <= i_f  < num_fields() is an index into the field_list vector, not the field index itself.*

----

.. glossary::      

   ``void FluxData::sum_neighbor_fluxes(FaceFluxes * ff, int axis, int face, int i_f)``
    *Accumulate (sum) the neighboring block's face fluxes associated with the given facet and field. Note 0 <= i_f  < num_fields() is an index into the field_list vector, not the field index itself.*

Flux Communication
==================

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

---- 

.. figure:: flux-refresh-1.png
    :width: 500px
    :align: center
    :alt: figure illustrating that fine block fluxes are coarsened before being communicated
    :figclass: align-center

    Communicating fluxes assuming constant time steps.               
      
---- 

.. figure:: flux-refresh-2.png
    :width: 500px
    :align: center
    :alt:  figure illustrating that fine block fluxes with shorter time steps are accumulated in the receiving coarse grid neighbor fluxes
    :figclass: align-center

    Communicating fluxes assuming adaptive time steps.               
      
Testing
=======

Multiple levels of testing are used, including unit tests for the lower level
:p:`Face`, :p:`FaceFluxes`, and :p:`FluxData` classes, and 
application testing for :p:`MethodFluxCorrect`.

Application tests include varying difficulties of meshes, physics,
boundary conditions, and floating-point precision.  Different levels
are summarized below:

*  **Mesh**

   * **M0** single block
   * **M1** unigrid
   * **M2** one additional refinement level
   * **M3** two additional refinement levels
     
* **Physics**

  * **PH** hydro
  * **PG** gravity
  * **PP** gravitating particles
  * **PC** cosmological expansion

* **Boundary**

  * **BP** periodic
  * **BR** reflecting

* **Floating-point precision**

  * **F4** single
  * **F8** double

* **Parallelism**
  
  * **p1**: single processor
  * **pc**: parallel across cores
  * **pn**: parallel across nodes

Documentation
=============

- **DD-1.** *Design*: add flux correction design to :p:`design/design-flux.rst` 
- **DD-2.** *Method*: add :p:`MethodFluxCorrect` method documentation to :p:`user/problem_method.rst`
- **DD-3.** *Testing*: add :p:`testing/testing_flux.rst` test documentation
- **DD-4.** *Parameters*: update :p:`doc/source/param/`  parameter documentation
   
==========
Milestones
==========

- **M-1.** :p:`MethodFluxCorrect` demonstrated  working for hydrodynamics
- **M-2.** :p:`MethodFluxCorrect` demonstrated  working for MHD

=====
Tasks
=====

- **T-1.** :p:`FaceFluxes` class design and implementation
- **T-2.** :p:`FluxData` class design and implementation
- **T-3.** :p:`MethodFluxCorrect` class design and implementation
