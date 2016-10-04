// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    [\ref Mesh] Declaration of the Hierarchy class

#ifndef MESH_HIERARCHY_HPP
#define MESH_HIERARCHY_HPP

class Factory;
class Simulation;

class CProxy_Block;

class Hierarchy {

  /// @class    Hierarchy
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Adaptive mesh refinement hierarchy

  friend class IoHierarchy;

public: // interface

  /// Empty constructor for Charm++ pup()
  Hierarchy() throw()
  : factory_(NULL),
    rank_(0),
    refinement_(0),
    max_level_(0),
    num_blocks_(0), 
    num_particles_(0), 
    num_zones_total_(0), 
    num_zones_real_(0), 
    block_array_(NULL),
    block_exists_(false)
  {
    for (int axis=0; axis<3; axis++) {
      root_size_[axis] = 0;
      lower_[axis] = 0.0;
      upper_[axis] = 0.0;
      blocking_[axis] = 0.0;
      periodicity_[axis][0] = false;
      periodicity_[axis][1] = false;
    }
  }
  
  /// Initialize a Hierarchy object
  Hierarchy (const Factory * factory, 
	     int rank,
	     int refinement,
	     int max_level) throw ();

  /// Delete the Hierarchy object
  virtual ~Hierarchy() throw ();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  //----------------------------------------------------------------------

  /// Set domain lower extent
  void set_lower(double x, double y, double z) throw ();

  /// Set domain upper extent
  void set_upper(double x, double y, double z) throw ();
  
  /// Set root-level grid size
  void set_root_size(int nx, int ny, int nz) throw ();

  /// Set root-level grid size
  void set_blocking(int nbx, int nby, int nbz) throw ();

  //----------------------------------------------------------------------

  /// Return the maximum refinement level (0 for unigrid)
  int max_level() const
  { return max_level_; }
  
  /// Return rank
  int rank() const throw ();

  /// Return domain lower extent
  void lower(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return root-level grid size
  void root_size(int * nx, int * ny = 0, int * nz = 0) const throw ();

  /// Set the periodicity of boundary conditions for domain faces
  inline void set_periodicity (int pxm,   int pxp, 
			       int pym=0, int pyp=0, 
			       int pzm=0, int pzp=0)
  {
    periodicity_[0][0] = pxm;
    periodicity_[0][1] = pxp;
    periodicity_[1][0] = pym;
    periodicity_[1][1] = pyp;
    periodicity_[2][0] = pzm;
    periodicity_[2][1] = pzp;
  }

  /// Return the periodicity of the boundary conditions for domain faces
  void periodicity (int * pxm,   int * pxp, 
		    int * pym=0, int * pyp=0, 
		    int * pzm=0, int * pzp=0)
  {
    if (pxm) (*pxm) = periodicity_[0][0];
    if (pxp) (*pxp) = periodicity_[0][1];
    if (pym) (*pym) = periodicity_[1][0];
    if (pyp) (*pyp) = periodicity_[1][1];
    if (pzm) (*pzm) = periodicity_[2][0];
    if (pzp) (*pzp) = periodicity_[2][1];
  }

  //----------------------------------------------------------------------

  /// Return whether Blocks have been allocated or not
  bool blocks_allocated() const throw()
  { 
    return block_exists_;
  }

  /// Return the number of Blocks
  size_t num_blocks(int * nbx, 
		    int * nby = 0,
		    int * nbz = 0) const throw();

  /// Deallocate local Blocks
  void deallocate_blocks() throw();

  /// Return pointer to the Block CHARM++ chare array
  CProxy_Block * block_array() const throw()
  { return block_array_;}

  /// Increment (decrement) number of mesh blocks
  void increment_block_count(int count)
  { num_blocks_ += count; }

  /// Increment (decrement) number of particles
  void increment_particle_count(int64_t count)
  { num_particles_ += count; }

  /// Increment (decrement) number of real_zones
  void increment_real_zone_count(int64_t count)
  { num_zones_real_ += count; }

  /// Increment (decrement) number of total_zones
  void increment_total_zone_count(int64_t count)
  { num_zones_total_ += count; }

  /// Return the number of blocks on this process
  size_t num_blocks() const throw()
  {  return num_blocks_;  }

  /// Return the number of particles on this process
  int64_t num_particles() const throw()
  {  return num_particles_;  }

  /// Return the number of real zones on this process
  int64_t num_zones_real() const throw()
  {  return num_zones_real_;  }

  /// Return the number of total zones on this process
  int64_t num_zones_total() const throw()
  {  return num_zones_total_;  }


  void create_forest (FieldDescr   * field_descr,
		      bool allocate_data,
		      bool testing          = false) throw();

  void create_subforest (FieldDescr   * field_descr,
			 bool allocate_data,
			 int min_level,
			 bool testing          = false) throw();


  /// Return the number of root-level Blocks along each rank
  void blocking (int * nbx, int * nby=0, int * nbz=0) const throw();

  /// Return the factory object associated with the Hierarchy
  const Factory * factory () const throw()
  { return factory_; }

protected: // attributes

  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  Factory * factory_;

  /// Rank of the hierarchy [ used for Charm++ pup() of Tree ]
  int rank_;
  
  /// Refinement of the hierarchy [ used for Charm++ pup() of Tree ]
  int refinement_;

  /// Maximum mesh level
  int max_level_;

  /// Maximum number of refinement levels

  /// Current number of blocks on this process
  int num_blocks_; 

  /// Current number of particles on this process
  int64_t num_particles_;

  /// Current number of total_zones on this process
  int64_t num_zones_total_; 

  /// Current number of real_zones on this process
  int64_t num_zones_real_; 

  /// Array of Blocks 
  CProxy_Block * block_array_;
  bool           block_exists_;

  /// Size of the root grid
  int root_size_[3];

  /// Lower extent of the hierarchy
  double lower_[3];

  /// Upper extent of the hierarchy
  double upper_[3];

  /// How the Forest is distributed into Blocks
  int blocking_[3];

  /// Periodicity of boundary conditions on faces
  bool periodicity_[3][2];
  
};


#endif /* MESH_HIERARCHY_HPP */

