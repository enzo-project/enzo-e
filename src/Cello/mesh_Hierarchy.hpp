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
    refinement_(0),
    min_level_(0),
    max_level_(0),
    num_blocks_(0),
    num_blocks_level_(),
    block_vec_(),
    num_particles_(0), 
    num_zones_total_(0), 
    num_zones_real_(0), 
    block_array_(),
    block_exists_(false)
  {
    for (int axis=0; axis<3; axis++) {
      root_size_[axis] = 0;
      lower_[axis] = 0.0;
      upper_[axis] = 0.0;
      blocking_[axis] = 0;
      periodicity_[axis] = 0;
    }
  }
  
  /// Initialize a Hierarchy object
  Hierarchy (const Factory * factory, 
	     int refinement,
	     int min_level,
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

  /// Return the minimum refinement level (0 for unigrid)
  int min_level() const
  { return min_level_; }

  /// Return the maximum refinement level (0 for unigrid)
  int max_level() const
  { return max_level_; }
  
  /// Return domain lower extent
  void lower(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return root-level grid size
  void root_size(int * nx, int * ny = 0, int * nz = 0) const throw ();

  /// Set the periodicity of boundary conditions for domain axes
  inline void set_periodicity (int px, int py=0, int pz=0)
  {
    periodicity_[0] = px;
    periodicity_[1] = py;
    periodicity_[2] = pz;
  }

  /// Return the periodicity of the boundary conditions for domain axes
  void get_periodicity (int * px, int * py=0, int * pz=0) const throw()
  {
    if (px) (*px) = periodicity_[0];
    if (py) (*py) = periodicity_[1];
    if (pz) (*pz) = periodicity_[2];
  }

  /// Sets the elements of `npi` to be the periodic image of position `x` that
  /// is nearest to position `y`
  ///
  /// If the domain is not periodic along a given axis `i`, then `npi[i] = x[i]`
  ///
  /// Each argument is expected to be an array of length `cello::rank()`
  void get_nearest_periodic_image(const double * x, const double *y,
				                          double * npi) const throw();

  /// For a position `x`, this sets `folded_x` to the periodic image of `x` 
  /// which is in the domain.
  ///
  /// If the domain is not periodic along a given axis `i`, 
  /// then `folded_x[i] = x[i]`
  ///
  /// Both arguments are expected to be an array of length `cello::rank()`.
  void get_folded_position(const double * x, double * folded_x) const throw();

  //----------------------------------------------------------------------

  /// Return whether Blocks have been allocated or not
  bool blocks_allocated() const throw()
  { 
    return block_exists_;
  }

  /// Deallocate local Blocks
  void deallocate_blocks() throw();

  /// Return pointer to the Block CHARM++ chare array
  CProxy_Block block_array() const throw()
  { return block_array_;}

  /// Return pointer to the Block CHARM++ chare array
  void set_block_array(CProxy_Block block_array) throw()
  { block_array_ = block_array;}

  /// Increment (decrement) number of mesh blocks
  void increment_block_count(int count, int level);

  /// Add Block to the list of blocks (block_vec_ and block_map_)
  void insert_block (Block * block)
  {
    block_vec_.push_back(block);
  }
  
  /// Remove Block from the list of blocks (block_vec_) and return
  /// true iff Block is found in the list
  bool delete_block (Block * block)
  {
    const int n = block_vec_.size();
    bool found = false;
    for (int i=0; i<n; i++) {
      if (found) block_vec_[i-1] = block_vec_[i];
      if (block_vec_[i] == block) found=true;
    }
    if (found) block_vec_.resize(n-1);
    return found;
  }
  
  /// Increment (decrement) number of particles
  void increment_particle_count(int64_t count);

  /// Increment (decrement) number of real_zones
  void increment_real_zone_count(int64_t count)
  { num_zones_real_ += count; }

  /// Increment (decrement) number of total_zones
  void increment_total_zone_count(int64_t count)
  { num_zones_total_ += count; }

  /// Return the number of blocks on this process
  size_t num_blocks() const throw()
  {  return num_blocks_;  }

  /// Return the number of blocks on this process for the given level
  size_t num_blocks(int level) const throw()
  {  return num_blocks_level_.at(level-min_level_);  }

  /// Return the ith block in this pe
  Block * block (int index_block)
  { return block_vec_.at(index_block); }

  /// Return the number of particles on this process
  int64_t num_particles() const throw()
  {  return num_particles_;  }

  /// Return the number of real zones on this process
  int64_t num_zones_real() const throw()
  {  return num_zones_real_;  }

  /// Return the number of total zones on this process
  int64_t num_zones_total() const throw()
  {  return num_zones_total_;  }

  CProxy_Block new_block_proxy (bool allocate_data) throw();
  
  void create_block_array (bool allocate_data) throw();

  void create_subblock_array (bool allocate_data) throw();


  /// Return the number of root-level Blocks along each rank
  void root_blocks (int * nbx, int * nby=0, int * nbz=0) const throw();

  /// Return the factory object associated with the Hierarchy
  const Factory * factory () const throw()
  { return factory_; }

protected: // attributes

  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  Factory * factory_;

  /// Refinement of the hierarchy [ used for Charm++ pup() of Tree ]
  int refinement_;

  /// Minimum mesh level (may be < 0, e.g. for multigrid)
  int min_level_;

  /// Maximum mesh level
  int max_level_;

  /// Maximum number of refinement levels

  /// Current number of blocks on this process
  int num_blocks_;

  /// Current number of blocks on this process per refinement level
  std::vector<int> num_blocks_level_;
  
  /// Pointers to Blocks on this process
  std::vector<Block *> block_vec_;

  /// Current number of particles on this process
  int64_t num_particles_;

  /// Current number of total_zones on this process
  int64_t num_zones_total_; 

  /// Current number of real_zones on this process
  int64_t num_zones_real_; 
  
  /// Array of Blocks 
  CProxy_Block block_array_;
  bool           block_exists_;

  /// Size of the root grid
  int root_size_[3];

  /// Lower extent of the hierarchy
  double lower_[3];

  /// Upper extent of the hierarchy
  double upper_[3];

  /// How the array-of-octrees is distributed into Blocks
  int blocking_[3];

  /// Periodicity of boundary conditions on faces
  int periodicity_[3];

public: // static attributes

  /// Current number of blocks on this node
  static int num_blocks_node;

  /// Current number of particles on this node
  static int64_t num_particles_node;

};


#endif /* MESH_HIERARCHY_HPP */

