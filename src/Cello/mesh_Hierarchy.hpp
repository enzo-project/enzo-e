// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    [\ref Mesh] Declaration of the Hierarchy class

#ifndef MESH_HIERARCHY_HPP
#define MESH_HIERARCHY_HPP


class Factory;
class Simulation;

#ifdef CONFIG_USE_CHARM
class CProxy_CommBlock;
#endif

class Hierarchy {

  /// @class    Hierarchy
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Adaptive mesh refinement hierarchy

  friend class IoHierarchy;

public: // interface

  /// Empty constructor for Charm++ pup()
  Hierarchy() throw() { }
  
  /// Initialize a Hierarchy object
  Hierarchy (
#ifndef CONFIG_USE_CHARM
	     Simulation * simulation,
#endif
	     const Factory * factory,
	      int dimension, int refinement,
	      int process_first, int process_last_plus
	      ) throw ();

  /// Delete the Hierarchy object
  virtual ~Hierarchy() throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
#endif

  //----------------------------------------------------------------------

  /// Set domain lower extent
  void set_lower(double x, double y, double z) throw ();

  /// Set domain upper extent
  void set_upper(double x, double y, double z) throw ();
  
  /// Set root-level grid size
  void set_root_size(int nx, int ny, int nz) throw ();

  //----------------------------------------------------------------------

  /// Return dimension
  int dimension() const throw ();

  /// Return domain lower extent
  void lower(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return root-level grid size
  void root_size(int * nx, int * ny = 0, int * nz = 0) const throw ();

  //----------------------------------------------------------------------

  /// Return whether CommBlocks have been allocated or not
  bool blocks_allocated() const throw()
  { 
#ifdef CONFIG_USE_CHARM
    return block_exists_;
#else
    return (block_.size() != 0);
#endif
  }

  /// Return the number of CommBlocks
  size_t num_blocks(int * nbx, 
		    int * nby = 0,
		    int * nbz = 0) const throw();

  /// Deallocate local CommBlocks
  void deallocate_blocks() throw();

# ifdef CONFIG_USE_CHARM
  /// Return pointer to the CommBlock CHARM++ chare array
  CProxy_CommBlock * block_array() const throw()
  { return block_array_;}

# else
  /// Return the total number of local blocks
  size_t num_local_blocks() const throw();

  /// Return the ith local CommBlock
  CommBlock * local_block(size_t i) const throw();
# endif

  /// Return the total number of blocks
  size_t num_blocks() const throw()
  { 
    WARNING("Hierarchy::num_blocks()",
	    "num_blocks_ initialization not implemented for AMR");
    return num_blocks_; 
  }

  void create_forest (FieldDescr   * field_descr,
		      int nx, int ny, int nz,
		      int nbx, int nby, int nbz,
		      bool allocate_blocks,
		      bool allocate_data,
		      bool testing          = false,
		      int process_first     = 0, 
		      int process_last_plus = -1) throw();


  /// Return the number of CommBlocks along each dimension
  void blocking (int * nbx, int * nby=0, int * nbz=0) const throw();

  /// Return the factory object associated with the Hierarchy
  const Factory * factory () const throw()
  { return factory_; }

  /// Return the layout of the patch, describing processes and blocking
  Layout * layout () throw();

  /// Return the layout of the patch, describing processes and blocking
  const Layout * layout () const throw();

  const GroupProcess * group_process()  const throw()
  { return group_process_; };

protected: // functions

  /// Allocate array, and optionally allocate element CommBlocks
  void allocate_array_
  (bool allocate_data = true,
   bool testing = false,
   const FieldDescr * field_descr = 0) throw ();

protected: // attributes

#ifndef CONFIG_USE_CHARM
  /// Simulation object (MPI only)
  Simulation * simulation_;
#endif

  /// 
  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  Factory * factory_;

  /// Dimension of the hierarchy [ used for Charm++ pup() of Tree ]
  int dimension_;
  
  /// Refinement of the hierarchy [ used for Charm++ pup() of Tree ]
  int refinement_;

  int num_blocks_; 

  /// Array of CommBlocks 
# ifdef CONFIG_USE_CHARM
  CProxy_CommBlock * block_array_;
  bool               block_exists_;
  Sync               block_sync_;
# else
  std::vector<CommBlock * > block_;
# endif

  /// Size of the root grid
  int root_size_[3];

  /// Lower extent of the hierarchy
  double lower_[3];

  /// Upper extent of the hierarchy
  double upper_[3];

  /// Parallel Group for distributing the Mesh across processors
  GroupProcess * group_process_;

  /// Layout: describes blocking, processor range, and processor mapping 
  Layout * layout_;

  /// How the Forest is distributed into CommBlocks
  int blocking_[3];

};


#endif /* MESH_HIERARCHY_HPP */

