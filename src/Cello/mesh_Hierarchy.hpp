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
  Hierarchy() throw() { }
  
  /// Initialize a Hierarchy object
  Hierarchy (
	     const Factory * factory,
	     int rank, int refinement,
	     int process_first, int process_last_plus
	     ) throw ();

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

  /// Return the total number of blocks
  size_t num_blocks() const throw()
  { 
    return num_blocks_; 
  }

  void increment_block_count(int increment)
  { num_blocks_ += increment; }

  void create_forest (FieldDescr   * field_descr,
		      bool allocate_data,
		      bool testing          = false,
		      int process_first     = 0, 
		      int process_last_plus = -1) throw();


  /// Return the number of Blocks along each rank
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

  int num_blocks_; 

  /// Array of Blocks 
  CProxy_Block * block_array_;
  bool           block_exists_;
  Sync           block_sync_;

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

