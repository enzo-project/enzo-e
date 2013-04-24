// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    [\ref Mesh] Declaration of the Block class
///

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class Block {

  friend class CommBlock;
  friend class IoBlock;

  /// @class    Block
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Block(int nx, int ny, int nz,
	int num_field_blocks,
	double xm, double xp,
	double ym, double yp,
	double zm, double zp) throw();

  /// Destructor
  ~Block() throw();

  /// Copy constructor
  Block(const Block & block) throw();

  /// Assignment operator
  Block & operator= (const Block & block) throw();

  /// Empty constructor
  Block() { }

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    bool up = p.isUnpacking();
    p | num_field_blocks_;
    // allocate field_block_[] vector first if unpacking
    if (up) field_block_.resize(num_field_blocks_);
    for (int i=0; i<num_field_blocks_; i++) {
      if (up) field_block_[i] = new FieldBlock;
      p | *field_block_[i];
    }
    PUParray(p,lower_,3);
    PUParray(p,upper_,3);
    // NOTE: change this function whenever attributes change
  }
#endif
  

  //----------------------------------------------------------------------

  /// Return the ith Field block
  const FieldBlock * field_block (int i=0) const throw();

  /// Return the ith Field block
  FieldBlock * field_block (int i=0) throw();

  /// Return domain lower extent
  /// B
  inline void lower(double * x, 
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = lower_[0];
    if (y) *y = lower_[1];
    if (z) *z = lower_[2];
  }

  //----------------------------------------------------------------------

  /// Return domain upper extent
  /// B
  inline void upper(double * x,
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = upper_[0];
    if (y) *y = upper_[1];
    if (z) *z = upper_[2];
  }


public: // virtual functions

  virtual void allocate (const FieldDescr * field_descr) throw();

private: // functions

  void copy_(const Block & block) throw();

private: // attributes

  /// Number of field blocks (required by CHARM++ PUP::er)
  int num_field_blocks_;

  /// Array of field blocks
  std::vector<FieldBlock *> field_block_;

  /// Lower extent of the box associated with the block [computable]
  /// B
  double lower_[3];

  /// Upper extent of the box associated with the block [computable]
  /// B
  double upper_[3];

  // NOTE: change pup() function whenever attributes change

};

#endif /* MESH_BLOCK_HPP */

