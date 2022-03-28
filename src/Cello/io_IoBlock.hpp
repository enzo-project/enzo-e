// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoBlock class

#ifndef IO_IO_BLOCK_HPP
#define IO_IO_BLOCK_HPP

class Block;

class IoBlock : public Io {

  /// @class    IoBlock
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between Block and Output classes

public: // interface

  /// Constructor
  IoBlock() throw();

  /// Destructor
  virtual ~IoBlock() throw()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoBlock);

  /// CHARM++ migration constructor for PUP::able
  IoBlock (CkMigrateMessage *m)
    : Io(m)
  {  }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;

    // NOTE: change this function whenever attributes change
    Io::pup(p);
    p | num_field_data_;
    PUParray(p,index_,3);
    PUParray(p,lower_,3);
    PUParray(p,upper_,3);
    p | cycle_;
    p | time_;
    p |  dt_;
    PUParray(p,array_,3);
  }

  /// Set block
  virtual void set_block (Block * block) throw();

  void lower (double lower[3]) {
    lower[0]=lower_[0];
    lower[1]=lower_[1];
    lower[2]=lower_[2];
  }
  void upper (double upper[3]) {
    upper[0]=upper_[0];
    upper[1]=upper_[1];
    upper[2]=upper_[2];
  }

  void index(int index3[3]) { index3[0]=index_[0]; index3[1]=index_[1]; index3[2]=index_[2]; }
  
  /// PACKING / UNPACKING
  
  /// Return the number of bytes required to serialize the data object
  virtual int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  virtual char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  virtual char * load_data (char * buffer);

  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Copy the values to the object
  virtual void save_to (void *); 

  virtual void print(const char * msg)
  {
    CkPrintf ("DEBUG_IO_BLOCK %d %s\n",CkMyPe(),msg);
    CkPrintf ("DEBUG_IO_BLOCK num_field_data_ %d\n", num_field_data_);
    CkPrintf ("DEBUG_IO_BLOCK index_          %d %d %d\n",
              index_[0], index_[1], index_[2]);
    CkPrintf ("DEBUG_IO_BLOCK lower_          %g %g %g\n",
              lower_[0], lower_[1], lower_[2]);
    CkPrintf ("DEBUG_IO_BLOCK upper_          %g %g %g\n",
              upper_[0], upper_[1], upper_[2]);
    CkPrintf ("DEBUG_IO_BLOCK cycle_          %d\n", cycle_);
    CkPrintf ("DEBUG_IO_BLOCK time_           %g\n", time_);
    CkPrintf ("DEBUG_IO_BLOCK dt_             %g\n", dt_);
    CkPrintf ("DEBUG_IO_BLOCK array_          %d %d %d\n",
              array_[0], array_[1], array_[2]);
  }
protected: // attributes

  int num_field_data_;
  int index_[3];
  double lower_[3];
  double upper_[3];
  int cycle_;
  double time_;
  double dt_;
  int array_[3];

};

#endif /* IO_IO_BLOCK_HPP */

