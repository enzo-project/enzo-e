// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12

/// @brief    [\ref Data] Declaration of the ParticleDescr class

#ifndef DATA_PARTICLE_DESCR_HPP
#define DATA_PARTICLE_DESCR_HPP

class ParticleDescr {

  /// @class    ParticleDescr
  /// @ingroup  Data
  /// @brief    [\ref Data] 

  //----------------------------------------------------------------------

  friend class Particle;

public: // interface

  /// Constructor
  ParticleDescr() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  //--------------------------------------------------
  // TYPES
  //--------------------------------------------------

  /// Create a new type and return its id

  int new_type(std::string type);

  /// Return the number of types of particles

  int num_types() const;

  /// Return the index for the given particle type

  int type_index (std::string type) const;

  /// Return the name of the given particle type given its index

  std::string type_name (int index) const;

  //--------------------------------------------------
  // CONSTANTS
  //--------------------------------------------------

  /// Create a new constant for the given type and return its id

  int new_constant (int it, std::string name, int type);

  /// Return the number of attributes of the given type.

  int num_constants(int it) const;

  /// Return the index for the given constant

  int constant_index (int it, std::string constant) const;

  /// Byte offsets of constants in constant array
  int constant_offset(int it, int ic) const;

  /// Return the index for the given constant
  std::string constant_name (int it, int ic) const;

  /// Return the number of bytes allocated for the given constant.
  int constant_bytes (int it,int ic) const;

  /// Return the constant array for the given particle type
  char * constant_array (int it);

   /// Return a pointer to the given constant for the given type
  char * constant_value (int it, int ic);

  //--------------------------------------------------
  // ATTRIBUTES
  //--------------------------------------------------

  /// Create a new attribute for the given type and return its id

  int new_attribute(int it, std::string name, int type);

  /// Return the number of attributes of the given type.

  int num_attributes(int it) const;

  /// Return the index for the given attribute

  int attribute_index (int it, std::string attribute) const;

  /// Return the index for the given attribute

  std::string attribute_name (int it, int ia) const;

  /// Byte offsets of attributes into block array.  Not including
  /// initial offset for 16-byte alignment.
  int attribute_offset(int it, int ia) const;

  /// Define which attributes represent position coordinates (-1 if not defined)
  void set_position (int it, int ix, int iy=-1, int iz=-1);

  /// Define which attributes represent velocity coordinates (-1 if not defined)
  void set_velocity (int it, int ix, int iy=-1, int iz=-1);

  //--------------------------------------------------
  // INTERLEAVING
  //--------------------------------------------------

  /// Set whether attributes are interleaved for the given type.

  void set_interleaved (int it, bool interleaved);

  /// Return whether attributes are interleaved or not

  bool interleaved (int it) const;

  /// Return the number of bytes use to represent a particle.
  int particle_bytes (int it) const;

  /// Return the data type of the given attribute.
  int attribute_type (int it,int ia) const;

  /// Return the number of bytes allocated for the given attribute.
  int attribute_bytes (int it,int ia) const;

  /// Return the attribute corresponding to the given position
  /// coordinate, -1 if none
  int attribute_position (int it, int axis)
  { return attribute_position_[it][axis]; }

  /// Return the attribute corresponding to the given velocity
  /// coordinate, -1 if none
  int attribute_velocity (int it, int axis)
  { return attribute_velocity_[it][axis]; }

  /// Return the stride of the given attribute if interleaved, otherwise 1.
  /// Computed as attribute\_bytes(it) / attribute\_bytes(it,ia).
  /// Must be evenly divisible.

  int stride(int it, int ia) const;

  //--------------------------------------------------
  // BATCHES
  //--------------------------------------------------

  /// Set batch size
  void set_batch_size(int batch_size)
  { batch_size_ = batch_size; }

  /// Return the current batch size.

  int batch_size() const;

  /// Return the batch and particle indices given a global particle
  /// index i.  This is useful e.g. for iterating over a range of
  /// particles, e.g. initializing new particles after insert().
  /// Basically just div / mod.  ASSUMES COMPRESSED.

  void index (int i, int * ib, int * ip) const;

  //--------------------------------------------------
  // GROUPING
  //--------------------------------------------------

  /// Return the Grouping object for the particle types

  Grouping * groups () { return & groups_; }

private: // functions

  /// Return true iff it and ia are in range
  bool check_(int it) const;
  bool check_ia_(int it, int ia) const;
  bool check_ic_(int it, int ic) const;

  /// increment value if needed so that it is a multiple of bytes
  int align_(int value, int bytes) const;
  //----------------------------------------------------------------------

private: // attributes

  //--------------------------------------------------
  // TYPES
  //--------------------------------------------------

  /// List of particle types
  std::vector<std::string> type_name_;

  /// Index of each particle type (inverse of type_)
  std::map<std::string,int> type_index_;

  /// Number of bytes used to store all attributes for each particle
  /// type, including any extra for alignment of stride.
  std::vector < int > particle_bytes_;

  //--------------------------------------------------
  // CONSTANTS
  //--------------------------------------------------

  /// List of particle constants
  std::vector < std::vector<std::string> > constant_name_;

  /// Index of each particle constant (inverse of constant_)
  std::vector < std::map<std::string,int> > constant_index_;

  /// Scalar type of each constant of each particle type.  Valid
  /// scalar types are in type_enum defined in cello.hpp
  std::vector < std::vector<int> > constant_type_;

  /// Arrays of byte offsets within constant_values_ for each type
  /// type.  Referenced as [it][ic]
  std::vector < std::vector <int> > constant_offset_;

  /// Number of bytes used by the given constant
  std::vector < std::vector<int> > constant_bytes_;

  /// Array of constant data [it][ic];
  std::vector< std::vector<char> > constant_array_;

  //--------------------------------------------------
  // ATTRIBUTES
  //--------------------------------------------------

  /// List of particle attributes
  std::vector < std::vector<std::string> > attribute_name_;

  /// Index of each particle attribute (inverse of attribute_)
  std::vector < std::map<std::string,int> > attribute_index_;

  /// Scalar type of each attribute of each particle type.  Valid
  /// scalar types are in type_enum defined in cello.hpp
  std::vector < std::vector<int> > attribute_type_;

  /// Attributes that define particle positions
  std::vector < std::vector <int> > attribute_position_;

  /// Attributes that define particle velocities
  std::vector < std::vector <int> > attribute_velocity_;

  /// Number of bytes used by the given attribute
  std::vector < std::vector<int> > attribute_bytes_;

  /// Whether attributes are interleaved.  (char since Charm++ pup
  /// doesn't recognize bool)

  std::vector < char > attribute_interleaved_;

  /// Arrays of byte offsets within a batch for attributes for each
  /// type.  Does not take into account byte alignment, since that's
  /// different for different blocks.  Referenced as [it][ia]

  std::vector < std::vector <int> > attribute_offset_;

  //--------------------------------------------------
  // GROUPING
  //--------------------------------------------------

  Grouping groups_;

  //--------------------------------------------------
  // BATCHES
  //--------------------------------------------------

  /// Number of particles per "batch".  Particles are allocated,
  /// deallocated, and operated on a batch at a time

  int batch_size_;
  
};

#endif /* DATA_PARTICLE_DESCR_HPP */

