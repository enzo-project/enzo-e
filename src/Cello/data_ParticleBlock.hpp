// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-08-16 14:57:21
/// @brief    [\ref Data] Declaration of the ParticleBlock class

#ifndef DATA_PARTICLE_BLOCK_HPP
#define DATA_PARTICLE_BLOCK_HPP

class ParticleBlock {

  /// @class    ParticleBlock
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  ParticleBlock(ParticleDescr * pd)
    : count_(),
      data_()
  {
    data_.resize(pd->num_types());
    count_.resize(pd->num_types());
    for (int i=0; i<pd->num_types(); i++) count_[i]=0;
  };

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | count_;
    p | data_;
  }

  void create (ParticleDescr * pd,
	       int index_type,
	       int count)
  {
    const int curr_size = data_[index_type].size();
    const int create_size  = pd->particle_size(index_type) * count;
    const int new_size = curr_size + create_size;

    count_[index_type] = count;
    data_[index_type].resize(new_size);

  };

  /// Return the encoded char vector of particle data for the given type
  void data (int index_type, std::vector<char> * data )
  { *data = data_[index_type]; }

  /// Return the number of particles of the given type
  int particle_count (int index_type) const
  { return count_[index_type]; }

  //----------------------------------------------------------------------

  /// Set local coordinates of given particle type using the (integer) values
  /// in the array 0 <= a[axis][value] < max(unsigned) such that 0
  /// is at the bottom of the containing block and max(unsigned) + 1 is at
  /// the top

//----------------------------------------------------------------------

  template <class T>
  void set_local_positions (ParticleDescr * pd,
			    int index_type,  const T ** a);

  /// Set global coordinates of given particle type to the (long double) values
  /// in the array: block.lower(axis) <= a[axis][value] < block.upper(axis)

  template <class T>
  void set_global_positions (ParticleDescr * pd,
			     int index_type,  const T ** a);

  /// Return local coordinates of particles of the given type using integer
  /// values relative to the block

  template <class T>
  void local_positions (ParticleDescr * pd,
			int index_type, T ** a);

  /// Return global coordinates of particles of the given type using
  /// long double values 

  template <class T>
  void global_positions (ParticleDescr * pd,
			 int index_type, T ** a);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Number of particles for each type
  std::vector< int > count_;

  /// Packed data for each type
  std::vector< std::vector<char> > data_;
};

//----------------------------------------------------------------------

template <class T>
void ParticleBlock::set_global_positions 
(ParticleDescr * pd,
 int index_type,
 const T ** a)
{
  // get number of bytes per particle
  int particle_size = pd->particle_size(index_type);

  // get offset of position (assumes only one position per particle)

  int position_offset = -1;
  for (int index_attribute=0; 
       index_attribute<pd->num_attributes(index_type);
       index_attribute++) {
    if (pd->attribute(index_type,index_attribute) == particle_attribute_position)
      position_offset = pd->attribute_offset(index_type,index_attribute);
  }

  ASSERT1 ("ParticleBlock::set_global_positions()",
	   "particle type %s does not have a position attribute!",
	   pd->type(index_type).c_str(),
	   position_offset != -1);

  char * position;

  for (int axis=0; axis<pd->rank(); axis++) {
    for (int particle_offset=0;
	 particle_offset < data_[index_type].size();
	 particle_offset += particle_size) {
      int index_position = particle_offset + position_offset;
      position = &(data_[axis][index_position]);
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void ParticleBlock::set_local_positions 
(ParticleDescr * pd,
 int index_type,
 const T ** u)
{
}

//----------------------------------------------------------------------

template <class T>
void ParticleBlock::local_positions 
(ParticleDescr * pd,
 int index_type,
 T ** a)
{
  for (int axis=0; axis<pd->rank(); axis++) {
    for (int index=0; index<data_[index_type].size(); index++) {
    }
  }
}

//----------------------------------------------------------------------

template <class T>
void ParticleBlock::global_positions 
(ParticleDescr * pd,
 int index_type,
 T ** a)
{
  for (int axis=0; axis<pd->rank(); axis++) {
    for (int index=0; index<data_[index_type].size(); index++) {
    }
  }
}

#endif /* DATA_PARTICLE_BLOCK_HPP */

