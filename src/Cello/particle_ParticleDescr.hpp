// See LICENSE_CELLO file for license and copyright information

/// @file     particle_ParticleDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 14 17:16:28 PDT 2014
/// @brief    [\ref Particle] Declaration of the ParticleDescr class

#ifndef PARTICLE_PARTICLE_DESCR_HPP
#define PARTICLE_PARTICLE_DESCR_HPP

class ParticleDescr {

  /// @class    ParticleDescr
  /// @ingroup  Particle
  /// @brief    [\ref Particle] 

  //----------------------------------------------------------------------

public: // interface

  /// Constructor
  ParticleDescr(int rank) throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Return the rank of the space the particles live in
  int rank () const { return rank_; }

  /// Return the number of bytes required to represent the given attribute
  int attribute_size (int particle_attribute) const throw(std::out_of_range)
  { return attribute_size_.at(particle_attribute); }

  /// Return the offset in bytes for the given attribute in the given type
  int attribute_offset (int index_type, int index_attribute) 
    const throw(std::out_of_range)
  { return attribute_offset_[index_type][index_attribute]; }

  /// Set the number of bytes used to store particle positions
  void set_attribute_size (int particle_attribute, int size) throw()
  { attribute_size_[particle_attribute] = size; }

  /// Add a particle type to the list of particle types used
  int new_type (std::string type) throw()
  {
    type_.push_back(type);
    return (type_.size()-1);
  }

  /// Add an attribute to the given type
  void add_attribute (int index_type, int attribute) throw ();

  /// Return the number of collections of particles
  int num_types() const throw()
  { return type_.size(); }

  /// Return type of the ith collection of particles
  std::string type(int index_type) const throw(std::out_of_range)
  { return type_.at(index_type); }

  /// Return name of the ith collection of particles
  int attribute(int index_type, int index_attribute) const
    throw(std::out_of_range)
  { return attribute_[index_type][index_attribute]; }

  /// Return number of attributes in the given particle type
  int num_attributes(int index_type) const throw(std::out_of_range)
  { return attribute_[index_type].size(); }
    
  /// Return number of bytes in the given particle type
  int particle_size(int index_type) const throw(std::out_of_range)
  { int bytes=1;
    for (int index_attribute=0;
	 index_attribute<num_attributes(index_type);
	 index_attribute++) {
      bytes += attribute_size_[attribute_[index_type][index_attribute]];
    }
    return bytes;
  }

  //----------------------------------------------------------------------

  Grouping * groups () throw () { return & groups_; }
  const Grouping * groups () const throw () { return & groups_; }

  /// Display this ParticleDescr for debugging
  void print () const;

private: // functions

  //----------------------------------------------------------------------

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Dimensionality of the particles
  int rank_;

  /// List of names of different groups of particles
  std::vector<std::string> type_;

  /// Number of bytes used to store each attribute
  std::vector<int> attribute_size_;

  /// Attributes for each particle type, e.g. type_[id_type][index_attribute] = particle_attribute_position
  std::vector < std::vector<int> > attribute_;  

  /// attribute_offset_[index_type][index_attribute] stores the offset
  /// of the attribute for the given particle type
  std::vector< std::vector<int> > attribute_offset_;
  
  /// String identifying each group
  Grouping groups_;

};

#endif /* PARTICLE_PARTICLE_DESCR_HPP */

