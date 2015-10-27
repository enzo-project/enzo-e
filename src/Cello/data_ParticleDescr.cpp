// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleDescr class
///
/// ParticleDescr is used to describe the content of ParticleData's
/// data Separate classes are used to avoid redundant storage between
/// ParticleBlocks, which are designed to be memory-efficient.
///
/// Particles have different types (e.g. tracer, dark matter, etc.),
/// and different types have different attributes (e.g. position,
/// velocity, mass, etc.)  ParticleDescr objects store which particle
/// types each ParticleBlock contains, and provides operations to
/// assist in "decoding" the data stored in ParticleBlocks.

#include "data.hpp"

//----------------------------------------------------------------------

ParticleDescr::ParticleDescr(int batch_size) throw()
  : type_name_(),
    type_index_(),
    attribute_name_(),
    attribute_index_(),
    attribute_bytes_(),
    particle_bytes_(),
    attribute_interleaved_(),
    attribute_offset_(),
    groups_(),
    batch_size_(batch_size)
{
}

//----------------------------------------------------------------------

void ParticleDescr::pup (PUP::er &p)
{
  p | type_name_;
  p | type_index_;
  p | attribute_name_;
  p | attribute_index_;
  p | attribute_bytes_;
  p | particle_bytes_;
  p | attribute_interleaved_;
  p | attribute_offset_;
  p | groups_;
  p | batch_size_;
}

//----------------------------------------------------------------------
// Types
//----------------------------------------------------------------------

int ParticleDescr::new_type(std::string type_name)
{
  const int nt = num_types();

#ifdef CELLO_CHECK
  for (int i=0; i<nt; i++) {
    ASSERT1("ParticleDescr::new_type()",
	    "Particle type %s already exists",
	    type_name.c_str(),
	    type_name_.at(i) != type_name);
  }
#endif

  type_name_.push_back(type_name);
  attribute_interleaved_.push_back(false);
  particle_bytes_.push_back(0);

  type_index_[type_name] = nt;

  attribute_name_. resize(nt + 1);
  attribute_index_.resize(nt + 1);
  attribute_bytes_.resize(nt + 1);
  attribute_offset_.resize(nt + 1);
  attribute_offset_[nt].push_back(0);

  return nt;
}

//----------------------------------------------------------------------

int ParticleDescr::num_types() const
{
  return type_name_.size();
}

//----------------------------------------------------------------------

int ParticleDescr::type_index (std::string type_name) const
{
  std::map<const std::string,int>::const_iterator it;
  it=type_index_.find(type_name);

  return (it != type_index_.end()) ? it->second : -1;
}

//----------------------------------------------------------------------

std::string ParticleDescr::type_name (int it) const
{
  ASSERT1("ParticleDescr::type_name",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));

  return type_name_[it];
}

//----------------------------------------------------------------------
// Attributes
//----------------------------------------------------------------------

int ParticleDescr::new_attribute
(int         it,
 std::string attribute_name,
 int         attribute_bytes)
{
#ifdef CELLO_CHECK
  
  ASSERT1("ParticleDescr::new_attribute",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));
  int b = attribute_bytes;
  ASSERT1("ParticleDescr::new_attribute",
	 "attribute_bytes %d must be a power of 2",
	 attribute_bytes,
	 (b&&!(b&(b-1))));
	 
#endif

  const int ia = num_attributes(it);

  attribute_name_[it].push_back(attribute_name);
  attribute_index_[it][attribute_name] = ia;
  attribute_bytes_[it].push_back(attribute_bytes);

  // compute offset of next attribute

  const int increment = attribute_interleaved_[it] ? 1 : batch_size_;

  attribute_offset_[it].push_back
    (attribute_offset_[it][ia] + increment * attribute_bytes_[it][ia]);

  // update particle bytes
  if (attribute_interleaved_[it]) {
    int max_bytes=0;
    int num_bytes=0;
    for (size_t ia=0; ia<attribute_bytes_[it].size(); ia++) {
      num_bytes += attribute_bytes_[it][ia];
      max_bytes = std::max(max_bytes,attribute_bytes_[it][ia]);
    }
    num_bytes = align_(num_bytes,max_bytes);
    particle_bytes_[it] = num_bytes;

  } else {
    particle_bytes_[it] += attribute_bytes;
  }

  return ia;
}

//----------------------------------------------------------------------

int ParticleDescr::num_attributes(int it) const
{
  ASSERT1("ParticleDescr::num_attributes",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));
  return attribute_name_[it].size();
}

//----------------------------------------------------------------------

int ParticleDescr::attribute_index (int it, std::string attribute_name) const
{
  ASSERT1("ParticleDescr::attribute_index",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));

  std::map<const std::string,int>::const_iterator iter;

  iter=attribute_index_[it].find(attribute_name);

  int index = (iter != attribute_index_[it].end()) ? iter->second : -1;

  if (index == -1) {
    WARNING2("ParticleDescr::attribute_index()",
	     "Trying to access unknown attribute %s in particle type \"%s\"",
	     attribute_name.c_str(),type_name(it).c_str());
  }

  return index;
}

//----------------------------------------------------------------------

std::string ParticleDescr::attribute_name (int it, int ia) const
{
  ASSERT2("ParticleDescr::attribute_name",
	  "Trying to access unknown particle attribute %d in type %d",
	  ia,it,
	  check_(it,ia));

  return attribute_name_[it][ia];
}

//----------------------------------------------------------------------

int ParticleDescr::attribute_bytes(int it,int ia) const
{
  ASSERT2("ParticleDescr::attribute_bytes",
	  "Trying to access unknown particle attribute %d in type %d",
	  ia,it,
	  check_(it,ia));

  return attribute_bytes_[it][ia];
}

//----------------------------------------------------------------------

int ParticleDescr::particle_bytes(int it) const
{
  ASSERT1("ParticleDescr::attribute_bytes",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));

  return particle_bytes_[it];

}

//----------------------------------------------------------------------

int ParticleDescr::stride(int it, int ia) const
{
  ASSERT2("ParticleDescr::stride",
	  "Trying to access unknown particle attribute %d in type %d",
	  ia,it,
	  check_(it,ia));

  return attribute_interleaved_[it] ? 
    particle_bytes_[it] / attribute_bytes(it,ia) : 1;
}

//----------------------------------------------------------------------

void ParticleDescr::set_interleaved (int it, bool interleaved)
{
  ASSERT1("ParticleDescr::set_interleaved",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));
  attribute_interleaved_[it] = interleaved;
}

//----------------------------------------------------------------------

bool ParticleDescr::interleaved (int it) const
{
  ASSERT1("ParticleDescr::interleaved",
	  "Trying to access unknown particle type %d",
	  it,
	  check_(it));
  return attribute_interleaved_.at(it);
}

//----------------------------------------------------------------------

int ParticleDescr::batch_size () const
{
  return batch_size_;
}

//----------------------------------------------------------------------

void ParticleDescr::index (int i, int * ib, int * ip) const
{
  *ib = i / batch_size_;
  *ip = i % batch_size_;
}

//----------------------------------------------------------------------

int ParticleDescr::attribute_offset (int it, int ia) const
{
  ASSERT2("ParticleDescr::attribute_offset",
	  "Trying to access unknown particle attribute %d in type %d",
	  ia,it,
	  check_(it,ia));

  return attribute_offset_[it][ia]; 
}

//======================================================================

bool ParticleDescr::check_(int it) const
{
  return (0 <= it && it < num_types());
}

bool ParticleDescr::check_(int it, int ia) const
{
  return ( (0 <= it && it < num_types()) &&
	   (0 <= ia && ia < num_attributes(it)));
}

//----------------------------------------------------------------------

int ParticleDescr::align_ (int value, int bytes) const
{
  int align = value % bytes;
  if (align != 0) {
    value += (bytes - align);
  }
  return value;
}
