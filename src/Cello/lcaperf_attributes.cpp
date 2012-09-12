// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_Attributes.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-20
/// @brief    Implemention of the Attributes class

#include "lcaperf.hpp"

//----------------------------------------------------------------------

Attributes::Attributes() throw ()
{
}

//----------------------------------------------------------------------

Attributes::~Attributes() throw ()
{
}

//----------------------------------------------------------------------

Attributes::Attributes(const Attributes & attributes) throw ()
/// @param     attributes  Object being copied
{
}

//----------------------------------------------------------------------

Attributes & Attributes::operator= (const Attributes & attributes) throw ()
/// @param     attributes  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================

void Attributes::create(std::string attribute)
{
  if (index_.find(attribute) != index_.end()) {
    fprintf (stderr,_ERROR "attribute '%s' already exists!\n",
	     __FILE__,__LINE__,attribute.c_str());
    return;
  }
  int size = index_.size();
  index_[attribute] = size;
  name_. push_back(attribute);
  value_.push_back(LCAP_STRING_NULL);
}

//----------------------------------------------------------------------

void Attributes::remove (std::string attribute)
{
  if (attribute == "*") {

    index_.clear();
    name_.clear();
    value_.clear();

  } else {

    if (index_.find(attribute) == index_.end()) {
      fprintf (stderr,_ERROR "attribute '%s' does not exist!\n",
	       __FILE__,__LINE__,attribute.c_str());
    } else {

      int index = index_[attribute];
      index_.erase(index_.find(attribute));

      name_.erase (name_.begin() + index);
      value_.erase(value_.begin() + index);
    
      // adjust indices of remaining attributes
      std::map<std::string,int>::iterator iter;
      for (iter = index_.begin();
	   iter != index_.end();
	   ++iter) {
	if (iter->second > index) --iter->second;
      }
    }
  }
}

//----------------------------------------------------------------------

void Attributes::assign(std::string attribute, std::string value)
{
  if (index_.find(attribute) == index_.end()) {
    printf (_ERROR "trying to access nonexistent attribute %s!\n",
	    __FILE__,__LINE__,attribute.c_str());
    return;
  }
  value_[index_[attribute]] = value;    
}

//----------------------------------------------------------------------

