// See LICENSE_CELLO file for license and copyright information

/// @file     lcaperf_Attributes.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-19
/// @brief    [\ref Lcaperf] Declaration of the Attributes class

#ifndef LCAPERF_ATTRIBUTES_HPP
#define LCAPERF_ATTRIBUTES_HPP

class Attributes {

  /// @class    Attributes
  /// @ingroup  lcaperf
  /// @brief    [\ref Lcaperf] Class for storing attribute values

public: // interface

  /// Constructor
  Attributes() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Attributes() throw();

  /// Copy constructor
  Attributes(const Attributes & Attributes) throw();

  /// Assignment operator
  Attributes & operator= (const Attributes & Attributes) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p | index_;
    p | name_;
    p | value_;
  }
#endif

  //----------------------------------------------------------------------

  /// Create a new attribute with the given name
  void create(std::string attribute);

  /// Delete an attribute
  void remove (std::string attribute);

  /// Assign a value to the named attribute
  void assign(std::string attribute, std::string value);

  /// Return the number of attributes
  int size () const
  { return name_.size(); }

  /// Return the name for the given attribute index
  std::string name (int index) const
  { return name_.at(index); };

  /// Return the index for the given attribute name
  int index (std::string name)
  { return index_[name]; };

  /// Return the current value of the given attribute
  std::string value (int i) const
  { return value_[i]; };

  /// Return the key for the current set of attribute values
  std::string get_key() const
  { std::string key = "";
    for (int i=0; i<size(); i++) {
      key = key + ":" + value(i);
    }
    return key;
  }
  
  /// Set current attribute values given a key
  void set_key(std::string key)
  { 
    int i=0,i_start,i_stop=-1;
    do {
      i_start = i_stop;
      i_stop = key.find(":",i_start+1);
      if (i < size()) value_[i++] = key.substr(i_start+1,i_stop-i_start-1);
    } while (i_stop != -1);
  }

  bool keys_match (std::string key_1, std::string key_2) const
  {
    // Obviously match if identical
    if (key_1 == key_2) return true;

    // If not identical, then check each field for matching wild cards '*'
    int i1_start,i1_stop=-1;
    int i2_start,i2_stop=-1;
    do {
      i1_start = i1_stop;
      i1_stop = key_1.find(":",i1_start+1);

      i2_start = i2_stop;
      i2_stop = key_2.find(":",i2_start+1);
      std::string field_1 = key_1.substr(i1_start+1,i1_stop-i1_start-1);
      std::string field_2 = key_2.substr(i2_start+1,i2_stop-i2_start-1);
      if (field_1 != field_2 && 
	  (field_1 != "*" && field_2 != "*")) return false;
				 
    } while (i1_stop != -1 && i2_stop != -1);
    // one key or other has more fields
    if (i1_stop != -1 || i2_stop != -1) return false;
    return true;
  }

protected: //functions

  //----------------------------------------------------------------------

protected: // attributes

  /// Map attribute name to attribute index
  std::map<std::string,int> index_;

  /// Map attribute index to attribute name
  std::vector<std::string>  name_;

  /// Map attribute index to attribute value
  std::vector<std::string>  value_;

};

#endif /* LCAPERF_ATTRIBUTES_HPP */

