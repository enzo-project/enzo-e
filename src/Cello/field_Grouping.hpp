// See LICENSE_CELLO file for license and copyright information

/// @file     field_Grouping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 16 15:27:42 PDT 2014
/// @brief    [\ref Field] Declaration of the Grouping class
///
/// This class serves to define groups of Fields, Particles, etc. into
/// named categories.  For example, one can define a "colour" group
/// for fields, and add all colour fields to the "colour" grouping.  The
/// API supports adding groups, adding "items" (field names, particle
/// set names, etc.) to groups, testing whether an item is included
/// in a group, returning the size of the group, and returning an
/// iterator for items in a group.
///
/// NOTE: This class is not named "Group" since it conflicts with a
/// Charm++ class by the same name

#ifndef FIELD_GROUPING_HPP
#define FIELD_GROUPING_HPP

class Grouping {

  /// @class    Grouping
  /// @ingroup  Field
  /// @brief    [\ref Field] 

public: // interface

  // /// Constructor
  // Grouping() throw();

  // /// Copy constructor
  // Grouping(const Grouping & Grouping) throw();

  // /// Assignment operator
  // Grouping & operator= (const Grouping & Grouping) throw();

  // /// Destructor
  // ~Grouping() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { p | groups_; }

  //----------------------------------------------------------------------

  /// Add an item to a group
  void add(std::string item, std::string group) 
    throw(std::out_of_range)
  {
    std::pair<std::string,std::string> value(item,group);
    groups_.insert(value);
  }

  /// Return whether the item is in the given group
  bool is_in(std::string item, std::string group) const
    throw(std::out_of_range)
  {
    std::pair<std::string,std::string> value(item,group);
    return groups_.find(value) != groups_.end();
  }

  /// Return the number of groups
  int size() const { return groups_.size(); }

  /// Return the number of items in the group
  int size(std::string item) const
  {
    int count = 0;
    std::set<std::pair<std::string,std::string> >::iterator it;
    for (it=groups_.begin(); it != groups_.end(); it++) {
      if (it->second == item) ++count;
    }
    return count;
  }

  std::string item (std::string group, int index)
  {
    int count = 0;
    std::set<std::pair<std::string,std::string> >::iterator it;
    for (it=groups_.begin(); it != groups_.end(); it++) {
      if (it->second == group) {
	if (count == index) return it->first;
	++count;
      }
    }
    return "";
  }
protected: // functions

  std::set<std::pair<std::string,std::string> > groups_;

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* FIELD_GROUPING_HPP */

