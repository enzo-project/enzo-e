// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    [\ref Field] Declaration for the FieldDescr class

#ifndef FIELD_FIELD_DESCR_HPP
#define FIELD_FIELD_DESCR_HPP

enum field_action_enum {
  field_action_unknown,  // Uninitialized action
  field_action_none,     // Do nothing if range exceeded
  field_action_assign,   // Assign field values to min / max if range exceeded
  field_action_warning,  // Issue warning if range exceeded
  field_action_error,    // Issue error if range exceeded
  field_action_timestep, // Retry with reduced timestep if range exceeded
  field_action_method    // Retry with alternate method if range exceeded
};

typedef int field_action_type;

class FieldDescr 
{

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    [\ref Field] Interface for the FieldDescr class

public: // functions

  /// Initialize a FieldDescr object
  FieldDescr() throw();

  /// Destructor
  ~FieldDescr() throw();

  /// Copy constructor
  FieldDescr(const FieldDescr & field_descr) throw();

  /// Assignment operator
  FieldDescr & operator= (const FieldDescr & field_descr) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    bool up = p.isUnpacking();

    // NOTE: change this function whenever attributes change
    p | field_name_;
    //    p | field_id_;
    WARNING("FieldDescr::pup",
	    "Skipping field_id_ [ p|field_id_ is an error in PGI due to const ]");
    p | group_name_;
    //    p | group_id_;
    WARNING("FieldDescr::pup",
	    "Skipping group_id_ [ p|group_id_ is an error in PGI due to const ]");
    p | alignment_;
    p | padding_;
    p | courant_;
    p | precision_;
    int num_fields = field_name_.size();
    if (up) centering_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) centering_[i] = new bool[3];
      PUParray(p,centering_[i],3);
    }
    if (up) ghosts_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) ghosts_[i] = new int[3];
      PUParray(p,ghosts_[i],3);
    }

    p | min_value_;
    p | max_value_;
    p | min_action_;
    p | max_action_;
  }
#endif

   // /// Set membership of a field in a group
   // void set_field_in_group(int id_field, int id_group) 
   //   throw(std::out_of_range);

  /// Set alignment
  void set_alignment(int alignment) throw();

  /// Set padding
  void set_padding(int padding) throw();

  /// Set courant
  void set_courant(double courant) throw();

  /// Set centering for a field
  void set_centering(int id_field, bool cx, bool cy=true, bool cz=true) 
    throw(std::out_of_range);

  /// Set ghosts for a field
  void set_ghosts(int id_field, int gx, int gy=0, int gz=0) 
    throw(std::out_of_range);

  /// Set precision for a field
  void set_precision(int id_field, precision_type precision) 
    throw(std::out_of_range);

  /// Set minimum bound and action
  void set_minimum (int id_field, double min_value,
		    field_action_type min_action) 
    throw(std::out_of_range);

  /// Set maximum bound and action
  void set_maximum (int id_field, double max_value, 
		    field_action_type max_action) 
    throw(std::out_of_range);

  /// Insert a new field
  int insert_field(const std::string & name_field) throw();

  /// Insert a new group
  void insert_group(const std::string & name_group) throw();

  //----------------------------------------------------------------------

  /// Return the number of fields
  int field_count() const throw();

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw(std::out_of_range);

  /// Return whether the field has been inserted
  bool is_field(const std::string & name) const throw();

  /// Return the integer handle for the named field
  int field_id(const std::string & name) const throw();

  /// Return the number of groups
  int group_count() const throw();

  /// Return name of the ith group
  std::string group_name(int id_group) const throw(std::out_of_range);

  /// Return whether the group has been inserted
  bool is_group(const std::string & name) const throw();

  /// Return the integer handle for the named group
  int group_id(const std::string & name) const throw();


   // /// Return whether the given field is in the given group
   // bool field_in_group(int id_field, int id_group) 
   //   const throw(std::out_of_range);


  /// alignment in bytes of fields in memory
  int alignment() const throw();

  /// padding in bytes between fields in memory
  int padding() const throw();

  /// courant number for fields
  double courant() const throw();

  /// centering of given field
  void centering(int id_field, bool * cx, bool * cy = 0, bool * cz = 0) const 
    throw(std::out_of_range);

  /// depth of ghost zones of given field
  void ghosts(int id_field, int * gx, int * gy = 0, int * gz = 0) const 
    throw(std::out_of_range);

  /// Return precision of given field
  precision_type precision(int id_field) const throw(std::out_of_range);

  /// Number of bytes per element required by the given field
  int bytes_per_element(int id_field) const throw();

  /// minimum value for the field
  double minimum_value(int id_field) const throw(std::out_of_range);

  /// minimum action for the field
  field_action_type minimum_action(int id_field) const
    throw(std::out_of_range);

  /// maximum value for the field
  double maximum_value(int id_field) const  throw(std::out_of_range);

  /// maximum action for the field
  field_action_type maximum_action(int id_field) const 
    throw(std::out_of_range);

private: // functions

  void copy_(const FieldDescr & field_descr) throw();

private: // attributes

  /// String identifying each field
  std::vector<std::string> field_name_;

  /// Index of each field in field_name_
  std::map<const std::string,int> field_id_;

  /// String identifying each group
  std::vector<std::string> group_name_;

  /// Index of each group in group_name_
  std::map<const std::string,int> group_id_;

  /// alignment of start of each field in bytes
  int alignment_;

  /// padding between fields in bytes
  int padding_;

  /// Courant number for fields
  double courant_;

  /// Precision of each field
  std::vector<precision_type> precision_;

  /// cell centering for each field
  std::vector<bool *> centering_;

  /// Ghost depth of each field
  std::vector<int *> ghosts_;

  /// minimum allowed value for each field
  std::vector<double> min_value_;

  /// maximum allowed value for each field
  std::vector<double> max_value_;

  /// what to do if a field violates its minimum value
  std::vector<field_action_type> min_action_;

  /// what to do if a field violates its maximum value
  std::vector<field_action_type> max_action_;

};

#endif /* FIELD_FIELD_DESCR_HPP */

