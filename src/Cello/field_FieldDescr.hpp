// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    [\ref Field] Declaration for the FieldDescr class

#ifndef FIELD_FIELD_DESCR_HPP
#define FIELD_FIELD_DESCR_HPP

class FieldDescr 
{

  /// @class    FieldDescr
  /// @ingroup  Field
  /// @brief    [\ref Field] Interface for the FieldDescr class
  ///
  /// This class is used to store information about Fields in general,
  /// including field names, how they are centered (cell-centered,
  /// face-centered, on corners, etc), number of ghost zones, padding,
  /// alignment, and precision.  There is one FieldDescr object per
  /// Simulation object.  Actual Field data are stored in FieldData
  /// objects, which are each associated with a unique Block of data.

public: // functions

  /// Initialize a FieldDescr object
  FieldDescr() throw();

  /// Destructor
  ~FieldDescr() throw();

  /// Copy constructor
  FieldDescr(const FieldDescr & field_descr) throw();

  /// Assignment operator
  FieldDescr & operator= (const FieldDescr & field_descr) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    bool up = p.isUnpacking();

    // NOTE: change this function whenever attributes change
    p | name_;
    p | num_permanent_;
    p | id_;
    p | groups_;
    p | alignment_;
    p | padding_;
    p | courant_;
    p | precision_;
    int num_fields = name_.size();
    if (up) centering_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) centering_[i] = new int[3];
      PUParray(p,centering_[i],3);
    }
    if (up) ghosts_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) ghosts_[i] = new int[3];
      PUParray(p,ghosts_[i],3);
    }
  }

  /// Set alignment
  void set_alignment(int alignment) throw();

  /// Set padding
  void set_padding(int padding) throw();

  /// Set courant
  void set_courant(double courant) throw();

  /// Set centering for a field
  void set_centering(int id_field, int cx, int cy=0, int cz=0) 
    throw(std::out_of_range);

  /// Set ghosts for a field
  void set_ghosts(int id_field, int gx, int gy=0, int gz=0) 
    throw(std::out_of_range);

  /// Set precision for a field
  void set_precision(int id_field, precision_type precision) 
    throw(std::out_of_range);

  /// Insert a new permanent field
  int insert_permanent(const std::string & name_field) throw();

  /// Insert a new temporary field
  int insert_temporary(const std::string & name_field) throw();

  /// Return the number of fields
  int field_count() const throw();

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw(std::out_of_range);

  /// Return whether the field has been inserted
  bool is_field(const std::string & name) const throw();

  /// Return the integer handle for the named field
  int field_id(const std::string & name) const throw();

  //----------------------------------------------------------------------
  // Properties
  //----------------------------------------------------------------------

  Grouping * groups () { return & groups_; }
  const Grouping * groups () const { return & groups_; }

  /// alignment in bytes of fields in memory
  int alignment() const throw();

  /// padding in bytes between fields in memory
  int padding() const throw();

  /// courant number for fields
  double courant() const throw();

  /// centering of given field
  void centering(int id_field, int * cx, int * cy = 0, int * cz = 0) const 
    throw(std::out_of_range);

  /// depth of ghost zones of given field
  void ghosts(int id_field, int * gx, int * gy = 0, int * gz = 0) const 
    throw(std::out_of_range);

  /// Return precision of given field
  int precision(int id_field) const throw(std::out_of_range);

  /// Number of bytes per element required by the given field
  int bytes_per_element(int id_field) const throw();

  /// Whether the field is permanent or temporary
  bool is_permanent (int id_field) const throw()
  { return (id_field < num_permanent_); }

  /// Return the number of permanent fields
  int num_permanent() const throw()
  { return num_permanent_; }

private: // functions

  void copy_(const FieldDescr & field_descr) throw();

  int insert_(const std::string & name_field) throw();

private: // attributes

  /// String identifying each field
  std::vector<std::string> name_;

  /// Number of permanent fields; remaining are temporary
  int num_permanent_;

  /// Index of each field in name_
  std::map<std::string,int> id_;

  /// String identifying each group
  Grouping groups_;

  /// alignment of start of each field in bytes
  int alignment_;

  /// padding between fields in bytes
  int padding_;

  /// Courant number for fields
  double courant_;

  /// Precision of each field
  std::vector<int> precision_;

  /// cell centering for each field
  std::vector<int *> centering_;

  /// Ghost depth of each field
  std::vector<int *> ghosts_;

};

#endif /* FIELD_FIELD_DESCR_HPP */

