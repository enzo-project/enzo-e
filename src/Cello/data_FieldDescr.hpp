// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    [\ref Data] Declaration for the FieldDescr class

#ifndef DATA_FIELD_DESCR_HPP
#define DATA_FIELD_DESCR_HPP

class FieldDescr 
{

  /// @class    FieldDescr
  /// @ingroup  Data
  /// @brief    [\ref Data] Interface for the FieldDescr class
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
    p | num_temporary_;
    p | id_;
    p | groups_;
    p | alignment_;
    p | padding_;
    p | precision_;
    int num_fields = name_.size();
    if (up) centering_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) centering_[i] = new int[3];
      PUParray(p,centering_[i],3);
    }
    if (up) ghost_depth_.resize(num_fields);
    for (int i=0; i<num_fields; i++) {
      if (up) ghost_depth_[i] = new int[3];
      PUParray(p,ghost_depth_[i],3);
    }
    p | conserved_;
  }

  /// Set alignment
  void set_alignment(int alignment) throw()
  { alignment_ = alignment; }

  /// Set padding
  void set_padding(int padding) throw()
  { padding_ = padding; }

  /// Set precision for a field
  void set_precision(int id_field, int precision) throw();

  /// Set centering for a field
  void set_centering(int id_field, int cx, int cy=0, int cz=0) throw();

  /// Set ghost_depth for a field
  void set_ghost_depth(int id_field, int gx, int gy=0, int gz=0) throw();

  /// Set whether a field is a conserved quantity
  void set_conserved(int id_field, bool conserved) throw();

  /// Insert a new permanent field
  int insert_permanent(const std::string & name_field) throw();

  /// Insert a new temporary field (does not require a name)
  int insert_temporary(const std::string & name_field = "") throw();

  /// Return the number of fields
  int field_count() const throw();

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw();

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
  int alignment() const throw()
  {  return alignment_;  }

  /// padding in bytes between fields in memory
  int padding() const throw()
  { return padding_; }

  /// Return precision of given field
  int precision(int id_field) const throw()
  {  return precision_.at(id_field); }

  /// centering of given field
  void centering(int id_field, int * cx, int * cy = 0, int * cz = 0) const 
    throw();

  /// return whether the field variable is centered in the cell
  bool is_centered(int id_field) const
  { return (centering_[0] == 0 && centering_[1] == 0 && centering_[2] == 0); }

  /// depth of ghost zones of given field
  void ghost_depth(int id_field, int * gx, int * gy = 0, int * gz = 0) const 
    throw();

  /// whether the field is a conserved quantity
  bool conserved(int id_field) const 
    throw()
  { return conserved_.at(id_field); }
  
  /// Number of bytes per element required by the given field
  int bytes_per_element(int id_field) const throw();

  /// Whether the field id refers to a valid permanent field
  bool is_permanent (int id_field) const throw()
  {
    return ( (0 <= id_field) &&
	     (id_field < num_permanent_)); }

    /// Whether the field id refers to a valid temporary field
  bool is_temporary (int id_field) const throw()
  { return ((num_permanent_ <= id_field) &&
	    (id_field < num_permanent_ + num_temporary_)); }

  /// Return the number of permanent fields
  int num_permanent() const throw()
  { return num_permanent_; }

private: // functions

  void copy_(const FieldDescr & field_descr) throw();

  int insert_(const std::string & name_field,
	      bool is_permanent = true) throw();

private: // attributes

  /// String identifying each field
  std::vector<std::string> name_;

  /// Number of permanent fields declared
  int num_permanent_;

    /// Number of temporary fields declared
  int num_temporary_;

  /// Index of each field in name_
  std::map<std::string,int> id_;

  /// Groupings of fields
  Grouping groups_;

  /// alignment of start of each field in bytes
  int alignment_;

  /// padding between fields in bytes
  int padding_;

  /// Precision of each field
  std::vector<int> precision_;

  /// cell centering for each field
  std::vector<int *> centering_;

  /// Ghost depth of each field
  std::vector<int *> ghost_depth_;

  /// Whether the field is conserved or not.  If so it is multiplied
  /// by density when interpolated or coarsened, then divided again.
  /// (int instead of bool otherwise charmc complains at compile time:
  /// "Error 1: taking address of temporary (-fpermissive)"
  std::vector<int> conserved_;

};

#endif /* DATA_FIELD_DESCR_HPP */
