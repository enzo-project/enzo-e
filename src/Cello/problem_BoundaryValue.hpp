// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-04-01
/// @brief    [\ref Problem] Declaration for the BoundaryValue component

#ifndef PROBLEM_BOUNDARY_VALUE_HPP
#define PROBLEM_BOUNDARY_VALUE_HPP

class Value;

class BoundaryValue : public Boundary
{

  /// @class    BoundaryValue
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Encapsulate a BoundaryValue conditions generator

  /// associates a value object with a list of one or more fields
  using ValueFListPair = std::pair<Value,std::vector<std::string>>;

public: // interface

  /// Create a new BoundaryValue
  ///
  /// @param[in] p is a reference to a parameter object
  /// @param[in] parameter_group is the string specifying the parameter_group
  ///     specifying the inflow boundary. Looks something like
  ///     ``"Boundary:<boundary-name>:"``
  /// @param[in] axis enum specifying the axis upon which the parameter applies
  /// @param[in] face enum specifying the face upon which the parameter applies
  ///
  /// @note
  /// axis and face are already parsed from the parameter group
  BoundaryValue(Parameters& p, const std::string& parameter_group,
                axis_enum axis, face_enum face) throw()
  : Boundary(axis,face,0),
    pairs_(BoundaryValue::construct_ValueFList_pairs_(p,parameter_group))
  { }

  /// Destructor
  virtual ~BoundaryValue() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(BoundaryValue);

  BoundaryValue(CkMigrateMessage *m)
    : Boundary (m),
      pairs_()
  { }

  /// returns a string that summarizes contents (for debugging)
  std::string debug_string() const throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;

public: // virtual functions

  /// Enforce BoundaryValue conditions

  virtual void enforce (Block   * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw();

protected: // functions

  template <class T>
  void copy_(T * field, double * value,
	     int ndx, int ndy, int ndz,
	     int nx,  int ny,  int nz,
	     int ix0, int iy0, int iz0) const throw ();

  /// Create a vector of ValueFListPair instances from the Parameters
  ///
  /// @param[in] p is a reference to a Parameter object
  /// @param[in] parameter_group is the string specifying the parameter_group
  ///     specifying the inflow boundary. Looks something like
  ///     ``"Boundary:<boundary-name>:"``
  static std::vector<ValueFListPair> construct_ValueFList_pairs_
  (Parameters& p, const std::string& parameter_group) throw();

protected: // attributes

  std::vector<ValueFListPair> pairs_;

};

#endif /* PROBLEM_BOUNDARY_VALUE_HPP */
