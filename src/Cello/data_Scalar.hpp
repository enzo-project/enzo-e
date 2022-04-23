// See LICENSE_CELLO file for license and copyright information

/// @file     data_Scalar.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-07-17
/// @brief    [\ref Data] Declaration of the Scalar class
///
/// The Scalar class is used to unify the interface of the global ScalarDescr
/// object and a given ScalarData object.

#ifndef DATA_SCALAR_HPP
#define DATA_SCALAR_HPP

template <class T>
class Scalar {

  /// @class    Scalar
  /// @ingroup  Data
  /// @brief    [\ref Data]

public: // interface

  /// Constructor
  Scalar(ScalarDescr * scalar_descr,
	 ScalarData<T> * scalar_data) throw()
    : scalar_descr_ (scalar_descr),
      scalar_data_ (scalar_data)
  {}

  /// Copy constructor
  Scalar(const Scalar<T> & scalar) throw()
  {
    scalar_descr_ = scalar.scalar_descr_;
    scalar_data_ = scalar.scalar_data_;
  }

  /// Assignment operator
  Scalar & operator= (const Scalar<T> & scalar) throw()
  {
    scalar_descr_ = scalar.scalar_descr_;
    scalar_data_ = scalar.scalar_data_;
    return *this;
  }

  /// Destructor
  ~Scalar() throw()
  {};

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    static bool warn[CONFIG_NODE_SIZE] = {false};
    const int in = cello::index_static();
    if (! warn[in]) {
      WARNING ("Scalar::pup()",
	       "Skipping since Scalar is intended as transient objects");
      warn[in] = true;
    }
  };

  /// Return the scalar descriptor for this scalar
  ScalarDescr * scalar_descr() { return scalar_descr_; }

  /// Return the scalar data for this scalar
  ScalarData<T> * scalar_data() { return scalar_data_; }

  void set_scalar_data(ScalarData<T> * scalar_data)
  { scalar_data_ = scalar_data; }

  void set_scalar_descr(ScalarDescr * scalar_descr)
  { scalar_descr_ = scalar_descr; }


  //--------------------------------------------------
  // ScalarDescr
  //--------------------------------------------------

  /// Reserve space for a new scalar (or array if n>1)
  int new_value (std::string name, int n=1)
  { return scalar_descr_->new_value (name,n); }

  /// Return the *first* index of the named scalar
  int index (std::string name) const
  { return scalar_descr_->index(name); }

  /// Return the name of the given scalar
  std::string name (int index) const
  { return scalar_descr_->name(index); }

  /// Return the length of the given scalar if array
  int length (int index) const
  { return scalar_descr_->length(index); }

  /// Return the number of variables currently stored (including scalar
  /// array elements)
  int size() const
  { return scalar_descr_->size(); }

  //--------------------------------------------------
  // ScalarData
  //--------------------------------------------------

  /// Allocate scalars
  void allocate()
  { return scalar_data_->allocate(scalar_descr_); }

  /// Return the specified scalar value
  T * value (int index)
  { return scalar_data_->value(scalar_descr_,index); }

private: // attributes

  /// Scalar descriptor for global scalar data
  ScalarDescr * scalar_descr_;

  /// Scalar data for the specific Block
  ScalarData<T> * scalar_data_;

  // NOTE: change pup() function whenever attributes change

};

#endif /* DATA_SCALAR_HPP */
