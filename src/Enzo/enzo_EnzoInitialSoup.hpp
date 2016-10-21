// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSoup.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-09-10
/// @brief    [\ref Enzo] Declaration of the EnzoInitialSoup class

#ifndef ENZO_SOUP_HPP
#define ENZO_SOUP_HPP

#define SOUP_IMAGE_NX  518
#define SOUP_IMAGE_NY 366
#define SOUP_CHAR_X0 10
#define SOUP_CHAR_Y0 29
#define SOUP_CHAR_DX 65
#define SOUP_CHAR_DY 92
#define SOUP_CHAR_NX 44
#define SOUP_CHAR_NY 44

class EnzoInitialSoup : public Initial {

  /// @class    EnzoInitialSoup
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the "alphabet soup" test

public: // interface

  /// Constructor
  EnzoInitialSoup
  (int cycle,
   double time,
   std::string filename,
   int rank,
   bool rotate,
   int     nx,     int ny,     int nz,
   double dpx, double dpy, double dpz,
   double dsx, double dsy, double dsz,
   double density,
   double pressure_in, double pressure_out) throw ();

  /// Constructor
  EnzoInitialSoup(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialSoup);

    /// CHARM++ migration constructor
  EnzoInitialSoup(CkMigrateMessage *m)
    : Initial (m),
      file_name_(""),
      rank_(0),
      rotate_(false),
      png_(NULL),
      density_(0.0),
      pressure_in_(0.0),
      pressure_out_(0.0),
      letter_(NULL)
  {
    for (int axis=0; axis<3; axis++) {
      array_[axis] = 1;
      d_pos_[axis] = 0.0;
      d_size_[axis] = 0.0;
    }
    for (int i=0; i<SOUP_IMAGE_NX*SOUP_IMAGE_NY; i++) {
      data_[i] = false;
    }
  }


  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Initial::pup(p);

    p | file_name_;
    p | rank_;
    PUParray(p,array_,3);
    PUParray(p,d_pos_,3);
    PUParray(p,d_size_,3);
    p | rotate_;
  }

  ~EnzoInitialSoup() throw();

public: // virtual functions

  /// Initialize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw();

  /// Return whether enforce() expects block != NULL
  virtual bool expects_blocks_allocated() const throw()
  { return true; }


private: // functions

private: // attributes

  // NOTE: change pup() function whenever attributes change

  // File name of PNG file containing font
  std::string file_name_;

  // dimensionality
  int rank_;

  // Whether to rotate letters
  bool rotate_;

  // size of alphabet array
  int array_[3];
  
  // size of random perturbations to positions (between 0 and 1)
  double d_pos_[3];

  // size of random perturbations to size
  double d_size_[3];

  /// Current pngwriter
  pngwriter * png_;

  /// Image data
  bool data_[SOUP_IMAGE_NX*SOUP_IMAGE_NY];

  /// initial density
  double density_;
  
  /// Internal and external pressure
  double pressure_in_;
  double pressure_out_;

  /// Letters in the array
  char * letter_;

  /// character position in font array
  static const int position_[26];
};

#endif /* ENZO_SOUP_HPP */

