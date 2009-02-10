//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Physical constants class

/**
 * 
 * @file      constants.hpp
 * @brief     Declaration and implementation of the Constants class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

class Constants 
{

private:

  static const Scalar pi_;
  static const Scalar G_;          // CGS +/- 0.0007e-8

public:

  Constants () throw () ;
  ~Constants () throw ();

  static const Scalar pi () throw () { return pi_; } ;
  static const Scalar G () throw ()  { return G_; } ;

};



