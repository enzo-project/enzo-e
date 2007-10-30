//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Physical constants class

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-04-10 
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



