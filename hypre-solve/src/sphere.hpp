//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Sphere class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-03-26
 *
 */

class Sphere
{

  //--------------------------------------------------------------------

 private:

  Scalar     *c_;   // Center
  Scalar      m_;   // Mass
  Scalar      r_;   // Radius

  static int  d_;   // Dimension

  //--------------------------------------------------------------------

 public:

  Sphere () throw ();
  Sphere (std::string parms) throw ();

  Sphere (Scalar mass,
	  Scalar radius, 
 	  Scalar p[3], 
	  Scalar v[3]) throw ();

  ~Sphere () throw ();

  Sphere (const Sphere & ) throw ();

  Sphere & operator = (const Sphere & ) throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  static void set_dim (int d) { d_ = d; }

  void check () throw () {assert (1 <= d_ && d_ <= 3);}
  void write (FILE *fp = 0) throw ();
  void print () throw ();
  void read (std::string parms) throw ();

  //--------------------------------------------------------------------
  // PRIVATE MEMBER FUNCTIONS
  //--------------------------------------------------------------------

private:

  void alloc_ (int d = 3) throw ();
  void dealloc_ () throw ();

  //--------------------------------------------------------------------
  
};
