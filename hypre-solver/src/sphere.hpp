//======================================================================
//
//        File: sphere.hpp
//
//     Summary: Sphere class header file
//
// Description:
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-03-26
//
//======================================================================

class Sphere
{

  //--------------------------------------------------------------------

 private:

  Scalar     *c_;   // Center
  Scalar      r_;   // Radius
  Scalar      m_;   // Mass

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
