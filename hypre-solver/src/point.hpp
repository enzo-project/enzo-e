
/// Point class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-04-10 
 *
 */

class Point 
{

public:

  Point (Scalar m=0.0, Scalar x1=0.0, Scalar x2=0.0, Scalar x3=0.0) throw () ;
  Point (std::string parms) throw ();
  ~Point () throw ();

  void print () throw ();
  void write (FILE *fp = 0) throw ();
  void read (std::string parms) throw ();

  static void set_dim (int d) { d_ = d; }

private:
  void alloc_ () throw ();
  void dealloc_ () throw ();
private:
  
  static int d_;

  Scalar *x_;
  Scalar m_;

};

