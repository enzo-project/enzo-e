
/// Domain class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-08-14 
 *
 */

class Domain 
{

private:

  int d_;
  Scalar xl_[3];
  Scalar xu_[3];

public:

  Domain () throw ();
  Domain (std::string parms) throw ();
  ~Domain () throw ();

  void print () throw ();
  void write (FILE *fp = 0) throw ();
  void read (std::string parms) throw ();

  void set_lower (Scalar xl_0, Scalar xl_1, Scalar xl_2) 
  { xl_[0] = xl_0;
    xl_[1] = xl_1;
    xl_[2] = xl_2; };

  void set_upper (Scalar xu_0, Scalar xu_1, Scalar xu_2) 
  { xu_[0] = xu_0;
    xu_[1] = xu_1;
    xu_[2] = xu_2; };

  /// Return the coordinates of the lower grid vertex.  No error checking on i.
  void lower(Scalar &xl_0, Scalar &xl_1, Scalar &xl_2) const throw ()
  { xl_0 = xl_[0];
    xl_1 = xl_[1];
    xl_2 = xl_[2]; };


  /// Return the coordinates of the upper grid vertex.  No error checking on i.
  void upper(Scalar &xu_0, Scalar &xu_1, Scalar &xu_2) const throw ()
  { xu_0 = xu_[0];
    xu_1 = xu_[1];
    xu_2 = xu_[2]; };



};

