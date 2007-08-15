
/// Problem class include file

/**
 * 
 * Everything required to define the problem to solve
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-04-02
 *
 */

#define BUFFER_LENGTH 255

class Problem {

  //----------------------------------------------------------------------


private:

  std::vector<Sphere *> spheres_;  // List of sphere masses
  std::vector<Point *>  points_;   // List of point masses
  Hierarchy hierarchy_;            // AMR mesh hierarchy
  Domain domain_;                  // Problem domain

  //----------------------------------------------------------------------

public:

  Problem () throw ();

  ~Problem () throw ();

  Problem (const Problem & ) throw ();

  Problem & operator = (const Problem & ) throw ();

  //--------------------------------------------------

  void print () throw ();
  void write (FILE *fp = 0) throw ();
  /// Read in a problem file two
  void read (std::string filename) throw ();

  int dimension ()          // Return the dimension
  { return hierarchy_.dimension(); } 

  Hierarchy & hierarchy ()   // Return the hierarchy
  { return hierarchy_; } 
  Domain & domain ()         // Return the domain
  { return domain_; } 
  Sphere & sphere (int i)   // Return a pointer to the ith Sphere
  { return * spheres_[i]; };
  int num_spheres ()        // Return the number of Spheres
  { return spheres_.size(); };

  Point & point (int i)     // Return a pointer to the ith Point
  { return * points_[i]; };
  int num_points ()         // Return the number of Points
  { return points_.size(); };

  Grid & grid (int i)       // Return a pointer to the ith Grid
  { return hierarchy_.grid(i); };
  int num_grids ()          // Return the number of Grids
  { return hierarchy_.num_grids(); };

  Level & level (int i)       // Return a pointer to the ith Level
  { return hierarchy_.level(i); };
  int num_levels ()           // Return the number of Levels
  { return hierarchy_.num_levels(); };

  //----------------------------------------------------------------------

private:

  int readline_ (FILE *, char * buffer, int n) throw ();
};
