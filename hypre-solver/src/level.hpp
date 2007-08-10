
/// Level class source file

/**
 * 
 * A Level object represents the set of grids in the same level within
 * a structured adaptive mesh refinement hierarchy.
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class Level
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  int n_;
  std::vector <Grid *>  grids_;

  //--------------------------------------------------------------------

 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  Level (int level) throw ();
  ~Level () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // Initialization

  void insert_grid (Grid & grid) throw ();

  // IO

  void print () throw ();
  void write (FILE * fp = 0) throw ();

  // Checking

  void check () throw () {  
    printf ("Level::check() is not implemented yet\n");
  }

  // List access

  Grid & grid (int i) { return * grids_.at(i); };
  int num_grids ()    { return grids_.size(); };

  //--------------------------------------------------------------------
  // PRIVATE MEMBER FUNCTIONS
  //--------------------------------------------------------------------

private:

};


