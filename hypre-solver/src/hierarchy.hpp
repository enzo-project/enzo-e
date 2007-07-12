//======================================================================
//
//        File: Hierarchy.hpp
//
//     Summary: Hierarchy class header file
//
// Description: A Hierarchy object represents a structured adaptive mesh
//              refinement hierarchy.  
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-05-02
//
//======================================================================

class Hierarchy
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  int d_;                           // Dimension
  std::vector <Level *> levels_;    // List of Levels
  std::vector <Grid *>  grids_;     // List of all grids

  //--------------------------------------------------------------------

 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  Hierarchy () throw ();

  ~Hierarchy () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // Initialization

  void insert_grid (Grid * grid) throw ();
  void init_levels () throw();
  void set_dim (int d) throw () { d_ = d; };

  // IO

  void print () throw ();
  void write (FILE * fp = 0) throw ();

  // Checking

  void check () throw () {  
    printf ("Hierarchy::check() is not implemented yet\n");
  }

  // Data access

  int dimension () throw () { return d_; } ;

  Level & level      (int i)            { return * levels_.at(i); };
  int     num_levels ()                 { return   levels_.size(); };

  Grid &  grid       (int i)            { return * grids_.at(i); };
  int     num_grids  ()                 { return   grids_.size(); };

  Grid &  grid       (int level, int j) { return   levels_[level]->grid(j); };
  int     num_grids  (int level)        { return   levels_[level]->num_grids(); };

  //--------------------------------------------------------------------
  // PRIVATE MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  void insert_in_level_ (int level, Grid & grid) throw ();

private:


};


