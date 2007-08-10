
/// Hierarchy class include file

/**
 * 
 * A Hierarchy object represents a structured adaptive mesh refinement
 * hierarchy.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */


class Hierarchy
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:

  int dimension_;                   // Dimension
  std::vector <Level *> levels_;    // List of Levels
  std::vector <Grid *>  grids_;     // List of all grids

  //  typedef std::vector <Grid *> grid_vector_type ;
  //  std::map<Grid *, grid_vector_type >  parent_; // List of each grid's parent
  std::map<Grid *, Grid * >  parent_; // List of each grid's parent

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
  void set_dim (int d) throw () { dimension_ = d; };

  // IO

  void print () throw ();
  void write (FILE * fp = 0) throw ();

  // Checking

  void check () throw () {  
    printf ("Hierarchy::check() is not implemented yet\n");
  }

  // Data access

  int dimension () throw () { return dimension_; } ;

  Level & level      (int i)            { return * levels_.at(i); };
  int     num_levels ()                 { return   levels_.size(); };

  Grid &  grid       (int i)            { return * grids_.at(i); };
  int     num_grids  ()                 { return   grids_.size(); };

  Grid &  grid       (int level, int j) { return   levels_[level]->grid(j); };
  int     num_grids  (int level)        { return   levels_[level]->num_grids(); };

  void    set_parent (Grid * grid, Grid * parent)  
  { parent_[grid]=parent; };
  Grid *  parent     (Grid * grid)      
  { return   parent_[grid]; };

  //--------------------------------------------------------------------
  // PRIVATE MEMBER FUNCTIONS
  //--------------------------------------------------------------------

private:

  void insert_in_level_ (int level, Grid & grid) throw ();



};


