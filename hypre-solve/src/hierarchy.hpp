
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

  friend class ItHierarchyLocalGrids;
  friend class ItHierarchyAllGrids;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  /// Dimension
  int                        dimension_;
  /// List of Levels
  std::vector <Level *>      levels_;
  /// List of each grid's parent
  std::map<Grid *, Grid * >  parent_;

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

protected:

  /// List of all grids, with extra sentinal 0-pointer at the end
  std::vector <Grid *>       grids0_;


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

  void init_grids () throw();

  void init_faces () throw();

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

  Grid &  grid       (int i)            { return * grids0_.at(i); };
  int     num_grids  ()                 { return   grids0_.size() - 1; };

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

  void init_grid_parents_() throw();
  void init_grid_levels_() throw();
  void init_grid_children_() throw();
  void init_grid_neighbors_() throw();

  void insert_in_level_ (int level, Grid & grid) throw ();



};


/// ItHierarchyLocalGrids class

/**
 * 
 * An ItHierarchyLocalGrids object allows iterating through all local
 * grids in a Hierarchy.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyLocalGrids
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:


  const Hierarchy * hierarchy_;
  int           curr_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyLocalGrids (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ~ItHierarchyLocalGrids () throw () {};
  
  /// Iterate through local Grids in the Hierarchy.
  Grid * operator++ (int) { 

    if (curr_ == hierarchy_->grids0_.size()) curr_ = 0;
    while (hierarchy_->grids0_[curr_] && ! hierarchy_->grids0_[curr_]->is_local()) curr_++;
    curr_ ++;
    return hierarchy_->grids0_[curr_-1];
  }

};

/// ItHierarchyAllGrids class

/**
 * 
 * An ItHierarchyAllGrids object allows iterating through all grids in
 * a Hierarchy, including non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyAllGrids
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:


  const Hierarchy * hierarchy_;
  int               curr_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyAllGrids (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ~ItHierarchyAllGrids () throw () {};
  
  /// Iterate through all Grids in the Hierarchy.
  Grid * operator++ (int) { 

    if (curr_ == hierarchy_->grids0_.size()) curr_ = 0;
    curr_ ++;
    return hierarchy_->grids0_[curr_-1];
  }

};
