
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
  friend class ItHierarchyLevels;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  /// Dimension
  int                        dimension_;
  /// List of each grid's parent
  std::map<Grid *, Grid * >         grid_parent_;

  HYPRE_SStructGrid   hypre_grid_;  // Struct for hypre grid
  HYPRE_SStructGraph  hypre_graph_; // Struct for hypre graph

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

protected:

  /// List of all grids, with extra sentinal 0-pointer at the end
  std::vector <Grid *>       grids0_;

  /// List of Levels,  with extra sentinal 0-pointer at the end
  std::vector <Level *>      levels0_;


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

  /// Set a Grid's parent
  void    set_parent (Grid * grid, Grid * parent) { grid_parent_[grid]=parent; };

  // -----------
  // Data access
  // -----------

  /// Return the dimension of the Hierarchy
  int dimension () throw () { return dimension_; } ;

  /// Return the ith Level of the Hierarchy.   No error checking.
  Level & level      (int i)            { return * levels0_.at(i); };

  /// Return the number of levels in the Hierarchy
  int     num_levels ()                 { return   levels0_.size() - 1; };

  /// Return the jth Grid.   No error checking.
  Grid &  grid       (int j)            { return * grids0_.at(j); };

  /// Return the number of grids
  int     num_grids  ()                 { return   grids0_.size() - 1; };

  /// Return the jth Grid in the ith Level of the Hierarchy.   No error checking.
  Grid &  grid       (int i, int j)     { return   levels0_[i]->grid(j); };

  /// Return the number of grids in the ith level.   No error checking.
  int     num_grids  (int i)            { return   levels0_[i]->num_grids(); };

  /// Return the Grid's parent Grid
  Grid *  parent     (Grid * grid)      { return   grid_parent_[grid]; };

  /// Return the HYPRE structure for the Hierarchy's hypre grid
  HYPRE_SStructGrid & hypre_grid ()     { return hypre_grid_; };

  /// Return the HYPRE structure for the Hierarchy's hypre graph
  HYPRE_SStructGraph & hypre_graph ()   { return hypre_graph_; };

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
 * An ItHierarchyLocalGrids object iterates through all local grids in
 * a Hierarchy.
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
 * An ItHierarchyAllGrids object iterates through all grids in a
 * Hierarchy, including non-local ones.
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

/// ItHierarchyLevels class

/**
 * 
 * An ItHierarchyLevels object iterates through each Level of the
 * Hierarchy in turn.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyLevels
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

  ItHierarchyLevels (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ~ItHierarchyLevels () throw () {};
  
  /// Iterate through all Levels in the Hierarchy.
  Level * operator++ (int) { 

    if (curr_ == hierarchy_->levels0_.size()) curr_ = 0;
    curr_ ++;
    return hierarchy_->levels0_[curr_-1];
  }

};
