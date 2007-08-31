
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

  friend class ItHierarchyGridsLocal;
  friend class ItHierarchyGridsAll;
  friend class ItHierarchyLevels;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  /// Dimension
  int                        dimension_;
  /// List of each grid's parent
  std::map<Grid *, Grid * >  grid_parent_;


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

  void init_faces (Domain & domain) throw();

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
  Grid *  parent     (Grid & grid)      { return   grid_parent_[&grid]; };

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


/// ItHierarchyGridsLocal class

/**
 * 
 * An ItHierarchyGridsLocal object iterates through all local grids in
 * a Hierarchy.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyGridsLocal
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  int               curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyGridsLocal (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ~ItHierarchyGridsLocal () throw () {};
  
  /// Iterate through local Grids in the Hierarchy.
  Grid * operator++ (int) { 

    if (curr_ == hierarchy_->grids0_.size()) curr_ = 0;
    while (hierarchy_->grids0_[curr_] && ! hierarchy_->grids0_[curr_]->is_local()) curr_++;
    curr_ ++;
    return hierarchy_->grids0_[curr_-1];
  }

};

/// ItHierarchyGridsAll class

/**
 * 
 * An ItHierarchyGridsAll object iterates through all grids in a
 * Hierarchy, including non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyGridsAll
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  int               curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyGridsAll (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ~ItHierarchyGridsAll () throw () {};
  
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

  int               curr_;
  const Hierarchy * hierarchy_;

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
