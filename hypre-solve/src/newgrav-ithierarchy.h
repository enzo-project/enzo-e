/// ItHierarchyGridsLocal class

/**
 * 
 * An ItHierarchyGridsLocal object iterates through all local grids in
 * a Hierarchy.
 * 
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

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyGridsLocal (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyGridsLocal (const ItHierarchyGridsLocal & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}
  void operator=(const ItHierarchyGridsLocal & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyGridsLocal () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

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

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyGridsAll (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyGridsAll (const ItHierarchyGridsAll & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}
  void operator=(const ItHierarchyGridsAll & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyGridsAll () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

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
 * Hierarchy in turn, from coarsest to finest.
 * 
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

  unsigned int      curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyLevels (Hierarchy & hierarchy) throw ()
    : curr_(0), hierarchy_(&hierarchy)
  { }

  ItHierarchyLevels (const ItHierarchyLevels & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}

  void operator=(const ItHierarchyLevels & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyLevels () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Iterate through all Levels in the Hierarchy.
  Level * operator++ (int) { 

    if (curr_ == hierarchy_->levels0_.size()) curr_ = 0;
    curr_ ++;
    return hierarchy_->levels0_[curr_-1];
  }

};

/// ItHierarchyLevelsReverse class

/**
 * 
 * An ItHierarchyLevelsReverse object iterates through each Level of the
 * Hierarchy in turn, from finest to coarsest.
 * 
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItHierarchyLevelsReverse
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  unsigned          curr_;
  const Hierarchy * hierarchy_;

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItHierarchyLevelsReverse (Hierarchy & hierarchy) throw ()
    : curr_(hierarchy.levels0_.size()-1), hierarchy_(&hierarchy)
  { }

  ItHierarchyLevelsReverse (const ItHierarchyLevelsReverse & it)
    : curr_(it.curr_), hierarchy_(it.hierarchy_)
  {}
  void operator=(const ItHierarchyLevelsReverse & it)
  { curr_ = it.curr_;  hierarchy_ = it.hierarchy_; }

  ~ItHierarchyLevelsReverse () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Iterate through all Levels in the Hierarchy from finest to coarsest
  Level * operator-- (int) { 
    printf ("DEBUG %s:%d levels0_.size()=%d\n",__FILE__,__LINE__,int(hierarchy_->levels0_.size()));
    printf ("DEBUG %s:%d levels0_ = %p\n",__FILE__,__LINE__,
	    hierarchy_->levels0_[0]);
//     if (curr_ == -1) curr_ = hierarchy_->levels0_.size()-1;
//     curr_ --;
     return hierarchy_->levels0_[curr_+1];
  }

};
