
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

  friend class ItLevelLocalGrids;
  friend class ItLevelAllGrids;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  int n_;

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

protected:

  /// Vector of Grid pointers, with extra sentinal 0-pointer at the end
  std::vector <Grid *>  grids0_;

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

  Grid & grid (int i) { return * grids0_.at(i); };
  int num_grids () const   { return grids0_.size() - 1; };

};


/// ItLevelLocalGrids class

/**
 * 
 * An ItLevelLocalGrids object allows iterating through all local
 * grids in a Level.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItLevelLocalGrids
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:


  const Level * level_;
  int           curr_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelLocalGrids (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ~ItLevelLocalGrids () throw () {};
  
  /// Iterate through local Grids in the Level.
  Grid * operator++ (int) { 

    if (curr_ == level_->grids0_.size()) curr_ = 0;
    while (level_->grids0_[curr_] && ! level_->grids0_[curr_]->is_local())  curr_++;
    curr_ ++;
    return level_->grids0_[curr_-1];
  }

};


/// ItLevelAllGrids class

/**
 * 
 * An ItLevelAllGrids object allows iterating through all grids in a
 * Level, even non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItLevelAllGrids
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:


  const Level * level_;
  int           curr_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelAllGrids (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ~ItLevelAllGrids () throw () {};
  
  /// Iterate through all Grids in the Level.
  Grid * operator++ (int) { 

    if (curr_ == level_->grids0_.size()) curr_ = 0;
    curr_ ++;
    return level_->grids0_[curr_-1];
  }

};


