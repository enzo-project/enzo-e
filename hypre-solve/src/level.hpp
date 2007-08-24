
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

  friend class ItLevelGridsLocal;
  friend class ItLevelGridsAll;

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  /// Level index, with root-level == 0
  int n_;
  /// Global index of lowest grid zone in the Level
  int il_[3];
  /// Global index of highest grid zone in the Level
  int iu_[3];

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

  /// Return the number of zones extended by all Grids in the Level
  int zones (int i) throw () { return iu_[i]-il_[i] + 1; };

  // List access

  /// Return the ith Grid in the Level
  Grid & grid (int i) 
  { return * grids0_.at(i); };

  /// Return the number of Grid's in the Level
  int num_grids () const   
  { return grids0_.size() - 1; };

};


/// ItLevelGridsLocal class

/**
 * 
 * An ItLevelGridsLocal object allows iterating through all local
 * grids in a Level.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItLevelGridsLocal
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  int           curr_;
  const Level * level_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelGridsLocal (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ~ItLevelGridsLocal () throw () {};
  
  /// Iterate through local Grids in the Level.
  Grid * operator++ (int) { 

    if (curr_ == level_->grids0_.size()) curr_ = 0;
    while (level_->grids0_[curr_] && ! level_->grids0_[curr_]->is_local())  curr_++;
    curr_ ++;
    return level_->grids0_[curr_-1];
  }

};


/// ItLevelGridsAll class

/**
 * 
 * An ItLevelGridsAll object allows iterating through all grids in a
 * Level, even non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItLevelGridsAll
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  int           curr_;
  const Level * level_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelGridsAll (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ~ItLevelGridsAll () throw () {};
  
  /// Iterate through all Grids in the Level.
  Grid * operator++ (int) { 

    if (curr_ == level_->grids0_.size()) curr_ = 0;
    curr_ ++;
    return level_->grids0_[curr_-1];
  }

};


