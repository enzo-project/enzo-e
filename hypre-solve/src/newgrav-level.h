//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Level class source file

/**
 * 
 * @file      level.hpp
 * @brief     Declaration of the Level class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
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

  /// Periodicity, or 0 if not periodic
  int period_[3];

  //--------------------------------------------------------------------
  // STATIC MEMBER DATA
  //--------------------------------------------------------------------

  static Domain domain_; // Only used by geomview stuff

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

protected:

  /// Vector of Grid pointers, with extra sentinal 0-pointer at the end
  std::vector <Grid *>  grids0_;

  //--------------------------------------------------------------------

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  Level (int level) throw ();
  ~Level () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  //--------------------------------------------------------------------
  // Initialization
  //--------------------------------------------------------------------

  void insert_grid (Grid & grid) throw ();

  //--------------------------------------------------------------------
  // IO
  //--------------------------------------------------------------------

  void print () throw ();

  /// Write the Level grids to the given open file in geomview format
  void geomview_grid  (FILE *fpr, bool full=true) throw ();

  /// Write the local Level grids to the given open file in geomview format
  void geomview_grid_local (FILE *fpr, bool full=true) throw ();

  /// Write the Level grid faces to the given open file in geomview format
  void geomview_face (FILE *fpr, bool full=true) throw ();

  /// Write the Level grid faces of the given type to the given open file in geomview format
  void geomview_face_types (FILE         * fpr, 
			    Faces::Label * types, 
			    int            num_types,
			    bool full=true) throw ();

  void write (FILE * fp = 0) throw ();

  //--------------------------------------------------------------------
  // Checking
  //--------------------------------------------------------------------

  void check () throw () {  
    NOT_IMPLEMENTED("Level::check()");
  }

  //--------------------------------------------------------------------
  // Data access
  //--------------------------------------------------------------------

  /// Return the number of zones extended by all Grids in the Level
  int zones (int i) throw () { return iu_[i]-il_[i] + 1; };

  /// Return the ith Grid in the Level
  Grid & return_grid (int i) 
  { return * grids0_.at(i); };

  /// Return the number of Grid's in the Level
  int num_grids () const   
  { return grids0_.size() - 1; };

  /// Return which Level this is in the Hierarchy.  Root is 0.
  int index () const
  { return n_; };

  //--------------------------------------------------------------------
  // STATIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Set static Domain object
  static void set_domain (Domain &domain)
  { domain_ = domain; };

};


/// ItLevelGridsLocal class

/**
 * 
 * An ItLevelGridsLocal object allows iterating through all local
 * grids in a Level.
 * 
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

  unsigned int  curr_;
  const Level * level_;


public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelGridsLocal (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ~ItLevelGridsLocal () throw () {};
  
  ItLevelGridsLocal (const ItLevelGridsLocal & it)
    : curr_(it.curr_), level_(it.level_)
  {}
  void operator=(const ItLevelGridsLocal & it)
  { curr_ = it.curr_;  level_ = it.level_; }

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

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

  unsigned      curr_;
  const Level * level_;

public:

  //--------------------------------------------------------------------
  // CONSTRUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItLevelGridsAll (Level & level) throw ()
    : curr_(0), level_(&level)
  { }

  ItLevelGridsAll (const ItLevelGridsAll & it)
    : curr_(it.curr_), level_(it.level_)
  {}
  void operator=(const ItLevelGridsAll & it)
  { curr_ = it.curr_;  level_ = it.level_; }

  ~ItLevelGridsAll () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Iterate through all Grids in the Level.
  Grid * operator++ (int) { 

    if (curr_ == level_->grids0_.size()) curr_ = 0;
    curr_ ++;
    return level_->grids0_[curr_-1];
  }

};


