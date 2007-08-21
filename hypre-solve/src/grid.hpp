
/// Grid class include file

/**
 * 
 *
 * A Grid is an orthogonal grid of zones in \f$ R^d, 1 \le d \le 3 \f$.  Zones
 * are congruent, and each zone contains a single double value at its
 * center
 *
 * @verbatim
         o---o---o---o---o---X
         |   |   |   |   |   |
         | * | * | * | * | * |
         |   |   |   |   |   |
         o---o---o---o---o---o
         |   |   |   |   |   |
         | * | * | * | * | * |
         |   |   |   |   |   |
         o---o---o---o---o---o
         |   |   |   |   |   |
         | * | * | * | * | * |
         |   |   |   |   |   |
         X---o---o---o---o---o
  
   o vertices
   X extreme vertices xl_ and xu_
   * unknowns
   @endverbatim
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Grid
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

 private:


  unsigned            id_;         /// Unique id > 0
  unsigned            id_parent_;  // Parent grid

  // data defined at grid creation (in Grid())

  int                 ip_;         // Owning processor

  Scalar              xl_[3];      // Position of lowest vertex
  Scalar              xu_[3];      // Position of highest vertex
  int                 il_[3];      // Global index of lowest vertex
  int                 n_[3];       // Number of zones

  Faces *             faces_;      // Faces class associated with grid

  // data computed at hierarchy creation (in read())

  std::vector<Grid *> neighbors_;  // Array of neighboring grids in this level
  std::vector<Grid *> children_;   // Array of child grids
  int                 level_;      // Containing Level number (0 = root)

  // data for unknowns (in allocate())

  Scalar *            u_;          // Unknowns (single cell-centered variable)
  ;                                // Stored as 3D fortran-style array
  ;                                // with ordering x-y-z

  //--------------------------------------------------------------------

 public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  /// Create a grid given the string "<id> <id_parent> <ip> <xl>[3] <xu>[3], <n>[3]
  Grid (std::string parms) throw ();

  ~Grid () throw ();

  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  // IO

  /// Write the grid definition to standard out in human-readable format

  void print () throw ();

  /// Write the grid to a file in compact format

  void write (FILE * fp = 0) throw ();

  /// Read the grid from the given string in compact format

  void read (std::string parms) throw ();

  //--------------------------------------------------------------------
  // Data access
  //--------------------------------------------------------------------

  /// Return the grid's integer id
  unsigned int id () throw () 
  { return id_; }; 

  /// Return the grid's parent's id
  unsigned int id_parent () throw () 
  { return id_parent_; }; 

  /// Set the grid to be a child.  Should only be called once per child.
  void set_child (Grid & child) throw () 
  { children_.push_back (&child); } ;

  /// Return the ith child.
  Grid & child (int i) 
  { return * children_.at(i); };

  /// Return the number of child grids
  int num_children () 
  { return children_.size(); };

  /// Set the grid to be a neighbor.  Should only be called once per neighbor.
  void set_neighbor (Grid & neighbor) throw () 
  { neighbors_.push_back (&neighbor); } ;

  /// Return the ith neighbor
  Grid & neighbor (int i) 
  { return * neighbors_.at(i); };

  /// Return the number of neighbors
  int num_neighbors () 
  { return neighbors_.size(); };

  /// Make g1 and g2 neighbors.  Should only be called once per pair of neighbors
  friend void assert_neighbors (Grid & g1, Grid & g2) throw ()
  { g1.set_neighbor(g2); 
    g2.set_neighbor(g1); };

  /// Set level for this grid.  0 = root
  void set_level (int level) throw () 
  { level_ = level; };

  /// Return the level for this grid
  int level () throw () 
  { return level_; };

  /// Return the global lower index of the lower grid unknown.  No error checking on i.
  void set_lower(int i0, int i1, int i2) throw ()
  { il_[0]=i0; il_[1]=i1; il_[2]=i2; };

  /// Return the number of unknowns along the ith coordinate.  No error checking on i.
  int num_unknowns(int i) throw ()
  { return n_[i]; };

  /// Return the coordinates of the lower grid vertex.  No error checking on i.
  Scalar x_lower(int i) throw ()
  { return xl_[i]; };

  /// Return the coordinates of the upper grid vertex.  No error checking on i.
  Scalar x_upper(int i) throw ()
  { return xu_[i]; };

  /// Return the global lower index of the lower grid unknown.  No error checking on i.
  int i_lower(int i) throw ()
  { return il_[i]; };

  /// Return the global upper index of the upper grid unknown.  No error checking on i.
  int i_upper(int i) throw ()
  { return il_[i] + n_[i] - 1; };

  /// Return the unknown u(i,j,k)
  Scalar & unknown(int i0, int i1, int i2) throw()
  { assert (u_);
    return u_[i0 + n_[0]*(i1 + n_[1]*i2)];
  }

  Faces & faces() throw()
  { return *faces_; }

  /// Processor owner
  int ip () throw () 
  { return ip_; };

  //--------------------------------------------------------------------
  // Query functions
  //--------------------------------------------------------------------

  /// Return true iff the two grids are next to each other, or are the same grid
  bool is_adjacent (Grid & grid) throw ();
  /// Return true iff the grid belongs to processor ip
  bool is_local () throw()
  { return mpi_.ip() == ip_; };

  //--------------------------------------------------------------------
  // Processing functions
  //--------------------------------------------------------------------

  /// Return index limits of ghost zones in the neighboring grid. Required by HYPRE.

  bool find_neighbor_indices (Grid & neighbor, int *gl, int *gu);

  //--------------------------------------------------------------------
  // STATIC MEMBER DATA
  //--------------------------------------------------------------------

private:

  static Mpi  mpi_;

public:

  static void set_mpi (Mpi &mpi)
  { mpi_ = mpi; };

  //--------------------------------------------------------------------
  // PRIVATE MEMBER FUNCTIONS
  //--------------------------------------------------------------------

private:

  //--------------------------------------------------------------------
 

};


