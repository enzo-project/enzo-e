//345678901234567890123456789012345678901234567890123456789012345678901234567890

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

  /// Unique id > 0
  int                 id_;
  /// Parent grid
  int                 id_parent_;

  // data defined at grid creation (in Grid())

  /// Owning processor
  int                 ip_;
  /// Position of lowest vertex
  Scalar              xl_[3];
  /// Position of highest vertex
  Scalar              xu_[3];
  /// Global index of lowest vertex
  int                 il_[3];
  /// Number of zones
  int                 n_[3];
  /// Faces class associated with grid
  Faces *             faces_;

  // data computed at hierarchy creation (in read())

  /// Containing Level number (0 = root)
  int                 level_;      

  /// Unknowns (single cell-centered variable) Stored as 3D
  /// fortran-style array
  //  Scalar *            u_;          

  /// Counters for nonstencil entries.  Required by
  /// hypre to maintain state between nonstencil and matrix initialization.
  int * counters_;

  //--------------------------------------------------------------------
  // STATIC MEMBER DATA
  //--------------------------------------------------------------------

  static Mpi    mpi_;
  static Domain domain_; // Only used by geomview stuff

  //--------------------------------------------------------------------
  // PROTECTED MEMBER DATA
  //--------------------------------------------------------------------

protected:

  friend class ItGridNeighbors;
  friend class ItGridChildren;

  /// Vector of neighboring Grids, with extra sentinal 0-pointer at the end
  std::vector<Grid *> neighbors0_; 
  /// Vector of child Grids
  std::vector<Grid *> children0_;   


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

  /// Write the grid outline to a geomview file 

  void geomview_grid (FILE *fpr, bool full=true) throw ();

  /// Write the grid's face data to a geomview file

  void geomview_face (FILE *fpr, bool full=true) throw ();

  /// Write the grid's face data of specified types to a geomview file

  void geomview_face_type (FILE         * fpr, 
			   Faces::Label * types, 
			   int            num_types, 
			   bool           full=true) throw ();

  /// Read the grid from the given string in compact format

  void read (std::string parms) throw ();

  //--------------------------------------------------------------------
  // Data access
  //--------------------------------------------------------------------

  /// Return the grid's integer id
  int id () throw () 
  { return id_; }; 

  /// Return the grid's parent's id
  int id_parent () throw () 
  { return id_parent_; }; 

  /// Set the grid to be a child.  Should only be called once per child.
  void set_child (Grid & child) throw () 
  { 
    children0_[children0_.size()-1] = & child;
    children0_.push_back (0); 
  } ;

  /// Return the ith child.
  Grid & child (int i) 
  { return * children0_.at(i); };

  /// Return the number of child grids
  int num_children () 
  { return children0_.size() - 1; };

  /// Set the grid to be a neighbor.  Should only be called once per neighbor.
  void set_neighbor (Grid & neighbor) throw () 
  { 
    neighbors0_[neighbors0_.size()-1] = & neighbor;
    neighbors0_.push_back (0); 
  } ;

  /// Determine the axis, face, and range of indices of zones adjacent 
  /// to a neighboring grid.  Returns false if the neighbor is not
  /// actually a neighbor.

  bool neighbor_shared_face (Grid & neighbor, int & axis, int & face, 
			     int & il0, int & il1, int & iu0, int & iu1) throw () ;

  /// Determine the "count"th axis (indexing from 0), face and
  /// corresponding range of coarse-grid indices of zones adjacent to
  /// the containing parent grid, and increment "count".  Returns true
  /// if the returned values are valid, or false if there is no
  /// "count"th face.

  bool parent_shared_face (Grid & parent, int & axis, int & face, 
			   int & il0, int & il1, int & iu0, int & iu1,
			   int & count) throw () ;

  /// Return the ith neighbor
  Grid & neighbor (int i) 
  { return * neighbors0_.at(i); };

  /// Return the number of neighbors
  int num_neighbors () const
  { return neighbors0_.size() - 1; };

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

  /// Set the global lower index of the lower grid unknown.  No error checking on i.
  void set_lower(int i0, int i1, int i2) throw ()
  { il_[0]=i0; il_[1]=i1; il_[2]=i2; };

  /// Return the number of unknowns along the ith coordinate.  No error checking on i.
  int num_unknowns(int i) throw ()
  { return n_[i]; };

  /// Return the total number of unknowns.
  int num_unknowns() throw ()
  { return n_[0]*n_[1]*n_[2]; };

  /// Return the coordinates of the lower grid vertex.  No error checking on i.
  void x_lower(Scalar &xl0, Scalar &xl1, Scalar &xl2) throw ()
  { 
    xl0 = xl_[0]; 
    xl1 = xl_[1]; 
    xl2 = xl_[2]; 
  };

  /// Return the coordinates of the upper grid vertex.  No error checking on i.
  void x_upper(Scalar &xu0, Scalar &xu1, Scalar &xu2) throw ()
  { 
    xu0 = xu_[0]; 
    xu1 = xu_[1]; 
    xu2 = xu_[2]; 
  };

  /// Return the global lower index of the lower grid unknown.  No error checking on i.
  int i_lower(int i) throw ()
  { return il_[i]; };

  /// Return the global lower indices of the lower grid unknown.  No error checking on il.
  void i_lower(int &il0, int &il1, int &il2) throw ()
  { 
    il0 = il_[0];
    il1 = il_[1];
    il2 = il_[2];
  };

  /// Return the global upper index of the upper grid unknown.  No error checking on i.
  int i_upper(int i) throw ()
  { return il_[i] + n_[i] - 1; };

  /// Return the global upper indices of the upper grid unknown.  No error checking on iu.
  void i_upper(int &iu0, int &iu1, int &iu2) throw ()
  { 
    iu0 = il_[0] + n_[0] - 1;
    iu1 = il_[1] + n_[1] - 1;
    iu2 = il_[2] + n_[2] - 1;
  };

  /// Return lower and upper global indices.  Upper indices are incremented by one.
  void indices(int ind[3][2]) throw ()
  { 
    ind[0][0] = il_[0];
    ind[1][0] = il_[1];
    ind[2][0] = il_[2];
    ind[0][1] = il_[0] + n_[0];
    ind[1][1] = il_[1] + n_[1];
    ind[2][1] = il_[2] + n_[2];
  };

  /// Return the unknown u(i,j,k)
  //  Scalar & unknown(int i0, int i1, int i2) throw()
  //  { assert (u_);
  //    return u_[i0 + n_[0]*(i1 + n_[1]*i2)];
  //  }

  /// Return the Faces object for this Grid.  If not allocated yet,
  /// create a new Faces object.

  Faces & faces() throw()
  { return *(faces_ ? faces_ : faces_=new Faces(n_)); }

  /// Return the mesh width along the given axis
  Scalar h(int i) throw ()
  { return (xu_[i] - xl_[i]) / n_[i]; };

  /// Return the mesh width along the given axis
  void h(Scalar &h0, Scalar &h1, Scalar &h2) throw ()
  { 
    h0 = (xu_[0] - xl_[0]) / n_[0];
    h1 = (xu_[1] - xl_[1]) / n_[1];
    h2 = (xu_[2] - xl_[2]) / n_[2];
  };    

  /// Return the grid size along the given axis
  int n(int i) throw ()
  { return n_[i]; };

  /// Return the grid size
  int n() throw ()
  {  return n_[0]*n_[1]*n_[2]; };    

  /// Return the center of the given zone
  void zone (int i0, int i1, int i2, Scalar &x0, Scalar &x1, Scalar &x2) throw ()
  { 
    x0 = xl_[0] + ((xu_[0] - xl_[0]) / n_[0]) * (i0 + 0.5); 
    x1 = xl_[1] + ((xu_[1] - xl_[1]) / n_[1]) * (i1 + 0.5); 
    x2 = xl_[2] + ((xu_[2] - xl_[2]) / n_[2]) * (i2 + 0.5); 
  };
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
  // Nonstencil / Matrix initialization functions
  //--------------------------------------------------------------------

  /// Return alias to the (i3[0],i3[1],i3[2])th counter
  int & counter (int i3[3])
  { 
    int i0=i3[0]-il_[0];
    int i1=i3[1]-il_[1];
    int i2=i3[2]-il_[2];
    assert (counters_);
    return counters_[i0 + n_[0]*(i1 + n_[1]*i2)];
  }

  /// Initialize the counters_ array to given value
  void init_counter (int value)
  {
    for (int i0=0; i0<n_[0]; i0++) {
      for (int i1=0; i1<n_[1]; i1++) {
	for (int i2=0; i2<n_[2]; i2++) {
	  counters_[i0 + n_[0]*(i1 + n_[1]*i2)] = value;
	}
      }
    }
  }

  //--------------------------------------------------------------------
  // Processing functions
  //--------------------------------------------------------------------

  /// Return index limits of ghost zones in the neighboring grid. Required by HYPRE.

  //  bool find_neighbor_indices (Grid & neighbor, int *gl, int *gu);

  //--------------------------------------------------------------------
  // STATIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Initialize static Mpi object
  static void set_mpi (Mpi &mpi)
  { mpi_ = mpi; };

  /// Set static Domain object
  static void set_domain (Domain &domain)
  { domain_ = domain; };

};


/// ItGridNeighbors class

/**
 * 
 * An ItGridNeighbors object allows iterating through all grids in a
 * Level, even non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItGridNeighbors
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  unsigned int curr_;
  const Grid * grid_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItGridNeighbors (Grid & grid) throw ()
    : curr_(0), grid_(&grid)
  { }

  ~ItGridNeighbors () throw () {};
  
  /// Iterate through all Grids in the Grid.
  Grid * operator++ (int) { 

    if (curr_ == grid_->neighbors0_.size()) curr_ = 0;
    curr_ ++;
    return grid_->neighbors0_[curr_-1];
  }

};

/// ItGridChildren class

/**
 * 
 * An ItGridChildren object allows iterating through all grids in a
 * Level, even non-local ones.
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 * @date 2007-05-02
 *
 */

class ItGridChildren
{

  //--------------------------------------------------------------------
  // PRIVATE MEMBER DATA
  //--------------------------------------------------------------------

private:

  unsigned int curr_;
  const Grid * grid_;

public:

  //--------------------------------------------------------------------
  // CONSTUCTORS AND DESTRUCTORS
  //--------------------------------------------------------------------

  ItGridChildren (Grid & grid) throw ()
    : curr_(0), grid_(&grid)
  { }

  ItGridChildren (const ItGridChildren & it)
    : curr_(it.curr_), grid_(it.grid_)
  {}
  void operator=(const ItGridChildren & it)
  { curr_ = it.curr_;  grid_ = it.grid_; }

  ~ItGridChildren () throw () {};
  
  //--------------------------------------------------------------------
  // PUBLIC MEMBER FUNCTIONS
  //--------------------------------------------------------------------

  /// Iterate through all Grids in the Grid.
  Grid * operator++ (int) { 

    if (curr_ == grid_->children0_.size()) curr_ = 0;
    curr_ ++;
    return grid_->children0_[curr_-1];
  }

};


