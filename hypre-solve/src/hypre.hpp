
/// Hypre class include file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

class Hypre {

public:

  Hypre ();

  void init_hierarchy (Hierarchy & hierarchy, 
		       Mpi       & mpi);
  void init_stencil   (Hierarchy & hierarchy);
  void init_graph     (Hierarchy & hierarchy);
  void init_matrix    (Hierarchy & hierarchy);
  void init_rhs       (Hierarchy & hierarchy);
  void init_solver    (Hierarchy & hierarchy);
  void solve          (Hierarchy & hierarchy);
  void evaluate       (Hierarchy & hierarchy);


private:

  HYPRE_SStructStencil stencil_;

private:

  void init_graph_children_               (Grid & grid);
  void init_graph_parent_                 (Hierarchy & hierarchy,
					   Grid & grid);
  void init_graph_neighbors_children_     (Grid & grid);
  void init_graph_parents_neighbor_       (Grid & grid);

};
