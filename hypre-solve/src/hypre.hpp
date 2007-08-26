
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

  HYPRE_SStructGrid    grid_;   // Struct for hypre grid
  HYPRE_SStructGraph   graph_;  // Struct for hypre graph
  HYPRE_SStructMatrix  matrix_; // Struct for hypre matrix
  HYPRE_SStructStencil stencil_;

private:

  // init_graph() functions

  void init_graph_children_               (Grid & grid);
  void init_graph_parent_                 (Hierarchy & hierarchy,
					   Grid & grid);
  void init_graph_neighbors_children_     (Grid & grid);
  void init_graph_parents_neighbor_       (Grid & grid);

  // init_matrix() functions

  void init_matrix_stencil_               (Grid & grid);
  void init_matrix_clear_                 (Grid & grid);
  void init_matrix_children_              (Grid & grid);
  void init_matrix_parent_                (Hierarchy & hierarchy,
					   Grid & grid);
  void init_matrix_neighbors_children_    (Grid & grid);
  void init_matrix_parents_neighbor_      (Grid & grid);

};
