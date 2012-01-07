// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Mesh] Declaration of the Node class
///

#ifndef MESH_NODE_HPP
#define MESH_NODE_HPP

class Node {

  /// @class    Node
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Node() throw();

  /// Destructor
  ~Node() throw();

  /// Copy constructor
  Node(const Node & node) throw();

  /// Assignment operator
  Node & operator= (const Node & node) throw();

  /// Refine the Node in its containing Tree
  void refine ();

  /// Create a copy of this Node
  void copy ();

  /// Balance the Node with respect to its neighbors.  Balancing may
  /// be with respect to some other level metric, e.g. containing
  /// Patch refinement level.
  void balance();

  /// Return the pointer to the contained data, e.g. Patch object.
  void data();

  /// Find the neighboring node in the specified direction and
  /// level. Consider including edge and corner neighbors, or sets of
  /// neighbors. May use NodeInfo to improve performance.
  void neighbor ();

  /// Return the parent node.
  void parent ();

  /// Return the child node.
  void child ();

  /// Return the node containing the given point
  Node * find_node (double x, double y, double z);


private: // functions


private: // attributes


};

#endif /* MESH_NODE_HPP */

