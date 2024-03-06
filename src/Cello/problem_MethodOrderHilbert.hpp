// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrderHilbert.hpp
/// @author   John Brennan (john.brennnan@mu.ie)
/// @date     2023-11-13
/// @brief    [\ref Problem] Declaration of the MethodOrderHilbert class for
///           generating the Hilbert ordering of blocks in the hierarchy. This
///           method uses the same communication pattern as MethodOrderMorton
///           and an iterative algorithm based on lookup tables to compute 
///           order indices for blocks. The iterative algorithm is based on
///           those presented in "Mengjuan Li et al 2023 (Efficient entry 
///           point encoding and decoding algorithms on 2D Hilbert space 
///           filling curve)" and "Lianyin Jia et al 2022 (Efficient 3D 
///           Hilbert curve encoding and decoding algorithms)".
///           Note: There are some typos in the lookup tables presented in
///           papers cited above. The tables used by this method should be
///           the corrected versions.

#ifndef PROBLEM_METHOD_ORDER_HILBERT_HPP
#define PROBLEM_METHOD_ORDER_HILBERT_HPP

class MethodOrderHilbert : public Method {

  /// @class    MethodOrderHilbert
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MethodOrderHilbert(int min_level) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodOrderHilbert);
  
  /// Charm++ PUP::able migration constructor
  MethodOrderHilbert (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    Method::pup(p);
    p | is_index_;
    p | is_count_;
    p | is_next_;
    p | is_weight_;
    p | is_weight_child_;
    p | is_sync_index_;
    p | is_sync_weight_;
    p | min_level_;
  }

  void compute_continue( Block * block);
  void compute_complete( Block * block);
  void send_weight(Block * block, int weight, bool self);
  void recv_weight(Block * block, int ic3[3], int weight, bool self);
  void send_index(Block * block, int index, int count, bool self);
  void recv_index(Block * block, int index, int count, bool self);

public: // virtual methods
  
  /// Apply the method to determine the Hilbert ordering of blocks
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "order_hilbert"; }

private: // methods

  /// Return the pointer to the Block's Hilbert ordering index 
  long long * pindex_(Block * block);

  /// Return the pointer to the number of Block indices
  long long * pcount_(Block * block);

  /// Return the pointer to the Index of the "next" block
  Index * pnext_(Block * block);

  /// Return the pointer to the Block's weight (including self)
  long long * pweight_(Block * block);

  /// Return the pointer to the given Block's child weight
  long long * pweight_child_(Block * block, int index);

  /// Return the pointer to the Block's Hilbert ordering index 
  Sync * psync_index_(Block * block);

  /// Return the pointer to the Block's weight (including self)
  Sync * psync_weight_(Block * block);

  /// Write child blocks of the given block to the children array in the hilbert order.
  void hilbert_children(Block * block, int* children);

  /// Return the index appearing after the given index in the hilber order.
  Index hilbert_next (Index index, int rank, bool is_leaf, int min_level);

  /// Write, to the given states array, the recursive states of the given index up to level m. 
  void hilbert_states(Index index, int m, int* states);

  /// Check if the given coordinate is the last in the hilbert order with state T. 
  bool is_last_child(int T, int coord);

  /// Convert a zyx coordinate to the corresponding hilbert index for a given state.
  int coord_to_hilbert_ind(int state, int coord);

  /// Convert a zyx coordinate to the next state for a given state.
  int coord_to_next_state(int state, int coord);

  /// Convert a hilbert index to the corresponding zyx coordinate for a given state.
  int hilbert_ind_to_coord(int state, int hilbert_ind);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Block Scalar<int> index
  int is_index_;
  /// Block Scalar<int> count
  int is_count_;
  /// Block Scalar<Index> next
  int is_next_;
  /// Block Scalar<int> weight (#decendent blocks + self)
  int is_weight_;
  /// Block Scalar<int> child weight (array of size cello::num_children())
  int is_weight_child_;
  /// Block Scalar<Sync> sync counter for index (coarse->fine)
  int is_sync_index_;
  /// Block Scalar<sync> sync counter for weight (fine->coarse)
  int is_sync_weight_;

  /// Minimum refinement level for ordering; may be < 0
  int min_level_;

  /// Look up tables for encoding/decoding Hilbert indices
  static int HPM[12][8];
  static int HNM[12][8];
  static int PHM[12][8];
  static int PNM[12][8];

  static int HCM[4][4];
  static int HSM[4][4];
  static int CHM[4][4];
  static int CSM[4][4];
};

#endif /* PROBLEM_METHOD_ORDER_HILBERT_HPP */

