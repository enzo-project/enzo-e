// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialRamsesRT.hpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     06 September 2022
/// @brief    [\ref Enzo] Initialization routine for EnzoMethodRamsesRT.
//                        This is needed to compute and store large-ish tables relevant
//                        to EnzoMethodRamsesRT that shouldn't be stored as attributes
//                        of the EnzoBlock class (e.g. Reading in HLL eigenvalue tables and cross sections
//                          for a variety of different stellar masses).
//
//                        These things can't be stored as attributes of EnzoMethodRamsesRT
//                        because attributes of the derived class get sliced away after
//                        the refresh steps.      

#ifndef ENZO_ENZO_INITIAL_RAMSES_RT_HPP
#define ENZO_ENZO_INITIAL_RAMSES_RT_HPP

class EnzoInitialRamsesRT : public Initial {

  /// @class    EnzoInitialRamsesRT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for Ramses RT

public: // interface

  /// CHARM++ constructor
  EnzoInitialRamsesRT(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialRamsesRT);

  /// CHARM++ migration constructor
  EnzoInitialRamsesRT(CkMigrateMessage *m)
     : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (   Block * block, const Hierarchy * hierarchy ) throw();

  int hll_table_f     (int i, int j) const throw() { return hll_table_f_[100*i+j]; }
  int hll_table_theta (int i, int j) const throw() { return hll_table_theta_[100*i+j]; }

  double hll_table_lambda_min (int i, int j) const throw() { return hll_table_lambda_min_[100*i+j]; }
  double hll_table_lambda_max (int i, int j) const throw() { return hll_table_lambda_max_[100*i+j]; }
  double hll_table_col3       (int i, int j) const throw() { return hll_table_col3_[100*i+j]; }
  double hll_table_col4       (int i, int j) const throw() { return hll_table_col4_[100*i+j]; }

private:
  void read_hll_eigenvalues(std::string hll_file) throw(); 
 
  std::vector<int> hll_table_f_, hll_table_theta_;
  std::vector<double> hll_table_lambda_min_, hll_table_lambda_max_;
  std::vector<double> hll_table_col3_, hll_table_col4_; // not sure what these columns are

};

#endif /* ENZO_ENZO_INITIAL_RAMSES_RT_HPP */
