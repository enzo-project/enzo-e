// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialM1Closure.cpp
/// @author   William Hicks (whicks@ucsd.edu)
/// @date     Tue September 6 2022
/// @brief    [\ref Enzo] EnzoMethodM1Closure initializer

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialM1Closure::EnzoInitialM1Closure
(const EnzoConfig * enzo_config) throw ()
  : Initial(enzo_config->initial_cycle, enzo_config->initial_time)
{
  return;
}

//----------------------------------------------------------------------

void EnzoInitialM1Closure::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
/*
  p | hll_table_f_;
  p | hll_table_lambda_min_;
  p | hll_table_lambda_max_;
  p | hll_table_col3_;
  p | hll_table_col4_;
*/
}

//----------------------------------------------------------------------

void EnzoInitialM1Closure::enforce_block
(Block * block, const Hierarchy  * hierarchy) throw()
{
  block->initial_done();

  const EnzoConfig * enzo_config = enzo::config();
  if (enzo_config->method_M1_closure_flux_function == "HLL") {
    read_hll_eigenvalues(enzo_config->method_M1_closure_hll_file);
  }

  return;
}

//---------------------------------------------------------------------

void EnzoInitialM1Closure::read_hll_eigenvalues(std::string hll_file) throw()
{
    std::fstream inFile;
    inFile.open(hll_file, std::ios::in);

    ASSERT("EnzoInitialM1Closure::read_hll_eigenvalues()", "hll_file failed to open!",
           inFile.is_open());

    // store table in vectors
    int line_count = 10201; // 101*101
    
    hll_table_f_.resize(line_count);
    hll_table_theta_.resize(line_count);
    hll_table_lambda_min_.resize(line_count);
    hll_table_lambda_max_.resize(line_count);
    hll_table_col3_.resize(line_count);
    hll_table_col4_.resize(line_count);


    int i = 0;
    while(inFile >> this->hll_table_f_[i] >> this->hll_table_theta_[i] >> 
                    this->hll_table_lambda_min_[i] >> 
                    this->hll_table_col3_[i] >> this->hll_table_col4_[i] >> 
                    this->hll_table_lambda_max_[i]) i++;

    inFile.close();


}
