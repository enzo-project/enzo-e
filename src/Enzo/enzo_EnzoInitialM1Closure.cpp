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
   Field field = block->data()->field();
   int mx,my,mz;  
   field.dimensions(0,&mx, &my, &mz); //field dimensions, including ghost zones
  
   const EnzoConfig * enzo_config = enzo::config();
         EnzoUnits * enzo_units = enzo::units();

   double inverse_Nunit = enzo_units->volume();
   double inverse_Funit = inverse_Nunit / enzo_units->velocity();

   for (int i=0; i<enzo_config->method_m1_closure_N_groups; i++) {
    std::string istring = std::to_string(i);
    enzo_float *  N_i = (enzo_float *) field.values("photon_density_" + istring);
    enzo_float * Fx_i = (enzo_float *) field.values("flux_x_" + istring);
    enzo_float * Fy_i = (enzo_float *) field.values("flux_y_" + istring);
    enzo_float * Fz_i = (enzo_float *) field.values("flux_z_" + istring);
    for (int j=0; j<mx*my*mz; j++)
    {
      N_i [j] = 1e-16 * inverse_Nunit;
      Fx_i[j] = 1e-16 * inverse_Funit;
      Fy_i[j] = 1e-16 * inverse_Funit;
      Fz_i[j] = 1e-16 * inverse_Funit; 
    }
  }

  if (enzo_config->method_m1_closure_flux_function == "HLL") {
    read_hll_eigenvalues(enzo_config->method_m1_closure_hll_file);
  }

  block->initial_done();
  
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
