// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovRandom.cpp
/// @author   Thomas Bolden (boldenth@gmail.com)
/// @date     Sat Mar  25 14:30:00 EST 2017
/// @brief    [\ref Enzo] Randomized Sedov Blasts initial conditions
///
/// This problem is designed to stress-test load balancing. Each block
/// has a variable number of supernovae with variable energy. 
/// Optionally, about half the blocks can be made empty to increase load
/// imbalance between blocks.

#include "cello.hpp"

#include "enzo.hpp"

#include <random>

//----------------------------------------------------------------------

EnzoInitialSedovRandom::EnzoInitialSedovRandom
(const EnzoConfig * config) throw ()
: Initial(config->initial_cycle, config->initial_time) 
{
  array_[0]        = config->initial_sedov_random_array[0];
  array_[1]        = config->initial_sedov_random_array[1];
  array_[2]        = config->initial_sedov_random_array[2];
  half_empty_      = config->initial_sedov_random_half_empty; // new
  grackle_cooling_ = config->initial_sedov_random_grackle_cooling; // new
  max_blasts_      = config->initial_sedov_random_max_blasts; // new
  radius_relative_ = config->initial_sedov_random_radius_relative;
  pressure_in_     = config->initial_sedov_random_pressure_in;
  pressure_out_    = config->initial_sedov_random_pressure_out;
  density_         = config->initial_sedov_random_density;
  te_multiplier_   = config->initial_sedov_random_te_multiplier; // new
}

//----------------------------------------------------------------------

void EnzoInitialSedovRandom::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | half_empty_;
  p | grackle_cooling_;
  p | max_blasts_;
  p | radius_relative_;
  p | pressure_in_;
  p | pressure_out_;
  p | density_;
  p | te_multiplier_;
  
}

//----------------------------------------------------------------------

void EnzoInitialSedovRandom::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()

{

  ASSERT("EnzoInitialSedovRandom",
   "Block does not exist",
   block != NULL);

  Field field = block->data()->field();

  ASSERT("EnzoInitialSedovRandom",
   "Insufficient number of fields",
   field.field_count() >= 4);

  enzo_float *  d = (enzo_float *) field.values
    (field.field_id("density"));
  
  enzo_float * te = (enzo_float *) field.values
    (field.field_id("total_energy"));

  // set up species fields
  // if(grackle_cooling_)

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double xmb,ymb,zmb;
  block->data()->lower(&xmb,&ymb,&zmb);

  double xpb,ypb,zpb;
  block->data()->upper(&xpb,&ypb,&zpb);

  double hx,hy,hz;
  field.cell_width(xmb,xpb,&hx,
                   ymb,ypb,&hy,
                   zmb,zpb,&hz);

  // Parameters

  const double sedov_radius = radius_relative_/array_[0];
  const double sedov_radius_2 = sedov_radius*sedov_radius;

  const int in = cello::index_static();

  const double sedov_te_in = 
    pressure_in_  / ((EnzoBlock::Gamma[in] - 1.0) * density_);
  const double sedov_te_out= 
    pressure_out_ / ((EnzoBlock::Gamma[in] - 1.0) * density_);

  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  int mx = nx + 2*gx;
  int my = ny + 2*gy;
  int mz = nz + 2*gz;

  // clear all fields

  for (int iv=0; iv<field.field_count(); iv++) {

    enzo_float * array = (enzo_float *) field.values (iv);

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<mx; ix++) {

          int i = INDEX(ix,iy,iz,mx,my);
          array[i] = 0.0;

        }
      }
    }

  }

  // background 

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<mx; ix++) {

        int i = INDEX(ix,iy,iz,mx,my);
        d[i]  = density_;
        te[i] = sedov_te_out;

        //if(grackle_cooling_){
        //  hi[i] = 0.75;
        //}

      }
    }
  }

  // array of explosions

  double xmd,ymd,zmd;
  hierarchy->lower(&xmd,&ymd,&zmd);

  double xpd,ypd,zpd;
  hierarchy->upper(&xpd,&ypd,&zpd);

  // bounds of possible explosions intersecting this Block

  double r = sedov_radius;
  
  int kxm = MAX((int)floor((xmb-xmd-r)/(xpd-xmd)*array_[0])-1,0);
  int kym = MAX((int)floor((ymb-ymd-r)/(ypd-ymd)*array_[1])-1,0);
  int kzm = MAX((int)floor((zmb-zmd-r)/(zpd-zmd)*array_[2])-1,0);
  int kxp = MIN( (int)ceil((xpb-xmd+r)/(xpd-xmd)*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil((ypb-ymd+r)/(ypd-ymd)*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil((zpb-zmd+r)/(zpd-zmd)*array_[2])+1,array_[2]);
  
  TRACE3 ("SEDOV: %d %d %d",kxp,kyp,kzp);

  double hxa = (xpd-xmd) / array_[0];
  double hya = (ypd-ymd) / array_[1];
  double hza = (zpd-zmd) / array_[2];

  // (kx,ky,kz) index bounds of explosions in domain

  // setup an array of random numbers to fill each block in array
  // add check for half to be zero
  // also need an array of multiple of sedov_te_in for each blast in this array

  // debugging
  // max_blasts_ = 1;

  int rn_array[kxp][kyp][kzp][max_blasts_+1];

  int alternate = 0;

  std::random_device rd;
  std::uniform_int_distribution<int> rand_blasts(1, max_blasts_);
  std::uniform_int_distribution<int> rand_te_mult(-te_multiplier_,
                                                  +te_multiplier_);
  std::uniform_real_distribution<double> rand_radius(-sedov_radius,
                                                     -sedov_radius+1.0);
  
  
  for(int iz = 0; iz<kzp; iz++){
    for(int iy = 0; iy<kyp; iy++){
      for(int ix = 0; ix<kxp; ix++){

        int rnum = rand_blasts(rd); // should give number in [1,max_blasts_]

        if(half_empty_ && alternate%2){ // every odd one should be empty
          rnum = 0;
        }

        rn_array[ix][iy][iz][0] = rnum;

        // debugging
        // std::cout << "rnum: " << rnum << "\n";

        for(int rnbeg = 1; rnbeg<=rnum; rnbeg++){
          rn_array[ix][iy][iz][rnbeg] = rand_te_mult(rd);  // adding te multiplier
        }

        alternate++;

      }
    }
  }

  // now rn_array is a 4 dimensional array such that:
  // rn_array[x][y][z][r]
  // x - x index of blast box
  // y - y index of blast box
  // z - z index of blast box
  // r - r=0 is number of blasts, rest are total energy multipliers for each blast

  for (int kz=kzm; kz<kzp; kz++) {

    double zc = hza*(0.5+kz);

    for (int ky=kym; ky<kyp; ky++) {

      double yc = hya*(0.5+ky);

      for (int kx=kxm; kx<kxp; kx++) {

        double xc = hxa*(0.5+kx);
        TRACE3("xc,yc,zc = %f %f %f",xc,yc,zc);

        // so now (xc, yc, zc) is set to center of block

        int blasts_current_box = rn_array[kx][ky][kz][0];

        for (int blast_num = 1; blast_num <= blasts_current_box; blast_num++){

          //std::cout << "Executing " << blast_num << " of " << blasts_current_box << std::endl;

          // loops over every cell including ghost cells
          // if cell is within sedov_radius, add total energy to cell

          double rz = rand_radius(rd);
          double ry = rand_radius(rd);
          double rx = rand_radius(rd);

          double snx = xc - 0.5 * hxa + sedov_radius + fabs(rx) * hxa; //xc + rx * hxa;
          double sny = yc - 0.5 * hya + sedov_radius + fabs(ry) * hya;
          double snz = zc - 0.5 * hza + sedov_radius + fabs(rz) * hza;

          // debugging
          // std::cout << "(snx, sny, snz) of blast " << blast_num << ": " 
          //          << "(" << snx << ", " << sny << ", " << snz << ")\n";

          // (snx, sny, snz) are coords of explosion center

          for (int iz=0; iz<mz; iz++) {
            
            double z = zmb + (iz - gz + 0.5)*hz;

            for (int iy=0; iy<my; iy++) {
              
              double y = ymb + (iy - gy + 0.5)*hy;

              for (int ix=0; ix<mx; ix++) {
                
                double x = xmb + (ix - gx + 0.5)*hx;

                double r2 = (snx - x) * (snx - x) + (sny - y) * (sny - y) + (snz - z) * (snz - z); //x*x + y*y + z*z;

                // cell center (x,y,z)

                int i = INDEX(ix,iy,iz,mx,my);

                if (r2 < sedov_radius_2) {
                  te[i] = pow(2,rn_array[kx][kz][ky][blast_num]) * sedov_te_in;
                }

              }
            }
          } // end for over all cells

        } // end for over blast_num
      }
    }
  }

  //delete [] rn_array;
  
}

