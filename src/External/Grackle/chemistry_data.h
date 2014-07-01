/***********************************************************************
/
/ Chemistry data structure
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __CHEMISTRY_DATA_H__
#define __CHEMISTRY_DATA_H__

typedef struct
{

  // Rank of dataset.
  gr_int grid_rank;

  // Dimension of dataset.
  gr_int *grid_dimension;

  // Dataset parameter values.
  gr_float **grid_parameters;

  // Heating values
  gr_float *heating_data;

  // Cooling values
  gr_float *cooling_data;

  // Length of 1D flattened data
  gr_int data_size;

#ifdef CONFIG_USE_CHARM
  void pup (PUP::er &p) {
    p | grid_rank;
    p | data_size;
    if (p.isUnpacking()) {
      grid_dimension = new gr_int[3];
      grid_parameters = new gr_float*[3];
      for (int i=0; i<3; i++) {
	grid_parameters[i] = new gr_float[grid_dimension[i]];
      }
      heating_data = new gr_float[data_size];
      cooling_data = new gr_float[data_size];
    }
    PUParray(p,grid_dimension,3);
    for (int i=0; i<3; i++) {
      PUParray(p,grid_parameters[i],grid_dimension[i]);
    }
    PUParray(p,heating_data,data_size);
    PUParray(p,cooling_data,data_size);

  }
#endif
} cloudy_data;

typedef struct 
{

  // adiabatic index
  gr_float Gamma;

  // HDF5 file containing Cloudy cooling table and UV background
  char *grackle_data_file;


  /****************************************
   *** chemistry and cooling parameters ***
   ****************************************/

  gr_int use_grackle;  // 0) grackle off, 1) grackle on
  gr_int with_radiative_cooling; // include cooling in chemistry solver
                                 // 0) no, 1) yes
  gr_int primordial_chemistry; // 1) HI, HII, HeI, HeII, HeIII, e
                               // 2) + H2, H2I+, H-
                               // 3) + D, D+, HD
  gr_int metal_cooling;        // 0) off, 1) on using Cloudy tables
  gr_int h2_on_dust;  // H2 formation on dust
                      // not well tested, should be left off for now

  // Use a CMB temperature floor.
  gr_int cmb_temperature_floor;

  /* additional H2 chemistry parameters
     best left unchanged. */
  
  gr_int three_body_rate;
  gr_int cie_cooling;
  gr_int h2_optical_depth_approximation;

  /* photo-electric heating from irradiated dust */

  gr_int photoelectric_heating;
  gr_float photoelectric_heating_rate; // in CGS

  /***************************************
   *** radiation background parameters ***
   ***************************************/

  gr_int UVbackground;

  struct UVBtable {
    gr_int Nz;

    gr_float zmin, zmax;    
    gr_float *z;

    gr_float *k24;
    gr_float *k25;
    gr_float *k26;
    gr_float *k27;
    gr_float *k28;
    gr_float *k29;
    gr_float *k30;
    gr_float *k31;

    gr_float *piHI;
    gr_float *piHeII;
    gr_float *piHeI;
  } UVbackground_table;

  gr_float UVbackground_redshift_on;
  gr_float UVbackground_redshift_off;
  gr_float UVbackground_redshift_fullon;
  gr_float UVbackground_redshift_drop;

  gr_int Compton_xray_heating;

  gr_float LWbackground_intensity;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  gr_int LWbackground_sawtooth_suppression;

  /**************************************
   *** primordial chemistry rate data ***
   **************************************/

  gr_int NumberOfTemperatureBins;   
  gr_int CaseBRecombination;
  gr_float TemperatureStart;        // range of temperature in K
  gr_float TemperatureEnd;

  /* 6 species rates */

  gr_float *k1;
  gr_float *k2;
  gr_float *k3;
  gr_float *k4;
  gr_float *k5;
  gr_float *k6;

  /* 9 species rates (including H2) */

  gr_float *k7;
  gr_float *k8;
  gr_float *k9;
  gr_float *k10;
  gr_float *k11;
  gr_float *k12;
  gr_float *k13;
  gr_float *k14;
  gr_float *k15;
  gr_float *k16;
  gr_float *k17;
  gr_float *k18;
  gr_float *k19;
  gr_float *k20;  /* currently not used */
  gr_float *k21;  /* currently not used */
  gr_float *k22;  /* 3-body H2 formation */
  gr_float *k23;  /* H2-H2 dissociation */

  gr_float *k13dd;  /* density dependent version of k13 (collisional H2
                    dissociation); actually 7 functions instead of 1. */

  /* Radiative rates for 6-species (for external field). */

  gr_float k24;
  gr_float k25;
  gr_float k26;

  /* Radiative rates for 9-species (for external field). */

  gr_float k27;
  gr_float k28;
  gr_float k29;
  gr_float k30;
  gr_float k31;

  /* 12 species rates (with Deuterium). */

  gr_float *k50;
  gr_float *k51;
  gr_float *k52;
  gr_float *k53;
  gr_float *k54;
  gr_float *k55;
  gr_float *k56;

  /* H2 formation on dust. */

  gr_int NumberOfDustTemperatureBins;   
  gr_float DustTemperatureStart;        // range of temperature in K
  gr_float DustTemperatureEnd;
  gr_float *h2dust;                     // function of Tgas and Tdust

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */

  gr_float *n_cr_n;
  gr_float *n_cr_d1;
  gr_float *n_cr_d2;

  /********************
   *** cooling data ***
   ********************/

  gr_int ih2co;                     // flag for H2 cooling (0-off/1-on)
  gr_int ipiht;                     // flag for photoionization cooling

  gr_float HydrogenFractionByMass;
  gr_float DeuteriumToHydrogenRatio;
  gr_float SolarMetalFractionByMass;

  /* 6 species rates */

  gr_float *ceHI;                   // collisional excitation rates
  gr_float *ceHeI;
  gr_float *ceHeII;
  gr_float *ciHI;                   // collisional ionization
  gr_float *ciHeI;
  gr_float *ciHeIS;
  gr_float *ciHeII;
  gr_float *reHII;                  // recombination
  gr_float *reHeII1;
  gr_float *reHeII2;
  gr_float *reHeIII;
  gr_float *brem;                   // free-free (Bremsstrahlung)
  gr_float comp;                    // Compton cooling
  gr_float comp_xray;               // X-ray compton heating coefficient
  gr_float temp_xray;               // X-ray compton heating temperature (K)
  gr_float gammah;                  // Photoelectric heating (code units)

  /* radiative rates (external field). */

  gr_float piHI;                    // photo-ionization cooling
  gr_float piHeI;                   //    (no temperature dependance)
  gr_float piHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi_g.F. */

  gr_float *hyd01k;
  gr_float *h2k01;
  gr_float *vibh;
  gr_float *roth;
  gr_float *rotl;

  gr_float *GP99LowDensityLimit;
  gr_float *GP99HighDensityLimit;

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  gr_float *GAHI;
  gr_float *GAH2;
  gr_float *GAHe;
  gr_float *GAHp;
  gr_float *GAel;

  /* 12 species rates (including HD) */

  gr_float *HDlte;
  gr_float *HDlow;
  gr_float *HDcool;

  /* CIE cooling */
  gr_float *cieco;

  /* Gas/grain energy transfer. */
  gr_float *gas_grain;

  // Tabulated mean molecular weight
  gr_float *mu;

  // Primordial cooling data

  cloudy_data cloudy_primordial;

  // Metal cooling data

  cloudy_data cloudy_metal;

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  gr_float cloudy_electron_fraction_factor;

#ifdef CONFIG_USE_CHARM
  void pup (PUP::er &p) {

  CkPrintf("INCOMPLETE: chemistry_data::pup()\n");

  p | Gamma;

  //  char *grackle_data_file;
  p | use_grackle;
  p | with_radiative_cooling;
  p | primordial_chemistry;
  p | metal_cooling;
  p | h2_on_dust;
  p | cmb_temperature_floor;

  /* additional H2 chemistry parameters
     best left unchanged. */
  
  p | three_body_rate;
  p | cie_cooling;
  p | h2_optical_depth_approximation;

  /* photo-electric heating from irradiated dust */

  p | photoelectric_heating;
  p | photoelectric_heating_rate; // in CGS

  /***************************************
   *** radiation background parameters ***
   ***************************************/

  p | UVbackground;

  /* struct UVBtable { */
  /*   p | Nz; */

  /*   p | zmin, zmax;     */
  /*   gr_float *z; */

  /*   gr_float *k24; */
  /*   gr_float *k25; */
  /*   gr_float *k26; */
  /*   gr_float *k27; */
  /*   gr_float *k28; */
  /*   gr_float *k29; */
  /*   gr_float *k30; */
  /*   gr_float *k31; */

  /*   gr_float *piHI; */
  /*   gr_float *piHeII; */
  /*   gr_float *piHeI; */
  /* } UVbackground_table; */

  p | UVbackground_redshift_on;
  p | UVbackground_redshift_off;
  p | UVbackground_redshift_fullon;
  p | UVbackground_redshift_drop;

  p | Compton_xray_heating;

  p | LWbackground_intensity;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  p | LWbackground_sawtooth_suppression;

  /**************************************
   *** primordial chemistry rate data ***
   **************************************/

  p | NumberOfTemperatureBins;   
  p | CaseBRecombination;
  p | TemperatureStart;        // range of temperature in K
  p | TemperatureEnd;

  /* 6 species rates */

  /* gr_float *k1; */
  /* gr_float *k2; */
  /* gr_float *k3; */
  /* gr_float *k4; */
  /* gr_float *k5; */
  /* gr_float *k6; */

  /* /\* 9 species rates (including H2) *\/ */

  /* gr_float *k7; */
  /* gr_float *k8; */
  /* gr_float *k9; */
  /* gr_float *k10; */
  /* gr_float *k11; */
  /* gr_float *k12; */
  /* gr_float *k13; */
  /* gr_float *k14; */
  /* gr_float *k15; */
  /* gr_float *k16; */
  /* gr_float *k17; */
  /* gr_float *k18; */
  /* gr_float *k19; */
  /* gr_float *k20;  /\* currently not used *\/ */
  /* gr_float *k21;  /\* currently not used *\/ */
  /* gr_float *k22;  /\* 3-body H2 formation *\/ */
  /* gr_float *k23;  /\* H2-H2 dissociation *\/ */

  /* gr_float *k13dd;  /\* density dependent version of k13 (collisional H2 */
  /*                   dissociation); actually 7 functions instead of 1. *\/ */

  /* Radiative rates for 6-species (for external field). */

  /* gr_float k24; */
  /* gr_float k25; */
  /* gr_float k26; */

  /* /\* Radiative rates for 9-species (for external field). *\/ */

  /* gr_float k27; */
  /* gr_float k28; */
  /* gr_float k29; */
  /* gr_float k30; */
  /* gr_float k31; */

  /* 12 species rates (with Deuterium). */

  /* gr_float *k50; */
  /* gr_float *k51; */
  /* gr_float *k52; */
  /* gr_float *k53; */
  /* gr_float *k54; */
  /* gr_float *k55; */
  /* gr_float *k56; */

  /* H2 formation on dust. */

  p | NumberOfDustTemperatureBins;   
  p | DustTemperatureStart;        // range of temperature in K
  p | DustTemperatureEnd;
  /* gr_float *h2dust;                     // function of Tgas and Tdust */

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */

  /* gr_float *n_cr_n; */
  /* gr_float *n_cr_d1; */
  /* gr_float *n_cr_d2; */

  /********************
   *** cooling data ***
   ********************/

  p | ih2co;                     // flag for H2 cooling (0-off/1-on)
  p | ipiht;                     // flag for photoionization cooling

  p | HydrogenFractionByMass;
  p | DeuteriumToHydrogenRatio;
  p | SolarMetalFractionByMass;

  /* 6 species rates */

  /* gr_float *ceHI;                   // collisional excitation rates */
  /* gr_float *ceHeI; */
  /* gr_float *ceHeII; */
  /* gr_float *ciHI;                   // collisional ionization */
  /* gr_float *ciHeI; */
  /* gr_float *ciHeIS; */
  /* gr_float *ciHeII; */
  /* gr_float *reHII;                  // recombination */
  /* gr_float *reHeII1; */
  /* gr_float *reHeII2; */
  /* gr_float *reHeIII; */
  /* gr_float *brem;                   // free-free (Bremsstrahlung) */
  p | comp;                    // Compton cooling
  p | comp_xray;               // X-ray compton heating coefficient
  p | temp_xray;               // X-ray compton heating temperature (K)
  p | gammah;                  // Photoelectric heating (code units)

  /* radiative rates (external field). */

  p | piHI;                    // photo-ionization cooling
  p | piHeI;                   //    (no temperature dependance)
  p | piHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi_g.F. */

  /* gr_float *hyd01k; */
  /* gr_float *h2k01; */
  /* gr_float *vibh; */
  /* gr_float *roth; */
  /* gr_float *rotl; */

  /* gr_float *GP99LowDensityLimit; */
  /* gr_float *GP99HighDensityLimit; */

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  /* gr_float *GAHI; */
  /* gr_float *GAH2; */
  /* gr_float *GAHe; */
  /* gr_float *GAHp; */
  /* gr_float *GAel; */

  /* 12 species rates (including HD) */

  /* gr_float *HDlte; */
  /* gr_float *HDlow; */
  /* gr_float *HDcool; */

  /* CIE cooling */
  /* gr_float *cieco; */

  /* Gas/grain energy transfer. */
  /* gr_float *gas_grain; */

  // Tabulated mean molecular weight
  /* gr_float *mu; */

  // Primordial cooling data

  p | cloudy_primordial;

  // Metal cooling data

  p | cloudy_metal;

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  p | cloudy_electron_fraction_factor;

  }
#endif
} chemistry_data;

#endif
