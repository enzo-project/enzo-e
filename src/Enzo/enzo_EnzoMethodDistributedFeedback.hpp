// See LICENSE_CELLO file for license and copyright information

/// @file	enzo_EnzoMethodStarMaker.hpp
/// @author     Andrew Emerick (aemerick11@gmail.com)
/// @date
/// @brief

#ifndef ENZO_ENZO_METHOD_DISTRIBUTED_FEEDBACK
#define ENZO_ENZO_METHOD_DISTRIBUTED_FEEDBACK

class EnzoMethodDistributedFeedback : public Method {

  /// @class   EnzoMethodDistributedFeedback
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Feedback Routines

public:

  EnzoMethodDistributedFeedback();

  /// Destructor
  virtual ~EnzoMethodDistributedFeedback() throw() {};

  /// CHarm++ Pup::able declarations
  PUPable_decl(EnzoMethodDistributedFeedback);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodDistributedFeedback (CkMigrateMessage *m)
    : Method (m)
    {  }

  /// Charm++ Pack / Unpack function
  void pup(PUP::er &p);

  /// Apply the method
  virtual void compute (Block * block) throw();

  void compute_ (Block * block);

  /// name
  virtual std::string name() throw()
  { return "feedback"; }

  // Compute the maximum timestep for this method
  virtual double timestep (Block * block) const throw();

  void add_ionization_feedback(Block * block,
                               double xpos, double ypos, double zpos,
                               const double & s49_tot, const int & willExplode);

  void inject_feedback(Block * block,
                       double xpos, double ypos, double zpos,
                       double m_eject, double E_51, double ke_f,
                       double metal_fraction,
                       enzo_float pvx = -9999.0, enzo_float pvy = -9999.0,
                       enzo_float pvz = -9999.0);

  void convert_momentum( enzo_float *vx, enzo_float *vy, enzo_float *vz, enzo_float *d,
                        const enzo_float &up, const enzo_float &vp, const enzo_float &wp,
                        const int &mx, const int &my, const int &mz,
                        const int &ix, const int &iy, const int &iz, int idir);

  void sum_mass_energy( enzo_float *px, enzo_float *py, enzo_float *pz, enzo_float * d,
                        enzo_float *ge, enzo_float *te,
                        const int &mx, const int& my, const int &mz,
                        const int &ix, const int& iy, const int &iz,
                        double &sum_mass, double &sum_energy, double & sum_ke);


  // AE NOTE: In final version, change metal field to something like
  //          double ** species, and have it contain all of the species fields
  //          that may exist in the simulation (metal, H, He, etc.) such that
  //          mass is deposited consistently over all of these fields. May
  //          need to also provide a M-dimensional array (M = num_species)
  //          containing the fraction of the ejecta mass that should go into
  //          each spcies. For example, if we have metal and multi species = 1:
  //              species   = {metal, HI, HII, HeI, HeII, HeIII, e};
  //              spec_frac = {ejecta_metal_frac,
  //                           0.0,
  //                           (1.0 - ejecta_metal_frac)*0.73,
  //                           0.0,
  //                           0.0,
  //                           (1.0 - ejecta_metal_frac)*(1.0-0.73),
  //                           electron_value};
  //            or soemthing like that to make sure we have proper conservation
  //            of species fields. Maybe just make these all vectors to
  //            make life easy.
  void add_feedback_to_grid( enzo_float * px, enzo_float * py, enzo_float *pz,
                             enzo_float * d, enzo_float *ge, enzo_float *te, enzo_float * metal,
                             const int &mx, const int &my, const int &mz,
                             const int &ix, const int &iy, const int &iz,
                             const double &dxc, const double &dyc, const double &dzc,
                             const double & mass_per_cell, const double & mom_per_cell,
                             const double & therm_per_cell, const double & metal_fraction);

  void compute_coefficients( enzo_float *px, enzo_float *py, enzo_float *pz, enzo_float *d,
                             enzo_float *ge, enzo_float* px_l, enzo_float* py_l, enzo_float *pz_l,
                             enzo_float *d_l, const int &mx, const int &my, const int &mz,
                             const int &ix, const int &iy, const int &iz,
                             double &A, double &B, double &C);


protected:

  double kinetic_fraction_;
  double time_first_sn_;

  int stencil_;
  int stencil_rad_;
  int number_of_feedback_cells_;
  int dual_energy_;

  bool shift_cell_center_;

  bool use_ionization_feedback_;

public:

  void set_shift_cell_center (bool val){ shift_cell_center_ = val;};
  bool get_shift_cell_center (void) {return shift_cell_center_;};

};


#endif /* EnzoMethodDistributedFeedback */
