// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHydro.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of EnzoMethodHydro class

#ifndef ENZO_ENZO_METHOD_HYDRO_HPP
#define ENZO_ENZO_METHOD_HYDRO_HPP

extern "C" void FORTRAN_NAME(woc_pgas2d_dual)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *geslice,
   enzo_float *pslice,
   enzo_float *uslice, enzo_float *vslice, enzo_float *wslice,
   enzo_float *eta1, enzo_float *eta2,
   int *idim, int *jdim,
   int *i1, int *i2, int *j1, int *j2,
   enzo_float *gamma, enzo_float *pmin);

extern "C" void FORTRAN_NAME(woc_pgas2d)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *pslice,
   enzo_float *uslice, enzo_float *vslice, enzo_float *wslice,
   int *idim, int *jdim,
   int *i1, int *i2, int *j1, int *j2,
   enzo_float *gamma, enzo_float *pmin);

extern "C" void FORTRAN_NAME(woc_calcdiss)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *uslice, enzo_float *v, 
   enzo_float *w, enzo_float *pslice,
   enzo_float *dx, enzo_float *dy, enzo_float *dz, 
   int *idim, int *jdim, int *kdim, int *i1, int *i2, int *j1, 
   int *j2, int *k, int *nzz, int *idir, int *dimx, int *dimy, 
   int *dimz, enzo_float *dt, enzo_float *gamma, int *idiff, int *iflatten, 
   enzo_float *diffcoef, enzo_float *flatten);

extern "C" void FORTRAN_NAME(woc_inteuler)
  (
   enzo_float *dslice, enzo_float *pslice, int *gravity, enzo_float *grslice, 
   enzo_float *geslice, enzo_float *uslice, enzo_float *vslice, enzo_float *wslice, 
   enzo_float *dxi, enzo_float *flatten, int *idim, int *jdim, int *i1, 
   int *i2, int *j1, int *j2, int *idual, enzo_float *eta1, 
   enzo_float *eta2, int *isteep, int *iflatten,
   int *iconsrec, int *iposrec,
   enzo_float *dt, enzo_float *gamma, int *ipresfree, enzo_float *dls, enzo_float *drs, 
   enzo_float *pls, enzo_float *prs, enzo_float *gels, enzo_float *gers, enzo_float *uls, 
   enzo_float *urs, enzo_float *vls, enzo_float *vrs, enzo_float *wls, enzo_float *wrs, 
   int *ncolor, enzo_float *colslice, enzo_float *colls, enzo_float *colrs);

extern "C" void FORTRAN_NAME(woc_twoshock)
  (
   enzo_float *dls, enzo_float *drs, enzo_float *pls, enzo_float *prs, 
   enzo_float *uls, enzo_float *urs, int *idim, int *jdim, int *i1,
   int *i2, int *j1, int *j2, enzo_float *dt, enzo_float *gamma, 
   enzo_float *pmin, int *ipresfree, enzo_float *pbar, enzo_float *ubar,
   int *gravity, enzo_float *grslice, int *idual, enzo_float *eta1);

extern "C" void FORTRAN_NAME(woc_flux_twoshock)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *geslice, enzo_float *uslice, 
   enzo_float *vslice, enzo_float *wslice, enzo_float *dx,
   enzo_float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
   int *j1, int *j2, enzo_float *dt, enzo_float *gamma, int *idiff,
   int *idual, enzo_float *eta1, int *ifallback,
   enzo_float *dls, enzo_float *drs, enzo_float *pls, enzo_float *prs, 
   enzo_float *gels, enzo_float *gers, 
   enzo_float *uls, enzo_float *urs, enzo_float *vls, enzo_float *vrs, 
   enzo_float *wls, enzo_float *wrs, enzo_float *pbar, enzo_float *ubar,
   enzo_float *df, enzo_float *ef, enzo_float *uf, enzo_float *vf, enzo_float *wf, enzo_float *gef, 
   enzo_float *ges,
   int *ncolor, enzo_float *colslice, enzo_float *colls, enzo_float *colrs, enzo_float *colf);

extern "C" void FORTRAN_NAME(woc_flux_hll)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *geslice, enzo_float *uslice, 
   enzo_float *vslice, enzo_float *wslice, enzo_float *dx,
   enzo_float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
   int *j1, int *j2, enzo_float *dt, enzo_float *gamma,
   int *idiff, int *idual, enzo_float *eta1, int *ifallback,
   enzo_float *dls, enzo_float *drs, enzo_float *pls, enzo_float *prs, 
   enzo_float *uls, enzo_float *urs, enzo_float *vls, enzo_float *vrs, 
   enzo_float *wls, enzo_float *wrs, enzo_float *gels, enzo_float *gers, 
   enzo_float *df, enzo_float *ef, enzo_float *uf, enzo_float *vf, enzo_float *wf, enzo_float *gef, 
   enzo_float *ges,
   int *ncolor, enzo_float *colslice, enzo_float *colls, enzo_float *colrs, enzo_float *colf);

extern "C" void FORTRAN_NAME(woc_flux_hllc)
  (
   enzo_float *dslice, enzo_float *eslice, enzo_float *geslice, enzo_float *uslice, 
   enzo_float *vslice, enzo_float *wslice, enzo_float *dx,
   enzo_float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
   int *j1, int *j2, enzo_float *dt, enzo_float *gamma,
   int *idiff, int *idual, enzo_float *eta1, int *ifallback,
   enzo_float *dls, enzo_float *drs, enzo_float *pls, enzo_float *prs, 
   enzo_float *uls, enzo_float *urs, enzo_float *vls, enzo_float *vrs, 
   enzo_float *wls, enzo_float *wrs, enzo_float *gels, enzo_float *gers, 
   enzo_float *df, enzo_float *ef, enzo_float *uf, enzo_float *vf, enzo_float *wf, enzo_float *gef, 
   enzo_float *ges,
   int *ncolor, enzo_float *colslice, enzo_float *colls, enzo_float *colrs, enzo_float *colf);

extern "C" void FORTRAN_NAME(woc_euler)
  (
   enzo_float *dslice, enzo_float *pslice, enzo_float *grslice, 
   enzo_float *geslice, enzo_float *uslice, enzo_float *vslice, enzo_float *wslice, 
   enzo_float *dx, enzo_float *diffcoef, int *idim, int *jdim, int *i1, 
   int *i2, int *j1, int *j2, enzo_float *dt, enzo_float *gamma, int *idiff, 
   int *gravity, int *idual, enzo_float *eta1, enzo_float *eta2, enzo_float *df, 
   enzo_float *ef, enzo_float *uf, enzo_float *vf, enzo_float *wf, enzo_float *gef,
   enzo_float *ges,
   int *ncolor, enzo_float *colslice, enzo_float *colf, enzo_float *dfloor);


class EnzoMethodHydro : public Method {

  /// @class    EnzoMethodHydro
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate ENZO's hydro methods

public: // interface

  /// Create a new EnzoMethodHydro object
  EnzoMethodHydro(const FieldDescr * field_descr,
		  std::string method,
		  enzo_float gamma,
		  bool gravity,
		  bool comoving_coordinates,
		  bool dual_energy,
		  enzo_float dual_energy_eta1,
		  enzo_float dual_energy_eta2,
		  std::string reconstruct_method,
		  bool reconstruct_conservative,
		  bool reconstruct_positive,
		  enzo_float ppm_density_floor,
		  enzo_float ppm_pressure_floor,
		  int ppm_pressure_free,
		  int ppm_diffusion,
		  int ppm_flattening,
		  int ppm_steepening,
		  std::string riemann_solver);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodHydro);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodHydro (CkMigrateMessage *m)
    : Method (m),
      method_(""),
      gamma_(0.0),
      gravity_(0),
      comoving_coordinates_(false),
      dual_energy_(0),
      dual_energy_eta1_(0.0),
      dual_energy_eta2_(0.0),
      reconstruct_method_(""),
      reconstruct_conservative_(0),
      reconstruct_positive_(0),
      ppm_density_floor_(0.0),
      ppm_pressure_floor_(0.0),
      ppm_pressure_free_(0),
      ppm_diffusion_(0),
      ppm_flattening_(0),
      ppm_steepening_(0),
      riemann_solver_("")
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
public: // virtual methods

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "hydro"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  void ppm_method_ (Block * block);
  void ppm_euler_x_ (Block * block, int iz);
  void ppm_euler_y_ (Block * block, int ix);
  void ppm_euler_z_ (Block * block, int iy);
  
protected: // attributes

  /// Which hydro method to call
  std::string method_;

  /// Gamma
  enzo_float gamma_;

  /// whether gravity is enabled (int for fortran call)
  int gravity_;

  /// comoving coordinates
  bool comoving_coordinates_;

  /// whether to use dual energy (int for fortran call)
  int dual_energy_;

  /// dual energy eta
  enzo_float dual_energy_eta1_;
  enzo_float dual_energy_eta2_;

  /// Flux reconstruction method
  std::string reconstruct_method_;

  /// Whether to enforce conservative flux reconstruction
  int reconstruct_conservative_;
  
  /// Whether to enforce positive flux reconstruction
  int reconstruct_positive_;

  /// minimum density
  enzo_float ppm_density_floor_;

  /// minimum pressure
  enzo_float ppm_pressure_floor_;

  /// Whether to assume pressure-free
  int ppm_pressure_free_;

  /// diffusion parameter
  int ppm_diffusion_;

  /// flattening parameter
  int ppm_flattening_;

  /// steepening parameter
  int ppm_steepening_;

  /// Riemann solver to use
  std::string riemann_solver_;

};
  
#endif /* ENZO_ENZO_METHOD_HYDRO_HPP */
