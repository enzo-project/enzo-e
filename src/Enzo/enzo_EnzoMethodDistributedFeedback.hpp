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

  void compute_ (Block * block) throw();

  /// name
  virtual std::string name() throw()
  { return "distributed_feedback"; }

  // Compute the maximum timestep for this method
  virtual double timestep (Block * block) const throw();

  void convert_momentum( double *vx, double *vy, double *vz,
                        const double &up, const double &vp, const double &wp,
                        const int &nx const int &ny, const int &nz,
                        const int &ix, const int &iy, const int &iz, int idir);

  void sum_mass_energy( double *px, double *py, double *pz, double * d,
                        double *ge, double *te,
                        const int &nx, const int& ny, const int &nz,
                        const int &ix, const int& iy, const int &iz,
                        double &sum_mass, double &sum_energy, double & sum_ke);

  void add_feedback_to_grid( double * px, double * py, double *pz,
                             double * d, double *ge, double *te,
                             const int &nx, const int &ny, const int &nz,
                             const int &ix, const int &iy, const int &iz,
                             const double &dxc, const double &dyc, const double &dzc,
                             const double & mass_per_cell, const double & mom_per_cell,
                             const double & therm_per_cell);

  void compute_coefficients( double *px, double *py, double *pz, double *d,
                             double *ge, double* px_l, double* py_l, double *pz_l,
                             double *d_l, const int &nx, const int &ny, const int &nz,
                             const int &ix, const int *iy, const int &iz,
                             double &A, double &B, double &C);

protected:

  double total_ejecta_mass_;
  double total_ejecta_energy_;
  double ejecta_metal_fraction_;

  int stencil_size_;
  int istencil_;
  int number_of_feedback_cells_;

  double m_eject_;

};


#endif /* EnzoMethodDistributedFeedback */
