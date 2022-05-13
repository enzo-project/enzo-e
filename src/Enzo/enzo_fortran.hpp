#ifndef ENZO_FORTRAN_HPP
#define ENZO_FORTRAN_HPP

extern "C" void FORTRAN_NAME(interpolate)
  (int *rank, enzo_float *pfield, int pdim[],
   int pis[], int pie[], int r[],
   enzo_float *field, int dim[], int is[], enzo_float *work,
   int *imethod, int *posflag,
   int *ierror);

extern "C" void FORTRAN_NAME(interp3d)
  (enzo_float * parent, enzo_float * work, 
   int *dim1, int *dim2, int *dim3, 
   int *start1, int *start2, int *start3,
   int *end1, int *end2, int *end3, 
   int *refine1, int *refine2, int *refine3, 
   enzo_float *grid,
   int *gdim1, int *gdim2, int *gdim3, 
   int *gstart1, int *gstart2, int *gstart3,
   int *wdim1, int *wdim2, int *wdim3,
   int *ierror);
  

extern "C" void FORTRAN_NAME(calc_dt)
  (int *rank, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz, 
   enzo_float *gamma, int *ipfree, enzo_float *aye,
   enzo_float *d, enzo_float *p, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dt);
 
extern "C" void FORTRAN_NAME(expand_terms)(
   int *rank, int *isize, int *idual, enzo_float *coef,
   int *imethod, enzo_float *gamma,
   enzo_float *p,  enzo_float *d, enzo_float *e, enzo_float *ge,
   enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dold, enzo_float *eold, enzo_float *geold,
   enzo_float *uold, enzo_float *vold, enzo_float *wold,
   int *icr, enzo_float *ecr, enzo_float *ecrold);

extern "C" void FORTRAN_NAME(turboinit)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v, enzo_float *w,
   int *in, int *jn, int *kn,
   int *ig, int *jg, int *kg);

extern "C" void FORTRAN_NAME(turboinit2d)
  (int *rank, int *nbox,
   enzo_float *u, enzo_float *v,
   int *in, int *jn,
   int *ig, int *jg);

extern "C" void FORTRAN_NAME(cello_init_turbulence_ou)
  (int * is_root,
   int * rank,
   double fomain_size[],
   double * gamma,
   int * apply_injection_rate,
   int * cooling_term,
   double * injection_rate,
   double * kmin,
   double * kmax,
   double * mach,
   int * read_sol,
   double * sol_weight);

extern "C" int FORTRAN_NAME(cello_turbou_state_size)
  (int * n_buffer_real, int * n_buffer_int);

extern "C" void FORTRAN_NAME(cello_get_turbou_state)
  (double * buffer_real, int * buffer_int);

extern "C" void FORTRAN_NAME(cello_put_turbou_state)
  (double * buffer_real, int * buffer_int);

extern "C" void FORTRAN_NAME(turbforceou)
  (int * mx, int * my, int * mz,
   int * ni, int * nj, int * nk,
   double * field_density, double * grid,
   double * wk, double * time, double * dt,
   int * cello_update_sol,
   int * cello_apply_cooling,
   int * cello_apply_forcing,
   int * cello_apply_injection_rate,
   int * cello_update_phases,
   int * cello_cooling_term,
   double * cello_gamma,
   double * cello_injection_rate,
   int * olap,
   double * r_gv
   );
 
extern "C" void FORTRAN_NAME(turbforceshift)
  (int * mx, int * my, int * mz,
   int * ni, int * nj, int * nk,
   double * field_density,
   double * field_momentum_x,
   double * field_momentum_y,
   double * field_momentum_z,
   double * field_jacobian,
   double * wk,
   int * cello_update_sol,
   int * cello_apply_injection_rate,
   int * cello_olap,
   double * cello_injection_rate,
   double * r_gv,
   double * r_av);
 
extern "C" void FORTRAN_NAME(turbforceupdate)
  (int * mx, int * my, int * mz,
   int * ni, int * nj, int * nk,
   double * field_density,
   double * field_momentum_x,
   double * field_momentum_y,
   double * field_momentum_z,
   double * field_energy,
   double * resid_density,
   double * resid_momentum_x,
   double * resid_momentum_y,
   double * resid_momentum_z,
   double * resid_energy,
   double * field_temperature,
   double * wk, double * dt,
   double * turbAcc,
   int * cello_update_sol,
   int * cello_apply_injection_rate,
   double * cello_injection_rate,
   int * cello_cooling_term,
   int * cello_apply_cooling,
   double * cello_gamma,
   double * cello_hc_alpha,
   double * cello_hc_sigma,
   double * cello_totemp,
   double * r_av );
 
extern "C" void FORTRAN_NAME(calc_dt_ppml)
  (int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz,
   enzo_float *dn, enzo_float *vx, enzo_float *vy, enzo_float *vz, 
   enzo_float *bx, enzo_float *by, enzo_float *bz, 
   enzo_float *dt);
 
extern "C" void FORTRAN_NAME(ppm_de)
  (enzo_float *d, enzo_float *E, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *ge,
   int *grav, enzo_float *gr_ax, enzo_float *gr_ay, enzo_float *gr_az,
   enzo_float *gamma, enzo_float *dt, int *cycle_number,
   enzo_float dx[], enzo_float dy[], enzo_float dz[],
   int *rank, int *in, int *jn, int *kn,
   int is[], int ie[],
   int *flatten, int *ipresfree,
   int * iconsrec, int *iposrec,
   int *diff, int *steepen, int *idual,
   enzo_float *eta1, enzo_float *eta2,
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   enzo_float *standard, int dindex[], int Eindex[],
   int uindex[], int vindex[], int windex[],
   int geindex[], enzo_float *temp,
   int *ncolor, enzo_float *colorpt, int *coloff,
   int colindex[], int *error,
   int *ie_error_x,
   int *ie_error_y,
   int *ie_error_z,
   int *num_ie_error
   );

extern "C" void FORTRAN_NAME(ppml)
  (enzo_float *dn,   enzo_float *vx,   enzo_float *vy,   enzo_float *vz,
   enzo_float *bx,   enzo_float *by,   enzo_float *bz,
   enzo_float *dnrx, enzo_float *vxrx, enzo_float *vyrx, enzo_float *vzrx,
   enzo_float *bxrx, enzo_float *byrx, enzo_float *bzrx,
   enzo_float *dnry, enzo_float *vxry, enzo_float *vyry, enzo_float *vzry,
   enzo_float *bxry, enzo_float *byry, enzo_float *bzry,
   enzo_float *dnrz, enzo_float *vxrz, enzo_float *vyrz, enzo_float *vzrz,
   enzo_float *bxrz, enzo_float *byrz, enzo_float *bzrz,
   enzo_float *dt, enzo_float dx[], enzo_float dy[], enzo_float dz[],
   int *in, int *jn, int *kn,
   int is[], int ie[],
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   enzo_float *standard, int dnindex[], 
   int vxindex[], int vyindex[], int vzindex[],
   int bxindex[], int byindex[], int bzindex[],
   enzo_float *f1,enzo_float *f2,enzo_float *f3,enzo_float *f4,
   enzo_float *f5,enzo_float *f6,enzo_float *f7,
   enzo_float *g1,enzo_float *g2,enzo_float *g3,enzo_float *g4,
   enzo_float *g5,enzo_float *g6,enzo_float *g7,
   enzo_float *h1,enzo_float *h2,enzo_float *h3,enzo_float *h4,
   enzo_float *h5,enzo_float *h6,enzo_float *h7,
   enzo_float *ex,enzo_float *ey,enzo_float *ez,
   enzo_float *qu1,enzo_float *qu2,enzo_float *qu3,enzo_float *qu4,
   enzo_float *qu5,enzo_float *qu6,enzo_float *qu7);

extern "C" void FORTRAN_NAME(OUpumpInit)
      ( enzo_float *gamma_, 
    	enzo_float *density_initial_,
    	enzo_float *pressure_initial_,
    	enzo_float *solenoidal_fraction_, 
    	enzo_float *mach_number_, 
    	enzo_float *kfmin_, 
    	enzo_float *kfmax_,
    	enzo_float *Lbox );

extern "C" void FORTRAN_NAME(OUpumpCompute)
  ( int *rank, int *mx,  int *my,  int *mz,         // rank and local block dimensions
    int *nig,  int *njg, int *nkg,              // root 
    int *gx,   int *gy,  int *gz,               // number of ghost zones
    enzo_float *vx,            // flow fields invilved
    enzo_float *vy,
    enzo_float *vz,
    enzo_float *density,
    int *ndx,int *ndy,int *ndz,                 // zone sizes
    int *ox,int *oy,int *oz,
    enzo_float *dt );                          // time step

extern "C" void FORTRAN_NAME(ppml_ig)
  (enzo_float *density,
   enzo_float *velox, enzo_float *veloy, enzo_float *veloz,
   enzo_float *bfieldx, enzo_float *bfieldy, enzo_float *bfieldz,
   enzo_float *pressure,
   enzo_float *dens_rx,
   enzo_float *velox_rx,enzo_float *veloy_rx,enzo_float *veloz_rx,
   enzo_float *bfieldx_rx,enzo_float *bfieldy_rx,enzo_float *bfieldz_rx,
   enzo_float *press_rx,
   enzo_float *dens_ry,
   enzo_float *velox_ry,enzo_float *veloy_ry,enzo_float *veloz_ry,
   enzo_float *bfieldx_ry,enzo_float *bfieldy_ry,enzo_float *bfieldz_ry,
   enzo_float *press_ry,
   enzo_float *dens_rz,
   enzo_float *velox_rz,enzo_float *veloy_rz,enzo_float *veloz_rz,
   enzo_float *bfieldx_rz,enzo_float *bfieldy_rz,enzo_float *bfieldz_rz,
   enzo_float *press_rz,
   enzo_float *b0, enzo_float *gamma,
   enzo_float *dt, enzo_float *hx, enzo_float *hy, enzo_float *hz,
   int *mx, int *my, int *mz,
   int *GridStartIndex, int *GridEndIndex,
   int *NumberOfSubgrids, int *leftface, int *rightface,
   int *istart, int *iend, int *jstart, int *jend,
   enzo_float *standard, int *dnindex,
   int *vxindex, int *vyindex, int *vzindex,
   int *bxindex, int *byindex, int *bzindex,
   int *pindex,
   enzo_float *f1,enzo_float *f2,enzo_float *f3,enzo_float *f4,
   enzo_float *f5,enzo_float *f6,enzo_float *f7,enzo_float *f8,
   enzo_float *g1,enzo_float *g2,enzo_float *g3,enzo_float *g4,
   enzo_float *g5,enzo_float *g6,enzo_float *g7,enzo_float *g8,
   enzo_float *h1,enzo_float *h2,enzo_float *h3,enzo_float *h4,
   enzo_float *h5,enzo_float *h6,enzo_float *h7,enzo_float *h8,
   enzo_float *ex,enzo_float *ey,enzo_float *ez,
   enzo_float *qu1,enzo_float *qu2,enzo_float *qu3,enzo_float *qu4,
   enzo_float *qu5,enzo_float *qu6,enzo_float *qu7,enzo_float *qu8);

extern "C" void FORTRAN_NAME(calc_dt_ppml_ig)
  (int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz,
   enzo_float *dn, enzo_float *vx, enzo_float *vy, enzo_float *vz, 
   enzo_float *bx, enzo_float *by, enzo_float *bz,
   enzo_float *pr, enzo_float *b0, enzo_float * gamma,
   enzo_float *dt);
#endif /* ENZO_FORTRAN_HPP */
