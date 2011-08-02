#ifndef ENZO_FORTRAN_HPP
#define ENZO_FORTRAN_HPP

extern "C" void FORTRAN_NAME(calc_dt)
  (int *rank, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz, 
   enzo_float *gamma, int *ipfree, enzo_float *aye,
   enzo_float *d, enzo_float *p, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dt);
 
 
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
   int *diff, int *steepen, int *idual,
   enzo_float *eta1, enzo_float *eta2,
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   enzo_float *standard, int dindex[], int Eindex[],
   int uindex[], int vindex[], int windex[],
   int geindex[], enzo_float *temp,
   int *ncolour, enzo_float *colourpt, int *coloff,
   int colindex[]);

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
   int bxindex[], int byindex[], int bzindex[]);

#endif /* ENZO_FORTRAN_HPP */
