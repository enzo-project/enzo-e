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

// extern "C" void FORTRAN_NAME(calc_dt_30)(
//                   int *rank, int *idim, int *jdim, int *kdim,
//                   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
// 			     hydro_method *ihydro, float *C2,
//                   FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
//                              float *vgz, float *gamma, int *ipfree, float *aye,
//                   float *d, float *p, float *u, float *v, float *w,
// 			     float *dt, float *dtviscous);
 
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

#endif /* ENZO_FORTRAN_HPP */
