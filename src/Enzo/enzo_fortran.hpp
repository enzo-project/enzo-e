extern "C" void FORTRAN_NAME(calc_dt)
  (int *rank, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   ENZO_FLOAT *dx, ENZO_FLOAT *dy, ENZO_FLOAT *dz, 
   float *gamma, int *ipfree, float *aye,
   float *d, float *p, float *u, float *v, float *w,
   float *dt, float *dtviscous);
 
 
extern "C" void FORTRAN_NAME(calc_dt_ppml)
  (int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   ENZO_FLOAT *dx, ENZO_FLOAT *dy, ENZO_FLOAT *dz,
   float *dn, float *vx, float *vy, float *vz, 
   float *bx, float *by, float *bz, 
   float *dt);
 
extern "C" void FORTRAN_NAME(ppm_de)
  (float *d, float *E, float *u, float *v, float *w,
   float *ge,
   int *grav, float *gr_ax, float *gr_ay, float *gr_az,
   float *gamma, float *dt, int *cycle_number,
   float dx[], float dy[], float dz[],
   int *rank, int *in, int *jn, int *kn,
   int is[], int ie[],
   int *flatten, int *ipresfree,
   int *diff, int *steepen, int *idual,
   float *eta1, float *eta2,
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   float *standard, int dindex[], int Eindex[],
   int uindex[], int vindex[], int windex[],
   int geindex[], float *temp,
   int *ncolour, float *colourpt, int *coloff,
   int colindex[]);

extern "C" void FORTRAN_NAME(ppml)
  (float *dn,   float *vx,   float *vy,   float *vz,
   float *bx,   float *by,   float *bz,
   float *dnrx, float *vxrx, float *vyrx, float *vzrx,
   float *bxrx, float *byrx, float *bzrx,
   float *dnry, float *vxry, float *vyry, float *vzry,
   float *bxry, float *byry, float *bzry,
   float *dnrz, float *vxrz, float *vyrz, float *vzrz,
   float *bxrz, float *byrz, float *bzrz,
   float *dt, float dx[], float dy[], float dz[],
   int *in, int *jn, int *kn,
   int is[], int ie[],
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   float *standard, int dnindex[], 
   int vxindex[], int vyindex[], int vzindex[],
   int bxindex[], int byindex[], int bzindex[]);
