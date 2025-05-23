c     See LICENSE_ENZO file for license and copyright information

#include "fortran.h"

c=======================================================================
c/////////////////////  SUBROUTINE ZEULER_SWEEP  \\\\\\\\\\\\\\\\\\\\\\\
c     
      subroutine zeuler_sweep(j, d, e, u, v, w, ge, in, jn, kn,
     &     gravity, gr_acc, idual, eta1, eta2,
     &     is, ie, js, je, ks, ke,
     &     gamma, pmin, dmin, dt, dx, dy, dz,
     &     idiff, iflatten, isteepen,
     &     iconsrec, iposrec,
     &     ipresfree,
     &     nsubgrids, lface, rface,
     &     fistart, fiend, fjstart, fjend,
     &     dindex, eindex, geindex,
     &     uindex, vindex, windex, array,
     &     ncolor, colorpt, coloff, colindex,
     &     dls, drs, flatten, pbar,
     &     pls, prs, pslice, ubar,
     &     uls, urs, vls, vrs, gels, gers,
     &     wls, wrs, diffcoef, dslice,
     &     eslice, uslice, vslice, wslice,
     &     df, ef, uf, vf,
     &     wf, grslice, geslice, gef,
     &     ges,
     &     colslice, colf, colls, colrs,
     &     ierror,
     &     ie_error_x,ie_error_y,ie_error_z,num_ie_error
     &     )
c     
c     CONTROL ROUTINE FOR Z-SWEEP DIRECT EULER VERSION OF PPM
c
c  written by: Greg Bryan
c  date:       June, 1994
c  modified1:  RH, May 2002
c
c  PURPOSE:  This routine servers as a wrapper for the eulerian version
c            of PPM that works on a two dimensional slice.  We extract
c            a slice, call INTEULR, R_SOLVER and then EULER.  Note
c            that this provides a natural way to do coarse-grain
c            parallelization on three dimension problems.
c
c  INPUTS:
c    d      - density field
c    dt     - timestep
c    dx,dy,dz - grid spacing
c    e      - total specific energy field
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    eta2   - (dual) selection parameter for total energy (typically ~0.1)
c    gamma  - ideal gas constant
c    ge     - gas specific energy field (used when idual = 1)
c    gravity - gravity flag (0 = off)
c    gr_acc - acceleration due to gravity in this direction
c    idiff  - diffusion flag (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    ie,je,ke - field active zone end index
c    iflatten - flattening flag (0 = off)
c    in,jn,kn - field dimensions
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    is,js,ks - field active zone start index
c    isteepen - steepening flag (0 = off)
c    j      - current slice position in y-direction
c    pmin   - minimum pressure
c    dmin   - minimum density
c    u      - x-velocity field
c    v      - y-velocity field
c    w      - z-velocity field
c
c    nsubgrids - number of subgrids (for fluxes)
c    l,rface  - position of left and right flux surfaces (for each dim)
c    fi,jstart,end - start and end indexes for the flux surfaces (for each dim)
c    d...index - indexes for each field into the flux 'array'
c    colorpt - the color 'array' (contains all color fields)
c    ncolor  - the number of color fields
c    coloff  - offsets into the color array (for each color field)
c
c  OUPUTS:
c    d      - density field
c    e      - total specific energy field
c    u      - x-velocity field
c    v      - y-velocity field
c    w      - z-velocity field
c    array  - array of subgrid fluxes
c    colorpt - color fields
c
c  LOCALS: (passed as temporaries)
c    diffcoef - diffusion coefficient in slice j
c    flatten - ammount of flattening (calculated in calcdiss)
c    dl,rs  - density at left and right edges of each cell
c    dslice - extracted 2d slice of the density   , d
c    el,rs  - total specific energy at left and right edges of each cell
c    eslice - extracted 2d slice of the energy    , e
c    geslice - extracted 2d slice of the gas energy, ge
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    pl,rs  - pressure at left and right edges of each cell
c    pslice - extracted 2d slice of the pressure  , p
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  EXTERNALS:
c    pgas2d - computes pressure from equation of state (on a slice)
c    inteuler - computes the Eulerian left and right states for a slice
c    R_SOLVER - Riemann solver (Lagrangean)
c    euler  - converts the lagrangean Riemann results to eulerian
c             coordinates and then computes the Eulerian fluxes for a slice
c    calcdiss - Calculate dissiptation and flattening coefficients
c
c-----------------------------------------------------------------------
      implicit NONE
#define FORTRAN
#include "fortran_types.h"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, idiff, idual, iflatten, isteepen, ipresfree,
     &        in, jn, kn, is, ie, j, js, je, ks, ke, nsubgrids,
     &        ncolor, coloff(ncolor)
      integer  iconsrec, iposrec
      ENZO_REAL dt, eta1, eta2, gamma, pmin, dmin
      ENZO_REAL d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn)
      ENZO_REAL dx(in),dy(jn),dz(kn)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      ENZO_REAL    array(1), colorpt(1)
      integer iriemannfallback
c
c  define local slices (passed as temps)
c
      ENZO_REAL dls(kn,in),    drs(kn,in),flatten(kn,in),
     &         pbar(kn,in),    pls(kn,in),
     &          prs(kn,in), pslice(kn,in),   ubar(kn,in),
     &          uls(kn,in),    urs(kn,in),    vls(kn,in),
     &          vrs(kn,in),    wls(kn,in),    wrs(kn,in),
     &     diffcoef(kn,in), dslice(kn,in), eslice(kn,in),
     &       uslice(kn,in), vslice(kn,in), wslice(kn,in),
     &     geslice(kn,in),    gef(kn,in)
      ENZO_REAL ges(kn,in)
      ENZO_REAL gels(kn,in),   gers(kn,in),
     &     df(kn,in),     ef(kn,in),grslice(kn,in),
     &     uf(kn,in),     vf(kn,in),     wf(kn,in)
      ENZO_REAL
     &     colslice(kn,in,ncolor),  colf(kn,in,ncolor),
     &     colls(kn,in,ncolor), colrs(kn,in,ncolor)
      integer ierror
      integer ie_error_x(*),ie_error_y(*),ie_error_z(*)
      integer num_ie_error
c
c  parameters
c
      integer zsweep
      parameter (zsweep = 3)
c
c  locals
c
      integer i, ic, idim, k, n, nxz, nyz, nzz, offset
      integer i1, i2
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Compute number of active zones
c
      nxz = ie - is + 1
      nyz = je - js + 1
      nzz = ke - ks + 1

      i1 = 1
      i2 = in
c
c Copy from field to slice
c    (note that we have permuted uvw here and not in the routine calls)
c
      do i=1, in
         do k=1, kn
            dslice(k,i) = d(i,j,k)
            eslice(k,i) = e(i,j,k)
            uslice(k,i) = w(i,j,k)
            vslice(k,i) = u(i,j,k)
            wslice(k,i) = v(i,j,k)
            if (dslice(k,i) < 0.0) then
               ierror = ENZO_ERROR_ZEULER_DSLICE
               print*, 'zeuler_sweep dslice ',i,j,k,dslice(k,i)
               if (num_ie_error .ge. 0) then
                  num_ie_error = num_ie_error + 1
                  ie_error_x(num_ie_error) = i
                  ie_error_y(num_ie_error) = j
                  ie_error_z(num_ie_error) = k
               endif
               return
            end if
            if (eslice(k,i) < 0.0) then
               ierror = ENZO_ERROR_ZEULER_ESLICE
               print*, 'zeuler_sweep eslice ',i,j,k,eslice(k,i)
               if (num_ie_error .ge. 0) then
                  num_ie_error = num_ie_error + 1
                  ie_error_x(num_ie_error) = i
                  ie_error_y(num_ie_error) = j
                  ie_error_z(num_ie_error) = k
               endif
               return
            end if
         enddo
         if (gravity .eq. 1) then
            do k=1, kn
               grslice(k,i) = gr_acc(i,j,k)
            enddo
         endif
         if (idual .eq. 1) then
            do k=1, kn
               geslice(k,i) = ge(i,j,k)
               if (geslice(k,i) < 0.0) then
                  ierror = ENZO_ERROR_ZEULER_GESLICE
                  print*, 'zeuler_sweep geslice ',i,j,k,geslice(k,i)
                  if (num_ie_error .ge. 0) then
                     num_ie_error = num_ie_error + 1
                     ie_error_x(num_ie_error) = i
                     ie_error_y(num_ie_error) = j
                     ie_error_z(num_ie_error) = k
                  endif
                  return
               end if
            enddo
         endif
         do ic=1, ncolor
            do k=1, kn
               colslice(k,i,ic) = 
     &              colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i)
            enddo
         enddo
      enddo
c
c  Compute the pressure on a slice
c
      if (idual .eq. 1) then
         call pgas2d_dual(
     &        dslice, eslice, geslice, pslice, 
     &        uslice, vslice, wslice, eta1, eta2,
     &        kn, in, ks-3, ke+3, i1, i2, gamma, pmin, ierror)
      else
         call pgas2d(
     &        dslice, eslice, pslice, uslice, vslice, wslice,
     &        kn, in, ks-3, ke+3, i1, i2, gamma, pmin, ierror)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &     call calcdiss(
     &     dslice, eslice, uslice, u, v, pslice, dz, dx, dy,
     &     kn, in, jn, ks, ke, i1, i2, j, nyz, zsweep,
     &     in, jn, kn, dt, gamma, idiff, iflatten,
     &     diffcoef, flatten
     &     )
c     
c  Compute Eulerian left and right states at zone edges via interpolation
c
      call inteuler(
     &     dslice, pslice, gravity, grslice, geslice,
     &     uslice, vslice, wslice, dz, flatten,
     &     kn, in, ks, ke, i1, i2, idual, eta1, eta2,
     &     isteepen, iflatten,
     &     iconsrec, iposrec,
     &     dt, gamma, ipresfree,
     &     dls, drs, pls, prs, gels, gers,
     &     uls, urs, vls, vrs, wls, wrs,
     &     ncolor, colslice, colls, colrs, ierror
     &     )
c
c  Compute (Lagrangian part of the) Riemann problem at each zone boundary
c
      call twoshock
     &     (
     &            dls, drs, pls, prs, uls, urs, 
     &            kn, in, ks, ke+1, i1, i2, 
     &            dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )

      iriemannfallback = 1
      
      CALL flux_twoshock (
     &     dslice, eslice, geslice, uslice, vslice, wslice,
     &     dz, diffcoef, 
     &     kn,in,
     &     ks, ke, i1, i2, dt, gamma,
     &     idiff, idual,
     &     eta1, iriemannfallback,
     &     dls, drs, pls, prs, gels, gers, uls, urs,
     &     vls, vrs, wls, wrs, pbar, ubar,
     &     df, ef, uf, vf, wf, gef, ges,
     &     ncolor, colslice, colls, colrs, colf,2)
c
c  Compute Eulerian fluxes and update zones-centered quantities
c
      call euler   (
     &     dslice, eslice, grslice, geslice,
     &     uslice, vslice, wslice, dz, diffcoef,
     &     kn, in, ks, ke, i1, i2, dt, 
     &     gamma, idiff, gravity, idual, eta1, eta2,
     &     df, ef, uf, vf, wf, gef, ges,
     &     ncolor, colslice, colf, dmin, ierror)
c     
c  If necessary, recompute the pressure to correctly set ge and e
c
      if (idual .eq. 1) then
         call pgas2d_dual
     &     (dslice, eslice, geslice, pslice, 
     &     uslice, vslice, wslice, eta1, eta2,
     &        kn, in, ks-3, ke+3, i1, i2, gamma, pmin, ierror)
      endif
c
c  Check this slice against the list of subgrids 
c     (all subgrid quantities are zero based)
c    Note that uf/vf/wf are switched to match u/v/wslice.
c
      do n=0, nsubgrids-1
        if (j .ge. fjstart(n*3+3)+1 .and. j .le. fjend(n*3+3)+1) then
          idim = fiend(n*3+3) - fistart(n*3+3) + 1
          do i=fistart(n*3+3)+1, fiend(n*3+3)+1
             offset = i-fistart(n*3+3) + (j-fjstart(n*3+3)-1)*idim
             array(dindex(n*6+5)+offset) = df(lface(n*3+3)+1, i)
             array(dindex(n*6+6)+offset) = df(rface(n*3+3)+2, i)
             array(eindex(n*6+5)+offset) = ef(lface(n*3+3)+1, i)
             array(eindex(n*6+6)+offset) = ef(rface(n*3+3)+2, i)
             if (nxz .gt. 1) then
                array(uindex(n*6+5)+offset) = vf(lface(n*3+3)+1, i)
                array(uindex(n*6+6)+offset) = vf(rface(n*3+3)+2, i)
             endif
             if (nyz .gt. 1) then
                array(vindex(n*6+5)+offset) = wf(lface(n*3+3)+1, i)
                array(vindex(n*6+6)+offset) = wf(rface(n*3+3)+2, i)
             endif
             array(windex(n*6+5)+offset) = uf(lface(n*3+3)+1, i)
             array(windex(n*6+6)+offset) = uf(rface(n*3+3)+2, i)
             if (idual .eq. 1) then
                array(geindex(n*6+5)+offset) = gef(lface(n*3+3)+1, i)
                array(geindex(n*6+6)+offset) = gef(rface(n*3+3)+2, i)
             endif
             do ic=1, ncolor
                array(colindex(n*6+5,ic)+offset) = 
     &                                       colf(lface(n*3+3)+1, i, ic)
                array(colindex(n*6+6,ic)+offset) = 
     &                                       colf(rface(n*3+3)+2, i, ic)
             enddo
          enddo
        endif
      enddo
c
c Copy from slice to field
c
      do i=i1, i2
         do k=1, kn
            d(i,j,k) = dslice(k,i)
            e(i,j,k) = eslice(k,i)
            w(i,j,k) = uslice(k,i)
            u(i,j,k) = vslice(k,i)
            v(i,j,k) = wslice(k,i)
         enddo
         if (idual .eq. 1) then
            do k=1, kn
               ge(i,j,k) = geslice(k,i)
            enddo
         endif
         do ic=1, ncolor
            do k=1, kn
               colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i) = 
     &               colslice(k,i,ic)
            enddo
         enddo
      enddo
c
      return
      end
