c     $Id$
c     See LICENSE_ENZO file for license and copyright information

#include "fortran.def"

c=======================================================================
c/////////////////////  SUBROUTINE XEULER_SWEEP  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine xeuler_sweep(k, d, e, u, v, w, ge, in, jn, kn,
     &                        gravity, gr_acc, idual, eta1, eta2,
     &                        is, ie, js, je, ks, ke, 
     &                        gamma, pmin, dt, dx, dy, dz,
     &                        idiff, iflatten, isteepen, ipresfree,
     &                        nsubgrids, lface, rface,
     &                        fistart, fiend, fjstart, fjend,
     &                        dindex, eindex, geindex,
     &                        uindex, vindex, windex, array,
     &                        ncolor, colorpt, coloff, colindex,
     &                        dls, drs, flatten, pbar,
     &                        pls, prs, pslice, ubar,
     &                        uls, urs, vls, vrs, gels, gers,
     &                        wls, wrs, diffcoef, dslice,
     &                        eslice, uslice, vslice, wslice,
     &                        df, ef, uf, vf,
     &                        wf, grslice, geslice, gef,
     &                        colslice, colf, colls, colrs
     &                       )
c
c  CONTROL ROUTINE FOR X-SWEEP DIRECT EULER VERSION OF PPM
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
c    gr_acc - acceleration due to gravity (in x-dimension)
c    idiff  - diffusion flag (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    ie,je,ke - field active zone end index
c    iflatten - flattening flag (0 = off)
c    in,jn,kn - field dimensions
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    is,js,ks - field active zone start index
c    isteepen - steepening flag (0 = off)
c    k        - current slice position in direction 3
c    pmin   - minimum pressure
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
c  LOCALS:  (passed as temporaries in argument list)
c    diffcoef - diffusion coefficient in slice k
c    df     - density flux
c    flatten - ammount of flattening (calculated in calcdiss)
c    dl,rs  - density at left and right edges of each cell
c    dslice - extracted 2d slice of the density   , d
c    dx,dy,dz - grid dimension
c    ef     - total energy flux
c    el,rs  - total specific energy at left and right edges of each cell
c    eslice - extracted 2d slice of the energy    , e
c    geslice - extracted 2d slice of the gas energy, ge
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    pl,rs  - pressure at left and right edges of each cell
c    pslice - extracted 2d slice of the pressure  , p
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    uf     - 1-momuntum flux
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vf     - 2-momentum flux
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wf     - 3-momentum flux
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c
c    colf   - color flux
c    coll,rs - color at left and right edges of each cell
c    colslice - extracted 2d slice of the color
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
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, idiff, idual, iflatten, isteepen, ipresfree,
     &        in, jn, kn, is, ie, js, je, ks, ke, k, nsubgrids,
     &        ncolor, coloff(ncolor)
      ENZO_REAL dt, eta1, eta2, gamma, pmin
      ENZO_REAL d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn),
     &        dx,dy,dz
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      ENZO_REAL array(1), colorpt(1)
c
c  define local slices (passed as temps)
c
      ENZO_REAL dls(in,jn),    drs(in,jn),flatten(in,jn),
     &       pbar(in,jn),    pls(in,jn),
     &        prs(in,jn), pslice(in,jn),   ubar(in,jn),
     &        uls(in,jn),    urs(in,jn),    vls(in,jn),
     &        vrs(in,jn),    wls(in,jn),    wrs(in,jn),
     &   diffcoef(in,jn),     df(in,jn),     ef(in,jn),
     &         uf(in,jn),     vf(in,jn),     wf(in,jn),
     &        gef(in,jn),   gels(in,jn),   gers(in,jn)
      ENZO_REAL dslice(in,jn), eslice(in,jn), grslice(in,jn),
     &     uslice(in,jn), vslice(in,jn),  wslice(in,jn), 
     &    geslice(in,jn)
      ENZO_REAL colslice(in,jn,ncolor),  colf(in,jn,ncolor),
     &        colls(in,jn,ncolor), colrs(in,jn,ncolor)
c
c  parameters
c
      integer xsweep
      parameter (xsweep = 1)
c
c  locals
c
      integer i, ic, idim, j, n, nxz, nyz, nzz, offset
      integer j1, j2
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Compute number of active zones
c
      nxz = ie - is + 1
      nyz = je - js + 1
      nzz = ke - ks + 1

      j1 = 1
      j2 = jn

!     j1 = js
!     j2 = je


c
c Copy from field to slice
c
      do j=j1, j2
         do i=1, in
            dslice(i,j) = d(i,j,k)
            eslice(i,j) = e(i,j,k)
            uslice(i,j) = u(i,j,k)
            vslice(i,j) = v(i,j,k)
            wslice(i,j) = w(i,j,k)
         enddo
         if (gravity .eq. 1) then
            do i=1, in
               grslice(i,j) = gr_acc(i,j,k)
            enddo
         endif
         if (idual .eq. 1) then
            do i=1, in
               geslice(i,j) = ge(i,j,k)
            enddo
         endif
         do ic=1, ncolor
            do i=1, in
               colslice(i,j,ic) = 
     &              colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i)
            enddo
         enddo
      enddo
c
c  Compute the pressure on a slice
c
      if (idual .eq. 1) then
         call pgas2d_dual(dslice, eslice, geslice, pslice,
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    in, jn, is-3, ie+3, j1, j2, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               in, jn, is-3, ie+3, j1, j2, gamma, pmin)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &   call calcdiss(
     &            dslice, eslice, uslice, v, w, pslice, dx, dy, dz,
     &            in, jn, kn, is, ie, j1, j2, k, nzz, xsweep,
c     &            in, jn, kn, is, ie, js, je, k, nzz, xsweep,
     &            in, jn, kn, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                )
c
c  Compute Eulerian left and right states at zone edges via interpolation
c
      call inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dx, flatten,
     &            in, jn, is, ie, j1, j2, idual, eta1, eta2,
     &            isteepen, iflatten, dt, gamma, ipresfree,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            ncolor, colslice, colls, colrs
     &             )
c
c  Compute (Lagrangian part of the) Riemann problem at each zone boundary
c
      call twoshock(
     &            dls, drs, pls, prs, uls, urs, 
     &            in, jn, is, ie+1, j1, j2, 
     &            dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )
c
c  Compute Eulerian fluxes and update zones-centered quantities
c
      call euler   (
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dx, diffcoef,
     &            in, jn, is, ie, j1, j2, dt, 
     &            gamma, idiff, gravity, idual, eta1, eta2,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            pbar, ubar,
     &            df, ef, uf, vf, wf, gef,
     &            ncolor, colslice, colls, colrs, colf
     &             )
c
c  If necessary, recompute the pressure to correctly set ge and e
c
      if (idual .eq. 1)
     &    call pgas2d_dual(dslice, eslice, geslice, pslice,
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    in, jn, is-3, ie+3, j1, j2, gamma, pmin)
c
c  Check this slice against the list of subgrids 
c     (all subgrid quantities are zero based)
c
      do n=0, nsubgrids-1
        if (k .ge. fjstart(n*3+1)+1 .and. k .le. fjend(n*3+1)+1) then
          idim = fiend(n*3+1) - fistart(n*3+1) + 1
          do j=fistart(n*3+1)+1, fiend(n*3+1)+1
             offset = j-fistart(n*3+1) + (k-fjstart(n*3+1)-1)*idim
             array(dindex(n*6+1)+offset) = df(lface(n*3+1)+1, j)
             array(dindex(n*6+2)+offset) = df(rface(n*3+1)+2, j)
             array(eindex(n*6+1)+offset) = ef(lface(n*3+1)+1, j)
             array(eindex(n*6+2)+offset) = ef(rface(n*3+1)+2, j)
             array(uindex(n*6+1)+offset) = uf(lface(n*3+1)+1, j)
             array(uindex(n*6+2)+offset) = uf(rface(n*3+1)+2, j)
             if (nyz .gt. 1) then
                array(vindex(n*6+1)+offset) = vf(lface(n*3+1)+1, j)
                array(vindex(n*6+2)+offset) = vf(rface(n*3+1)+2, j)
             endif
             if (nzz .gt. 1) then
                array(windex(n*6+1)+offset) = wf(lface(n*3+1)+1, j)
                array(windex(n*6+2)+offset) = wf(rface(n*3+1)+2, j)
             endif
             if (idual .eq. 1) then
                array(geindex(n*6+1)+offset) = gef(lface(n*3+1)+1, j)
                array(geindex(n*6+2)+offset) = gef(rface(n*3+1)+2, j)
             endif
             do ic=1, ncolor
                array(colindex(n*6+1,ic)+offset) = 
     &                                       colf(lface(n*3+1)+1, j, ic)
                array(colindex(n*6+2,ic)+offset) = 
     &                                       colf(rface(n*3+1)+2, j, ic)
             enddo
          enddo
        endif
      enddo
c
c Copy from slice to field
c
      do j=j1, j2
         do i=1, in
            d(i,j,k) = dslice(i,j)
            e(i,j,k) = eslice(i,j)
            u(i,j,k) = uslice(i,j)
            v(i,j,k) = vslice(i,j)
            w(i,j,k) = wslice(i,j)
         enddo
         if (idual .eq. 1) then
            do i=1, in
               ge(i,j,k) = geslice(i,j)
            enddo
         endif
         do ic=1, ncolor
            do i=1, in
               colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i) = 
     &               colslice(i,j,ic)
            enddo
         enddo
      enddo
c
      return
      end
