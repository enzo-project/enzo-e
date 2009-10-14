*****************************************************************************
*                                                                           *
* Copyright 2004 Greg Bryan                                                 *
* Copyright 2004 Laboratory for Computational Astrophysics                  *
* Copyright 2004 Board of Trustees of the University of Illinois            *
* Copyright 2004 Regents of the University of California                    *
*                                                                           *
* This software is released under the terms of the "Enzo Public License"    *
* in the accompanying LICENSE file.                                         *
*                                                                           *
*****************************************************************************

#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE CALCDIFF  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calcdiss(
     &            dslice, eslice, uslice, v, w, pslice, dx, dy, dz,
     &            idim, jdim, kdim, i1, i2, j1, j2, k, nzz, idir,
     &            dimx, dimy, dimz, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                    )
c
c  COMPUTES THE DIFFUSION COEFFICIENTS AND FLATTENING FOR PPM (CW84)
c
c  written by: Greg Bryan
c  date:       July, 1994
c  modified1:  Robert Harkness
c  date:       July, 2003
c              Mods to prevent out-of-bounds memory references with
c              diffusion.  Works with new euler sweep range mod by
c              Alexei Kritsuk to at least allow for flattening.
c              Calculation for idiff = 2 was incorrect in any case
c  modified2:  James Bordner
c  date:       2003-11-21
c              Fixed minor index bug; moved some loop invariants outside
c              inner loops; general housecleaning
c  modified3:  Alexei Kritsuk
c  date:       2004-08-10
c              Fixed iflatten = 1 procedure (A.2), shock capturing (A.1);
c              cleaned and enabled iflatten = 3 (A.7-10)
c
c  PURPOSE:  Given the values of idiff and iflatten, we calculate
c     the diffusion coefficients and amount of local flattening using
c     one of a number of routines.  They are:
c     a) idiff = 0  -> no difussion
c     b) idiff = 1  -> simple diffusion
c       Computes the diffusion coefficient given in Colella &
c       Woodward eq. A.3 for use later when calculating the diffusive
c       flux.  The essential logic is given in 4.4:
c
c            diffcoef =  K max|-Div u(j+1/2), 0|
c
c       where Div is the discrete undivided difference approximation to
c       the multidimensional divergence of u.  This version calculates
c       the diffusion coefficient for a single slice (k).
c     c) idiff = 2  -> alternate diffusion
c       Here we use the scheme described in CW84, eq. A11.  Although,
c       this is intended for local grid motion, here we adapt it for
c       the local amount of explicit diffusion.
c     d) iflatten = 0 -> no flattening
c     e) iflatten = 1 -> The simple flattening scheme described in eqs
c        A1 & A2 of CW84.
c     f) iflatten = 2 -> Flattening needed for Lagrangean hydrodynamics
c        (note: not L+remap).  Stuck in for completeness, but not tested.
c        This is described in equations. A.4 - A.6.
c     g) iflatten = 3 -> Flattening recomended for multidimensional
c        calculations.  Described in CW84 A7-A10.
c     NOTE: this routine is intended to be called for x,y and z directions,
c        cyclically permuting the arguments to obtain the correct result
c        for each direction (i.e. call with dx,dy,dz,... for x-dir, dy,dz,dx
c        for y-dir, and dz,dx,dy for z-dir).  Also, insure idir is set
c        correctly, since a few calculation must know which sweep
c        direction this is.  
c           permute: v,w,dx,dy,dz,i1,i2,j1,j2,k,nzz,idim,jdim,kdim
c           do not permute: dimx,dimy,dimz
c
c  INPUT ARGUMENTS:
c     d      - density field
c     dimx,y,z - dimensions of v and w field arrays
c     dslice - single slice of the density in (sweep, orthodim) direction
c     dx     - distance between Eulerian zone edges in direction 1
c     dy,dz  - distance between Eulerian zone edges in 2,3 directions
c     eslice - single slice of the total energy in (sweep, orthodim) direction
c     i1,i1  - active range in direction-1
c     idiff  - type of explicit diffusion to be used (see above)
c     idir   - direction of sweep (1=x, 2=y, 3=z)
c     idim   - first dimension of slices (direction-1)
c     iflatten - type of local flattening to be used (see above)
c     j2,j2  - active range in direction-2
c     jdim   - second dimension of slices (direction-2)
c     k      - current slice index
c     kdim   - dimension of direction-3
c     nzz    - number of active zones in the z-direction
c     pslice - single slice of the pressure in (sweep, orthodim) direction
c     uslice - single slice of the pressure in (sweep, orthodim) direction
c     v      - 2-velocity field
c     w      - 3-velocity field
c
c  OUTPUT ARGUMENTS:
c     diffcoef - diffusion coefficient in slice k
c     flatten  - amount of local flattening (0-1)
c
c  LOCALS:
c     di     - temporary for 1/dslice
c     wflag  - CW84 eq. A1
c
c  PARAMETERS
c     K      - free parameter in diffusion coefficient calculation [0.1]
c     (for other parameters see CW84)
c
c-----------------------------------------------------------------------
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer dimx, dimy, dimz, i1, i2, idiff, idim, idir, iflatten, 
     &        j1, j2, jdim, k, kdim, nzz
      real    gamma, dt, dx(idim), dy(jdim), dz(kdim)
      real    dslice(idim,jdim), eslice(idim,jdim), uslice(idim,jdim),
     &        diffcoef(idim,jdim), flatten(idim,jdim), pslice(idim,jdim)
      real    v(dimx,dimy,dimz), w(dimx,dimy,dimz)
c
c  locals
c
      integer i, j, nyz, is, isp, ism
      real    cj2s, di(ijkn), difftemp(ijkn), flattemp(ijkn),
     &        kappa(ijkn), kappa_tilde, omega(ijkn),
     &        qa, qb, qc, s, sigma(ijkn),
     &        sigma_tilde, vdiff1(ijkn), vdiff2(ijkn), wdiff1(ijkn),
     &        wdiff2(ijkn), wflag(ijkn), Z, ZE
      logical ljx1,ljy1,ljz1,ljx2,ljy2,ljz2,lj1,lj2
      logical lkx1,lky1,lkz1,lkx2,lky2,lkz2,lk1,lk2
      real dp1,dp2,de1,de2,dee,dpp,omega_tilde
c
c  parameters
c
      real    epsilon, kappa1, kappa2, Kparam, nu1, nu2, nu3,
     &        omega1, omega2, sigma1, sigma2
      parameter (epsilon = 0.33  , kappa1 = 2.0 , kappa2 =  0.01,
     &            Kparam = 0.1   ,    nu1 = 2.0 , nu2    =  0.1 ,
     &               nu3 = 0.2333, omega1 = 0.75, omega2 = 10.0 ,
     &            sigma1 = 0.5   , sigma2 = 1.0)
c
c Notice, for 1D calculations nu2 = 0 and nu3 = 0.3333, see CW84, p. 200.
c
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c  Compute number of active zones in transverse direction
c

      nyz = j2 - j1 + 1
c
c  Precompute logical masks for k loop bounds
c
      lk1  = (1.lt.k .and. k.lt.kdim)
      lkx1 = (1.lt.k .and. k.lt.dimx)
      lky1 = (1.lt.k .and. k.lt.dimy)
      lkz1 = (1.lt.k .and. k.lt.dimz)

      lk2  = (2.lt.k .and. k.lt.kdim-1)
      lkx2 = (2.lt.k .and. k.lt.dimx-1)
      lky2 = (2.lt.k .and. k.lt.dimy-1)
      lkz2 = (2.lt.k .and. k.lt.dimz-1)
c
c  Loop across the k slice
c
      do j=j1, j2
c
c  Precompute logical masks for j loop bounds
c
         lj1  = (1.lt.j .and. j.lt.jdim)
         ljx1 = (1.lt.j .and. j.lt.dimx)
         ljy1 = (1.lt.j .and. j.lt.dimy)
         ljz1 = (1.lt.j .and. j.lt.dimz)

         lj2  = (2.lt.j .and. j.lt.jdim-1)
         ljx2 = (2.lt.j .and. j.lt.dimx-1)
         ljy2 = (2.lt.j .and. j.lt.dimy-1)
         ljz2 = (2.lt.j .and. j.lt.dimz-1)
c
c  All (well, almost all) routines below need this quantity (CW84, eq. a1)
c
         do i=i1-2, i2+2
            qb =   abs(pslice(i+1,j) - pslice(i-1,j))  ! AK abs()
     &           / min(pslice(i+1,j),  pslice(i-1,j))
            if ( qb.gt.epsilon .and.
     &          uslice(i-1,j).gt.uslice(i+1,j) ) then
               wflag(i) = 1.0
            else
               wflag(i) = 0.0
            endif
         end do
c
c  Diffusion routine (B) needs this calculation, which depends on v & w,
c     so it depends on the sweep direction.
c     
         if (idiff .eq. 1) then
c     
            if (idir .eq. 1) then

               if (dimy.gt.1 .and. ljy1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(i,j-1,k) + v(i-1,j-1,k)) 
     &                    -      (v(i,j+1,k) + v(i-1,j+1,k))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               end if
c     
               if (dimz.gt.1 .and. lkz1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(i,j,k-1) + w(i-1,j,k-1))
     &                    -      (w(i,j,k+1) + w(i-1,j,k+1))
                  end do
               else 
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               end if
c     
c     
            elseif (idir .eq. 2) then
c     
c     
               if (dimz.gt.1 .and. ljz1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(k,i,j-1) + v(k,i-1,j-1))
     &                    -      (v(k,i,j+1) + v(k,i-1,j+1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               endif
c     
               if (dimx.gt.1 .and. lkx1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(k-1,i,j) + w(k-1,i-1,j))
     &                    -      (w(k+1,i,j) + w(k+1,i-1,j))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               endif
c
c     
            elseif (idir .eq. 3) then
c
c     
c     
               if (dimx.gt.1 .and. ljx1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(j-1,k,i) + v(j-1,k,i-1)) 
     &                    -      (v(j+1,k,i) + v(j+1,k,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               end if
c     
               if (dimy.gt.1 .and. lky1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(j,k-1,i) + w(j,k-1,i-1)) 
     &                    -      (w(j,k+1,i) + w(j,k+1,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               endif
c     
            end if
c
         endif
c
c  Diffusion routine (C) needs this calculation, which depends on v & w,
c    so it depends on the sweep direction.
c     
         if (idiff .eq. 2) then
c     
            if (idir .eq. 1) then
c     
               if (dimy.gt.1 .and. ljy2) then
                  do i=i1-2, i2+2
                     vdiff2(i) = (v(i,j-2,k) + v(i-1,j-2,k))
     &                    -      (v(i,j+2,k) + v(i-1,j+2,k))
                  end do
               else
                  do i=i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               endif
c     
               if (dimz.gt.1 .and. lkz2) then
                  do i=i1-2, i2+2
                     wdiff2(i) = (w(i,j,k-2) + w(i-1,j,k-2)) 
     &                    -      (w(i,j,k+2) + w(i-1,j,k+2))
                  end do
               else
                  do i=i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               endif
c     
c     
            elseif (idir .eq. 2) then
c     
               if (dimz.gt.1 .and. ljz2) then
                  do i=i1-2, i2+2
                     vdiff2(i) = (v(k,i,j-2) + v(k,i-1,j-2)) 
     &                    -      (v(k,i,j+2) + v(k,i-1,j+2))
                  end do
               else
                  do i=i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               endif
c     
               if (dimx.gt.1 .and. lkx2) then
                  do i=i1-2, i2+2
                     wdiff2(i) = (w(k-2,i,j) + w(k-2,i-1,j)) 
     &                    -      (w(k+2,i,j) + w(k+2,i-1,j))
                  end do
               else
                  do i=i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               endif
c     
c
            elseif (idir .eq. 3) then
c     
               if (dimx.gt.1 .and. ljx2) then
                  do i = i1-2, i2+2
                     vdiff2(i) = (v(j-2,k,i)+v(j-2,k,i-1)) 
     &                    -      (v(j+2,k,i)+v(j+2,k,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               end if
c     
               if (dimy.gt.1 .and. lky2) then
                  do i = i1-2, i2+2
                     wdiff2(i) = (w(j,k-2,i)+w(j,k-2,i-1)) 
     &                    -      (w(j,k+2,i)+w(j,k+2,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               end if
c     
c                
            endif
c
         endif
c  
c-----------------------------------------------------------------------
c (B)  Compute diffusion coefficient
c
         if (idiff .eq. 1) then
c
            do i=i1, i2+1
               diffcoef(i,j) = uslice(i-1,j) - uslice(i,j)
            end do
c     
c
            if (nyz.gt.1 .and. lj1) then
               do i=i1, i2+1
                  diffcoef(i,j) = diffcoef(i,j) + (0.25*(dx(i)+dx(i-1))
     &                 /   (0.5*(dy(j+1)+dy(j-1)) + dy(j)))*vdiff1(i)
               end do
            endif
c     
            if (nzz.gt.1 .and. lk1) then
               do i=i1, i2+1
                  diffcoef(i,j) = diffcoef(i,j) + (0.25*(dx(i)+dx(i-1))
     &                 /   (0.5*(dz(k+1)+dz(k-1)) + dz(k)))*wdiff1(i)
               end do
            endif
c     
            do i=i1, i2+1
               diffcoef(i,j) = Kparam*max(0.0, diffcoef(i,j))
            end do
c
         endif
c
c-----------------------------------------------------------------------
c (E)  Construct flattening parameter (eqns A1 and A2), also 68-74
c
         if (iflatten .eq. 1) then

            do i=i1-1, i2+1
               if (abs(pslice(i+2,j) - pslice(i-2,j))/
     &            min(pslice(i+2,j),pslice(i-2,j)) .lt. epsilon) then
                 qa = 1.0
               else
                 qa =     ( pslice(i+1,j) - pslice(i-1,j) )
     &                /   ( pslice(i+2,j) - pslice(i-2,j) )
               endif
               flattemp(i) = min( 1.0, (qa-omega1)*omega2*wflag(i) )
c               flattemp(i) = wflag(i) !AK: Svetsov 2001, ShockWaves 11, 229
               flattemp(i) = max(0.0, flattemp(i))
c              flattemp(i) = 1.0-flattemp(i)  !AK: typo in CW84 eq. (A.2)
             end do

            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)

c     Now, choose the maximum (eq. A2, first part)

            do i=i1-1, i2+1
               if (pslice(i+1,j) - pslice(i-1,j) .lt. 0.0) then
                  flatten(i,j) = max(flattemp(i),flattemp(i+1))
               else
                  flatten(i,j) = max(flattemp(i),flattemp(i-1))
               endif
c        There are some additional recipes used by others  !AK
c               flatten(i,j) = 0.025 !Runacres+1 astro-ph/0405315v1, p.7;
c               flatten(i,j) = 1.0 !Godunov interpolation (constant)
c               flatten(i,j) = max(flatten(i,j), 0.1)
            end do

         endif

c-----------------------------------------------------------------------
c (F) Construct second type flattening parameter (eq. A4-A6)
c
         if (iflatten .eq. 2) then

            do i=i1-3, i2+3
               di(i) = 1.0/dslice(i,j)
            end do

            do i=i1-1, i2+1

               is = i + sign(2.0, pslice(i+1,j) - pslice(i-1,j))

               omega(i) = max(0.0, omega1 * (omega2 
     &              - (pslice(i+1,j) - pslice(i-1,j))
     &              / (pslice(i+2,j) - pslice(i-2,j)) ))

               Z = sqrt(( max(pslice(i+2,j),pslice(i-2,j)) +
     &              0.5*(pslice(i+2,j)+pslice(i-2,j))*(gamma-1.0))
     &              / max(di(i+2),di(i-2)))

               kappa_tilde = 
     &              (Z + sqrt(gamma*pslice(is,j)*dslice(is,j))) / Z

               kappa(i)   = max(0.0, (kappa_tilde - kappa1)
     &              /                (kappa_tilde + kappa2))

               flattemp(i) = min(wflag(i)*omega(i), kappa(i))

            end do

            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)

c     Choose minimum (eq. A6)

            do i=i1-1, i2+1
               flatten(i,j) = 
     &              max(flattemp(i-1),flattemp(i),flattemp(i+1))
            end do

         endif
c-----------------------------------------------------------------------
c (G) Construct third type flattening parameter (eq. A7-A9)
c
         if (iflatten .eq. 3 .or. idiff .eq. 2) then
c     
            do i=i1-3, i2+3
               di(i) = 1.0/dslice(i,j)
            end do
c     
            do i=i1-1, i2+1
c     
               dp1 = pslice(i+1,j) - pslice(i-1,j)
               dp2 = pslice(i+2,j) - pslice(i-2,j)
               de1 = eslice(i+1,j) - eslice(i-1,j)
               de2 = eslice(i+2,j) - eslice(i-2,j)
c
               dpp = 0.0
               dee = 0.0
               if (dp2 .ne. 0.0) dpp = dp1/dp2
               if (de2 .ne. 0.0) dee = de1/de2
               omega_tilde = max( dpp , dee )
c
               ism = i + sign(2.0, dp1)        ! post-shock
               isp = i - sign(2.0, dp1)        ! upstream
               s   =   - sign(1.0, dp1)        ! i+s - zone upstream from i
               if (dp1 .eq. 0.0) s = 0.0
c
c              Compute sigma(i): strength of a shock near zone i
c     
               sigma_tilde = 
     &              wflag(i)*abs(dp2)/min(pslice(i+2,j),pslice(i-2,j))
               sigma(i) = max( 0.0, (sigma_tilde - sigma1)
     &              /               (sigma_tilde + sigma2) )
c
c              Compute omega(i): steepness of a shock near zone i 
c
               omega(i) = 
     &              max(0.0, omega2*(omega_tilde - omega1)) ! as in Prometheus
c     &              max(0.0, omega1*(omega2 - omega_tilde)) !AK: a CW84 typo?
c     
c              Compute an estimate of the Lagrangean shock speed W_j
c
               Z = sqrt(( max(pslice(i+2,j),pslice(i-2,j)) +
     &              0.5*(pslice(i+2,j)+pslice(i-2,j))*(gamma-1.0))
     &              / max(di(i+2),di(i-2)))
c
c              Compute W_j^E; Possible typos in kappa_tilda definition in 
c              CW84 (A.9). !AK
c
               ZE    = s*Z/dslice(ism,j) + uslice(ism,j) + tiny
               cj2s = sqrt(gamma*pslice(isp,j)/dslice(isp,j))
c
c              Compute kappa(i): the wavelength of the noise
c
               kappa_tilde = abs((ZE - uslice(isp,j) + s*cj2s)/ZE)
               kappa(i) = max(0.0, (kappa_tilde - kappa1)
     &              /              (kappa_tilde + kappa2))
c     
               flattemp(i) =
     &              min(kappa(i),wflag(i)*omega(i),wflag(i)*sigma(i))
c     
            end do
c     
         endif
c     
         if (iflatten .eq. 2 .or. iflatten .eq. 3) then
c     
            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)
c     
c     Now, choose the maximum (eq. A.2, first part; also A.10)
c     
            do i=i1-1, i2+1
               if (pslice(i+1,j) - pslice(i-1,j) .lt. 0.0) then
                  flatten(i,j) = max(flattemp(i),flattemp(i+1))
               else
                  flatten(i,j) = max(flattemp(i),flattemp(i-1))
               endif
            end do
c     
         endif
c     
c-----------------------------------------------------------------------
c (C)  Compute complex diffusion coefficient
c     

         if (idiff .eq. 2) then

            do i=i1-1, i2+1
               difftemp(i) = uslice(i-2,j) - uslice(i+1,j)
            end do
            
            if (nyz.ne.1 .and. lj2) then
               do i=i1-1, i2+1
                  qa = 0.5*(0.5*(dx(i+1)+dx(i-2))+dx(i)+dx(i-1))
     &                 / (0.5*(dy(j+2)+dy(j-2)) + dy(j+1)+dy(j)+dy(j-1))
                  difftemp(i) = difftemp(i) + qa*vdiff2(i)
               end do
            end if

            if (nzz.ne.1 .and. lk2) then
               do i=i1-1, i2+1
                  qb = 0.5*(0.5*(dx(i+1)+dx(i-2))+dx(i)+dx(i-1))
     &                 / (0.5*(dz(k+2)+dz(k-2)) + dz(k+1)+dz(k)+dz(k-1))
                  difftemp(i) = difftemp(i) + qb*wdiff2(i)
               end do
            end if

            do i=i1-1, i2+1
               difftemp(i) = max(0.0, nu1*difftemp(i))
            end do

            difftemp(i1-2) = difftemp(i1-1)
            difftemp(i2+2) = difftemp(i2+1)

            do i=i1-1, i2+1

               s  = - sign(1.0, pslice(i-1,j) - pslice(i,j))
               diffcoef(i,j) = s*min(
     &              max(difftemp(i-1),difftemp(i),difftemp(i+1)),
     &              (nu2 + nu3*max(sigma(i-1)*kappa(i-1)**2
     &              ,              sigma(i  )*kappa(i  )**2 ))
     &              * dx(i+1)*s/dt)

            end do

         endif
c
c  Done
c
      end do
c
      return
      end

c=======================================================================
c////////////////////////  SUBROUTINE CALC_DT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
c
      subroutine calc_dt(rank, idim, jdim, kdim, 
     &                   i1, i2, j1, j2, k1, k2, 
     &                   dx, dy, dz, gamma, ipfree, aye,
     &                   d, p, u, v, w, dt, dtviscous)
c
c  COMPUTES TIME STEP FOR NEXT CYCLE
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:  Alexei Kritsuk, Jan. 2001; changed the expression for
c                 hydro time step for PPM_DE. Now it follows the linear 
c                 stability condition for multidimensional Godunov scheme 
c                 (Godunov 1959). It is now safe to use CourantSafetyNumber =
c                 0.8 for PPM_DE runs.
c
c  PURPOSE:  Computes the new timestep using the Courant condition.
c            (For rank < 3, the unused fields and cell widths may be
c             null)
c
c  INPUTS:
c    rank    - rank of fields
c    i,j,dim - declared dimensions of fields
c    i,j,k1  - start index of active region in fields (0 based)
c    i,j,k2  - end index of active region in fields (0 based)
c    dx,y,z  - cell widths along each dimension
c    gamma   - ratio of specific heats
c    ipfree  - pressure free flag (1 = on, 0 = off)
c    aye     - expansion factor (or 1 if not using comoving coordinates)
c    d,p     - density and pressure fields
c    u,v,w   - velocity fields (x,y,z)
c    dtviscous - viscous time for stability (if used)
c
c  OUTPUTS:
c    dt      - minimum allowed dt (without Courant safety factor)
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer idim, jdim, kdim, i1, i2, j1, j2, k1, k2, rank, ipfree
      REALSUB dx(idim), dy(jdim), dz(kdim)
      real    dt, gamma, aye, dtviscous
      real    d(idim,jdim,kdim), p(idim,jdim,kdim), u(idim,jdim,kdim),
     &        v(idim,jdim,kdim), w(idim,jdim,kdim)
c
c     Locals
c
      integer i,j,k
      real    cs, dt1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Set initial timestep to a large number
c
      dt = huge
c
c     one-dimensional version
c
      if (rank .eq. 1) then
         do i = i1+1, i2+1
            cs = max(sqrt(gamma*p(i,1,1)/d(i,1,1)), tiny)
            if (ipfree .eq. 1) cs = tiny
            dt = min(dt, real(dx(i)*aye/(cs + abs(u(i,1,1)))))
         enddo
      endif
c
c     two-dimensional version
c
      if (rank .eq. 2) then
         do j = j1+1, j2+1
            do i = i1+1, i2+1
               cs = max(sqrt(gamma*p(i,j,1)/d(i,j,1)), tiny)
               if (ipfree .eq. 1) cs = tiny

c
c 	Here is Godunov's formula (to make sure ppm works with 0.8)
c
 	       dt = min(dt, real(aye/((cs + abs(u(i,j,1)))/dx(i)+
     &                       (cs + abs(v(i,j,1)))/dy(j)) ))
c
c       The way it was originally in enzo
c
c               dt = min(dt, real(dx(i)*aye/(cs + abs(u(i,j,1)))),
c     &                      real(dy(j)*aye/(cs + abs(v(i,j,1)))))
            enddo
c
         enddo
      endif
c
c     three-dimensional version
c
      if (rank .eq. 3) then
         do k = k1+1, k2+1
            do j = j1+1, j2+1
               do i = i1+1, i2+1
                  if (d(i,j,k) .ne. d(i,j,k) .or.
     &                p(i,j,k) .ne. p(i,j,k))
     &               write(6,*) 'calc_dt',d(i,j,k),p(i,j,k),i,j,k
                  cs = max(sqrt(gamma*p(i,j,k)/d(i,j,k)), tiny)
                  if (ipfree .eq. 1) cs = tiny

c     
c     Godunov's formula.
c     
                  dt1 = aye/((cs + abs(u(i,j,k)))/dx(i) +
     &                 (cs + abs(v(i,j,k)))/dy(j) +
     &                 (cs + abs(w(i,j,k)))/dz(k))
c     
c     The way it was originally in enzo
c     
c     dt1 = min(dx(i)*aye/(cs + abs(u(i,j,k))),
c     &                    dy(j)*aye/(cs + abs(v(i,j,k))),
c     &                    dz(k)*aye/(cs + abs(w(i,j,k))))

                  dt = min(dt, dt1)
!                  if (dt1 .lt. 1.0e-5) write(6,1000) dt1,d(i,j,k),
!    &                  p(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
 1000             format('calc_dt (small dt): dt,d,p,uvw=',1p,6e12.3)
               enddo
c
            enddo
         enddo
      endif
c
      return
      end


c=======================================================================
c//////////////////////////  SUBROUTINE EULER  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine euler(
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dx, diffcoef,
     &            idim, jdim, i1, i2, j1, j2, dt, 
     &            gamma, idiff, gravity, idual, eta1, eta2,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            pbar, ubar,
     &            df, ef, uf, vf, wf, gef,
     &            ncolor, colslice, colls, colrs, colf
     &                 )
c
c  SOLVES THE EULERIAN CONSERVATION LAWS USING FLUXES FROM THE RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  June, 1993 by Greg Bryan (Lagrange+remap --> Eulerian)
c  modified2:  July, 1994 by GB; switched to slicewise and put all the
c                                information into the argument line
c  modified3:  July, 1994 by GB; moved diffusion coefficient to own routine
c  modified4:  Febuary, 2005 by Alexei Kritsuk; fixed a bug in RAREFACTION1
c                            and a bug in update eq. 3.1 with idiff != 0.
c  modified5:  Sept. 2005 by DC, fixed the flux array to work with cosmology and AMR
c
c  PURPOSE:  Updates the conservation laws in Eulerian form using
c    fluxes in the sweep-direction computed by the Riemann solver.  This
c    versions works on a single two dimensional slice.  It also adds
c    diffusive fluxes, if requested.
c
c  INPUT:
c    diffcoef - diffusion coefficient in slice k
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dl,rs  - density at left and right edges of each cell
c    dx     - distance between Eulerian zone edges in sweep direction
c    eslice - extracted 2d slice of the energy, e
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    eta2   - (dual) selection parameter for total energy (typically ~0.1)
c    gamma  - parameter in ideal gas law
c    geslice - extracted 2d slice of the gas energy, ge
c    gravity - gravity flag (0 = off)
c    grslice - acceleration in this dimension in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idiff  - integer flag for standard artificial diffusion (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  OUTPUT:
c    dslice - extracted 2d slice of the density, d
c    geslice - extracted 2d slice of the gas energy, ge
c    eslice - extracted 2d slice of the energy, e
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  LOCALS:
c
#define NO_GRAVITY_SECOND_ORDER_CORRECTION
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, i1, i2, idiff, idim, idual, j1, j2, jdim, ncolor
      real    dt, eta1, eta2, gamma, velold
      real diffcoef(idim,jdim),  dslice(idim,jdim),      dx(idim     ),
     &       eslice(idim,jdim), grslice(idim,jdim), geslice(idim,jdim),
     &       uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim)
      real      dls(idim,jdim),     drs(idim,jdim),
     &         pbar(idim,jdim),     pls(idim,jdim),    gels(idim,jdim),
     &          prs(idim,jdim),    ubar(idim,jdim),    gers(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim),     vls(idim,jdim),
     &          vrs(idim,jdim),     wls(idim,jdim),     wrs(idim,jdim),
     &           df(idim,jdim),      ef(idim,jdim),      uf(idim,jdim),
     &           vf(idim,jdim),      wf(idim,jdim),     gef(idim,jdim)
      real colslice(idim,jdim,ncolor),   colls(idim,jdim,ncolor),
     &        colrs(idim,jdim,ncolor),    colf(idim,jdim,ncolor)
c
c  local declarations
c
      integer i, j, n
      real      alpha      ,    c0(ijkn),    cb(ijkn),  cbar(ijkn),
     &             d0(ijkn),    db(ijkn),  dbar(ijkn),   dnu(ijkn),
     &            dub(ijkn),  duub(ijkn),  duvb(ijkn),  duwb(ijkn),
     &           dueb(ijkn),    eb(ijkn),    l0(ijkn),  lbar(ijkn),
     &             p0(ijkn),   geb(ijkn), dugeb(ijkn),    pb(ijkn),
     &             qa      ,    qb      ,    qc      ,    sn(ijkn),
     &             u0(ijkn),    ub(ijkn),   upb(ijkn),    vb(ijkn),
     &             wb(ijkn),    z0(ijkn),
     &           uold(ijkn), pcent(ijkn),eratio(ijkn),  frac      ,
     &           dadx(ijkn),  dddx(ijkn), colb(ijkn,MAX_COLOR)
      real dele, delu, uleft, uright, uadvect, eold, tempu
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     write(6,*) 'EULER: dt =',dt,' idiff =',idiff
c     write(6,*) 'EULER: idim =',idim,' jdim =',jdim
c     write(6,*) 'EULER: i1   =',i1,  ' i2   =',i2
c     write(6,*) 'EULER: j1   =',j1,  ' j2   =',j2
c
c  Set constants
c
      qa = (gamma + 1.0)/(2.0*gamma)
      qb = (gamma - 1.0)/(gamma + 1.0)
c
c  Loop over sweep lines (in this slice)
c
      do j=j1, j2
c
c  Evaluate time-averaged quantities 
c   (see Colella, Siam J Sci Stat Comput 1982, 3, 77.  Appendix)
c
       do i=i1, i2+1
c
          sn(i)    = sign(1.0, -ubar(i,j))
c
c  Collect values of interest depending on which way fluid is flowing
c
          if (sn(i) .lt. 0.0) then
             u0(i) = uls(i,j)
             p0(i) = pls(i,j)
             d0(i) = dls(i,j)
c             vb(i) = vls(i,j)
c             wb(i) = wls(i,j)
c             geb(i) = gels(i,j)
          else
             u0(i) = urs(i,j)
             p0(i) = prs(i,j)
             d0(i) = drs(i,j)
c             vb(i) = vrs(i,j)
c             wb(i) = wrs(i,j)
c             geb(i) = gers(i,j)
          endif
          c0(i)    = sqrt(max(gamma*p0(i)/d0(i), tiny))
          z0(i)    = c0(i)*d0(i)*sqrt(max(
     &               1.0 + qa*(pbar(i,j)/p0(i)-1.0), tiny))
c          write(24,*) i, u0(i), uls(i,j), urs(i,j), c0(i)
       enddo
c
c      Repeat for colour variables (now moved below)
c
#ifdef UNUSED
       do n=1,ncolor
          do i=i1, i2+1
             if (sn(i) .lt. 0.0) then
                colb(i,n) = colls(i,j,n)
             else
                colb(i,n) = colrs(i,j,n)
             endif
          enddo
       enddo
#endif
c
c  Compute equivalent bar (inside shock & rarefaction) values for density
c    and sound speed
c
       do i=i1, i2+1
          dbar(i)  = 1.0/(1.0/d0(i) - (pbar(i,j)-p0(i))/
     &                                max(z0(i)**2, tiny))
          cbar(i)  = sqrt(max(gamma*pbar(i,j)/dbar(i), tiny))
       enddo
c
c
c  Find lambda values for the shock and rarefaction
c
       do i=i1, i2+1
         if (pbar(i,j) .lt. p0(i)) then
            l0(i)   = u0(i)*sn(i) + c0(i)
            lbar(i) = sn(i)*ubar(i,j)+cbar(i)
         else
            l0(i)   = u0(i)*sn(i) + z0(i)/d0(i)
c            l0(i)   = ubar(i,j)*sn(i) + z0(i)/d0(i)
            lbar(i) = l0(i)
         endif
c
c  Compute values for inside a rarefaction fan
c     (method described in Colella, 1982)
c
#define RAREFACTION2 //AK
c
#ifdef RAREFACTION0
c
c  (quasi-)isothermal hydrodynamics => use RAREFACTION0
c  sound speed is assumed to be 1.0 in expression for db(i)
c  (see eqs (3.11) - (3.16) in Balsara (1994) ApJ 420, p. 200)
c
         ub(i)    = -cbar(i)*sn(i)
         db(i)    = d0(i)*exp(-sn(i)*u0(i) - 1.0)
         pb(i)    = db(i)*cbar(i)**2
c
#endif /* RAREFACTION0 */
c
#ifdef RAREFACTION1
c
c  Compute values for inside a rarefaction fan
c  (method described in Colella (1982) SIAM J Sci Stat Comp 3, p. 107)
c
         cb(i)    = (-sn(i)*ubar(i,j) + 2.0*cbar(i)/(gamma - 1.0))*qb
         ub(i)    = -cb(i)*sn(i)
c         alpha    = (c0(i)**2 - cb(i)**2)/(c0(i)**2 + cb(i)**2)/
c     &              (gamma - 1.0)
c         db(i)    = d0(i)*(1.0 - alpha)/(1.0 + alpha)
         db(i)    = d0(i)*(cb(i)/c0(i))**(2.0/(gamma - 1.0)) !AK
         pb(i)    = db(i)*cb(i)**2/gamma
c
#endif /* RAREFACTION1 */
c
#ifdef RAREFACTION2
c
c  Compute values for inside a rarefaction fan
c     (linear interpolation between end states, as suggested in PPM ref)
c
         frac     = l0(i) - lbar(i)
         if (frac .lt. tiny) frac = tiny
         if (frac .gt. 1.0)  frac = 1.0
         frac     = (0.0 - lbar(i))/frac
         frac     = min(max(frac, 0.0), 1.0)
         pb(i)    = p0(i)*frac + pbar(i,j)*(1.0 - frac)
         db(i)    = d0(i)*frac + dbar(i  )*(1.0 - frac)
         ub(i)    = u0(i)*frac + ubar(i,j)*(1.0 - frac)
c
#endif /* RAREFACTION2 */
c
       enddo
c
c  Cull appropriate states depending on where eulerian position is in solution
c    (lbar >= 0 --> inside post-shock region,
c     l0   <  0 --> outside shock/rarefaction wave,
c     otherwise --> inside rarefaction wave).
c
       do i=i1, i2+1
         if (lbar(i) .ge. 0.0) then !AK
            pb(i) = pbar(i,j)
            db(i) = dbar(i  )
            ub(i) = ubar(i,j)
         endif
         if (l0  (i) .lt. 0.0) then !AK
            pb(i) = p0(i)
            db(i) = d0(i)
            ub(i) = u0(i)
         endif
c         write(22,*) i, ub(i), u0(i), ubar(i,j), sn(i)
c
c         if (db(i) .gt. dbar(i) .and. db(i) .gt. d0(i))
c     &      write(6,1000) i, j, db(i), dbar(i), d0(i)
c 1000  format('euler: ', 2(i4), 3(e12.4))
c
       enddo
c
#define NO_CAVITATE_AVOID
#ifdef CAVITATE_AVOID
       do i=i1, i2+1
         if (ubar(i,j)*ub(i) .lt. 0.0 .and. 
     &       -db(i)*ub(i)*dt/dx(i) .gt. dslice(i,j) .or.
     &        db(i)*ub(i)*dt/dx(i) .gt. dslice(i-1,j)     ) then
            write(6,1010) i,j,ubar(i,j),ub(i),db(i),
     &              dslice(i,j),dslice(i-1,j)
 1010       format('cav:',2i6,8g12.4)
            ub(i) = 0.0
         endif
      enddo
#endif
c
c  Collect values of interest depending on which way fluid is flowing
c     (note: new placement for this - uses ub instead of ubar for consistency)
c
       do i=i1, i2+1
          if (ub(i) .gt. 0.0) then !AK
             vb(i) = vls(i,j)
             wb(i) = wls(i,j)
             geb(i) = gels(i,j)
          else
             vb(i) = vrs(i,j)
             wb(i) = wrs(i,j)
             geb(i) = gers(i,j)
          endif
       enddo
c
c      Repeat for colour variables
c
       do n=1,ncolor
          do i=i1, i2+1
             if (ub(i) .gt. 0.0) then
                colb(i,n) = colls(i,j,n)
             else
                colb(i,n) = colrs(i,j,n)
             endif
          enddo
       enddo
c
c  Dual energy formalism: if sound speed squared is less than eta1*v^2 
c    then discard pbar,dbar,ubar in favour of p0,d0,u0.  This amounts
c     to assuming that we are outside the shocked region but the flow is
c     hypersonic so this should be true.  This is inserted because the
c     mechanism above occasionally fails in hypersonic flows.
c
#ifdef UNUSED
       if (idual .eq. 1) then
          do i=i1, i2+1
             if (gamma*pbar(i,j)/min(dls(i,j),drs(i,j)) .lt. 
     &                                         eta1*ubar(i,j)**2) then
                 pb(i) = p0(i)
                 db(i) = d0(i)
                 ub(i) = u0(i)
              endif
           enddo
        endif
#endif /* UNUSED */
c
c  Calculate total specific energy corresponding to this state
c     (and the specific gas energy).
c
       do i=i1, i2+1
         eb(i) = pb(i)/((gamma-1.0)*db(i)) + 
     &           0.5*(ub(i)**2 + vb(i)**2 + wb(i)**2)
c         geb(i) = pb(i)/((gamma-1.0)*db(i))
       enddo
c
c  Compute terms in differenced hydro equations (eq. 3.1)
c
       if (idiff .ne. 0) then
c
c     ...with diffusion
c
         do i=i1,i2+1
            upb(i)  =  pb(i)*ub(i)
c
            dub(i)  =  ub(i)*db(i)  !AK
c
            duub(i) = dub(i)*ub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*uslice(i-1,j) - dslice(i,j)*uslice(i,j))
c
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
c
c     (should we add diffusion to the cross velocities?  I doubt it)
c     I don't. This doubt kills Noh test problem in 2D at high resolution.
c     Diffusion has to be added to cross velocities. !AK        May 2005.
c
            duvb(i) = dub(i)*vb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*vslice(i-1,j) - dslice(i,j)*vslice(i,j))
            duwb(i) = dub(i)*wb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*wslice(i-1,j) - dslice(i,j)*wslice(i,j))
c
            dueb(i) = dub(i)*eb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*eslice(i-1,j) - dslice(i,j)*eslice(i,j))
c
c	This update must be the last !AK
c
            dub(i)  =  dub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)               - dslice(i,j)            )
         enddo
c
c        If using dual energy formalism, compute dugeb
c
         if (idual .eq. 1) then
            do i=i1,i2+1
              dugeb(i) = dub(i)*geb(i) + diffcoef(i,j)*
     &         (dslice(i-1,j)*geslice(i-1,j) - dslice(i,j)*geslice(i,j))
            enddo
         endif
c
       else
c
c     ...and without
c
         do i=i1, i2+1
            upb(i)  =  pb(i)*ub(i)
            dub(i)  =  ub(i)*db(i)
            duub(i) = dub(i)*ub(i)
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
            dueb(i) = dub(i)*eb(i)
         enddo
c
         if (idual .eq. 1) then
            do i=i1, i2+1
               dugeb(i) = dub(i)*geb(i)
            enddo
         endif
c
       endif
c
c  Copy into flux slices (to return to caller)
c

!      do i=i1, i2+1
!         df(i,j) = dt*dub(i)
!         ef(i,j) = dt*(dueb(i) + upb(i))
!         uf(i,j) = dt*(duub(i) + pb(i))
!         vf(i,j) = dt*duvb(i)
!         wf(i,j) = dt*duwb(i)
!      enddo

       do i=i1, i2+1
          qc = dt/dx(i)
          df(i,j) = qc*dub(i)
          ef(i,j) = qc*(dueb(i) + upb(i))
          uf(i,j) = qc*(duub(i) + pb(i))
          vf(i,j) = qc*duvb(i)
          wf(i,j) = qc*duwb(i)
       enddo

c
       do n=1,ncolor
          do i=i1, i2+1
c            colf(i,j,n) = (dt/dx(i))*dub(i)*colb(i,n) ! color*dens conserved
             colf(i,j,n) = (dt/dx(i))*ub(i)*colb(i,n) ! color conserved
          enddo
       enddo
c
c      Do the same for the gas energy if using the dual energy formalism
c         (note that we do not include the source term)
c
       if (idual .eq. 1) then
          do i=i1, i2+1
             gef(i,j) = dt*dugeb(i)
          enddo
       endif
c     
c  Update conservation laws  (eq. 3.1)
c
       do i=i1, i2
          qc = dt/dx(i)
          dnu(i)      =              dslice(i,j) +
     &                     qc*( dub(i)- dub(i+1))
c          write(21,*) i, dnu(i), dub(i), db(i), ub(i)
c          uadvect = (dnu(i)-dslice(i,j))/(qc*(db(i+1)-db(i)))
c          if (uslice(i,j)/(uadvect) .gt. 1.8 .and. uadvect .ne. 0.0)
c     &        write(6,*) i,j,uadvect,uslice(i,j),dt
c          if (dnu(i) .le. 0.0 .or. eslice(i,j)*dslice(i,j) +
c     &      qc*(dueb(i)-dueb(i+1)+upb(i)-upb(i+1)) .lt. 0.0) then
          if (dnu(i) .le. 0.0 .or. eslice(i,j) .lt. 0.0) then
             if (dnu(i) .le. 0.0) then
                write (6,*) 'euler: dnu <= 0:', i, j, dub(i), dub(i+1)
             else
                write (6,*) 'euler: eslice < 0: ', i, j
                write (6,*) qc*(dueb(i)-dueb(i+1)),
     &                     qc*(upb(i)-upb(i+1)), eslice(i,j)*dslice(i,j)
                write (6,*) eb(i-1), eb(i), eb(i+1)
                write (6,*) dueb(i-1), dueb(i), dueb(i+1)
                write (6,*) upb(i-1), upb(i), upb(i+1)
             endif
             write (6,*) 'e1',dnu(i), dslice(i,j), qc, dt, 
     &                      lbar(i), l0(i)
             write (6,*) 'e2',db(i-1), db(i), db(i+1)
             write (6,*) 'e3',ub(i-1), ub(i), ub(i+1)
             write (6,*) 'e4',pb(i-1), pb(i), pb(i+1)
             write (6,*) 'e5',d0(i-1), d0(i), d0(i+1)
             write (6,*) 'e6',u0(i-1), u0(i), u0(i+1)
             write (6,*) 'e7', p0(i-1), p0(i), p0(i+1)
             write (6,*) 'e8', ubar(i-1,j), ubar(i,j), ubar(i+1,j)
             write (6,*) 'e9', pbar(i-1,j), pbar(i,j), pbar(i+1,j)
             write (6,*) 'e0', dls(i-1,j), dls(i,j), dls(i+1,j)
             write (6,*) 'e1', drs(i-1,j), drs(i,j), drs(i+1,j)
             write (6,*) 'e2', uls(i-1,j), uls(i,j), uls(i+1,j)
             write (6,*) 'e3', urs(i-1,j), urs(i,j), urs(i+1,j)
             write (6,*) 'e4', pls(i-1,j), pls(i,j), pls(i+1,j)
             write (6,*) 'e5', prs(i-1,j), prs(i,j), prs(i+1,j)
             write (6,*) 'e6', dslice(i-1,j), dslice(i,j), dslice(i+1,j)
             write (6,*) 'e7', uslice(i-1,j), uslice(i,j), uslice(i+1,j)
             write (6,*) 'e8', geslice(i-1,j), geslice(i,j), 
     &                         geslice(i+1,j)
             write (6,*) 'e9', eslice(i-1,j), eslice(i,j), eslice(i+1,j)
             write (6,*) 'e0', geslice(i-1,j)*dslice(i-1,j)*(gamma-1.0), 
     &                   geslice(i  ,j)*dslice(i  ,j)*(gamma-1.0), 
     &                   geslice(i+1,j)*dslice(i+1,j)*(gamma-1.0)
             write (6,*) 'e1', ub(i),db(i),diffcoef(i,j)
             if (gravity .eq. 1) 
     &        write (6,*) 'e2',grslice(i-1,j), grslice(i,j), 
     &                   grslice(i+1,j)
             write(6,*) 'e3', dub(i), dub(i+1)
             ERROR_MESSAGE
c
          endif
c           A small cheat: if idual is on, assume a cosmo sim, and institute
c              a minimum density.  This should be a parameter and passed in.
          if (idual .eq. 1 .and. gravity .eq. 1) 
     &          dnu(i) = max(dnu(i),1.0e-3)
          uold(i)     = uslice(i,j)

c     Reformulation to avoid overflows [BWO v1_1_6]

          uslice(i,j) = uslice(i,j)*(dslice(i,j)/dnu(i)) +
     &         (qc/dnu(i))*(duub(i)-duub(i+1) +  pb(i) -  pb(i+1))
          vslice(i,j) = vslice(i,j)*(dslice(i,j)/dnu(i)) +
     &         (qc/dnu(i))*(duvb(i)-duvb(i+1))
          velold = wslice(i,j)
          wslice(i,j) = wslice(i,j)*(dslice(i,j)/dnu(i)) +
     &         (qc/dnu(i))*(duwb(i)-duwb(i+1))

          eold = eslice(i,j)
          eslice(i,j) = max(0.1*eslice(i,j),
     &                  (eslice(i,j)*dslice(i,j) +
     &         qc*(dueb(i)-dueb(i+1) + upb(i) - upb(i+1)))/dnu(i) )
       enddo
c
c      Colour variables (note: colf already multiplied by dt)
c
       do n=1,ncolor
          do i=i1, i2
c     Removed dx from this, and placed earlier during colf calculation [BWO]
             if (colslice(i,j,n)              +
     &         (colf(i,j,n)-colf(i+1,j,n)) .lt. 0.0) then
                write(6,*)'euler_c0:',i,j,n,dx(i)
                write(6,*)'euler_c1:',colf(i-1,j,n),colf(i,j,n),
     &                      colf(i+1,j,n)
                write(6,*)'euler_c2:',colslice(i-1,j,n),colslice(i,j,n),
     &                     colslice(i+1,j,n)
                write(6,*)'euler_c3:',colb(i-1,n),colb(i,n),colb(i+1,n)
                write(6,*)'euler_c4:',colls(i-1,j,n),colls(i,j,n),
     &                     colls(i+1,j,n)
                write(6,*)'euler_c5:',colrs(i-1,j,n),colrs(i,j,n),
     &                     colrs(i+1,j,n)
                write(6,*)'euler_c6:',ub(i-1),ub(i),ub(i+1)
                write(6,*)'euler_c7:',u0(i-1),u0(i),u0(i+1)
                write(6,*)'euler_c8:',ubar(i-1,j),ubar(i,j),ubar(i+1,j)
c                stop
             endif
             colslice(i,j,n) =  colslice(i,j,n)              +
     &            (colf(i,j,n)-colf(i+1,j,n)) ! c conserved
             colslice(i,j,n) = max(colslice(i,j,n),tiny)
          enddo
       enddo
c
c      Conservation law for gas energy, if using the dual energy formalism
c           (this includes both the flux term and a source term - yuck).
c         Here, we compute the ratio of thermal energies derived the
c            two different ways and then use that ratio to choose how to
c            compute the pressure at the center of zone i.  This is needed
c            for the source term in the gas energy equation.
c
       if (idual .eq. 1) then
          do i=i1-1, i2+1
             eratio(i) = (eslice(i,j) - 
     &         0.5*(uslice(i,j)**2 + vslice(i,j)**2 + wslice(i,j)**2)
     &                    )/eslice(i,j)
          enddo
          do i=i1, i2
c             if (max(eratio(i-1),eratio(i),eratio(i+1)) .gt. eta2) then
                pcent(i) = max((gamma-1.0)*geslice(i,j)*dslice(i,j), 
     &                         tiny)
c             else
c                pcent(i) = 0.5*(pb(i)+pb(i+1))
c             endif
          enddo
          do i=i1, i2
             if (geslice(i,j) .lt. 0.0) write(6,*) i,j,geslice(i,j)
             qc = dt/dx(i)
             geslice(i,j) = max((geslice(i,j)*dslice(i,j) +
     &            qc*(dugeb(i)-dugeb(i+1) + pcent(i)*(ub(i)-ub(i+1))  
     &                                                       ))/dnu(i)
     &                          ,0.5*geslice(i,j))
c             if (geslice(i,j) .lt. tiny) geslice(i,j) = tiny
             if (geslice(i,j) .lt. 0.0) then
                write(6,*) i,j,dslice(i,j),dnu(i),pcent(i),eslice(i,j),
     &               qc*(dugeb(i)-dugeb(i+1)),
     &               qc*(pcent(i)*(ub(i)-ub(i+1))),
     &               geslice(i,j)*dslice(i,j),qc,
     &               uslice(i,j),uslice(i+1,j)
                write(6,*) pb(i),pb(i+1),p0(i),p0(i+1),
     &                     pbar(i,j),pbar(i+1,j),lbar(i),l0(i),
     &                     pls(i,j),prs(i,j),pls(i+1,j),prs(i+1,j)
                write(6,*) ub(i),ub(i+1),geb(i),geb(i+1),
     &                     db(i),db(i+1),dub(i),dub(i+1)
                ERROR_MESSAGE
             endif
          enddo
       endif
c
c  If there is gravity, the compute the second order correction to the
c   acceleration due a slope in both density and acceleration.
c
#ifdef GRAVITY_SECOND_ORDER_CORRECTION
c
c      Compute slopes and enforce limited monotonoctiy on dddx
c
       if (gravity .eq. 1) then
          do i=i1,i2
             dadx(i) = grslice(i+1,j) - grslice(i-1,j)
             dddx(i) =  dslice(i+1,j) -  dslice(i-1,j)
c
             dddx(i) = 2.0*( dslice(i,j) - max(dslice(i,j) -
     &                   0.5*dddx(i), min(dslice(i,j), dslice(i-1,j))))
             dddx(i) = 2.0*(-dslice(i,j) + max(dslice(i,j) +
     &                   0.5*dddx(i), min(dslice(i,j), dslice(i+1,j))))
          enddo
c
          do i=i1, i2
             grslice(i,j) = grslice(i,j) + 
     &                      0.5*dadx(i)*dddx(i)/(12.0*dslice(i,j))
          enddo
       endif
c
#endif /* GRAVITY_SECOND_ORDER_CORRECTION */
c
c  If there is gravity, add the gravity terms here (eq. 3.1 or 3.8).
c    (Note: the acceleration is already time-centered).
c
       if (gravity .eq. 1) then
          do i=i1, i2

#define GRAVITY_METHOD1

#ifdef GRAVITY_METHOD1 /* (right) */
             tempu = uslice(i,j)
             uslice(i,j) = uslice(i,j) + 
     &                     dt*grslice(i,j)*0.5*(dslice(i,j)/dnu(i)+1.0)
             eold = eslice(i,j)
             eslice(i,j) = eslice(i,j) + (dt*grslice(i,j))*0.5*
     &                     (uslice(i,j) + uold(i)*dslice(i,j)/dnu(i))
            if (eslice(i,j) .le. 0) then
            write(6,*) 'eu1a',i,j,eslice(i,j),
     &                     dslice(i,j),dnu(i),uslice(i,j),grslice(i,j),
     &                     dt*grslice(i,j)*0.5*(dslice(i,j)/dnu(i)+1.0),
     &                     dt*grslice(i,j)*0.5*
     &                     (uslice(i,j) + uold(i)*dslice(i,j)/dnu(i))
            write(6,*) 'eu1b', eold, dt, grslice(i,j)
            write(6,*) 'eu1c', uslice(i,j), uold(i) 
            write(6,*) 'eu1d', dslice(i,j), dnu(i), tempu
         endif

         eslice(i,j) = max(eslice(i,j), tiny)

#endif /* GRAVITY_METHOD1 */

#ifdef GRAVITY_METHOD2 /* (wrong) */
             uslice(i,j) = uslice(i,j) + dt*grslice(i,j)
             eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*0.5*
     &                                      (uslice(i,j) + uold(i))
#endif /* GRAVITY_METHOD2 */
          enddo
       endif
c
c  Update the new density
c
       do i=i1, i2
          dslice(i,j) = dnu(i)
       enddo
c
      enddo
c
      return
      end

*****************************************************************************
*                                                                           *
* Copyright 2004 Greg Bryan                                                 *
* Copyright 2004 Laboratory for Computational Astrophysics                  *
* Copyright 2004 Board of Trustees of the University of Illinois            *
* Copyright 2004 Regents of the University of California                    *
*                                                                           *
* This software is released under the terms of the "Enzo Public License"    *
* in the accompanying LICENSE file.                                         *
*                                                                           *
*****************************************************************************

c=======================================================================
c/////////////////////  SUBROUTINE F_ERROR  \\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine f_error (sourcefile, linenumber)
c
c     PRINT ERROR MESSAGE AND EXIT PROGRAM
c=======================================================================


      implicit none

      CHARACTER sourcefile*(*)
      INTEGER linenumber

      CALL fc_error (sourcefile // char(0), linenumber)

      return
      end

c=======================================================================
c/////////////////////  SUBROUTINE F_WARNING  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine f_warning (sourcefile, linenumber)
c
c     PRINT WARNING MESSAGE AND CONTINUE
c=======================================================================


      implicit none

      CHARACTER sourcefile*(*)
      INTEGER linenumber

      CALL fc_warning (sourcefile // char(0), linenumber)
      return
      end

c=======================================================================
c///////////////////////  SUBROUTINE INTEULER  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dxi, flatten,
     &            idim, jdim, i1, i2, j1, j2, idual, eta1, eta2,
     &            isteep, iflatten, dt, gamma, ipresfree,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            ncolor, colslice, colls, colrs
     &                    )
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January,1991
c  modified1:  June, 1993 by Greg Bryan (changed to eulerian)
c  modified2:  July, 1994 by Greg Bryan (changed to slicewise)
c  modified3:  April, 1995 by GB (added gravity)
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  This version computes the Eulerian corrections
c    to the left and right states described in section three of Colella &
c    Woodward (1984), JCP.  The routine works on one two dimensional
c    slice at a time.
c
c  INPUT:
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dxi    - distance between Eulerian zone edges in sweep direction
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    flatten - ammount of flattening (calculated in calcdiss)
c    gamma  - ideal gas law constant
c    gravity - gravity flag (0 = off)
c    grslice - acceleration in this direction in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idual  - dual energy formalism flag (0 = off)
c    iflatten - integer flag for flattener (eq. A1, A2) (0 = off)
c    isteep - integer flag for steepener (eq. 1.14,1.17,3.2) (0 = off)
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pslice - extracted 2d slice of the pressure, p
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    wslice - extracted 2d slice of the 3-velocity, w
c    
c  OUTPUT:
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c    vl,rs  - 2-velocity at left and right edges of each cell
c    wl,rs  - 3-velocity at left and right edges of each cell
c
c  EXTERNALS:
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, i1, i2, idim, idual, iflatten, ipresfree, isteep,
     &        j1, j2, jdim, ncolor
      real    dt, eta1, eta2, gamma, tempw
      real   dslice(idim,jdim),     dxi(idim     ),  pslice(idim,jdim),
     &       uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim),
     &      grslice(idim,jdim), geslice(idim,jdim)
      real      dls(idim,jdim),     drs(idim,jdim), flatten(idim,jdim),
     &          pls(idim,jdim),
     &          prs(idim,jdim),    gels(idim,jdim),    gers(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim),     vls(idim,jdim),
     &          vrs(idim,jdim),     wls(idim,jdim),     wrs(idim,jdim)
      real colslice(idim,jdim,ncolor),   colls(idim,jdim,ncolor),
     &        colrs(idim,jdim,ncolor)
c
c  local declarations
c
      integer i, j, ic
      real steepen(ijkn),tmp1(ijkn),tmp2(ijkn),tmp3(ijkn),tmp4(ijkn)
      real qa,qb,qc,qd,qe,s1,s2      
      real  c1(ijkn), c2(ijkn), c3(ijkn), c4(ijkn), c5(ijkn), c6(ijkn)
     &    ,dp(ijkn),pl(ijkn),pr(ijkn),p6(ijkn)
     &    ,du(ijkn),ul(ijkn),ur(ijkn),u6(ijkn)
     &    ,dla(ijkn),dra(ijkn),pla(ijkn),pra(ijkn),ula(ijkn),ura(ijkn)
     &    ,vla(ijkn),vra(ijkn),wla(ijkn),wra(ijkn)
     &    ,plm(ijkn),prm(ijkn),ulm(ijkn),urm(ijkn)
     &    ,dl0(ijkn),dr0(ijkn),pl0(ijkn),pr0(ijkn)
     &    ,plp(ijkn),prp(ijkn),ulp(ijkn),urp(ijkn),ul0(ijkn),ur0(ijkn)
     &    ,vl0(ijkn),vr0(ijkn),wl0(ijkn),wr0(ijkn)
     &    , cs(ijkn),d2d(ijkn),dxb(ijkn)
     &    , cm(ijkn), c0(ijkn), cp(ijkn),char1(ijkn),char2(ijkn)
     &    ,betalp(ijkn),betalm(ijkn),betal0(ijkn),cla(ijkn)
     &    ,betarp(ijkn),betarm(ijkn),betar0(ijkn),cra(ijkn)
     &    ,gela(ijkn),gera(ijkn),gel0(ijkn),ger0(ijkn)
     &    ,colla(ijkn,MAX_COLOR),colra(ijkn,MAX_COLOR)
     &    ,coll0(ijkn,MAX_COLOR),colr0(ijkn,MAX_COLOR)
c
c  parameters
c
      real ft
      parameter(ft = 4.0/3.0)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     write(6,*) 'INTEULER: dt =',dt,' isteep =',isteep
c     write(6,*) 'INTEULER: iflatten =',iflatten
c     write(6,*) 'INTEULER: idim =',idim,' jdim =',jdim
c     write(6,*) 'INTEULER: i1   =',i1,  ' i2   =',i2
c     write(6,*) 'INTEULER: j1   =',j1,  ' j2   =',j2
c
c     Compute coefficients used in interpolation formulae (from eq. 1.6)
c
#ifdef OPTIMIZE_CELL_WIDTH
#  define DXI   dxi(i1)
#  define DXIMM dxi(i1)
#  define DXIM  dxi(i1)
#  define DXIP  dxi(i1)
#  define C1    c1(i1)
#  define C2    c2(i1)
#  define C3    c3(i1)
#  define C4    c4(i1)
#  define C5    c5(i1)
#  define C6    c6(i1)
#else
#  define DXI   dxi(i)
#  define DXIMM dxi(i-2)
#  define DXIM  dxi(i-1)
#  define DXIP  dxi(i+1)
#  define C1    c1(i)
#  define C2    c2(i)
#  define C3    c3(i)
#  define C4    c4(i)
#  define C5    c5(i)
#  define C6    c6(i)
#endif

#ifndef OPTIMIZE_CELL_WIDTH
      do i=i1-2,i2+2
#endif
         qa    = DXI/(DXIM + DXI + DXIP)
         C1 = qa*(2.0*DXIM + DXI)/(DXIP + DXI)
         C2 = qa*(2.0*DXIP + DXI)/(DXIM + DXI)
#ifndef OPTIMIZE_CELL_WIDTH
      enddo
#endif
c
#ifndef OPTIMIZE_CELL_WIDTH
      do i=i1-1,i2+2
#endif
         qa    = DXIMM + DXIM + DXI + DXIP
         qb    = DXIM/(DXIM + DXI)
         qc    = (DXIMM + DXIM)/(2.0*DXIM + DXI)
         qd    = (DXIP + DXI)/(2.0*DXI + DXIM)
         qb    = qb + 2.0*DXI*qb/qa*(qc-qd)
         C3 = 1.0 - qb
         C4 = qb
         C5 =  DXI/qa*qd
         C6 = -DXIM/qa*qc
#ifndef OPTIMIZE_CELL_WIDTH
      enddo
#endif
c
c    Loop over sweep lines (in this slice)
c
      do 400 j=j1, j2
c
c     Precompute steepening coefficients if needed (eqns 1.14-1.17, plus 3.2)
c
      if (isteep .ne. 0) then
         do i=i1-2,i2+2
            qa     = DXIM + DXI + DXIP
            d2d(i) = (dslice(i+1,j) - dslice(i,j))/(DXIP + DXI)
            d2d(i) = (d2d(i) - (dslice(i,j)-dslice(i-1,j))
     &               /(DXI+DXIM))/qa
            dxb(i) = 0.5*(DXI + DXIP)
         enddo
         do i=i1-1,i2+1
            qc = abs(dslice(i+1,j) - dslice(i-1,j))
     &           - 0.01*min(abs(dslice(i+1,j)),abs(dslice(i-1,j)))
            s1 = (d2d(i-1) - d2d(i+1))*(dxb(i-1)**3 + dxb(i)**3)
     &           /((dxb(i) + dxb(i-1))*
     &           (dslice(i+1,j) - dslice(i-1,j) + tiny))
            if (d2d(i+1)*d2d(i-1) .gt. 0.0) s1 = 0.0
            if (qc .le. 0.0) s1 = 0.0
            s2 = max(0.0, min(20.0*(s1-0.05), 1.0))
            qa = abs(dslice(i+1,j) - dslice(i-1,j))/
     &           min(dslice(i+1,j),  dslice(i-1,j))
            qb = abs(pslice(i+1,j) - pslice(i-1,j))/
     &           min(pslice(i+1,j),  pslice(i-1,j))
            if (gamma*0.1*qa .ge. qb) then
               steepen(i) = s2
            else
               steepen(i) = 0.0
            endif
         enddo
      endif
c
c     Precompute left and right characteristic distances
c
      do i=i1-1,i2+1
         cs(i) = sqrt(gamma*pslice(i,j)/dslice(i,j))
         if (ipresfree .eq. 1) cs(i) = tiny
         char1(i) = max(0.0, dt*(uslice(i,j)+cs(i)))/(2.0*DXI)
         char2(i) = max(0.0,-dt*(uslice(i,j)-cs(i)))/(2.0*DXI)
      enddo
c
c Brought dt and dx terms here to avoid possible overflows [BWO v1_1_6]
c
      do i=i1-1,i2+1
         cm(i) = (uslice(i,j)-cs(i))*(dt/(2.0*DXI))
         c0(i) = (uslice(i,j)      )*(dt/(2.0*DXI))
         cp(i) = (uslice(i,j)+cs(i))*(dt/(2.0*DXI))
      enddo
c
c     Compute left and right states for each variable
c       (steepening, if requested, is only applied to density)
c
      call intvar(dslice(1,j), idim, i1, i2, isteep, steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            tmp1, tmp2, tmp3, tmp4, dla, dra, dl0, dr0)
c
      call intvar(pslice(1,j), idim, i1, i2, 0     , steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            dp, pl, pr, p6, pla, pra, pl0, pr0)
c
      call intvar(uslice(1,j), idim, i1, i2, 0     , steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            du, ul, ur, u6, ula, ura, ul0, ur0)
c
      call intvar(vslice(1,j), idim, i1, i2, 0     , steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            tmp1, tmp2, tmp3, tmp4, vla, vra, vl0, vr0)
c
      call intvar(wslice(1,j), idim, i1, i2, 0     , steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            tmp1, tmp2, tmp3, tmp4, wla, wra, wl0, wr0)
c
      if (idual .eq. 1)
     & call intvar(geslice(1,j), idim, i1, i2, 0   , steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            tmp1, tmp2, tmp3, tmp4, gela, gera, gel0, ger0)
c
      do ic=1,ncolor
        call intvar(colslice(1,j,ic), idim, i1,i2, 0, steepen, iflatten, 
     &            flatten, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &            tmp1, tmp2, tmp3, tmp4, colla(1,ic), colra(1,ic),
     &            coll0(1,ic), colr0(1,ic))
      enddo
c
#ifdef UNUSED
      do i=i1,i2+1
         if (dla(i)/dslice(i-1,j) .gt. 4.0) then
            write(6,*) 'inteuler left:',i,j
            write(6,*) dla(i-1),dla(i),dla(i+1)
            write(6,*) dslice(i-1,j),dslice(i,j),dslice(i+1,j)
            write(6,*) dra(i-1),dra(i),dra(i+1)
            write(6,*) ula(i-1),ula(i),ula(i+1)
            write(6,*) pla(i-1),pla(i),pla(i+1)
            write(6,*) uslice(i-1,j),uslice(i,j),uslice(i+1,j)
         endif
      enddo
#endif /* UNUSED */
c
c
c Correct the initial guess from the linearized gas equations
c
c     First, compute averge over characteristic domain of dependance (3.5)
c
c
      do i=i1,i2+1
        plm(i)= pr(i-1)-cm(i-1)*(dp(i-1)-(1.0-ft*cm(i-1))*p6(i-1))
        prm(i)= pl(i  )-cm(i  )*(dp(i  )+(1.0+ft*cm(i  ))*p6(i  ))
        plp(i)= pr(i-1)-cp(i-1)*(dp(i-1)-(1.0-ft*cp(i-1))*p6(i-1))
        prp(i)= pl(i  )-cp(i  )*(dp(i  )+(1.0+ft*cp(i  ))*p6(i  ))
      enddo
c
      do i=i1,i2+1
        ulm(i)= ur(i-1)-cm(i-1)*(du(i-1)-(1.0-ft*cm(i-1))*u6(i-1))
        urm(i)= ul(i  )-cm(i  )*(du(i  )+(1.0+ft*cm(i  ))*u6(i  ))
        ulp(i)= ur(i-1)-cp(i-1)*(du(i-1)-(1.0-ft*cp(i-1))*u6(i-1))
        urp(i)= ul(i  )-cp(i  )*(du(i  )+(1.0+ft*cp(i  ))*u6(i  ))
      enddo
c
c     Compute correction terms (3.7)
c
c     Changed the order of operations here to prevent overflow [BWO v1_1_6]

      do i = i1, i2+1
c     cla(i) = sqrt(max(gamma*pla(i)*dla(i), 0.0))
c     cra(i) = sqrt(max(gamma*pra(i)*dra(i), 0.0))
         cla(i) = sqrt(max(gamma*pla(i),0.0)) * sqrt(max(dla(i),0.0))
         cra(i) = sqrt(max(gamma*pra(i),0.0)) * sqrt(max(dra(i),0.0))
         if ( (cla(i) .gt. huge) .or.
     &        (cla(i) .eq.0)) then
            write(6,*) 'inteuler: cla', cla(i), pla(i),dla(i)
            stop
         endif
         if ( (cra(i) .gt. huge) .or.
     &        (cra(i).eq.0)) then
            write(6,*) 'inteuler: cra', cra(i), pra(i),dra(i)
            stop
         endif
      enddo
c
c     a) left side
c
      do i = i1, i2+1
         betalp(i) = (ula(i)-ulp(i)) + (pla(i)-plp(i))/cla(i)
         betalm(i) = (ula(i)-ulm(i)) - (pla(i)-plm(i))/cla(i)
         betal0(i) = ((pla(i)-pl0(i))/cla(i))/cla(i) + 
     &               1.0/dla(i) - 1.0/dl0(i)
      enddo
c
c     Add gravity component
c      
      if (gravity .eq. 1) then
        do i = i1, i2+1
c          if (cla(i) .gt. 0.3*ula(i)) then
c          if (gamma*pla(i)/dla(i) .gt. eta2*ula(i)**2) then
           betalp(i) = betalp(i) - 0.25*dt*(grslice(i-1,j)+grslice(i,j))
           betalm(i) = betalm(i) - 0.25*dt*(grslice(i-1,j)+grslice(i,j))
c          endif
        enddo
      endif
c
      do i = i1, i2+1
         betalp(i) = -betalp(i)/(2.0*cla(i))
         betalm(i) = +betalm(i)/(2.0*cla(i))
      enddo
c
      do i = i1, i2+1
         if (cp(i-1) .le. 0) betalp(i) = 0
         if (cm(i-1) .le. 0) betalm(i) = 0
         if (c0(i-1) .le. 0) betal0(i) = 0
      enddo
c
c     b) right side
c
      do i = i1, i2+1
         betarp(i) = (ura(i)-urp(i)) + (pra(i)-prp(i))/cra(i)
         betarm(i) = (ura(i)-urm(i)) - (pra(i)-prm(i))/cra(i)
         betar0(i) = ((pra(i)-pr0(i))/cra(i))/cra(i) + 
     &               1.0/dra(i) - 1.0/dr0(i)
      enddo
c
c     Add gravity component
c      
      if (gravity .eq. 1) then
        do i = i1, i2+1
c          if (cra(i) .gt. 0.3*ura(i)) then
c           if (gamma*pra(i)/dra(i) .gt. eta2*ura(i)**2) then
           betarp(i) = betarp(i) - 0.25*dt*(grslice(i-1,j)+grslice(i,j))
           betarm(i) = betarm(i) - 0.25*dt*(grslice(i-1,j)+grslice(i,j))
c          endif
        enddo
      endif
c
      do i = i1, i2+1
         betarp(i) = -betarp(i)/(2.0*cra(i))
         betarm(i) = +betarm(i)/(2.0*cra(i))
      enddo
c
      do i = i1, i2+1
         if (cp(i) .ge. 0) betarp(i) = 0
         if (cm(i) .ge. 0) betarm(i) = 0
         if (c0(i) .ge. 0) betar0(i) = 0
      enddo
c
c     Finally, combine to create corrected left/right states (eq. 3.6)
c
      do i=i1,i2+1
c   Changed order of operations to avoid overflow in 32 bits [BWO v1_1_6]
         pls(i,j) = pla(i) + ((betalp(i)+betalm(i))*cla(i))*cla(i)
         prs(i,j) = pra(i) + ((betarp(i)+betarm(i))*cra(i))*cra(i)
c
         uls(i,j) = ula(i) + (betalp(i)-betalm(i))*cla(i)
         urs(i,j) = ura(i) + (betarp(i)-betarm(i))*cra(i)
c
         dls(i,j) = 1.0/(1.0/dla(i) - (betal0(i)+betalp(i)+betalm(i)))
         drs(i,j) = 1.0/(1.0/dra(i) - (betar0(i)+betarp(i)+betarm(i)))
      enddo
c
c     Take the appropriate state from the advected variables
c
      do i=i1,i2+1
         if (uslice(i-1,j) .le. 0.0) then
            vls(i,j)  = vla(i)
            wls(i,j)  = wla(i)
            gels(i,j) = gela(i)
         else
            vls(i,j)  = vl0(i)
            wls(i,j)  = wl0(i)
            gels(i,j) = gel0(i)
         endif
c
         if (uslice(i  ,j) .ge. 0.0) then
            vrs(i,j)  = vra(i)
            wrs(i,j)  = wra(i)
            gers(i,j) = gera(i)
         else
            vrs(i,j)  = vr0(i)
            wrs(i,j)  = wr0(i)
            gers(i,j) = ger0(i)
         endif
      enddo
c
      do ic=1,ncolor
         do i=i1,i2+1
            if (uslice(i-1,j) .le. 0.0) then
               colls(i,j,ic) = colla(i,ic)
            else
               colls(i,j,ic) = coll0(i,ic)
            endif
c
            if (uslice(i  ,j) .ge. 0.0) then
               colrs(i,j,ic) = colra(i,ic)
            else
               colrs(i,j,ic) = colr0(i,ic)
            endif
         enddo
      enddo
c
c  Dual energy formalism: if sound speed squared is less than eta1*v^2 
c    then discard the corrections and use pla, ula, dla.  This amounts
c     to assuming that we are outside the shocked region but the flow is
c     hypersonic so this should be true.  This is inserted because the
c     corrections are inaccurate for hypersonic flows.
c
      if (idual .eq. 1) then
         do i=i1, i2+1
            if (gamma*pla(i)/dla(i) .lt. (eta2*ula(i))*ula(i) .or.
     &          max(abs(cm(i-1)),abs(c0(i-1)),abs(cp(i-1))) .lt. 
     &          1.0e-3
     &          .or. dls(i,j)/dla(i) .gt. 5.0) then
               pls(i,j) = pla(i)
               uls(i,j) = ula(i)
               dls(i,j) = dla(i)
            endif
            if (gamma*pra(i)/dra(i) .lt. (eta2*ura(i))*ura(i) .or.
     &          max(abs(cm(i)),abs(c0(i)),abs(cp(i))) .lt. 1.0e-3
     &          .or. drs(i,j)/dra(i) .gt. 5.0) then
               prs(i,j) = pra(i)
               urs(i,j) = ura(i)
               drs(i,j) = dra(i)
            endif
         enddo
      endif
c
c     Testing code
c
#define CORRECTION
c
#ifdef NO_CORRECTION
      do i=i1, i2+1
         pls(i,j) = pla(i)
         uls(i,j) = ula(i)
         dls(i,j) = dla(i)
         prs(i,j) = pra(i)
         urs(i,j) = ura(i)
         drs(i,j) = dra(i)
      enddo
#endif /* NO_CORRECTION */
c
c     Enforce minimum values.
c
      do i=i1, i2+1
         pls(i,j) = max(pls(i,j), tiny)
         prs(i,j) = max(prs(i,j), tiny)
         dls(i,j) = max(dls(i,j), tiny)
         drs(i,j) = max(drs(i,j), tiny)
      enddo
c
c     If approximating pressure free conditions, then the density should be
c       reset to the pre-corrected state.
c
      if (ipresfree .eq. 1) then
         do i=i1, i2+1
            dls(i,j) = dla(i)
            drs(i,j) = dra(i)
         enddo
      endif
c
400   continue
c
      return
      end

c=======================================================================
c///////////////////////  SUBROUTINE INTVAR  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intvar(qslice, idim, i1, i2, isteep, steepen, iflatten, 
     &                  flatten, c1, c2, c3, c4, c5, c6, char1, char2,
     &                  c0, dq, ql, qr, q6, qla, qra, ql0, qr0)
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Greg Bryan
c  date:       March, 1996
c  modified1:
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  This version computes the Eulerian corrections
c    to the left and right states described in section three of Colella &
c    Woodward (1984), JCP.  The routine works on a single variable in
c    one dimension.
c
c  INPUT:
c    qslice   - one dimensional field of quantity q (one of d,e,u,v...)
c    idim     - declared dimension of 1D fields
c    i1, i2   - start and end indexes of active region
c    isteep   - steepening flag (1 = on, 0 = off); only apply to density!
c    steepen    - steepening coefficients
c    iflatten - flattening flag (1 = on, 0 = off)
c    flatten  - flattening coefficients
c    c1-6     - precomputed grid coefficients
c    char1,2  - characteristic distances for +/- waves (for average)
c    c0       - characteristic distance (for lagrangean cell face)
c    dq, ql, qr, q6 - 1D field temporaries
c    
c  OUTPUT:
c    qla, qra - left and right state values (from char1,2)
c    ql0, qr0 - left and right state values (from c0)
c
c  EXTERNALS:
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  argument declarations
c
      integer idim, i1, i2, iflatten, isteep
      real c1(idim), c2(idim), c3(idim), c4(idim), c5(idim), c6(idim),
     &     char1(idim), char2(idim), c0(idim),
     &     qla(idim), qra(idim), ql0(idim), qr0(idim)
      real qslice(idim), steepen(idim), flatten(idim)
c
c  parameters
c
      real ft
      parameter(ft = 4.0/3.0)
c
c  local declarations (arrays passed as temps)
c
      integer i
      real qplus, qmnus

c  note that in Matt Turk's version of intvar temp1, temp2, temp3
c  are all promoted to real*8 to avoid overflow.  I avoid this here,
c  but note that this may be necessary.  If so, uncomment the next 
c  two lines and then comment out the following line. [BWO v1_1_6]
c <OLD>
c      real*8 temp1, temp2, temp3
c </OLD>
c <NEW>
      real temp1, temp2, temp3
c </NEW>

      real dq(idim), ql(idim), qr(idim), q6(idim)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     Compute average linear slopes (eqn 1.7)
c      Monotonize (eqn 1.8)
c
#ifdef OPTIMIZE_CELL_WIDTH
#  define C1 c1(i1)
#  define C2 c2(i1)
#  define C3 c3(i1)
#  define C4 c4(i1)
#  define C5 c5(i1)
#  define C6 c6(i1)
#else
#  define C1 c1(i)
#  define C2 c2(i)
#  define C3 c3(i)
#  define C4 c4(i)
#  define C5 c5(i)
#  define C6 c6(i)

#endif
      do i = i1-2, i2+2
         qplus = qslice(i+1)-qslice(i  )
         qmnus = qslice(i  )-qslice(i-1)
         dq(i) = C1*qplus + C2*qmnus
         temp1 = min(abs(dq(i)), 2.0*abs(qmnus), 2.0*abs(qplus))
         if (qplus*qmnus .gt. 0) then
            dq(i) = temp1*sign(1.0, dq(i))
         else
            dq(i) = 0.0
         endif
      enddo
c     
c     Construct left and right values (eqn 1.6)
c
      do i = i1-1, i2+2
         ql(i) = C3*qslice(i-1) + C4*qslice(i) +
     &           C5*    dq(i-1) + C6*dq(i)
         qr(i-1) = ql(i)
      enddo
c
c     Steepen if asked for (use precomputed steepening parameter)
c
      if (isteep .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = (1.0-steepen(i))*ql(i) + 
     &              steepen(i)*(qslice(i-1)+0.5*dq(i-1))
            qr(i) = (1.0-steepen(i))*qr(i) + 
     &              steepen(i)*(qslice(i+1)-0.5*dq(i+1))
         enddo
      endif
c
c     Monotonize again (eqn 1.10)
c
      do i=i1-1,i2+1
         temp1 = (qr(i)-qslice(i))*(qslice(i)-ql(i))
         temp2 = qr(i)-ql(i)
         temp3 = 6.0*(qslice(i)-0.5*(qr(i)+ql(i)))
         if (temp1 .le. 0.0) then
            ql(i) = qslice(i)
            qr(i) = qslice(i)
         endif
c     Factored out temp2 [JOB v1_1_6]
         if (temp2*(temp2-temp3) .lt. 0.0) then
            ql(i) = 3.0*qslice(i) - 2.0*qr(i)
         endif	
         if (temp2*(temp2+temp3) .lt. 0.0) then
            qr(i) = 3.0*qslice(i) - 2.0*ql(i)
         endif
      enddo
c
c     If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
c
      if (iflatten .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = qslice(i)*flatten(i) + ql(i)*(1.0-flatten(i))
            qr(i) = qslice(i)*flatten(i) + qr(i)*(1.0-flatten(i))
         enddo
      endif
c
c    Now construct left and right interface values (eqn 1.12 and 3.3)
c
      do i = i1-1, i2+1
         q6(i) = 6.0*(qslice(i)-0.5*(ql(i)+qr(i)))
         dq(i) = qr(i) - ql(i)
         if(q6(i).lt.-1e30)then
           write(6,*) 'intvar: q60a',q6(i),qslice(i)
           write(6,*) 'intvar: q60b',ql(i),qr(i)
         endif
      enddo
c
      do i = i1, i2+1
        qla(i)= qr(i-1)-char1(i-1)*(dq(i-1)-(1.0-ft*char1(i-1))*q6(i-1))
        qra(i)= ql(i  )+char2(i  )*(dq(i  )+(1.0-ft*char2(i  ))*q6(i  ))
        if (qra(i).lt.-1e30)then
          write(6,*)  'intvar: qra0a',qra(i),ql(i),char2(i)
          write(6,*)  'intvar: qra0b',dq(i),ft,q6(i)
        endif
      enddo
c
      do i=i1,i2+1
         ql0(i) = qr(i-1)-c0(i-1)*(dq(i-1)-(1.0-ft*c0(i-1))*q6(i-1))
         qr0(i) = ql(i  )-c0(i  )*(dq(i  )+(1.0+ft*c0(i  ))*q6(i  ))
      enddo
c
c
      return
      end

c=======================================================================
c//////////////////////  SUBROUTINE PGAS2D_DUAL  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pgas2d_dual(dslice, eslice, geslice, pslice, 
     &                       uslice, vslice, wslice, eta1, eta2,
     &                       idim, jdim, i1, i2, j1, j2, gamma, pmin)
c
c  COMPUTES GAS PRESSURE ON A SLICE (DUAL ENERGY VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1994
c  modified1:
c
c  PURPOSE: Computes gas pressure needed for interpolations for PPM
c           solver.  EOS used is p=(gamma-1)*d*[E-(u**2+v**2+w**2)/2],
c           (eqn 2.1) Only works for gamma law gas with gamma > 0.
c           This particular version uses the dual energy formalism
c           (see PPM gravity method paper) to determine the pressure
c           from a combination of the total energy and the gas energy.
c
c  INPUTS:
c     dslice  - density slice
c     eslice  - total specific energy slice
c     eta1    - selection parameter for gas energy (typically ~0.001)
c     eta2    - selection parameter for total energy (typically ~0.1)
c     gamma   - gamma law parameter
c     geslice - specific gas energy slice
c     i1,i2   - starting and ending addresses for dimension 1
c     idim    - declared leading dimension of slices
c     j1,j2   - starting and ending addresses for dimension 2
c     jdim    - declared second dimension of slices
c     pmin    - minimum allowed pressure
c     uslice  - velocity-1 slice
c     vslice  - velocity-2 slice
c     wslice  - velocity-3 slice
c
c  OUTPUT ARGUMENTS: 
c     pslice  - pressure slice
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer i1, i2, idim, j1, j2, jdim
      real    eta1, eta2, gamma, pmin
      real dslice(idim,jdim), eslice(idim,jdim), pslice(idim,jdim),
     &     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim),
     &    geslice(idim,jdim)
c
c  locals
c
      integer i, j, im1, ip1
      real    demax, ge1, ge2, ke
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c  compute the pressure
c
      do j=j1,j2
         do i=i1,i2
c
c           Compute the specific energy from the total energy
c
            ke = 0.5*(uslice(i,j)**2 + vslice(i,j)**2 + wslice(i,j)**2)
            ge1 = eslice(i,j) - ke
c
c           Find the maximum nearby total energy (not specific)
c
            im1 = max(i-1,i1)
            ip1 = min(i+1,i2)
            demax = max(dslice(i  ,j)*eslice(i  ,j),
     &                  dslice(im1,j)*eslice(im1,j),
     &                  dslice(ip1,j)*eslice(ip1,j))
c
c           If the ratio of the gas energy to the max nearby total energy is
c             is > eta2 then use the gas energy computed from the total energy
c             to update the specific energy (total energy is better).
c
            if (geslice(i,j) .le. 0.0) write(6,*) 
     &            'pga:',i,j,geslice(i,j),demax
            if (ge1*dslice(i,j)/demax .gt. eta2) geslice(i,j) = ge1
            if (geslice(i,j) .le. 0.0) write(6,*) 
     &            'pgb:',i,j,geslice(i,j),dslice(i,j),demax,eslice(i,j)
c
c           If the ratio of the specific gas energy to specific total energy 
c             (in this cell) is < eta1 then use the gas energy to update
c             the specific energy (gas energy is better).
c
            if (ge1/eslice(i,j) .gt. eta1) then
               ge2 = ge1
            else
               ge2 = geslice(i,j)
            endif
c
c           If pressure is below the minimum, set it to the minimum
c
            ge2 = max(ge2, pmin/((gamma - 1.0)*dslice(i,j)))
c
c           Update the total energy
c
            eslice(i,j) = eslice(i,j) - ge1 + ge2
c
c           Compute the pressure with the total energy derived gas energy
c
            pslice(i,j) = (gamma - 1.0)*dslice(i,j)*ge2
c
         enddo
      enddo
c
      return
      end
c=======================================================================
c////////////////////////  SUBROUTINE PGAS2D  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine pgas2d(dslice, eslice, pslice, uslice, vslice,
     &                  wslice, idim, jdim, i1, i2, j1, j2,
     &                  gamma, pmin)
c
c  COMPUTES GAS PRESSURE ON A SLICE
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  July, 1994; by Greg Bryan; changed to a slice
c
c  PURPOSE: Computes gas pressure needed for interpolations for PPM
c           solver.  EOS used is p=(gamma-1)*d*[E-(u**2+v**2+w**2)/2],
c          (eqn 2.1)
c          Only works for gamma law gas with gamma.gt.0
c          All quantities are zone centered (defined over nxz,nyz,nzz).
c
c  INPUTS:
c     dslice  - density slice
c     eslice  - total specific energy slice
c     gamma   - gamma law parameter
c     i1,i2   - starting and ending addresses for dimension 1
c     idim    - declared leading dimension of slices
c     j1,j2   - starting and ending addresses for dimension 2
c     jdim    - declared second dimension of slices
c     pmin    - minimum allowed pressure
c     uslice  - velocity-1 slice
c     vslice  - velocity-2 slice
c     wslice  - velocity-3 slice
c
c  OUTPUT ARGUMENTS: 
c     pslice  - pressure slice
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer i1, i2, idim, j1, j2, jdim
      real    gamma, pmin
      real dslice(idim,jdim), eslice(idim,jdim), pslice(idim,jdim),
     &     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim)
c
c  locals
c
      integer i, j
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c  compute the pressure
c
      do j=j1,j2
         do i=i1,i2
            pslice(i,j) =
     &           (gamma-1.0)*dslice(i,j)*(
     &              eslice(i,j) - 0.5*(uslice(i,j)**2 +
     &                                 vslice(i,j)**2 +
     &                                 wslice(i,j)**2  )
     &                                   )
            if (pslice(i,j) .lt. pmin) pslice(i,j) = pmin
         enddo
      enddo
c
      return
      end


c=======================================================================
c/////////////////////////  SUBROUTINE PPM_DE  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine ppm_de(d, e, u, v, w, ge,
     &                  gravity, gr_xacc, gr_yacc, gr_zacc,
     &                  gamma, dt, nhy, dx, dy, dz,
     &                  rank, in, jn, kn, start, pend,
     &                  iflatten, ipresfree,
     &                  idiff, isteepen, idual, eta1, eta2,
     &                  nsubgrids, lface, rface,
     &                  fistart, fiend, fjstart, fjend,
     &                  array, dindex, eindex,
     &                  uindex, vindex, windex, geindex, tmp,
     &                  ncolor, colorpt, coloff, colindex)
c
c  PERFORMS PPM (DIRECT EULERIAN) UPDATE FOR ONE TIMESTEP (WRAPPER)
c
c  written by: Greg Bryan
c       (some of this is based on code originally written by Jim Stone)
c  date:       June, 1993
c  modified1:  November, 1994; by GB
c              Modified to work as an entirely callable routine within
c              an AMR package.  There is now only one exterior parameter:
c              MAX_ANY_SINGLE_DIRECTION.
c  modified2:  Summer(?), 2002; by Robert Harkness
c              Extended sweeps into ghost zones to update boundary conditions
c              within the directionally split scheme.
c  modified3:  July, 2003; by Alexei Kritsuk
c              Corrected permutation pattern for 2D runs. 
c
c  PURPOSE:  Advances the fluid equations by one timestep using the
c    PPM Direct Eulerian algorithm in 3-D.  This routine is a callable
c    wrapper for the 1d sweeps.  It assumes (in classic fortran style)
c    that the arguments are defined to dimension three; however, it does
c    the right thing if called with one or more collapsed dimensions
c    (i.e. with the width set to 1).
c
c  EXTERNALS:
c    x,y,zeuler_sweep - routines to compute the Eulerian step in
c                       one dimension
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     pend     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (unused until modification3)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer gravity, idiff, idual, iflatten, isteepen, nhy, rank,
     &        ipresfree, pend(3), in, jn, kn, nsubgrids, start(3),
     &        ncolor, coloff(ncolor)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real d(in,jn,kn), e(in,jn,kn), u(in,jn,kn),
     &     v(in,jn,kn), w(in,jn,kn),ge(in,jn,jn),
     &     gr_xacc(in,jn,kn), gr_yacc(in,jn,kn), gr_zacc(in,jn,kn),
     &     dx(in), dy(jn), dz(kn)
      real dt, eta1, eta2, gamma, pmin
      real array(1), colorpt(1)
c
c  Parameters
c
c
c  Locals
c
      integer i, ie, is, ixyz, j, je, js, k, ke, ks,
     &        n, nxz, nyz, nzz, ms
      integer i1, i2, j1, j2, k1, k2
      real tmp(1)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error check
c
      if (max(in,jn,kn) .gt. MAX_ANY_SINGLE_DIRECTION) then
         write(6,*) 'PPM_DE: A grid dimension is too long.'
         write(6,*) '   (increase max_any_single_direction.)'
         ERROR_MESSAGE
      endif
c
c  Convert arguments to usable form
c
      is = start(1) + 1
      js = start(2) + 1
      ks = start(3) + 1
      ie = pend(1) + 1
      je = pend(2) + 1
      ke = pend(3) + 1
      ms = max(in*jn, jn*kn, kn*in)
      nxz = ie-is+1
      nyz = je-js+1
      nzz = ke-ks+1

      i1 = 1
      i2 = in
      j1 = 1
      j2 = jn
      k1 = 1
      k2 = kn

!     i1 = is
!     i2 = ie
!     j1 = js
!     j2 = je
!     k1 = ks
!     k2 = ke



!     write(*,'("AAARRGS",9i5)') is,ie,in,js,je,jn,ks,ke,kn
c
c  Set minimum pressure (better if it were a parameter)
c
      pmin = tiny
c
c  Loop over directions, using a Strang-type splitting
c
!      ixyz = mod(nhy,3)
!      do n=ixyz,ixyz+2
      ixyz = mod(nhy,rank)
      do n=ixyz,ixyz+rank-1
c
c  Update in x-direction
c
!        if (mod(n,3) .eq. 0 .and. nxz .gt. 1) then
        if (mod(n,rank) .eq. 0 .and. nxz .gt. 1) then
c
c         write (6,*) 'SOLVER: starting x-sweep'
c
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(k)
!         do k=ks, ke
          do k=k1, k2
            call xeuler_sweep(k, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_xacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7),
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
c  Update in y-direction
c
!        if (mod(n,3) .eq. 1 .and. nyz .gt. 1) then
        if (mod(n,rank) .eq. 1 .and. nyz .gt. 1) then
c
c         write (6,*) 'SOLVER: starting y-sweep'
c
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(i)
!         do i=is, ie
          do i=i1, i2
            call yeuler_sweep(i, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_yacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7),
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
c  Update in z-direction
c
!        if (mod(n,3) .eq. 2 .and. nzz .gt. 1) then
        if (mod(n,rank) .eq. 2 .and. nzz .gt. 1) then
c
c         write (6,*) 'SOLVER: starting z-sweep'
c
c*$* ASSERT CONCURRENT CALL
c$DOACROSS LOCAL(j)
!         do j=js, je
          do j=j1, j2
            call zeuler_sweep(j, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_zacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7),
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
      enddo
c
      return
      end

c=======================================================================
c///////////////////////  SUBROUTINE TWOSHOCK  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine twoshock(
     &           dls, drs, pls, prs, uls, urs, 
     &           idim, jdim, i1, i2, j1, j2, dt, gamma, pmin, ipresfree,
     &           pbar, ubar, gravity, grslice, idual, eta1
     &                   )

c  RIEMANN SOLVER USING TWO-SHOCK APPROXIMATION
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  July, 1994 by Greg Bryan; made slicewise version
c  modified2:  March, 1995 by Greg Bryan; added gravity
c  modified3:  June 2003 by Alexei Kritsuk - wrong sign in PPM!
c  modified4:  February 2004 by Robert Harkness 
c              Iteration closed to tolerance
c
c  PURPOSE:
c     This routine performs an approximate Riemann solution for each
c     zone edge described by dls, drs (etc) where dls(i) is the density
c     on the left side of zone interface (i), and drs(i) is the density
c     on the right side of zone interface (i).  Note that this differs
c     from the zone-centered way of indexing things.
c
c  INPUT:
c    dt     - timestep in problem time
c    eta1   - (dual) selection parameter for total energy (typically ~0.001)
c    gamma  - ideal gas law constant
c    gravity - gravity flag (0 = off, 1 = on)
c    grslice - acceleration in this direction in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idual  - dual energy formalism flag (0 = off)
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    pmin   - minimum allowed pressure
c    ul,rs  - 1-velocity at left and right edges of each cell
c    
c  OUTPUT:
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    
c  EXTERNALS:
c
c  LOCALS:
c    tolerance - close iteration to tolerance
c
c  PARAMETERS:
c    numiter - max number of Newton iterations to perform [8]
c
c-----------------------------------------------------------------------
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, i1, i2, idim, idual, ipresfree, j1, j2, jdim
      real    dt, eta1, gamma, pmin
      real      dls(idim,jdim),     drs(idim,jdim),    pbar(idim,jdim),
     &          pls(idim,jdim),     prs(idim,jdim),    ubar(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim), grslice(idim,jdim)
c
c  local declarations
c
      integer i, j, n

c  NOTE: mturk promotes the following declarations explicitly to real*8 in
c  order to prevent overflows.  I do not do this here because it is not in 
c  accordance with our general good programming practices.  However, if you
c  are getting PPM over/underflow problems and are compiling with 32-bit 
c  internal precision, it may be necessary to explicitly promote this to
c  real*8.  [BWO v1_1_6]

      real     cl(ijkn),    cr(ijkn), dpdul(ijkn), dpdur(ijkn),
     &         ps(ijkn),    ubl(ijkn),   ubr(ijkn),
     &         qa      ,     zl(ijkn),    zr(ijkn)

      integer itok, jjj
      integer mask(ijkn), jfail(ijkn)
      real    old_ps(ijkn), delta_ps(ijkn)
      real    fails(ijkn)

!     To reduce this to the old fixed iteration count method
!     simply set tolerance = 1.0e-50

!     Tolerance can be relaxed to reduce the cost of this routine

#ifdef r4
      real tolerance
      parameter (tolerance = 1.0e-07)
#endif

#ifdef r8
      real tolerance
      parameter (tolerance = 1.0e-14)
#endif

c
c  parameters
c
      integer numiter
      parameter( numiter = 8 )
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
c=======================================================================

c     write(6,*) 'TWOSHOCK: dt =',dt,'  idim =',idim,'  jdim=',jdim
c     write(6,*) 'TWOSHOCK: i1 =',i1,'  i2   =',i2
c     write(6,*) 'TWOSHOCK: j1 =',j1,'  j2   =',j2

      qa = (gamma + 1.0)/(2.0*gamma)
c
c  Loop through this slice, doing 1 sweep at a time
c
      do j=j1,j2
c
c  If pressfree conditions are needed, set pbar to zero and ubar to the
c     average of left and right velocity states.
c
       if (ipresfree .eq. 1) then
        do i=i1, i2
           pbar(i,j) = pmin
           ubar(i,j) = 0.5*(uls(i,j)+urs(i,j))
           pls(i,j)  = pmin
           prs(i,j)  = pmin
        enddo
c
c  Otherwise, solve Riemann problem
c
       else
c
c  First guess at pbar and left- and right-ubar
c     (or is it a plus? :    &     + cr(i)*cl(i)*(uls(i,j)-urs(i,j))
c  it is a plus: see van Leer 1979, JCP 32, 101 (eq. 60 on page 109).
c
c       if (gravity .eq. 1) then
       if (gravity .eq. 99) then
        do i=i1,i2
          cl  (i) = sqrt(gamma*pls(i,j)*dls(i,j))
          cr  (i) = sqrt(gamma*prs(i,j)*drs(i,j))
          ps(i) = (cr(i)*pls(i,j) + cl(i)*prs(i,j)
     &          +  cr(i)*cl(i)*(uls(i,j) - urs(i,j)
     &                        + 0.5*dt*(grslice(i-1,j)-grslice(i,j)) 
     &                                           ) )/(cr(i)+cl(i))

!               ^
!               changed to + June 2003

          if (ps(i) .lt. pmin) ps(i) = pmin
        enddo
       else
        do i=i1,i2

c     Changed order of operations to avoid overflow [BWO v1_1_6]

          cl  (i) = sqrt(gamma*pls(i,j)) * sqrt(dls(i,j))
          cr  (i) = sqrt(gamma*prs(i,j)) * sqrt(drs(i,j))

          ps(i) = (cr(i)*pls(i,j) + cl(i)*prs(i,j)
     &          +  cr(i)*cl(i)*(uls(i,j) - urs(i,j) 
     &                                           ) )/(cr(i)+cl(i))

c               ^
c               changed to + June 2003

            if (ps(i) .lt. pmin) ps(i) = pmin
c          if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))
c          if (ps(i) .lt. pls(i,j) .and. ps(i) .lt. prs(i,j))
c     &         ps(i) = min(pls(i,j), prs(i,j))
        enddo
       endif
c
c  Newton iterations to compute succesive guesses for pbar, l and r ubar
c    (there are some gravity changes in here but they have not really
c     been worked out yet).
c

       do i=i1,i2
          old_ps(i) = ps(i)
          delta_ps(i) = ps(i)
          mask(i) = 1
       end do

       do n=2, numiter

                                !       zl, zr depend on ps only
                                !       ubl, ubr depend on ps only
                                !       dpdul, dpdur depend on ps only

          do i=i1,i2
             if ( mask(i) .gt. 0 ) then
                zl (i) = cl(i)*sqrt((1.0+qa*(ps(i)/pls(i,j)-1.0)))
                zr (i) = cr(i)*sqrt((1.0+qa*(ps(i)/prs(i,j)-1.0)))
                ubl(i) = uls(i,j) - (ps(i)-pls(i,j))/zl(i)
                ubr(i) = urs(i,j) + (ps(i)-prs(i,j))/zr(i)
             end if
          enddo

          do i=i1,i2
             if ( mask(i) .gt. 0 ) then
                dpdul(i) = ((-4.0*zl(i)/dls(i,j))*zl(i)**2)
     &               /      ( 4.0*zl(i)*(zl(i)/dls(i,j))
     &               -   (gamma+1.0)*(ps(i)-pls(i,j)))
                dpdur(i) = (( 4.0*zr(i)/drs(i,j))*zr(i)**2)
     &               /      ( 4.0*zr(i)*(zr(i)/drs(i,j))
     &               -   (gamma+1.0)*(ps(i)-prs(i,j)))
                ps (i) = ps(i) + (ubr(i)-ubl(i))*(dpdur(i)*(dpdul(i)
     &               /(dpdur(i)-dpdul(i))))

                if (dpdur(i).gt.1e30)then
                   write(6,*) 'twoshock: dpdur:',zr(i),drs(i,j),
     &                  ps(i),prs(i,j)
                endif

c     if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))

                if (ps(i) .lt. pmin) ps(i) = pmin

                delta_ps(i) = ps(i) - old_ps(i)
                old_ps(i) = ps(i)
                if ( abs(delta_ps(i) / ps(i)) .lt. tolerance ) then
                   mask(i) = 0
                end if
             end if
          enddo

       enddo

!       itok = 0

!       do i=i1,i2
!       if ( mask(i) .gt. 0 ) then
!         if ( abs(delta_ps(i) / ps(i)) .gt. 10.0 * tolerance ) then
!           itok = itok+1
!         end if
!       end if
!       end do

!       jjj = 0

!       if ( itok /= 0 ) then
!         write(0,'("Riemann failed ",i4," times")') itok
!         do i=i1,i2
!           if( mask(i) .gt. 0 ) then
!             jjj = jjj + 1
!             fails(jjj) = delta_ps(i) / ps(i)
!             jfail(jjj) = i
!           end if
!         end do
!         write(0,'((8(1x,i4,1pe11.3)))') (jfail(i), fails(i), i=1,jjj)
!       end if

c
c  Compute final values of resolved state
c
        do i=i1,i2
          if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))
          pbar(i,j) = ps(i)
          ubar(i,j) = ubl(i) + (ubr(i)-ubl(i))*dpdur(i)/
     &                         (dpdur(i)-dpdul(i))
        enddo
c
c  Dual energy formalism: if sound speed squared from pbar is less than 
c     eta1*v^2 then discard pbar and use average of pls/prs.
c
#ifdef UNUSED
        if (idual .eq. 1) then
           do i=i1, i2
              if (gamma*ps(i)/min(dls(i,j),drs(i,j)) .lt.
     &             eta1*ubar(i,j)**2) then
                 pbar(i,j) = 0.5*(pls(i,j) + prs(i,j))
                 ubar(i,j) = 0.5*(uls(i,j) + urs(i,j))
              endif
           enddo
        endif
#endif /* UNUSED */
c
       endif
c
      enddo
c
      return
      end


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
      real    dt, eta1, eta2, gamma, pmin
      real    d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn),
     &        dx(in), dy(jn), dz(kn)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real    array(1), colorpt(1)
c
c  define local slices (passed as temps)
c
      real    dls(in,jn),    drs(in,jn),flatten(in,jn),
     &       pbar(in,jn),    pls(in,jn),
     &        prs(in,jn), pslice(in,jn),   ubar(in,jn),
     &        uls(in,jn),    urs(in,jn),    vls(in,jn),
     &        vrs(in,jn),    wls(in,jn),    wrs(in,jn),
     &   diffcoef(in,jn),     df(in,jn),     ef(in,jn),
     &         uf(in,jn),     vf(in,jn),     wf(in,jn),
     &        gef(in,jn),   gels(in,jn),   gers(in,jn)
      real dslice(in,jn), eslice(in,jn), grslice(in,jn),
     &     uslice(in,jn), vslice(in,jn),  wslice(in,jn), 
     &    geslice(in,jn)
      real colslice(in,jn,ncolor),  colf(in,jn,ncolor),
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


c=======================================================================
c/////////////////////  SUBROUTINE YEULER_SWEEP  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine yeuler_sweep(i, d, e, u, v, w, ge, in, jn, kn,
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
c  CONTROL ROUTINE FOR Y-SWEEP DIRECT EULER VERSION OF PPM
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
c    i      - current slice position in x-direction
c    idiff  - diffusion flag (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    ie,je,ke - field active zone end index
c    iflatten - flattening flag (0 = off)
c    in,jn,kn - field dimensions
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    is,js,ks - field active zone start index
c    isteepen - steepening flag (0 = off)
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
c  LOCALS: (passed as temporaries)
c    diffcoef - diffusion coefficient in slice i
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
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, idiff, idual, iflatten, ipresfree, isteepen, 
     &        i, in, jn, kn, is, ie, js, je, ks, ke, nsubgrids,
     &        ncolor, coloff(ncolor)
      real    dt, eta1, eta2, gamma, pmin
      real    d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn),
     &        dx(in), dy(jn), dz(kn)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real    array(1), colorpt(1)
c
c  define local slices (passed as temps)
c
      real     dls(jn,kn),    drs(jn,kn),flatten(jn,kn),
     &        pbar(jn,kn),    pls(jn,kn),
     &         prs(jn,kn), pslice(jn,kn),   ubar(jn,kn),
     &         uls(jn,kn),    urs(jn,kn),    vls(jn,kn),
     &         vrs(jn,kn),    wls(jn,kn),    wrs(jn,kn),
     &    diffcoef(jn,kn),     df(jn,kn),     ef(jn,kn),
     &          uf(jn,kn),     vf(jn,kn),     wf(jn,kn),
     &         gef(jn,kn),   gels(jn,kn),   gers(jn,kn)
      real  dslice(jn,kn), eslice(jn,kn),grslice(jn,kn),
     &      uslice(jn,kn), vslice(jn,kn), wslice(jn,kn),
     &     geslice(jn,kn)
      real colslice(jn,kn,ncolor),  colf(jn,kn,ncolor),
     &        colls(jn,kn,ncolor), colrs(jn,kn,ncolor)
c
c  parameters
c
      integer ysweep
      parameter (ysweep = 2)
c
c  locals
c
      integer idim, ic, j, k, n, nxz, nyz, nzz, offset
      integer k1, k2
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Compute number of active zones
c
      nxz = ie - is + 1
      nyz = je - js + 1
      nzz = ke - ks + 1

      k1 = 1
      k2 = kn

!     k1 = ks
!     j2 = ke

c
c Copy from field to slice
c    (note that we have permuted uvw here and not in the routine calls)
c
      do k=k1, k2
         do j=1, jn
            dslice(j,k) = d(i,j,k)
            eslice(j,k) = e(i,j,k)
            uslice(j,k) = v(i,j,k)
            vslice(j,k) = w(i,j,k)
            wslice(j,k) = u(i,j,k)
         enddo
         if (gravity .eq. 1) then
            do j=1, jn
               grslice(j,k) = gr_acc(i,j,k)
            enddo
         endif
         if (idual .eq. 1) then
            do j=1, jn
               geslice(j,k) = ge(i,j,k)
            enddo
         endif
         do ic=1, ncolor
            do j=1, jn
               colslice(j,k,ic) = 
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
     &                    jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &   call calcdiss(
     &            dslice, eslice, uslice, w, u, pslice, dy, dz, dx,
     &            jn, kn, in, js, je, k1, k2, i, nxz, ysweep,
c     &            jn, kn, in, js, je, ks, ke, i, nxz, ysweep,
     &            in, jn, kn, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                )
c
c  Compute Eulerian left and right states at zone edges via interpolation
c
      call inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dy, flatten,
     &            jn, kn, js, je, k1, k2, idual, eta1, eta2,
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
     &            jn, kn, js, je+1, k1, k2, 
     &            dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )
c
c  Compute Eulerian fluxes and update zones-centered quantities
c
      call euler   (
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dy, diffcoef,
     &            jn, kn, js, je, k1, k2, dt, 
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
     &                    jn, kn, js-3, je+3, k1, k2, gamma, pmin)
c
c  Check this slice against the list of subgrids 
c     (all subgrid quantities are zero based)
c    Note that uf/vf/wf are switched to match u/v/wslice.
c
      do n=0, nsubgrids-1
        if (i .ge. fistart(n*3+2)+1 .and. i .le. fiend(n*3+2)+1) then
          idim = fiend(n*3+2) - fistart(n*3+2) + 1
          do k=fjstart(n*3+2)+1, fjend(n*3+2)+1
             offset = i-fistart(n*3+2) + (k-fjstart(n*3+2)-1)*idim
             array(dindex(n*6+3)+offset) = df(lface(n*3+2)+1, k)
             array(dindex(n*6+4)+offset) = df(rface(n*3+2)+2, k)
             array(eindex(n*6+3)+offset) = ef(lface(n*3+2)+1, k)
             array(eindex(n*6+4)+offset) = ef(rface(n*3+2)+2, k)
             if (nxz .gt. 1) then
                array(uindex(n*6+3)+offset) = wf(lface(n*3+2)+1, k)
                array(uindex(n*6+4)+offset) = wf(rface(n*3+2)+2, k)
             endif
             array(vindex(n*6+3)+offset) = uf(lface(n*3+2)+1, k)
             array(vindex(n*6+4)+offset) = uf(rface(n*3+2)+2, k)
             if (nzz .gt. 1) then
                array(windex(n*6+3)+offset) = vf(lface(n*3+2)+1, k)
                array(windex(n*6+4)+offset) = vf(rface(n*3+2)+2, k)
             endif
             if (idual .eq. 1) then
                array(geindex(n*6+3)+offset) = gef(lface(n*3+2)+1, k)
                array(geindex(n*6+4)+offset) = gef(rface(n*3+2)+2, k)
             endif
             do ic=1, ncolor
                array(colindex(n*6+3,ic)+offset) = 
     &                                       colf(lface(n*3+2)+1, k, ic)
                array(colindex(n*6+4,ic)+offset) = 
     &                                       colf(rface(n*3+2)+2, k, ic)
             enddo
          enddo
        endif
      enddo
c
c Copy from slice to field
c
      do k=k1, k2
         do j=1, jn
            d(i,j,k) = dslice(j,k)
            e(i,j,k) = eslice(j,k)
            v(i,j,k) = uslice(j,k)
            w(i,j,k) = vslice(j,k)
            u(i,j,k) = wslice(j,k)
         enddo
         if (idual .eq. 1) then
            do j=1, jn
               ge(i,j,k) = geslice(j,k)
            enddo
         endif
         do ic=1, ncolor
            do j=1, jn
               colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i) = 
     &               colslice(j,k,ic)
            enddo
         enddo
      enddo
c
      return
      end


c=======================================================================
c/////////////////////  SUBROUTINE ZEULER_SWEEP  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine zeuler_sweep(j, d, e, u, v, w, ge, in, jn, kn,
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
c  CONTROL ROUTINE FOR Z-SWEEP DIRECT EULER VERSION OF PPM
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
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, idiff, idual, iflatten, isteepen, ipresfree,
     &        in, jn, kn, is, ie, j, js, je, ks, ke, nsubgrids,
     &        ncolor, coloff(ncolor)
      real    dt, eta1, eta2, gamma, pmin
      real    d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn),
     &        dx(in), dy(jn), dz(kn)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real    array(1), colorpt(1)
c
c  define local slices (passed as temps)
c
      real      dls(kn,in),    drs(kn,in),flatten(kn,in),
     &         pbar(kn,in),    pls(kn,in),
     &          prs(kn,in), pslice(kn,in),   ubar(kn,in),
     &          uls(kn,in),    urs(kn,in),    vls(kn,in),
     &          vrs(kn,in),    wls(kn,in),    wrs(kn,in),
     &     diffcoef(kn,in), dslice(kn,in), eslice(kn,in),
     &       uslice(kn,in), vslice(kn,in), wslice(kn,in),
     &      geslice(kn,in),    gef(kn,in),
     &         gels(kn,in),   gers(kn,in),
     &           df(kn,in),     ef(kn,in),grslice(kn,in),
     &           uf(kn,in),     vf(kn,in),     wf(kn,in)
      real colslice(kn,in,ncolor),  colf(kn,in,ncolor),
     &        colls(kn,in,ncolor), colrs(kn,in,ncolor)
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

!     i1 = is
!     i2 = ie

c
c Copy from field to slice
c    (note that we have permuted uvw here and not in the routine calls)
c
      do i=i1, i2
         do k=1, kn
            dslice(k,i) = d(i,j,k)
            eslice(k,i) = e(i,j,k)
            uslice(k,i) = w(i,j,k)
            vslice(k,i) = u(i,j,k)
            wslice(k,i) = v(i,j,k)
         enddo
         if (gravity .eq. 1) then
            do k=1, kn
               grslice(k,i) = gr_acc(i,j,k)
            enddo
         endif
         if (idual .eq. 1) then
            do k=1, kn
               geslice(k,i) = ge(i,j,k)
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
         call pgas2d_dual(dslice, eslice, geslice, pslice, 
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    kn, in, ks-3, ke+3, i1, i2, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               kn, in, ks-3, ke+3, i1, i2, gamma, pmin)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &   call calcdiss(
     &            dslice, eslice, uslice, u, v, pslice, dz, dx, dy,
     &            kn, in, jn, ks, ke, i1, i2, j, nyz, zsweep,
c     &            kn, in, jn, ks, ke, is, ie, j, nyz, zsweep,
     &            in, jn, kn, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                )
c
c  Compute Eulerian left and right states at zone edges via interpolation
c
      call inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dz, flatten,
     &            kn, in, ks, ke, i1, i2, idual, eta1, eta2,
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
     &            kn, in, ks, ke+1, i1, i2, 
     &            dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )
c
c  Compute Eulerian fluxes and update zones-centered quantities
c
      call euler   (
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dz, diffcoef,
     &            kn, in, ks, ke, i1, i2, dt, 
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
     &                    kn, in, ks-3, ke+3, i1, i2, gamma, pmin)
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
