c     See LICENSE_ENZO file for license and copyright information

#include "fortran.h"
c=======================================================================
c/////////////////////////  SUBROUTINE TURBOINIT \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine turboinit2d(rank, nbox, u, v, in, jn,  ig, jg)
c
c  PROVIDES QUASI-ISOTROPIC 2D ZERO-MEAN SOLENOIDAL VELOCITIES FOR LARGE-SCALE  
c  TURBULENCE DRIVING (AND/OR INITIAL CONDITIONS) IN A PERIODIC BOX
c
c  written by: Alexei Kritsuk
c
c  date:       August 2014
c
c  PURPOSE: 
c
c For definitions of Multidimensional DFT see
c http://en.wikipedia.org/wiki/Discrete_Fourier_transform
c
c  EXTERNALS:
c
c  INPUTS:
c     in,jn     - dimensions of field arrays
c     ig,jg     - global zone number (for the left-most zone) for each dimension
c     rank      - dimension of problem
c     u         - x-velocity field
c     v         - y-velocity field
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c  Arguments
c
      ENZO_REAL u(in,jn), v(in,jn)
      integer rank, nbox, in, jn, ig, jg
c
c  Locals
c
      integer i, j, imo, m, nmode
      integer(kind=4) seed
      integer ky, kx, kfi, kfa
      ENZO_REAL, dimension (: , :), allocatable :: mode
      ENZO_REAL x, y, pi, km, amp, pha, vel, twopi, delta
      parameter (pi=3.14159265)

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error checks
c
      if (rank .ne. 2) then
         write(6,*) 'TURBOINIT2D: Only 2D grids are supported.'
         write(0,*) 'stop_all_cpus in turboinit2d'
         stop
      endif
c
      if (max(in,jn) .gt. MAX_ANY_SINGLE_DIRECTION) then
         write(6,*) 'TURBOINIT2D: A grid dimension is too long.'
         write(6,*) '   (increase max_any_single_direction.)'
         write(0,*) 'stop_all_cpus in turboinit'
         stop
      endif
c

c      write(6,*) 'TURBOINIT2D: ig =',ig,  ' jg =',jg

c
c     Set min and max wavenumbers (should be a parameter)
c
c      kfi  = nbox/32.0
c      kfa  = nbox/30.0
      kfi  = nint(nbox/24.0)
      kfa  = nint(nbox/23.7)
c
c     Define modes involved
c
      nmode = 0
      do kx=-kfa, kfa 
         do ky=-kfa, kfa 
            km = sqrt(real(kx*kx+ky*ky))
            if (km .lt. kfa .and. km .gt. kfi) nmode = nmode+1
         enddo
      enddo

      if (nmode .lt. 16) then
         write(6,*) 'TURBOINIT2D: Too few modes.',nmode,kfi,kfa,km
         stop
      else
         write(6,*) 'TURBOINIT2D: nmodes =',nmode
      endif

      allocate(mode(rank, nmode))
      nmode = 1
      do kx=-kfa, kfa 
         do ky=-kfa, kfa
            km = sqrt(real(kx*kx+ky*ky))
            if (km .lt. kfa .and. km .gt. kfi) then 
               mode(1,nmode) = kx
               mode(2,nmode) = ky
               nmode = nmode+1
            endif
         enddo
      enddo
c 
c initialize velocity arrays with zeroes
c
      do j=1, jn
         do i=1, in
            u(i,j) = 0.0
            v(i,j) = 0.0
         enddo
      enddo
c
c fill-in velocity arrays with values
c
      seed = 12345
      call srand(seed)
      twopi = pi*2.0
      delta = 1.0/nbox
c          
c use random phases "pha" uniformly sampled from [0, 2*pi)
c and random amplitudes "amp" sampled from a normal distribution (ideally)
c
      do m=1, nmode
         amp = rand()
         pha = rand()*twopi
         km  = sqrt(real(mode(1,m)**2 + mode(2,m)**2))
         if (km .gt. 0) amp = amp/km
         do j=1, jn
            y = delta*(j+jg-1)
            do i=1, in
	       x = delta*(i+ig-1)
               vel = amp*sin(twopi*(mode(1,m)*x + mode(2,m)*y) + pha)
               if (km .gt. 0) then
                  u(i,j) = u(i,j) + vel*mode(2,m)
                  v(i,j) = v(i,j) - vel*mode(1,m)
               endif
            enddo
         enddo
      enddo
c
c normalize to get rms 2D Mach = 1.0 for seed 12345 and given k_f
c
      do j=1, jn
         do i=1, in
	    u(i,j) = u(i,j)/4.2846841
            v(i,j) = v(i,j)/4.2846841
         enddo
      enddo
c
      return
      end
c-------------------------------------------------------------------------
