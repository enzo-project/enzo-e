c#define STATIONARY_FORCING
c***********************************************************************
      subroutine calcTurbAcc2d(rank, nbox, grid, ni, nj, seed, kk, uv)
c
c  PROVIDES QUASI-ISOTROPIC ZERO-MEAN SOLENOIDAL ACCELERATION FOR 2D
c  TURBULENCE DRIVING (AND/OR INITIAL CONDITIONS) IN A PERIODIC BOX
c  IF STATIONARY_FORCING IS DEFINED, THE ACCELERATION IS STATIONARY
c  (I.E. ONLY AMPLITUDE CHANGES WITH TIME SLIGHTLY TO KEEP THE ENERGY
c  INJECTION RATE CONSTANT).
c  IF NOT, THIS SUBROUTINE PROVIDES DELTA-CORRELATED IN TIME FORCING.
c
c  written by: Alexei Kritsuk (adapted for ADPDIS3D by Dmitry Kotov)
c
c  date:       August 2014
c
c  PURPOSE: 
c
c  EXTERNALS:
c
c  INPUT:
c     rank      - dimension of problem
c     nbox      - grid resolution
c     grid      - grid coordinates
c     ni, nj    - dimensions of field array
c     seed      - seed for random numbers generator
c     kk        - ratios to compute kfi = nint(nbox/kk(1)) and kfa = nint(nbox/kk(2))
c  OUTPUT:
c     uv        - x and y-velocity/acceleration field
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c  Arguments
c
      real*8 uv(2,ni,nj), grid(3,ni,nj)
      integer rank, nbox, ni, nj
      integer(kind=8) seed
c
c  Locals
c
      integer i, j, imo, m, nmode, p, h
      integer kx, ky, kfi, kfa
      integer, allocatable :: new_seed (:)
      real*8, dimension (: , :), allocatable :: mode
      real*8 x, y, pi, km, amp, pha, vel, twopi, delta, kk(2)
      parameter (pi=3.14159265)
#ifdef STATIONARY_FORCING
      integer(kind=4) seed1
#endif
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
c      if (max(in,jn) .gt. 1000000) then
c         write(6,*) 'TURBOINIT2D: A grid dimension is too long.'
c         write(6,*) '   (increase max_any_single_direction.)'
c         write(0,*) 'stop_all_cpus in turboinit'
c         stop
c      endif
c

c      write(6,*) 'TURBOINIT2D: ig =',ig,  ' jg =',jg

c
c     Set min and max wavenumbers (should be a parameter)
c
      kfi  = nint(nbox/kk(1))
      kfa  = nint(nbox/kk(2))
c      kfi  = nbox/32.0
c      kfa  = nbox/30.0
c scaling tests
c      kfi  = nint(nbox/14.0)
c      kfa  = nint(nbox/13.3)
c 512
c      kfi  = nint(nbox/24.0)
c      kfa  = nint(nbox/23.7)
c 2048b, 4096a
c      kfi  = nint(nbox/48.0)
c      kfa  = nint(nbox/47.0)
c 2048a and 8192
c      kfi  = nint(nbox/96.0)
c      kfa  = nint(nbox/95.0)
c
c     Define modes involved
c
      nmode = 0
      do kx=-kfa, kfa 
         do ky=-kfa, kfa 
            km = sqrt(real(kx*kx+ky*ky,8))
            if (km .lt. kfa .and. km .gt. kfi) nmode = nmode+1
         enddo
      enddo

      if (nmode .lt. 16) then
      write(6,*) 'TURBOINIT2D: Too few modes.',nmode,kfi,kfa,km
         stop
      endif

      allocate(mode(rank, nmode))
      nmode = 1
      do kx=-kfa, kfa 
         do ky=-kfa, kfa
            km = sqrt(real(kx*kx+ky*ky,8))
	    if (km .lt. kfa .and. km .gt. kfi) then 
               mode(1,nmode) = kx
               mode(2,nmode) = ky
               nmode = nmode+1
            endif
         enddo
      enddo
      nmode = nmode-1
c
c      write(6,*) 'TURBOINIT2D: nmodes =',nmode
c 
c initialize velocity arrays with zeroes
c
      uv = 0d0
c      do j=1, jn
c         do i=1, in
c            u(i,j) = 0.0
c            v(i,j) = 0.0
c         enddo
c      enddo
c
c fill-in velocity arrays with values; if seed is 12345 => steady force
c
#ifdef STATIONARY_FORCING
      seed1 = 12346
      call srand(seed1)
#else
      seed = seed + 1
c
c	p is set to the size of the seed array; usually p=6
c
      call random_seed(size=p)
      allocate (new_seed(p))
c
c	assign values to new_seed(p)
c
      do h=1, p
        new_seed(h) = h*seed	
      enddo 		  
      call random_seed(put=new_seed(1:p))
#endif
      twopi = pi*2.0
      delta = 1.0/nbox
c          
c use random phases "pha" uniformly sampled from [0, 2*pi)
c use random amplitudes "amp" sampled from a normal distribution (ideally)
c
      do m=1, nmode
#ifdef STATIONARY_FORCING
         amp = rand()
         pha = rand()*twopi
#else
         call random_number(amp)
         amp = amp + 0.5
         call random_number(pha)
         pha = pha*twopi
#endif
         km  = sqrt(real(mode(1,m)**2 + mode(2,m)**2,8))
         if (km .gt. 0) amp = amp/km
         do j=1, nj
!            y = delta*(j+jg-1)
            do i=1, ni
	           x = grid(1,i,j) !delta*(i+ig-1)
               y = grid(2,i,j)
               vel = amp*sin(twopi*(mode(1,m)*x + mode(2,m)*y) + pha)
               if (km .gt. 0) then
                  uv(1,i,j) = uv(1,i,j) + vel*mode(2,m)
                  uv(2,i,j) = uv(2,i,j) - vel*mode(1,m)
               endif
            enddo
         enddo
      enddo
c
c normalize to get rms 2D Mach = 1.0 for seed 12345 and given k_f=512/24
c
#ifdef STATIONARY_FORCING
    uv = uv / 4.2846841d0 
c      do j=1, jn
c         do i=1, in
c	    u(i,j) = u(i,j)/4.2846841
c            v(i,j) = v(i,j)/4.2846841
c         enddo
c      enddo
#endif
c
c     clean up.
c
      deallocate(mode)
      deallocate(new_seed)
c
      return
      end
c***********************************************************************
c     2D Turbulent forcing source term
c     Updates res:    res = res + source*dt
c     by DK
c INPUT:
c       nc          - number of components
c       ni, nj, nk  - domain dimensions (nk should be 1)
c       w           - flow field
c       grid        - grid coordinates
c       jac         - jacobian
c       wk          - work array, size 2*ni*nj
c       dt          - time step
c       update_sol  - flag to update w
c OUTPUT:
c       res         - residual (updated)
c***********************************************************************
      subroutine turbForce2D(nc, ni, nj, nk, nig, njg, nkg, 
     *                          w, grid, jac, wk, turbAcc,
     *                          dt, seed, res, update_sol)
      implicit none
#ifdef USING_MPI
      include 'mpif.h'
      integer err
      real*8 gvmpi(5)
#endif      
      integer :: nc, ni, nj, nk, argi(1), nbox, i, j, ol2, nig, njg, nkg
      integer(kind=8) :: seed
      real*8 :: grid(3,ni,nj), w(nc,ni,nj), wk(2,ni,nj), res(nc,ni,nj)
      real*8 :: rho, dt, dumx, ss, kk(2), de, turbAcc(4,ni,nj)
      real*8 :: jac(ni,nj), norm, eps, dru, drv, ruavr, rvavr, gv(5)
      logical :: update_sol
      
      if(nk.ne.1) stop "turbForceNorm2d is for 2D only"
      if(nig .ne. njg) stop "turbForceNorm2d assume ni=nj"
      
      if(argi(1).eq.0) return 
      
      nbox = nig - argi(1) ! nig - olap
      ol2 = argi(1)/2
      
      call calcTurbAcc2d(2, nbox, grid, ni, nj, seed, kk, wk)
 
      gv = 0d0
      do j=1+ol2,nj-ol2
      do i=1+ol2,ni-ol2
c        ss = jac(i,j)
        gv(1) = gv(1) + w(2,i,j)*wk(1,i,j)+w(3,i,j)*wk(2,i,j)  ! <rho*(xacceleration*xvelocity + yacceleration*yvelocity)>
        gv(2) = gv(2) + w(1,i,j)*(wk(1,i,j)**2 + wk(2,i,j)**2)    ! <rho*(xacceleration**2 + yacceleration**2)>
        gv(3) = gv(3) + w(1,i,j)*wk(1,i,j)  ! <rho*xacceleration>
        gv(4) = gv(4) + w(1,i,j)*wk(2,i,j)  ! <rho*yacceleration>
        gv(5) = gv(5) + w(1,i,j)
      enddo
      enddo
      
#ifdef USING_MPI
      call MPI_AllReduce(gv, gvmpi, 5, MPI_DOUBLE_PRECISION, MPI_SUM,
     *                  MPI_COMM_WORLD, err )
      gv = gvmpi
#endif     

      gv(3:4) = gv(3:4)/gv(5) ! average momentum
      
      norm = ( sqrt(gv(1)**2+gv(2)*dt*eps*2.d0*nbox**2) - gv(1) )/gv(2)
            
      do j=1,nj
      do i=1,ni
        rho = w(1,i,j)
        dru = (wk(1,i,j)-gv(3))*norm
        drv = (wk(2,i,j)-gv(4))*norm
        de = (dru**2+drv**2)*0.5d0*rho + w(2,i,j)*dru + w(3,i,j)*drv
        res(2,i,j) = res(2,i,j) + dru*rho
        res(3,i,j) = res(3,i,j) + drv*rho
        res(5,i,j) = res(5,i,j) + de
        turbAcc(4,i,j) = turbAcc(4,i,j) + de
      enddo
      enddo

      if(update_sol) then
      do j=1,nj
      do i=1,ni
        w(1:nc,i,j) = w(1:nc,i,j) + res(1:nc,i,j)
      enddo
      enddo
      endif
      
      end subroutine



