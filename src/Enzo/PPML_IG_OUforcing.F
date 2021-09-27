c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

C ************************************************************************
C *** Ornstein-Uhlenbeck forcing source terms for 3D turbulence simulation
C *** Adopted by Dmitry Kotov from AREPO code developed by Volker Springel
C ************************************************************************

      subroutine OUpumpInit(gamma,rho0,P0,SolWeight,Mach,kmin,kmax,argx)
      implicit none

      include 'PPML_IG_OUforcing.h'

      integer :: ikmax, ikx, iky, ikz, N, i, j
      integer :: myid
      ENZO_REAL :: gamma, rho0, P0, Mach, Eturb, c_s
      ENZO_REAL :: kmax, kmin, boxSize, kx, ky, kz, kk, kc, amp
      ENZO_REAL :: RandomPhase, argx(3)
      ENZO_REAL, parameter :: PI_2 = 3.14159265d0*2.d0 
      
c      call GETDATA("gamma", gamma, argi)
c      call GETDATA("Lx",argx(1),argi)
c      call GETDATA("Ly",argx(2),argi)
c      call GETDATA("Lz",argx(3),argi)

c      if ((argx(1) .ne. argx(2)).or.(argx(2).ne.argx(3))) 
c     &  stop "Box sizes must be equal"

      boxSize = argx(1)                           ! default 1.0

c      call GETDATA("SolWeight", SolWeight, argi)  ! default 1.0
c      call GETDATA("Mach", Mach, argi)            ! default 0.2
c      call GETDATA("kfi", kmin, argi)             ! default 6.27
c      call GETDATA("kfa", kmax, argi)             ! default 12.57
      
      call initrandom()
      
!AK what's P0? mean pressure? same Q for rho0
c      P0      = 1.d0
c      rho0    = 1.d0

!AK totTime is time counter used to control the addition of new random phases (should be kept for restarts).
      totTime = 0d0
      c_s     = sqrt(gamma*P0/rho0)
!AK If this is the so-called dynamical time, then it corresponds to a force acting on Lbox/2.
!AK A more general expression would use 2pi/kc as the spatial scale determining the time correlation.
! TDecay is the correlation time (= 1/2 of eddy-turnover time scale [at what spatial scale?]) 
! Autocorrelation time of the force TDecay [same as the integral time scale (T=Linjection/Velo_injection)]
      TDecay  = boxSize/(2.d0*Mach*c_s)
!AK DtFreq time interval between the updates of phases (10 updates per one correlation time).
      DtFreq  = 0.1d0*TDecay
!AK ideally, kc should be a weighted average of all wave numbers from the spherical shell [kmin,kmax].
!AK But in this case a centered parabola is used as a weight, so this should be ok.
      kc      = 0.5d0*(kmin + kmax)
!AK Why division by kc instead of multiplication? This makes no sense to me.
!AK should be: Eturb = 0.5*(Mach*c_s)**3/(2*pi/kc)=0.159*0.5*(Mach*c_s)**3*kc~0.006 for Mach=0.2 and kc~10
      Eturb   = 0.119d0 * 0.5d0*(Mach*c_s)**3/kc     ! default 5.d-5 (?)
!AK The variance of the OU process =sigma^2 (sigma=kc*(Mach*c_s)**2/sqrt(2*pi)~0.2 for Mach=0.2).
      OUVar   = sqrt(Eturb/TDecay)

!8/31/16 12:18 PM, Andreas Bauer wrote:
!ST_decay        2.5
!ST_energy       0.00005
!ST_DtFreq       0.25
!ST_Kmin         6.27
!ST_Kmax         12.57
!ST_SolWeight    1.
!ST_AmplFac      1.
!ST_Seed         42
!ST_SpectForm    1
!sigma is given by StOUVar = sqrt(All.StEnergy / All.StDecay)
!The used spectral form is parameterization 1 in our implementation, i.e. a paraboloid from k_min to k_max.

c      if(myid.eq.0) then
c        print *, 'c_s = ',    c_s
c        print *, 'TDecay = ', TDecay
c        print *, 'DtFreq = ', DtFreq
c        print *, 'OUVar = ',  OUVar
c        print *, 'Eturb = ',  Eturb
c        print *, 'kc = ',     kc
c      endif

      ikmax = boxSize*kmax/PI_2
      
      NModes = 0
      do ikx = 0, ikmax
         kx = PI_2*ikx/boxSize
         do iky = 0, ikmax
            ky = PI_2*iky/boxSize
            do ikz = 0, ikmax
               kz = PI_2*ikz/boxSize
        
               kk = sqrt(kx*kx + ky*ky + kz*kz)
        
               if((kk.ge.kmin).and.(kk.le.kmax)) then
                  NModes = NModes + 4
               endif
            
            enddo
         enddo
      enddo
      
      allocate(Mode(3,NModes))
      allocate(OUPhase(2,3,NModes), Ampl(NModes), Phase(2,3,NModes))

!AK correction for projecting out non-solenoidal part, see Schmidt, Hillebrandt & Niemeyer [2005, p. 706, eq. (40)]
      WeightNorm =sqrt(3.d0)/sqrt(1.d0-2.d0*SolWeight+3.d0*SolWeight**2)
      
      N = 0
      do ikx = 0, ikmax
         kx = PI_2*ikx/boxSize
         do iky = 0, ikmax
            ky = PI_2*iky/boxSize
            do ikz = 0, ikmax
               kz = PI_2*ikz/boxSize
        
               kk = sqrt(kx*kx + ky*ky + kz*kz)
        
               if((kk.ge.kmin).and.(kk.le.kmax)) then
                  amp = 1.d0 - 4.d0*(kk-kc)**2/(kmax-kmin)**2
          
                  N = N + 1
                  Ampl(N) = amp
                  Mode(1,N) = kx
                  Mode(2,N) = ky
                  Mode(3,N) = kz
          
                  N = N + 1
                  Ampl(N) = amp
                  Mode(1,N) = kx
                  Mode(2,N) = -ky
                  Mode(3,N) = kz

                  N = N + 1
                  Ampl(N) = amp
                  Mode(1,N) = kx
                  Mode(2,N) = ky
                  Mode(3,N) = -kz

                  N = N + 1
                  Ampl(N) = amp
                  Mode(1,N) = kx
                  Mode(2,N) = -ky
                  Mode(3,N) = -kz
               endif
               
            enddo
         enddo
      enddo

!AK if solution needs to be written, write random seeds as well
c      call GETDATA('readinsol', argx, argi)

      if (argi(1) .ne. 0) then
        open (42, file='OUturb.bin', form='UNFORMATTED', 
     &      access='STREAM', action='read')
        read (42) argi(1)
        if (argi(1) .ne. NModes) then
          stop "wrong number of modes in OUturb.bin"
        endif
        print *, 'reading phases for ', NModes, ' modes'
        read (42) OUPhase, Ampl
        close (42)
      else
        do j = 1, NModes
           do i = 1, 3
              OUPhase(1,i,j) = RandomPhase()*OUVar
              OUPhase(2,i,j) = RandomPhase()*OUVar
           enddo
        enddo
        
        call calcPhases()
      endif
      
      end subroutine

C **********************************************************************      
      subroutine cleanupTurbForcing3D(myid)
      implicit none
      include 'PPML_IG_OUforcing.h'
      integer :: myid

      if (myid.eq.0) then
        open (42, file='OUturb.bin', form='UNFORMATTED',
     &        access='STREAM')
        write (42) NModes, OUPhase, Ampl
        close (42)
      endif

      deallocate(Mode, OUPhase, Ampl, Phase)
      end subroutine

C **********************************************************************
      ENZO_REAL function RandomPhase()
      implicit none
      ENZO_REAL r0, r1, randomnr
      ENZO_REAL, parameter :: PI_2 = 3.14159265d0*2.d0 
      
      r0 = randomnr()
      r1 = randomnr()
      
      RandomPhase = sqrt(-2.d0*log(r0))*cos(PI_2*r1)
      end function 

C **********************************************************************
      subroutine calcPhases()   
      implicit none
      include 'PPML_IG_OUforcing.h'
      integer :: i, j
      ENZO_REAL :: ka, kb, kk, kk1, diva, divb, curla, curlb
      
      do j = 1, NModes
        ka = 0d0
        kb = 0d0
        kk = 0d0
        
        do i = 1, 3
          kk = kk + Mode(i,j)**2
          ka = ka + Mode(i,j)*OUPhase(2,i,j)
          kb = kb + Mode(i,j)*OUPhase(1,i,j)
        enddo
        
        kk1 = 1.d0/kk
        
!AK This is the projection operator
!AK Compare to eq. (39) from Schmidt et al. (2005, Comb. Theory & Modelling, vol. 9, pp. 693-720)
        do i = 1, 3
          diva = Mode(i,j)*ka*kk1
          divb = Mode(i,j)*kb*kk1
          curla = OUPhase(1,i,j) - divb
          curlb = OUPhase(2,i,j) - diva
          Phase(1,i,j) = SolWeight*curla + (1.d0 - SolWeight)*divb
          Phase(2,i,j) = SolWeight*curlb + (1.d0 - SolWeight)*diva
        enddo
        
      enddo      
      
      end subroutine

C **********************************************************************
      subroutine updatePhases()   
      implicit none
      include 'PPML_IG_OUforcing.h'
      ENZO_REAL :: damping, dfact, RandomPhase
      integer :: i, j, n

!AK damping is the decay factor; TDecay is the correlation length;
!AK DtFreq is the time interval for updating the phases.
      damping = exp(-DtFreq/TDecay)
      dfact = sqrt(1.d0 - damping**2)
      
      do j = 1, NModes
      do i = 1, 3
      do n = 1, 2
        OUPhase(n,i,j) = OUPhase(n,i,j)*damping +
     *                   OUVar*dfact*RandomPhase()
      enddo
      enddo
      enddo
            
      end subroutine      

C **********************************************************************
c     3D Turbulent forcing source term a la Ornstein-Uhlenbeck 
c     Updates res:    res = res + source*dt
c     by DK
c INPUT:
c       nc             - number of components
c       ni, nj, nk     - local domain dimensions
c       nig, njg, nkg  - global domain dimensions
c       w              - flow field
c       grid           - grid coordinates
c       dt             - time step
c       update_sol     - flag to update w
c OUTPUT:
c       res            - residual (updated)
c***********************************************************************
      subroutine turbForce3D(nc, ni, nj, nk, nig, njg, nkg, 
     *                       w, grid, dt, res, update_sol)
      implicit none
      include 'PPML_IG_OUforcing.h'
#ifdef USING_MPI
      include 'mpif.h'
      integer err
      ENZO_REAL gvmpi(4)
#endif            
      integer :: nc, ni, nj, nk, argi(1), ol2, i, j, k, m, dir
      integer :: nig, njg, nkg
      ENZO_REAL :: grid(3,ni,nj,nk), w(nc,ni,nj,nk), res(nc,ni,nj,nk)
      ENZO_REAL :: rho, dt, dumx, acc(3), du(3), gv(4), kr
      logical :: update_sol
      
      if(nk.eq.1) stop "turbForce3D is for 3D only"
      if((nig.ne.njg).or.(nig.ne.nkg)) then
        stop "turbForceNorm2d assume ni=nj"
      endif
      
! boolean switch: return if no need to pump
c      call GETDATA("apply_forcing", dumx, argi)
c      if(argi(1).eq.0) return 
      
!AK overlap (number of ghost zones on one side)
      call GETDATA("olap", dumx, argi)
      ol2 = argi(1)/2

      if(dt*0.5d0 .gt. DtFreq) then ! require at least two dt in DtFreq
        print *, "dt = ", dt, " > DtFreq = ", DtFreq
        stop "dt too big"
      endif
      
      totTime = totTime + dt
      if(totTime .gt. DtFreq) then
        call updatePhases()
        call calcPhases()
        totTime = totTime - DtFreq
      endif
 
      gv = 0d0
      do k=1+ol2,nk-ol2
      do j=1+ol2,nj-ol2
      do i=1+ol2,ni-ol2
        
        acc = 0d0        
        do m=1, NModes
          kr = 0d0
          do dir=1,3
            kr = kr + Mode(dir,m)*grid(dir,i,j,k)
          enddo
!AK why do we need both sin() and cos() and why Phase is not in the arguments?
          acc(1:3) = acc(1:3) + Ampl(m)*( Phase(1,1:3,m)*cos(kr) +
     &                                    Phase(2,1:3,m)*sin(kr) )
          
        enddo
        
        du = acc*(2.d0*WeightNorm)*dt
                
        rho = w(1,i,j,k)
        res(2:4,i,j,k) = res(2:4,i,j,k) + du(1:3)*rho 
        res(5,i,j,k) = res(5,i,j,k) + dot_product(du,du)*0.5d0*rho + 
     *               dot_product(w(2:4,i,j,k), du(1:3))

        gv(1:3) = gv(1:3) + du(1:3)*rho
        gv(4) = gv(4) + rho      
      enddo
      enddo
      enddo
      
#ifdef USING_MPI
! gv is a 4-vector which includes sums for 3 momentum components (1-3) and sum density (4)
      call MPI_AllReduce(gv, gvmpi, 4, MPI_DOUBLE_PRECISION, MPI_SUM,
     *                  MPI_COMM_WORLD, err )
      gv = gvmpi
#endif 

! force total impulse (momentum) increment to zero (for active zones only)
      du = -gv(1:3)/gv(4)
      do k=1+ol2,nk-ol2
      do j=1+ol2,nj-ol2
      do i=1+ol2,ni-ol2
        rho = w(1,i,j,k)
        res(2:4,i,j,k) = res(2:4,i,j,k) + du(1:3)*rho 
        res(5,i,j,k) = res(5,i,j,k) + dot_product(du,du)*0.5d0*rho + 
     *               dot_product(w(2:4,i,j,k), du(1:3))
      enddo
      enddo
      enddo
      
! update solution
      if(update_sol) then
!AK do we need to initialize res(1,...) and res(6:8,...) with zeroes?
        w = w + res
      endif
      
      end subroutine
