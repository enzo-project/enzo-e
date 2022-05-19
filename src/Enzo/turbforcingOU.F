C     **********************************************************************
C     *** Forcing terms for OU turbulence (2D or 3D)
C     *** Adapted from AREPO code developed by Volker Springel
C     **********************************************************************
      subroutine cello_init_turbulence_OU
     &     (cello_is_root,
     &     cello_rank,
     &     cello_domain_size,
     &     cello_gamma,
     &     cello_apply_injection_rate,
     &     cello_cooling_term,
     &     cello_injection_rate,
     &     cello_kfi,
     &     cello_kfa,
     &     cello_mach,
     &     cello_read_sol,
     &     cello_sol_weight)

      implicit none
      include 'turbforcingOU.h'
C
      integer :: cello_is_root
      integer :: cello_rank
      real*8 :: cello_domain_size(3)
      real*8 :: cello_gamma;
      integer :: cello_apply_injection_rate
      integer :: cello_cooling_term
      real*8 :: cello_injection_rate
      real*8 :: cello_kfi
      real*8 :: cello_kfa
      real*8 :: cello_mach
      integer :: cello_read_sol
      real*8 :: cello_sol_weight
      integer :: modes_input
C
      real*8 :: Mach, gamma, P0, rho0, c_s, Edot
      real*8 :: kmax, kmin, boxSize(3), kx, ky, kz, kk, kc, amp
      integer :: ikmax(3), ikx, iky, ikz, N, i, j
      real*8 :: RandomPhase, argx(3)
      real*8, parameter :: PI_2 = 3.14159265358979323846d0*2.d0
      integer :: myid, seed, modesPerK
      logical :: apply_injection_rate

      apply_injection_rate = (cello_apply_injection_rate .ne. 0)

      Ndims = cello_rank
      gamma = cello_gamma
      boxSize(1) = cello_domain_size(1)
      boxSize(2) = cello_domain_size(2)
      boxSize(3) = cello_domain_size(3)
      if (Ndims.eq.3) then
         modesPerK = 4
         boxVolume = boxSize(1) * boxSize(2) * boxSize(3)
      elseif(Ndims.eq.2) then
         modesPerK = 2
         boxVolume = boxSize(1) * boxSize(2)
      else
         stop "only dim=2 or dim=3 is supported"
      endif

      Edot = cello_injection_rate
      SolWeight = cello_sol_weight
      Mach = cello_mach
      kmin = cello_kfi
      kmax = cello_kfa

      P0 = 1.d0
      totTime = 0d0
      rho0 = 1.d0
      c_s = sqrt(gamma*P0/rho0)
      kc = 0.5d0*(kmin + kmax)
      if (apply_injection_rate) then
         Mach = (PI_2 * Edot / (rho0 * kc))**(1.d0/3.d0) / c_s
      end if
      TDecay = PI_2/(Mach*c_s*kc)
      DtFreq = 0.1d0*TDecay
      OUVar = (Mach*c_s)/TDecay

      ikmax(1:3) = nint(boxSize(1:3)*(kmax/PI_2))
      if (Ndims.eq.2) ikmax(3) = 0

      NModes = 0
      do ikx = 0, ikmax(1)
         kx = PI_2*dble(ikx)/boxSize(1)
         do iky = 0, ikmax(2)
            ky = PI_2*dble(iky)/boxSize(2)
            do ikz = 0, ikmax(3)
               kz = PI_2*dble(ikz)/boxSize(3)

               kk = sqrt(kx*kx + ky*ky + kz*kz)

               if((kk.ge.kmin).and.(kk.le.kmax)) then
                  NModes = NModes + modesPerK
               endif

            enddo
         enddo
      enddo

      allocate(Mode(3,NModes))
      allocate(OUPhase(2,3,NModes), Ampl(NModes), Phase(2,3,NModes))

      if (Ndims.eq.3) then
         WeightNorm =
     &        1.d0/sqrt(1.d0-2.d0*SolWeight+3.d0*SolWeight**2)
      else
         WeightNorm =
     &        1.d0/sqrt(1.d0-2.d0*SolWeight+2.d0*SolWeight**2)
      endif

      WeightNorm = WeightNorm * 2.d0 / sqrt(dble(NModes))

      N = 0
      do ikx = 0, ikmax(1)
         kx = PI_2*dble(ikx)/boxSize(1)
         do iky = 0, ikmax(2)
            ky = PI_2*dble(iky)/boxSize(2)
            do ikz = 0, ikmax(3)
               kz = PI_2*dble(ikz)/boxSize(3)

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

                  if (Ndims.eq.3) then
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
               endif

            enddo
         enddo
      enddo

C      if (cello_is_root .ne. 0) then
C         print *, "accDims = ", Ndims
C         print '(A,3(G10.4))', " Box = ", boxSize
C         print *, "NModes = ", NModes
C         print *, "Mach = ", Mach
C         print *, 'c_s = ', c_s
C         print *, 'TDecay = ', TDecay
C         print *, 'DtFreq = ', DtFreq
C         print *, 'OUVar = ', OUVar
C         print *, 'kc = ', kc
C         print *, 'cello_read_sol = ',cello_read_sol
C      endif

      if (cello_read_sol .ne. 0) then
         open (42, file='OUturb.bin', form='UNFORMATTED',
     &        access='STREAM', action='read')
         read (42) modes_input
         if (modes_input .ne. NModes) then
            stop "wrong number of modes in OUturb.bin"
         endif
         if (cello_is_root .ne. 0) print *, 'reading phases for ',
     &        NModes, ' modes'
         read (42) seed, totTime, OUPhase, Phase
         close (42)
         call initrandomx(seed)
      else
         call initrandomx(100)

         do j = 1, NModes
            do i = 1, 3
               OUPhase(1,i,j) = RandomPhase()*OUVar
               OUPhase(2,i,j) = RandomPhase()*OUVar
            enddo
         enddo

         call calcPhases()
      endif

      end subroutine


C     **********************************************************************
      subroutine saveTurbForcingOUState(myid, nameSuffix)
      implicit none
      include 'turbforcingOU.h'
      integer :: myid, seed, crandseed
      character(*) :: nameSuffix

      if (myid.eq.0) then
         seed = crandseed()
         open (42, file='OUturb.bin'//trim(nameSuffix),
     &        form='UNFORMATTED', access='STREAM')
         write (42) NModes, seed, totTime, OUPhase, Phase
         close (42)
      endif

      end subroutine
C     **********************************************************************
      subroutine cleanupTurbForcingOU(myid)
      implicit none
      include 'turbforcingOU.h'
      integer :: myid

      call saveTurbForcingOUState(myid, char(0))
      deallocate(Mode, OUPhase, Ampl, Phase)

      end subroutine
C     **********************************************************************
      real*8 function RandomPhase()
      implicit none
      real*8 r0, r1, randomnr
      real*8, parameter :: PI_2 = 3.14159265358979323846d0*2.d0

      r0 = randomnr()
      r1 = randomnr()

      RandomPhase = sqrt(-2.d0*log(r0))*cos(PI_2*r1)
      end function
C     **********************************************************************
      subroutine calcPhases()
      implicit none
      include 'turbforcingOU.h'
      integer :: i, j
      real*8 :: ka, kb, kk, kk1, diva, divb, curla, curlb

      do j = 1, NModes
         ka = 0d0
         kb = 0d0
         kk = 0d0

         do i = 1, Ndims
            kk = kk + Mode(i,j)**2
            ka = ka + Mode(i,j)*OUPhase(2,i,j)
            kb = kb + Mode(i,j)*OUPhase(1,i,j)
         enddo

         kk1 = 1.d0/kk

         do i = 1, Ndims
            diva = Mode(i,j)*ka*kk1
            divb = Mode(i,j)*kb*kk1
            curla = OUPhase(1,i,j) - divb
            curlb = OUPhase(2,i,j) - diva
            Phase(1,i,j) = SolWeight*curla + (1.d0 - SolWeight)*divb
            Phase(2,i,j) = SolWeight*curlb + (1.d0 - SolWeight)*diva
         enddo

      enddo

      end subroutine
C     **********************************************************************
      subroutine updatePhases()
      implicit none
      include 'turbforcingOU.h'
      real*8 :: damping, dfact, RandomPhase
      integer :: i, j, n,crandseed

      damping = exp(-DtFreq/TDecay)
      dfact = sqrt(1.d0 - damping**2)
      do j = 1, NModes
         do i = 1, 3
            do n = 1, 2
               OUPhase(n,i,j) = OUPhase(n,i,j)*damping +
     &              OUVar*dfact*RandomPhase()
            enddo
         enddo
      enddo
      end subroutine
C     **********************************************************************
      subroutine cello_turbou_state_size(n_real_state,n_integer_state)
      implicit none
      include 'turbforcingOU.h'
      integer :: n_real_state, n_integer_state
      n_real_state = 1+2*(NModes*2*3)
      n_integer_state = 2

      end subroutine
C     **********************************************************************
      subroutine cello_get_turbou_state(real_state,integer_state)
      implicit none
      integer num_modes
      real*8 :: real_state(2*3*num_modes)
      integer :: integer_state(*),crandseed
      include 'turbforcingOU.h'
      integer i,j,k,n

      integer_state(1) = NModes
      integer_state(2) = crandseed()

      k=1
      real_state(k) = totTime
      k=k+1
      do j = 1, NModes
         do i = 1, 3
            do n = 1, 2
               real_state(k) = OUPhase(n,i,j)
               real_state(k+1) = Phase(n,i,j)
               k=k+2
            enddo
         enddo
      enddo
      end subroutine
C     **********************************************************************
      subroutine cello_put_turbou_state(real_state,integer_state)
      implicit none
      integer num_modes
      real*8 :: real_state(2*3*num_modes)
      integer :: integer_state(*)
      include 'turbforcingOU.h'
      integer i,j,k,n

      NModes = integer_state(1)
      call setcrandseed (integer_state(2))

      k=1
      TotTime = real_state(k)
      k=k+1
      do j = 1, NModes
         do i = 1, 3
            do n = 1, 2
               OUPhase(n,i,j) = real_state(k)
               Phase  (n,i,j) = real_state(k+1)
               k=k+2
            enddo
         enddo
      enddo
      end subroutine
C     **********************************************************************
c     3D Turbulent forcing source term
c     Updates res:    res = res + source*dt
c     by DK
c     INPUT:
c     ni, nj, nk  - local domain dimensions
c     w           - flow field
c     grid        - grid coordinates
c     temperature - temperature field
c     wk          - work array 3*ni*nj*nk
c     dt          - time step
c     update_sol  - flag to update w
c     OUTPUT:
c     res         - residual (updated)
c     turbAcc     - injected turbulent accelerations and energy (updated)
c***********************************************************************
      subroutine turbForceOU
     &     (mx,my,mz,
     &     ni, nj, nk,
     &     field_density,
     &     grid,
     &     wk, time, dt,
     &     update_sol,
     &     cello_apply_cooling,
     &     cello_apply_forcing,
     &     cello_apply_injection_rate,
     &     cello_update_phases,
     &     cello_cooling_term,
     &     cello_gamma,
     &     cello_injection_rate,
     &     cello_olap,
     &     r_gv
     &     )
      implicit none
      include 'turbforcingOU.h'

      integer :: cello_apply_forcing
      integer :: cello_apply_injection_rate
      integer :: cello_cooling_term
      integer :: cello_apply_cooling
      real*8 :: cello_gamma
      real*8 :: cello_injection_rate
      integer :: cello_olap
      real*8 :: r_gv(4)
      integer :: cello_update_phases

      integer :: mx,my,mz
      integer :: ni, nj, nk, i, j, k, m, dir
      integer :: ol2, ol2k, cooling_term
      real*8 :: grid(3,mx,my,mz)
      real*8 :: field_density(mx,my,mz)
      real*8 :: wk(3,mx,my,mz)
      real*8 :: gamma
      real*8 :: deTurb
      real*8 :: rho, time, dt, dumx, acc(3), da(3), kr, de
      logical :: update_sol, apply_injection_rate
      real*8 :: timeUpdate

      if(cello_apply_forcing .eq. 0) return
      apply_injection_rate = (cello_apply_injection_rate .ne. 0)

      ol2 = cello_olap/2
      ol2k = ol2
      if (nk.eq.1) ol2k = 0

      timeUpdate = (floor(time / DtFreq)+1)*DtFreq

      if ( cello_update_phases.ne.0. ) then
         totTime = totTime + dt
         if ( (time       .le. timeUpdate) .and.
     &        (timeUpdate .lt. time + dt)) then
            call updatePhases()
            call calcPhases()
            totTime = totTime - DtFreq
         endif
      endif

      r_gv = 0d0
      do k=1,nk
         do j=1,nj
            do i=1,ni

               acc = 0d0
               do m=1, NModes
                  kr = 0d0
                  do dir=1,3
                     kr = kr + Mode(dir,m)*grid(dir,i,j,k)
                  enddo

                  acc(1:3) = acc(1:3) + Ampl(m)*( Phase(1,1:3,m)*cos(kr) -
     &                 Phase(2,1:3,m)*sin(kr) )
               enddo

               acc = acc*WeightNorm

               if ((k.ge.1+ol2k).and.(k.le.nk-ol2k).and.
     &              (j.ge.1+ol2).and.(j.le.nj-ol2).and.
     &              (i.ge.1+ol2).and.(i.le.ni-ol2)) then
                  rho = field_density(i,j,k)
                  r_gv(1:3) = r_gv(1:3) + acc*rho
                  r_gv(4) = r_gv(4) + rho
               endif

               wk(1:3,i,j,k) = acc
            enddo
         enddo
      enddo
      end subroutine
C======================================================================
c     jac         - jacobian
      subroutine turbForceShift
     &     (mx,my,mz,
     &     ni, nj, nk,
     &     field_density,
     &     field_momentum_x,
     &     field_momentum_y,
     &     field_momentum_z,
     &     field_jacobian, wk,
     &     update_sol,
     &     cello_apply_injection_rate,
     &     cello_olap,
     &     cello_injection_rate,
     &     r_gv,
     &     r_av
     &     )
      implicit none
      include 'turbforcingOU.h'

      integer :: ni, nj, nk, i, j, k
      integer :: cello_apply_injection_rate
      integer :: cello_olap
      integer :: mx,my,mz
      real*8 :: field_density(mx,my,mz)
      real*8 :: field_momentum_x(mx,my,mz)
      real*8 :: field_momentum_y(mx,my,mz)
      real*8 :: field_momentum_z(mx,my,mz)
      real*8 :: field_jacobian(mx,my,mz)
      real*8 :: wk(3,mx,my,mz)
      real*8 :: cello_injection_rate
      logical :: update_sol, apply_injection_rate
      real*8 :: r_gv(4), r_av(2)

      real*8 :: da(3), acc(3), Edot
      integer ol2, ol2k
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     force total impulse to zero using reduction r_gv
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      da = r_gv(1:3)/r_gv(4)
      do k=1,nk
         do j=1,nj
            do i=1,ni
               wk(1:3,i,j,k) = wk(1:3,i,j,k) - da
            enddo
         enddo
      enddo

      apply_injection_rate = (cello_apply_injection_rate .ne. 0)

      if (apply_injection_rate) then

C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     compute r_av reduction
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         Edot = cello_injection_rate * boxVolume

         ol2 = cello_olap/2
         ol2k = ol2
         if (nk.eq.1) ol2k = 0

         r_av = 0d0
         do k=1+ol2k,nk-ol2k
            do j=1+ol2,nj-ol2
               do i=1+ol2,ni-ol2
                  acc = wk(1:3,i,j,k)
                  r_av(1) = r_av(1) + 0.5d0 * field_jacobian(i,j,k) *
     &                 field_density(i,j,k) * dot_product(acc, acc)
                  if (Ndims.eq.2) then
                     r_av(2) = r_av(2) +
     &                    field_jacobian(i,j,k) *
     &                    (field_momentum_x(i,j,k)*wk(1,i,j,k) +
     &                    field_momentum_y(i,j,k)*wk(2,i,j,k))
                  else if (Ndims.eq.3) then
                     r_av(2) = r_av(2) +
     &                    field_jacobian(i,j,k) *
     &                    (field_momentum_x(i,j,k)*wk(1,i,j,k) +
     &                    field_momentum_y(i,j,k)*wk(2,i,j,k) +
     &                    field_momentum_z(i,j,k)*wk(3,i,j,k))
                  endif
               enddo
            enddo
         enddo
      endif
      end subroutine
C======================================================================
      subroutine turbForceUpdate
     &     (mx, my, mz,
     &     ni, nj, nk,
     &     field_density,
     &     field_momentum_x,
     &     field_momentum_y,
     &     field_momentum_z,
     &     field_energy,
     &     resid_density,
     &     resid_momentum_x,
     &     resid_momentum_y,
     &     resid_momentum_z,
     &     resid_energy,
     &     field_temperature, wk, dt,
     &     turbAcc, update_sol,
     &     cello_apply_injection_rate,
     &     cello_injection_rate,
     &     cello_cooling_term,
     &     cello_apply_cooling,
     &     cello_gamma,
     &     cello_hc_alpha,
     &     cello_hc_sigma,
     &     cello_totemp,
     &     r_av
     &     )
      implicit none
      include 'turbforcingOU.h'

      integer :: mx,my,mz
      integer :: ni, nj, nk, i, j, k
      real*8 :: field_density(mx,my,mz)
      real*8 :: field_momentum_x(mx,my,mz)
      real*8 :: field_momentum_y(mx,my,mz)
      real*8 :: field_momentum_z(mx,my,mz)
      real*8 :: field_energy(mx,my,mz)
      integer :: cello_apply_injection_rate
      integer :: cello_cooling_term
      integer :: cello_apply_cooling
      real*8 :: cello_gamma
      real*8 :: cello_hc_alpha
      real*8 :: cello_hc_sigma
      real*8 :: resid_density(mx,my,mz)
      real*8 :: resid_momentum_x(mx,my,mz)
      real*8 :: resid_momentum_y(mx,my,mz)
      real*8 :: resid_momentum_z(mx,my,mz)
      real*8 :: resid_energy(mx,my,mz)
      real*8 :: dt
      real*8 :: field_temperature(mx,my,mz), wk(3,mx,my,mz)
      real*8 :: turbAcc(4, mx,my,mz)
      real*8 :: cello_injection_rate
      real*8 :: cello_totemp
      logical :: update_sol, apply_injection_rate
      real*8 :: r_av(2)

      real*8 :: deTurb, rho, du(3), Acorr, acc(3), de, Edot, gamma, E0
      real*8 :: gm1, hc_alpha, hc_sigma, totemp
      integer nd1, cooling_term

      apply_injection_rate = (cello_apply_injection_rate .ne. 0)

      if (apply_injection_rate) then
         Edot = cello_injection_rate * boxVolume
         Acorr = solveQuadEquation(r_av(1)*dt, r_av(2), -Edot)
      else
         Acorr = 1.d0
      endif
      gamma    = cello_gamma

      E0 = 1.d0 / (gamma - 1.d0)
      gm1 = gamma-1
      hc_alpha = cello_hc_alpha
      hc_sigma = cello_hc_sigma

C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     use computed Acorr
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      cooling_term = cello_cooling_term

      if (cooling_term.eq.-1) then
! backwards compatibility for "apply_cooling" option
         if (cello_apply_cooling .eq. 0) then
            cooling_term = 0
         else
            cooling_term = 1
         endif
      endif

      nd1 = Ndims + 1
      do k=1,nk
         do j=1,nj
            do i=1,ni
               acc = Acorr * wk(1:3,i,j,k)
               du = acc*dt

               rho = field_density(i,j,k)
               if (Ndims.eq.2) then
                  deTurb = dot_product(du(1:Ndims),du(1:Ndims))*
     &                 0.5d0*rho +
     &                 (field_momentum_x(i,j,k) * du(1) +
     &                 field_momentum_y(i,j,k) * du(2))
               else if (Ndims.eq.3) then
                  deTurb = dot_product(du(1:Ndims),du(1:Ndims))*
     &                 0.5d0*rho +
     &                 (field_momentum_x(i,j,k) * du(1) +
     &                 field_momentum_y(i,j,k) * du(2) +
     &                 field_momentum_z(i,j,k) * du(3))
               end if
               de = addCooling (deTurb, rho,
     &              field_energy(i,j,k), field_temperature(i,j,k))
               resid_momentum_x(i,j,k) = resid_momentum_x(i,j,k) +
     &              du(1)*rho
               resid_momentum_y(i,j,k) = resid_momentum_y(i,j,k) +
     &              du(2)*rho
               if (Ndims.eq.3) then
                  resid_momentum_z(i,j,k) = resid_momentum_z(i,j,k) +
     &                 du(3)*rho
               end if
               resid_energy(i,j,k) = resid_energy(i,j,k) + de

               turbAcc(1:Ndims,i,j,k) =
     &              turbAcc(1:Ndims,i,j,k) + acc(1:Ndims)
               turbAcc(4,i,j,k) = turbAcc(4,i,j,k) + de
            enddo
         enddo
      enddo

      if(update_sol) then
         field_density = field_density + resid_density
         field_momentum_x = field_momentum_x + resid_momentum_x
         field_momentum_y = field_momentum_y + resid_momentum_y
         if (Ndims.eq.3) then
            field_momentum_z = field_momentum_z + resid_momentum_z
         endif
         field_energy = field_energy + resid_energy
      endif

      contains
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real*8 function solveQuadEquation(a, b, c)
      real*8 :: a, b, c, sd
      sd = sqrt(b * b - 4 * a * c)
      solveQuadEquation = (sd - b) / (2.d0 * a)
      end function
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      real*8 function addCooling(deTurb, rho, E, temp)
      external cooling_F2, cooling_F6, cooling_F40
      integer, parameter :: iTask = 1, iOpt = 0, NRW = 32, NIW = 21
      real*8, parameter :: ngas = 1.d0
      integer :: iState, IW(NIW)
      real*8 :: deTurb, rho, e, temp, T, T2, X
      real*8 :: tm0, tm1, Y(5), RW(NRW), dummy

      select case (cooling_term)
      case (0)
c     No cooling
         addCooling = deTurb
      case (1)
c     Empiric cooling formula
         T = temp*2.d3          ! Convert to SI

         addCooling = deTurb - dt * ngas * rho*
     &        (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
     &        1.d-5             ! Convert from SI
      case (2)
c     Stiff version of case 1 => use ODE solver
         X = gm1*totemp/rho
         Y(1) = E
         Y(2) = deTurb/dt
         Y(3) = temp - X*E
         Y(4) = X
         Y(5) = rho
         tm0 = 0d0
         tm1 = dt
         iState = 1
 002     CALL DLSODE (cooling_F2, 1, Y, tm0, tm1, 1, 1.D-8, 1.D-10,
     &        iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

         IF(iState .EQ. -1) THEN
            iState = 2
            GOTO 002
         ENDIF
         addCooling = Y(1) - E
!     T = temp*2.d3 ! Convert to SI
!     X = deTurb - dt * ngas * rho*
!     &   (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
!     &   1.d-5      ! Convert from SI
!     print *, IW(11), addCooling, X
      case (3)
c     Linear cooling (as in W.Schmidt & P.Grete 2019)
         addCooling = deTurb - dt * hc_alpha  * (E - E0)
      case (4)
c     Stiff version of case 3 => use analytical solution
         X = deTurb/(dt * hc_alpha)
         addCooling = (E0 + X - E)*(1.0 - exp(-hc_alpha * dt))
      case (40)
c     Stiff version of case 3 => use ODE solver (for testing)
         Y(1) = E
         Y(2) = deTurb/dt
         Y(3) = E0
         Y(4) = hc_alpha
         tm0 = 0d0
         tm1 = dt
         iState = 1
 040     CALL DLSODE (cooling_F40, 1, Y, tm0, tm1, 1, 1.D-8, 1.D-10,
     &        iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

         IF(iState .EQ. -1) THEN
            iState = 2
            GOTO 040
         ENDIF
         addCooling = Y(1) - E
!     print *, IW(11), addCooling,
!     &  (E0 +  deTurb/(dt * hc_alpha) - E)*(1.0 - exp(-hc_alpha * dt))
      case (5)
c     Stefan's law (as in D.Porter 2002)
         T2 = temp * temp
         addCooling = deTurb - dt * hc_sigma * T2 * T2
      case (6)
c     Stiff version of case 5 => use ODE solver
         X = gm1*totemp/rho
         Y(1) = E
         Y(2) = deTurb/dt
         Y(3) = temp - X*E
         Y(4) = X
         Y(5) = hc_sigma
         tm0 = 0d0
         tm1 = dt
         iState = 1
 006     CALL DLSODE (cooling_F6, 1, Y, tm0, tm1, 1, 1.D-8, 1.D-10,
     &        iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

         IF(iState .EQ. -1) THEN
            iState = 2
            GOTO 006
         ENDIF
         addCooling = Y(1) - E
!     print *, IW(11), addCooling
!     &  deTurb - dt * hc_sigma * (temp ** 4)
      case default
         stop "unsupported cooling term"
      end select

      end function
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end subroutine
C     **********************************************************************
      subroutine cooling_F2(NEQ, tm, Y, Ydot)
      implicit none
      integer :: NEQ
      real*8, parameter :: ngas = 1.d0
      real*8 :: tm, Y(*), Ydot(NEQ)
      real*8 :: deTurbDot, T, temp, rho

      deTurbDot = Y(2)
      temp = Y(3) + Y(4) * Y(1) ! T = (T0 - k*E0) + k*E
      rho = Y(5)

      T = temp*2.d3             ! Convert to SI

      Ydot(1) = deTurbDot - ngas * rho*
     &     (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
     &     1.d-5                ! Convert from SI

      end subroutine
C     **********************************************************************
      subroutine cooling_F6(NEQ, tm, Y, Ydot)
      implicit none
      integer :: NEQ
      real*8, parameter :: ngas = 1.d0
      real*8 :: tm, Y(*), Ydot(NEQ)
      real*8 :: deTurbDot, T2, temp, hc_sigma

      deTurbDot = Y(2)
      temp = Y(3) + Y(4) * Y(1) ! T = (T0 - k*E0) + k*E
      hc_sigma = Y(5)

      T2 = temp * temp
      Ydot(1) = deTurbDot - hc_sigma * T2 * T2

      end subroutine
C     **********************************************************************
      subroutine cooling_F40(NEQ, tm, Y, Ydot)
      implicit none
      integer :: NEQ
      real*8, parameter :: ngas = 1.d0
      real*8 :: tm, Y(*), Ydot(NEQ)
      real*8 :: deTurbDot, E0, hc_alpha

      deTurbDot = Y(2)
      E0 = Y(3)
      hc_alpha = Y(4)

      Ydot(1) = deTurbDot - hc_alpha  * (Y(1) - E0)
      end subroutine
