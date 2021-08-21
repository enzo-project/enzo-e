C **********************************************************************
C *** Forcing terms for OU turbulence (2D or 3D)
C *** Adapted from AREPO code developed by Volker Springel
C **********************************************************************
      subroutine initTurbForcingOU(myid)
      implicit none
      include 'turbforcingOU.h'
      real*8 :: Mach, gamma, P0, rho0, c_s, Edot
      real*8 :: kmax, kmin, boxSize(3), kx, ky, kz, kk, kc, amp
      integer :: argi(1), ikmax(3), ikx, iky, ikz, N, i, j
      real*8 :: RandomPhase, argx(3)
      real*8, parameter :: PI_2 = 3.14159265358979323846d0*2.d0 
      integer :: myid, seed, modesPerK
      logical :: apply_injection_rate

      call GETDATA('apply_injection_rate', argx, argi)
      apply_injection_rate = (argi(1) .ne. 0)

      call GETDATA('accDims', argx, argi)
      Ndims = argi(1)

      call GETDATA("gamma", gamma, argi)
      call GETDATA("Lx",boxSize(1),argi)
      call GETDATA("Ly",boxSize(2),argi)
      call GETDATA("Lz",boxSize(3),argi)
      if (Ndims.eq.3) then
        modesPerK = 4
        boxVolume = boxSize(1) * boxSize(2) * boxSize(3)
      elseif(Ndims.eq.2) then
        modesPerK = 2
        boxVolume = boxSize(1) * boxSize(2)
      else
        stop "only dim=2 or dim=3 is supported"
      endif        

      call GETDATA("InjectionRate", Edot, argi)   ! default 0.006
      call GETDATA("SolWeight", SolWeight, argi)  ! default 1.0
      call GETDATA("Mach", Mach, argi)            ! default 0.2
      call GETDATA("kfi", kmin, argi)             ! default 6.27
      call GETDATA("kfa", kmax, argi)             ! default 12.57
            
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
     *    1.d0/sqrt(1.d0-2.d0*SolWeight+3.d0*SolWeight**2)
      else
        WeightNorm = 
     *    1.d0/sqrt(1.d0-2.d0*SolWeight+2.d0*SolWeight**2)
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

      if(myid.eq.0) then
        print *, "accDims = ", Ndims
        print '(A,3(G10.4))', " Box = ", boxSize
        print *, "NModes = ", NModes
        print *, "Mach = ", Mach
        print *, 'c_s = ', c_s
        print *, 'TDecay = ', TDecay
        print *, 'DtFreq = ', DtFreq
        print *, 'OUVar = ', OUVar
        print *, 'kc = ', kc
      endif

      call GETDATA('readinsol', argx, argi)
      if (argi(1) .ne. 0) then
        open (42, file='OUturb.bin', form='UNFORMATTED', 
     &      access='STREAM', action='read')
        read (42) argi(1)
        if (argi(1) .ne. NModes) then
          stop "wrong number of modes in OUturb.bin"
        endif
        if (myid.eq.0) print *, 'reading phases for ', NModes, ' modes'
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
C **********************************************************************      
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
C **********************************************************************
      subroutine cleanupTurbForcingOU(myid)
      implicit none
      include 'turbforcingOU.h'
      integer :: myid

      call saveTurbForcingOUState(myid, char(0))
      deallocate(Mode, OUPhase, Ampl, Phase)

      end subroutine
C **********************************************************************
      real*8 function RandomPhase()
      implicit none
      real*8 r0, r1, randomnr
      real*8, parameter :: PI_2 = 3.14159265358979323846d0*2.d0 
      
      r0 = randomnr()
      r1 = randomnr()
      
      RandomPhase = sqrt(-2.d0*log(r0))*cos(PI_2*r1)
      end function 
C **********************************************************************
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
C **********************************************************************
      subroutine updatePhases()   
      implicit none
      include 'turbforcingOU.h'
      real*8 :: damping, dfact, RandomPhase
      integer :: i, j, n
      
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
c     3D Turbulent forcing source term
c     Updates res:    res = res + source*dt
c     by DK
c INPUT:
c       nc          - number of components
c       ni, nj, nk  - local domain dimensions
c       w           - flow field
c       grid        - grid coordinates
c       jac         - jacobian
c       temperature - temperature field
c       wk          - work array 3*ni*nj*nk
c       dt          - time step
c       update_sol  - flag to update w
c OUTPUT:
c       res         - residual (updated)
c       turbAcc     - injected turbulent accelerations and energy (updated)
c***********************************************************************
      subroutine turbForceOU(nc, ni, nj, nk,
     *                       w, grid, jac, temperature, wk, dt, res,
     *                       turbAcc, update_sol)
      implicit none
      include 'turbforcingOU.h'
#ifdef USING_MPI
      include 'mpif.h'
      integer :: err
      real*8 :: gvmpi(4), avmpi(2)
#endif            
      integer :: nc, ni, nj, nk, argi(1), i, j, k, m, dir
      integer :: nd1, ol2, ol2k, cooling_term
      real*8 :: grid(3,ni,nj,nk), w(nc,ni,nj,nk), res(nc,ni,nj,nk)
      real*8 :: jac(ni,nj,nk), temperature(ni,nj,nk), wk(3,ni,nj,nk) 
      real*8 :: turbAcc(4, ni,nj,nk), hc_alpha, hc_sigma, gamma
      real*8 :: totemp, gm1, deTurb, E0, Edot, Acorr, av(2)
      real*8 :: rho, dt, dumx, acc(3), da(3), du(3), gv(4), kr, de
      logical :: update_sol, apply_injection_rate
            
      call GETDATA("apply_forcing", dumx, argi)
      if(argi(1).eq.0) return 

      call GETDATA('apply_injection_rate', dumx, argi)
      apply_injection_rate = (argi(1) .ne. 0)
      
      call GETDATA("cooling_term", dumx, argi)
      cooling_term = argi(1)
      
      if (cooling_term.eq.-1) then 
        ! backwards compatibility for "apply_cooling" option
        call GETDATA("apply_cooling", dumx, argi)
        if (argi(1).eq.0) then
          cooling_term = 0
        else
          cooling_term = 1
        endif
      endif

      call GETDATA("hc_alpha", hc_alpha, argi)
      call GETDATA("hc_sigma", hc_sigma, argi)
      call GETDATA("totemp", totemp, argi )
      call GETDATA("gamma", gamma, argi)
      gm1 = gamma-1

      E0 = 1.d0 / (gamma - 1.d0)

      call GETDATA('olap', dumx, argi)
      ol2 = argi(1)/2
      ol2k = ol2
      if (nk.eq.1) ol2k = 0

      if(dt*2.d0 .gt. DtFreq) then ! require at least two dt in DtFreq
        print *, "dt = ", dt, " > DtFreq = ", DtFreq
        stop "dt too big"
      endif
      
      totTime = totTime + dt
      if(totTime .gt. DtFreq) then
        call updatePhases()
        call calcPhases()
        totTime = totTime - DtFreq
      endif

      nd1 = Ndims + 1
      gv = 0d0
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
     &                                    Phase(2,1:3,m)*sin(kr) )
        enddo

        acc = acc*WeightNorm

        if ((k.ge.1+ol2k).and.(k.le.nk-ol2k).and.
     *      (j.ge.1+ol2).and.(j.le.nj-ol2).and.
     *      (i.ge.1+ol2).and.(i.le.ni-ol2)) then
          rho = w(1,i,j,k) 
          gv(1:3) = gv(1:3) + acc*rho
          gv(4) = gv(4) + rho
        endif

        wk(1:3,i,j,k) = acc
      enddo
      enddo
      enddo
      
#ifdef USING_MPI
      call MPI_AllReduce(gv, gvmpi, 4, MPI_DOUBLE_PRECISION, MPI_SUM,
     *                  MPI_COMM_WORLD, err )
      gv = gvmpi
#endif 

      ! force total impulse to zero
      da = gv(1:3)/gv(4)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        wk(1:3,i,j,k) = wk(1:3,i,j,k) - da
      enddo
      enddo
      enddo

      if (apply_injection_rate) then
        call GETDATA("InjectionRate", dumx, argi)
        Edot = dumx * boxVolume

        av = 0d0
        do k=1+ol2k,nk-ol2k
        do j=1+ol2,nj-ol2
        do i=1+ol2,ni-ol2
          acc = wk(1:3,i,j,k)
          av(1) = av(1) + 0.5d0 * 
     &            jac(i,j,k) * w(1,i,j,k) * dot_product(acc, acc)
          av(2) = av(2) +
     &            jac(i,j,k) * dot_product(w(2:4,i,j,k), acc)
        enddo
        enddo
        enddo

#ifdef USING_MPI
        call MPI_AllReduce(av, avmpi, 2, MPI_DOUBLE_PRECISION, MPI_SUM,
     *                     MPI_COMM_WORLD, err )
        av = avmpi
#endif
        Acorr = solveQuadEquation(av(1)*dt, av(2), -Edot)
      else
        Acorr = 1.d0
      endif

      do k=1,nk
      do j=1,nj
      do i=1,ni
        acc = Acorr * wk(1:3,i,j,k)
        du = acc*dt

        rho = w(1,i,j,k)
        deTurb = dot_product(du(1:Ndims),du(1:Ndims))*0.5d0*rho + 
     *           dot_product(w(2:nd1,i,j,k), du(1:Ndims))
        de = addCooling(deTurb, rho, w(5,i,j,k), temperature(i,j,k))
        res(2:nd1,i,j,k) = res(2:nd1,i,j,k) + du(1:Ndims)*rho
        res(5,i,j,k) = res(5,i,j,k) + de

        turbAcc(1:Ndims,i,j,k) = turbAcc(1:Ndims,i,j,k) + acc(1:Ndims)
        turbAcc(4,i,j,k) = turbAcc(4,i,j,k) + de
      enddo
      enddo
      enddo
      
      if(update_sol) then
        w = w + res
      endif

      contains
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      real*8 function solveQuadEquation(a, b, c)
        real*8 :: a, b, c, sd
        sd = sqrt(b * b - 4 * a * c)
        solveQuadEquation = (sd - b) / (2.d0 * a)
      end function
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      real*8 function addCooling(deTurb, rho, E, temp)
      external cooling_F2, cooling_F6, cooling_F40
      integer, parameter :: iTask = 1, iOpt = 0, NRW = 32, NIW = 21
      real*8, parameter :: ngas = 1.d0
      integer :: iState, IW(NIW)
      real*8 :: deTurb, rho, e, temp, T, T2, X
      real*8 :: tm0, tm1, Y(5), RW(NRW), dummy

      select case (cooling_term)
      case (0)
c No cooling        
        addCooling = deTurb
      case (1)
c Empiric cooling formula        
        T = temp*2.d3 ! Convert to SI

        addCooling = deTurb - dt * ngas * rho*
     *   (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
     *   1.d-5      ! Convert from SI
      case (2)
c Stiff version of case 1 => use ODE solver
        X = gm1*totemp/rho
        Y(1) = E
        Y(2) = deTurb/dt
        Y(3) = temp - X*E
        Y(4) = X
        Y(5) = rho
        tm0 = 0d0
        tm1 = dt
        iState = 1
002       CALL DLSODE (cooling_F2, 1, Y, tm0, tm1, 1, 1.D-8, 1.D-10,
     *         iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

        IF(iState .EQ. -1) THEN
          iState = 2
          GOTO 002
        ENDIF
        addCooling = Y(1) - E
!        T = temp*2.d3 ! Convert to SI
!        X = deTurb - dt * ngas * rho*
!     *   (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
!     *   1.d-5      ! Convert from SI
!        print *, IW(11), addCooling, X
      case (3)
c Linear cooling (as in W.Schmidt & P.Grete 2019)
        addCooling = deTurb - dt * hc_alpha  * (E - E0)
      case (4)
c Stiff version of case 3 => use analytical solution
        X = deTurb/(dt * hc_alpha)
        addCooling = (E0 + X - E)*(1.0 - exp(-hc_alpha * dt))
      case (40)
c Stiff version of case 3 => use ODE solver (for testing)
        Y(1) = E
        Y(2) = deTurb/dt
        Y(3) = E0
        Y(4) = hc_alpha
        tm0 = 0d0
        tm1 = dt
        iState = 1
040     CALL DLSODE (cooling_F40, 1, Y, tm0, tm1, 1, 1.D-8, 1.D-10,
     *       iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

        IF(iState .EQ. -1) THEN
          iState = 2
          GOTO 040
        ENDIF
        addCooling = Y(1) - E 
!        print *, IW(11), addCooling,
!     *  (E0 +  deTurb/(dt * hc_alpha) - E)*(1.0 - exp(-hc_alpha * dt))
      case (5)
c Stefan's law (as in D.Porter 2002)
        T2 = temp * temp
        addCooling = deTurb - dt * hc_sigma * T2 * T2
      case (6)
c Stiff version of case 5 => use ODE solver
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
     *       iTask, iState, iOpt, RW, NRW, IW, NIW, dummy, 22) ! 21 TO USE JAC; 22 - NO JAC

        IF(iState .EQ. -1) THEN
          iState = 2
          GOTO 006
        ENDIF
        addCooling = Y(1) - E
!        print *, IW(11), addCooling
!     *  deTurb - dt * hc_sigma * (temp ** 4)
      case default
        stop "unsupported cooling term"
      end select

      end function
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      end subroutine
C **********************************************************************
      subroutine cooling_F2(NEQ, tm, Y, Ydot)
        implicit none
        integer :: NEQ
        real*8, parameter :: ngas = 1.d0
        real*8 :: tm, Y(*), Ydot(NEQ)
        real*8 :: deTurbDot, T, temp, rho

        deTurbDot = Y(2)
        temp = Y(3) + Y(4) * Y(1) ! T = (T0 - k*E0) + k*E
        rho = Y(5)

        T = temp*2.d3 ! Convert to SI

        Ydot(1) = deTurbDot - ngas * rho*
     *   (1.d7*exp(-1.184d5/(T+1.d3)) + 1.4d-2*sqrt(T)*exp(-92/T))*
     *   1.d-5      ! Convert from SI

      end subroutine
C **********************************************************************
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
C **********************************************************************
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
