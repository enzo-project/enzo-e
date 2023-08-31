      common /TURB3D/ OUPhase, Ampl, Phase, Mode, totTime, DtFreq
      common /TURB3D/ OUVar, TDecay, SolWeight, WeightNorm, NModes
      real*8 :: totTime, TDecay, DtFreq, SolWeight, WeightNorm
      real*8 :: OUVar ! Ornstein-Uhlenbeck var
      integer :: NModes
      real*8, pointer :: OUPhase(:,:,:) ! Ornstein-Uhlenbeck phases
      real*8, pointer :: Ampl(:)
      real*8, pointer :: Phase(:,:,:)
      real*8, pointer :: Mode(:,:)