c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

#ifdef INCLUDE_GRAVITY
      Subroutine POTOKy_ig(Nx,gamma,Dx,Dt,Bp,Qlx,Qrx,Qpm,Qdx,Qdz,Sy,
     &     Qrp,Fm)
#else
      Subroutine POTOKy_ig(Nx,gamma,Dx,Dt,Bp,Qlx,Qrx,Qpm,Qdx,Qdz,
     &     Qrp,Fm)
#endif

      Implicit NONE
      Integer nx,i,j,k,m,l

      ENZO_REAL Qlx(Nx,8),Qrx(Nx,8),Qpm(Nx,8),Qdx(Nx,8),Qdz(Nx,8)
      ENZO_REAL Qrp(Nx,8),Fm(Nx,8),Bp(Nx,3)

      ENZO_REAL qp(8),px(8),qc(8),qvl(8),qvr(8),bs(3)
      ENZO_REAL vm(8),ql(8,8),qr(8,8),uv(8),vl(8)
      ENZO_REAL sl(8),qul(8),qur(8),tl(8)
      ENZO_REAL tr(8),fr(8),f(8),ur(8),qv(8)
      ENZO_REAL udx(8),udz(8),ax(8,8),az(8,8),av(8)
      ENZO_REAL tvl(8),tvr(8),tp(8),tc(8),tdx(8),tdz(8)
      ENZO_REAL DT,DX,DTX,DXI,TL6,gamma
      ENZO_REAL bl,br,fl1,fr1,ul_v,ur_v,cfl,qsr(8)
      ENZO_REAL Syl(8),Syr(8)
#ifdef INCLUDE_GRAVITY
      ENZO_REAL Sy(Nx,8)
#endif

      Do i=1,Nx
         Do m=1,8
            qrp(i,m)=0.
            fm(i,m)=0.
         Enddo
      Enddo

      Do i=3,Nx-3

         Do m=1,8

            Qvl(m)=Qlx(i,m)
            Qvr(m)=Qrx(i,m)
            Qp(m)=Qpm(i,m)
            Udx(m)=Qdx(i,m)
            Udz(m)=Qdz(i,m)
#ifdef INCLUDE_GRAVITY
            Syl(m)=Sy(i,m)
#else
            Syl(m)=0.0d0
#endif
            Tvl(m)=Qlx(i+1,m)
            Tvr(m)=Qrx(i+1,m)
            Tp(m)=Qpm(i+1,m)
            Tdx(m)=Qdx(i+1,m)
            Tdz(m)=Qdz(i+1,m)
#ifdef INCLUDE_GRAVITY
            Syr(m)=Sy(i+1,m)
#else
            Syr(m)=0.0d0
#endif
         Enddo

         Do m=1,3
            Bs(m)=Bp(i,m)
         Enddo


         DTX=DT/DX

         CALL VECTEGy_ig(gamma,Qp,Bs,Qc)
         CALL VECTEGy_ig(gamma,Tp,Bs,Tc)

         IF(QC(7).GT.0.)THEN

            DXI=QC(7)*DTX

            CALL FPML_ig(DXI,QVL,QVR,QP,VM)

            CALL VECTLRY_ig(gamma,QP,BS,QL,QR,UV)
            CALL MATR_AX_ig(gamma,QP,BS,AX)
            CALL MATR_AZ_ig(gamma,QP,BS,AZ)

            DO L=1,8
               VL(L)=0.
            ENDDO

            DO M=1,8
               AV(M)=0.
               DO L=1,8
                  AV(M)=AV(M)+AX(M,L)*UDX(L)+AZ(M,L)*UDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.-SYL(M)*DT/2.
            ENDDO

            DO M=1,8

               IF(QC(M).GT.0.)THEN

                  DXI=QC(M)*DTX

                  CALL FPML_ig(DXI,QVL,QVR,QP,QV)

                  SL(M)=0.
                  DO L=1,8
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO

                  DO L=1,8
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF

            ENDDO

            DO M=1,8
               QUL(M)=VM(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,8
               QUL(M)=QVR(M)
            ENDDO

         ENDIF


         IF(TC(1).LT.0.)THEN

            DXI=-TC(1)*DTX

            CALL FPMR_ig(DXI,TVL,TVR,TP,VM)

            CALL VECTLRY_ig(gamma,TP,BS,QL,QR,UV)
            CALL MATR_AX_ig(gamma,TP,BS,AX)
            CALL MATR_AZ_ig(gamma,TP,BS,AZ)

            DO L=1,8
               VL(L)=0.
            ENDDO

            DO M=1,8
               AV(M)=0.
               DO L=1,8
                  AV(M)=AV(M)+AX(M,L)*TDX(L)+AZ(M,L)*TDZ(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.-SYR(M)*DT/2.
            ENDDO

            DO M=1,8

               IF(TC(M).LT.0.)THEN

                  DXI=-TC(M)*DTX

                  CALL FPMR_ig(DXI,TVL,TVR,TP,QV)

                  SL(M)=0.
                  DO L=1,8
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO

                  DO L=1,8
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF

            ENDDO

            DO M=1,8
               QUR(M)=VM(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,8
               QUR(M)=TVL(M)
            ENDDO

         ENDIF

         DO M=1,8
            TL(M)=QUL(M)
            TR(M)=QUR(M)
         ENDDO

         TL(6)=QVR(6)
         TR(6)=TVL(6)

         CALL CONS_ig(gamma,TL,PX)
         CALL CONS_ig(gamma,TR,VL)
         CALL HLLDy_ig(gamma,PX,VL,BS,FR,UR) !HLLD method

         DO M=1,8
            FM(I,M)=FR(M)
         ENDDO





         IF(QC(7).GT.0.)THEN

            DXI=1.-QC(7)*DTX

            CALL FQML_ig(DXI,QVL,QVR,QP,VM)

            CALL VECTLRY_ig(gamma,QP,BS,QL,QR,UV)
            CALL MATR_AX_ig(gamma,QP,BS,AX)
            CALL MATR_AZ_ig(gamma,QP,BS,AZ)

            DO L=1,8
               VL(L)=0.
            ENDDO

            DO M=1,8
               AV(M)=0.
               DO L=1,8
                  AV(M)=AV(M)+AX(M,L)*UDX(L)+AZ(M,L)*UDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.-SYL(M)*DT/2.
            ENDDO

            DO M=1,8

               IF(QC(M).GT.0.)THEN

                  DXI=1.-QC(M)*DTX

                  CALL FQML_ig(DXI,QVL,QVR,QP,QV)

                  SL(M)=0.
                  DO L=1,8
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO

                  DO L=1,8
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF

            ENDDO

            DO M=1,8
               QUL(M)=VM(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,8
               QUL(M)=QVR(M)
            ENDDO

         ENDIF

         IF(TC(1).LT.0.)THEN

            DXI=-TC(1)*DTX

            CALL FQML_ig(DXI,TVL,TVR,TP,VM)

            CALL VECTLRY_ig(gamma,TP,BS,QL,QR,UV)
            CALL MATR_AX_ig(gamma,TP,BS,AX)
            CALL MATR_AZ_ig(gamma,TP,BS,AZ)

            DO L=1,8
               VL(L)=0.
            ENDDO

            DO M=1,8
               AV(M)=0.
               DO L=1,8
                  AV(M)=AV(M)+AX(M,L)*TDX(L)+AZ(M,L)*TDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.
            ENDDO

            DO M=1,8

               IF(TC(M).LT.0.)THEN

                  DXI=-TC(M)*DTX

                  CALL FQML_ig(DXI,TVL,TVR,TP,QV)

                  SL(M)=0.
                  DO L=1,8
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO

                  DO L=1,8
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF

            ENDDO

            DO M=1,8
               QUR(M)=VM(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,8
               QUR(M)=TVL(M)
            ENDDO

         ENDIF

         DO M=1,8
            TL(M)=QUL(M)
            TR(M)=QUR(M)
         ENDDO

         TL(6)=QVR(6)
         TR(6)=TVL(6)

         CALL CONS_ig(gamma,TL,PX)
         CALL CONS_ig(gamma,TR,VL)
         CALL HLLDy_ig(gamma,PX,VL,BS,FR,UR) !HLLD method


         CALL PRIM_ig(gamma,UR,PX)

         DO M=1,8
            QRP(I,M)=PX(M)
         ENDDO

      ENDDO

      return
      end
