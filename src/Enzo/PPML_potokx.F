c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      Subroutine POTOKx(Nx,Dx,Dt,Qlx,Qrx,Qpm,Qdy,Qdz,Qrp,Fm) 
      Implicit NONE

      Integer nx,i,m,l

      ENZO_REAL Qlx(Nx,7),Qrx(Nx,7),Qpm(Nx,7),Qdy(Nx,7),Qdz(Nx,7)
      ENZO_REAL Qrp(Nx,7),Fm(Nx,7) 
      
      ENZO_REAL qp(7),px(7),qc(7),qvl(7),qvr(7)
      ENZO_REAL vm(7),ql(7,7),qr(7,7),uv(7),vl(7)
      ENZO_REAL sl(7),qul(7),qur(7),tl(7)
      ENZO_REAL tr(7),fr(7),ur(7),qv(7)
      ENZO_REAL udy(7),udz(7),ay(7,7),az(7,7),av(7)
      ENZO_REAL tvl(7),tvr(7),tp(7),tc(7),tdy(7),tdz(7)
      ENZO_REAL DT,DX,DTX,DXI,TL5

      Do i=1,Nx  
         Do m=1,7
            qrp(i,m)=0.
            fm(i,m)=0.
         Enddo                  
      Enddo                  

      Do i=3,Nx-3
         
         Do m=1,7
            
            Qvl(m)=Qlx(i,m)
            Qvr(m)=Qrx(i,m)
            Qp(m)=Qpm(i,m)
            Udy(m)=Qdy(i,m)
            Udz(m)=Qdz(i,m)
            
            Tvl(m)=Qlx(i+1,m)
            Tvr(m)=Qrx(i+1,m)
            Tp(m)=Qpm(i+1,m)
            Tdy(m)=Qdy(i+1,m)
            Tdz(m)=Qdz(i+1,m)

         Enddo

         DTX=DT/DX
         
         CALL VECTEGx(QP,QC)
         CALL VECTEGx(TP,TC)

         IF(QC(7).GT.0.)THEN                  
            
            DXI=QC(7)*DTX

            CALL FPML(DXI,QVL,QVR,QP,VM) 

            CALL VECTLRX(QP,QL,QR,UV)
            CALL MATR_AY(QP,AY)
            CALL MATR_AZ(QP,AZ)
            
            DO L=1,7
               VL(L)=0.
            ENDDO

            DO M=1,7
               AV(M)=0.
               DO L=1,7
                  AV(M)=AV(M)+AY(M,L)*UDy(L)+AZ(M,L)*UDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.
            ENDDO

            DO M=1,7  
               
               IF(QC(M).GT.0.)THEN                  

                  DXI=QC(M)*DTX

                  CALL FPML(DXI,QVL,QVR,QP,QV) 

                  SL(M)=0.   
                  DO L=1,7
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO

                  DO L=1,7
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF
               
            ENDDO
            
            DO M=1,7 
               QUL(M)=VM(M)-VL(M)
            ENDDO
            
         ELSE

            DO M=1,7 
               QUL(M)=QVR(M)
            ENDDO

         ENDIF 

         IF(TC(1).LT.0.)THEN                  
            
            DXI=-TC(1)*DTX
            
            CALL FPMR(DXI,TVL,TVR,TP,VM) 

            CALL VECTLRX(TP,QL,QR,UV)
            CALL MATR_AY(TP,AY)
            CALL MATR_AZ(TP,AZ)

            DO L=1,7
               VL(L)=0.
            ENDDO

            DO M=1,7
               AV(M)=0.
               DO L=1,7
                  AV(M)=AV(M)+AY(M,L)*TDy(L)+AZ(M,L)*TDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.
            ENDDO

            DO M=1,7  
               
               IF(TC(M).LT.0.)THEN                  

                  DXI=-TC(M)*DTX

                  CALL FPMR(DXI,TVL,TVR,TP,QV) 

                  SL(M)=0.   
                  DO L=1,7
                     SL(M)=SL(M)+QL(M,L)*(VM(L)-QV(L)+AV(L))
                  ENDDO
                  
                  DO L=1,7
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO
                  
               ENDIF

            ENDDO
            
            DO M=1,7
               QUR(M)=VM(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,7 
               QUR(M)=TVL(M)
            ENDDO

         ENDIF 

         do M=1,7
            TL(M)=QUL(M)
            TR(M)=QUR(M)
         enddo

         TL5=0.5d0*(TL(5)+TR(5))
         
         TL(5)=TL5
         TR(5)=TL5

         CALL CONS(TL,PX)
         CALL CONS(TR,VL)
         CALL HLLDx(PX,VL,FR,UR) !HLLD method

         Do M=1,7  
            Fm(I,M)=FR(M)      
         Enddo 







         IF(QC(7).GT.0.)THEN                  

            DXI=1.-QC(7)*DTX

            CALL FQML(DXI,QVL,QVR,QP,VM) 

            CALL VECTLRX(QP,QL,QR,UV)
            CALL MATR_AY(QP,AY)
            CALL MATR_AZ(QP,AZ)

            DO L=1,7
               VL(L)=0.
            ENDDO
            
            DO M=1,7
               AV(M)=0.
               DO L=1,7
                  AV(M)=AV(M)+AY(M,L)*UDy(L)+AZ(M,L)*UDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.
            ENDDO

            DO M=1,7  
               
               IF(QC(M).GT.0.)THEN                  

                  DXI=1.-QC(M)*DTX

                  CALL FQML(DXI,QVL,QVR,QP,QV) 

                  SL(M)=0.   
                  DO L=1,7
                     SL(M)=SL(M)+QL(M,L)*(QP(L)-QV(L)+AV(L))
                  ENDDO
                  
                  DO L=1,7
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO

               ENDIF
               
            ENDDO
            
            DO M=1,7 
               QUL(M)=QP(M)-VL(M)
            ENDDO
            
         ELSE

            DO M=1,7 
               QUL(M)=QVR(M)
            ENDDO

         ENDIF 

         IF(TC(1).LT.0.)THEN                  
            
            DXI=-TC(1)*DTX

            CALL FQML(DXI,TVL,TVR,TP,VM) 

            CALL VECTLRX(TP,QL,QR,UV)
            CALL MATR_AY(TP,AY)
            CALL MATR_AZ(TP,AZ)

            DO L=1,7
               VL(L)=0.
            ENDDO

            DO M=1,7
               AV(M)=0.
               DO L=1,7
                  AV(M)=AV(M)+AY(M,L)*TDy(L)+AZ(M,L)*TDz(L)
               ENDDO
               AV(M)=AV(M)*DTX/2.
            ENDDO

            DO M=1,7  
               
               IF(TC(M).LT.0.)THEN                  

                  DXI=-TC(M)*DTX

                  CALL FQML(DXI,TVL,TVR,TP,QV) 

                  SL(M)=0.   
                  DO L=1,7
                     SL(M)=SL(M)+QL(M,L)*(TP(L)-QV(L)+AV(L))
                  ENDDO
                  
                  DO L=1,7
                     VL(L)=VL(L)+QR(L,M)*SL(M)
                  ENDDO
                  
               ENDIF

            ENDDO
            
            DO M=1,7 
               QUR(M)=TP(M)-VL(M)
            ENDDO

         ELSE

            DO M=1,7 
               QUR(M)=TVL(M)
            ENDDO

         ENDIF 


         DO M=1,7
            TL(M)=QUL(M)
            TR(M)=QUR(M)
         ENDDO
         
         TL5=0.5d0*(TL(5)+TR(5))
         
         TL(5)=TL5
         TR(5)=TL5

         CALL CONS(TL,PX)
         CALL CONS(TR,VL)
         CALL HLLDx(PX,VL,FR,UR) !HLLD method


         CALL PRIM(UR,PX)

         DO M=1,7  
            QRP(I,M)=PX(M)      
         ENDDO 

      Enddo
      
      Return
      End
