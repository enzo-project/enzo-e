c     See LICENSE_PPML file for license and copyright information

#include "fortran.h"

      Subroutine Monot_ig(md,nx,gamma,dv1,dv2,bpx,qp,qm,qxl,qxr)
      Implicit NONE

      External Xmedian
	
      Integer mk,imd,md,nx,i,j,k,m

      ENZO_REAL q1(8),q2(8),q3(8),q4(8),q5(8)
      ENZO_REAL qpl(8),qpr(8),vl(8),vr(8),bs(3)
      ENZO_REAL a1(8),a2(8),a3(8),a4(8),a5(8)
      ENZO_REAL u1(8),u2(8),u3(8),u4(8),u5(8)
      ENZO_REAL b1(8),b2(8),c1(8),c2(8)
      ENZO_REAL ql(8,8),qr(8,8),uv(8)
      ENZO_REAL qp(nx,8),qm(nx,8),qxl(nx,8),qxr(nx,8)
      ENZO_REAL dv1(nx),dv2(nx),bpx(nx,3)
      
      ENZO_REAL uwl,uwr,vx1,vx2,vx3,vx4,vx5,vl1,vr1,dml,gamma
      ENZO_REAL dmr,dwl,dwr,sjn,sjm,sjx,sjy,uvr,uvl,Xmedian 
      ENZO_REAL o1

      o1 = 1.0d0
      
      Do i=1,Nx  
         Do m=1,8
            qxl(i,m)=0.
            qxr(i,m)=0.
         Enddo		  
      Enddo		  

      Do i=3,Nx-2  

         Do m=1,8
            Q1(m)=Qp(i-2,m)	 
            Q2(m)=Qp(i-1,m)	 
            Q3(m)=Qp(i,m)	 
            Q4(m)=Qp(i+1,m)	 
            Q5(m)=Qp(i+2,m)	 
         Enddo

         Do m=1,8
            Qpl(m)=qm(i-1,m)	 
            Qpr(m)=qm(i,m)	 
         Enddo
		 
         Do m=1,3
            bs(m)=Bpx(i,m)
         Enddo		  
		 
         if(md.eq.1)CALL VECTLRx_ig(gamma,q3,bs,QL,QR,UV)
         if(md.eq.2)CALL VECTLRy_ig(gamma,q3,bs,QL,QR,UV)
         if(md.eq.3)CALL VECTLRz_ig(gamma,q3,bs,QL,QR,UV)
        
         CALL AMPLTD_ig(QL,q1,u1)
         CALL AMPLTD_ig(QL,q2,u2)
         CALL AMPLTD_ig(QL,q3,u3)
         CALL AMPLTD_ig(QL,q4,u4)
         CALL AMPLTD_ig(QL,q5,u5)
         CALL AMPLTD_ig(QL,qpl,c1)
         CALL AMPLTD_ig(QL,qpr,c2)

         Do m=1,8  

            uwl = Xmedian(u3(m),c1(m),u2(m))
            uwr = Xmedian(u3(m),c2(m),u4(m))

            b1(m) = Xmedian(u3(m),uwl,3.d0*u3(m)-2.d0*uwr)
            b2(m) = Xmedian(u3(m),uwr,3.d0*u3(m)-2.d0*uwl)

         Enddo 

         imd=0

         Do m=1,8  

            if((b1(m)-c1(m))**2.gt.1e-12.or.(b2(m)-c2(m))**2.
     & 	          gt.1e-12) imd=1
	  
         enddo
	 
         If(imd.eq.0)then 	 

            CALL AMPLTD_ig(QR,c2,b2)
	 
            Do m=1,8  
               qxr(i,m)=b2(m)
            enddo

            CALL AMPLTD_ig(QR,c1,b1)
	 
            Do m=1,8  
               qxl(i,m)=b1(m)
            enddo

         Else 

            do m=1,8
	  
               vx1=u1(m)
               vx2=u2(m)
               vx3=u3(m)
               vx4=u4(m)
               vx5=u5(m)
               
               CALL W5RECM_ig(vx1,vx2,vx3,vx4,vx5,VL1,VR1)
               
               vl(m)=vl1 
               vr(m)=vr1 
               
            enddo

            Do m=1,8

               If((b1(m)-u3(m))**2.lt.1e-12.and.(b2(m)-u3(m))**2.
     &	  lt.1.e-12)then

                  dml = Xmedian(u3(m),vl(m),c1(m))
                  dmr = Xmedian(u3(m),vr(m),c2(m))
                  
                  dwl = Xmedian(u3(m),dml,u2(m))
                  dwr = Xmedian(u3(m),dmr,u4(m))

                  uwl = Xmedian(u3(m),3.d0*u3(m)-2.d0*dwr,dml)         
                  uwr = Xmedian(u3(m),3.d0*u3(m)-2.d0*dwl,dmr)
		
                  c1(m) = Xmedian(uwl,dml,c1(m))         
                  c2(m) = Xmedian(uwr,dmr,c2(m))

               Else
	  
                  sjn = 0.5d0*(u4(m)-u2(m))
                  sjm = 2.d0*Xmedian(0.,u4(m)-u3(m),u3(m)-u2(m))
                  sjx = Xmedian(0.,sjn,sjm)

                  sjn = 0.5d0*(u5(m)-u3(m))
                  sjm = 2.d0*Xmedian(0.,u5(m)-u4(m),u4(m)-u3(m))
                  sjy = Xmedian(0.,sjn,sjm)

                  dwr = 0.5d0*(u3(m)+u4(m))-(sjy-sjx)/6.d0

                  sjn = 0.5d0*(u2(m)-u4(m))
                  sjm = 2.d0*Xmedian(0.,u2(m)-u3(m),u3(m)-u4(m))
                  sjx = Xmedian(0.,sjn,sjm)

                  sjn = 0.5d0*(u1(m)-u3(m))
                  sjm = 2.d0*Xmedian(0.,u1(m)-u2(m),u2(m)-u3(m))
                  sjy = Xmedian(0.,sjn,sjm)

                  dwl = 0.5d0*(u3(m)+u2(m))-(sjy-sjx)/6.d0

                  uwr = Xmedian(dwr,vr(m),c2(m))
                  uwl = Xmedian(dwl,vl(m),c1(m))

                  dmr = Xmedian(u3(m),uwr,u4(m))
                  dml = Xmedian(u3(m),uwl,u2(m))

                  uvr = Xmedian(u3(m),3.d0*u3(m)-2.d0*dml,dmr)         
                  uvl = Xmedian(u3(m),3.d0*u3(m)-2.d0*dmr,dml)

                  c1(m) = Xmedian(uvl,vl(m),c1(m))         
                  c2(m) = Xmedian(uvr,vr(m),c2(m))

               Endif

            Enddo

            CALL AMPLTD_ig(QR,c2,b2)
	 
            Do m=1,8  
               qxr(i,m)=b2(m)
            enddo

            CALL AMPLTD_ig(QR,c1,b1)

            Do m=1,8  
               qxl(i,m)=b1(m)
            enddo

         Endif 

         if(b1(1).lt.1.e-5) then 

            Do m=1,8 

               uwl=q4(m)-q3(m)
               uwr=q3(m)-q2(m)
               
               u1(m)=0.5d0*(SIGN(o1,uwl)+SIGN(o1,uwr))*
     &              min(abs(uwl),abs(uwr))

            Enddo
	 
            Do m=1,8  
               qxl(i,m)=q3(m)-u1(m)/2.d0
            Enddo

            if(abs(qxl(i,1)-q3(1)).ge.0.8d0*q3(1))then

               Do m=1,8  
                  qxl(i,m)=q3(m)
               Enddo

            Endif

         Endif

         If(b2(1).lt.1.0e-5) then 

            Do m=1,8  

               uwl=q4(m)-q3(m)
               uwr=q3(m)-q2(m)

               u1(m)=0.5d0*(SIGN(o1,uwl)+SIGN(o1,uwr))*
     &              min(abs(uwl),abs(uwr))

            Enddo
            
            Do m=1,8  
               qxr(i,m)=q3(m)+u1(m)/2.d0
            Enddo

            If(abs(qxr(i,1)-q3(1)).ge.0.8d0*q3(1))then

               Do m=1,8  
                  qxr(i,m)=q3(m)
               Enddo

            Endif

         Endif

         if(abs(q4(1)-q3(1))/max(q4(1),q3(1)).gt.0.8d0.
     &        and.dv2(i).lt.0.)then

            do m=1,8
               qxr(i,m)=q3(m)
            enddo

         endif

         if(abs(q2(1)-q3(1))/max(q2(1),q3(1)).gt.0.8d0.
     &        and.dv1(i).lt.0.)then

            do m=1,8
               qxl(i,m)=q3(m)
            enddo

         endif
		 
         
      Enddo
		 
      Return
      End
