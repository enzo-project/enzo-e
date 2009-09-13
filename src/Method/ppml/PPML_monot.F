      Subroutine Monot(md,nx,dv1,dv2,qp,qm,qxl,qxr)
	 Implicit NONE

	     External Xmedian
	
	     Integer mk,imd,md,nx,i,j,k,m

	     REAL q1(7),q2(7),q3(7),q4(7),q5(7)
	     REAL qpl(7),qpr(7),vl(7),vr(7)
         REAL a1(7),a2(7),a3(7),a4(7),a5(7)
         REAL u1(7),u2(7),u3(7),u4(7),u5(7)
	     REAL b1(7),b2(7),c1(7),c2(7)
	     REAL ql(7,7),qr(7,7),uv(7)
         REAL qp(nx,7),qm(nx,7),qxl(nx,7),qxr(nx,7)
         REAL dv1(nx),dv2(nx)

         REAL uwl,uwr,vx1,vx2,vx3,vx4,vx5,vl1,vr1,dml
         REAL dmr,dwl,dwr,sjn,sjm,sjx,sjy,uvr,uvl,Xmedian 

         Do i=1,Nx  
           Do m=1,7
             qxl(i,m)=0.
		     qxr(i,m)=0.
           Enddo		  
         Enddo		  

         Do i=3,Nx-2  

         Do m=1,7
          Q1(m)=Qp(i-2,m)	 
          Q2(m)=Qp(i-1,m)	 
          Q3(m)=Qp(i,m)	 
          Q4(m)=Qp(i+1,m)	 
          Q5(m)=Qp(i+2,m)	 
         Enddo

         Do m=1,7
          Qpl(m)=qm(i-1,m)	 
          Qpr(m)=qm(i,m)	 
         Enddo
		 
         if(md.eq.1)CALL VECTLRx(q3,QL,QR,UV)
         if(md.eq.2)CALL VECTLRy(q3,QL,QR,UV)
         if(md.eq.3)CALL VECTLRz(q3,QL,QR,UV)
        
         CALL AMPLTD(QL,q1,u1)
         CALL AMPLTD(QL,q2,u2)
         CALL AMPLTD(QL,q3,u3)
         CALL AMPLTD(QL,q4,u4)
         CALL AMPLTD(QL,q5,u5)
         CALL AMPLTD(QL,qpl,c1)
         CALL AMPLTD(QL,qpr,c2)

c	  Do m=1,7  
c         c1(m)=(7.e0*(u2(m)+u3(m))-(u1(m)+u4(m)))/12.e0
c         c2(m)=(7.e0*(u3(m)+u4(m))-(u2(m)+u5(m)))/12.e0
c	  enddo

		     Do m=1,7  

	          uwl = Xmedian(u3(m),c1(m),u2(m))
	          uwr = Xmedian(u3(m),c2(m),u4(m))

	          b1(m) = Xmedian(u3(m),uwl,3.e0*u3(m)-2.e0*uwr)
	          b2(m) = Xmedian(u3(m),uwr,3.e0*u3(m)-2.e0*uwl)

	        Enddo 

			 imd=0

            Do m=1,7  

			   if((b1(m)-c1(m))**2.gt.1e-12.or.(b2(m)-c2(m))**2.
     1 	          gt.1e-12) imd=1
	  
            enddo
	 
      If(imd.eq.0)then 	 

        CALL AMPLTD(QR,c2,b2)
	 
         Do m=1,7  
	      qxr(i,m)=b2(m)
         enddo

        CALL AMPLTD(QR,c1,b1)
	 
         Do m=1,7  
          qxl(i,m)=b1(m)
         enddo

      Else 

        do m=1,7
	  
			vx1=u1(m)
			vx2=u2(m)
			vx3=u3(m)
			vx4=u4(m)
			vx5=u5(m)

			CALL W5RECM(vx1,vx2,vx3,vx4,vx5,VL1,VR1)

			vl(m)=vl1 
			vr(m)=vr1 

        enddo

        Do m=1,7

            If((b1(m)-u3(m))**2.lt.1e-12.and.(b2(m)-u3(m))**2.
     1	  lt.1.e-12)then

				dml = Xmedian(u3(m),vl(m),c1(m))
				dmr = Xmedian(u3(m),vr(m),c2(m))

				dwl = Xmedian(u3(m),dml,u2(m))
				dwr = Xmedian(u3(m),dmr,u4(m))

				uwl = Xmedian(u3(m),3.e0*u3(m)-2.e0*dwr,dml)         
				uwr = Xmedian(u3(m),3.e0*u3(m)-2.e0*dwl,dmr)
		
				c1(m) = Xmedian(uwl,dml,c1(m))         
				c2(m) = Xmedian(uwr,dmr,c2(m))

            Else
	  
				sjn = 0.5e0*(u4(m)-u2(m))
				sjm = 2.e0*Xmedian(0.,u4(m)-u3(m),u3(m)-u2(m))
				sjx = Xmedian(0.,sjn,sjm)

				sjn = 0.5e0*(u5(m)-u3(m))
				sjm = 2.e0*Xmedian(0.,u5(m)-u4(m),u4(m)-u3(m))
				sjy = Xmedian(0.,sjn,sjm)

				dwr = 0.5e0*(u3(m)+u4(m))-(sjy-sjx)/6.e0

				sjn = 0.5e0*(u2(m)-u4(m))
				sjm = 2.e0*Xmedian(0.,u2(m)-u3(m),u3(m)-u4(m))
				sjx = Xmedian(0.,sjn,sjm)

				sjn = 0.5e0*(u1(m)-u3(m))
				sjm = 2.e0*Xmedian(0.,u1(m)-u2(m),u2(m)-u3(m))
				sjy = Xmedian(0.,sjn,sjm)

				dwl = 0.5e0*(u3(m)+u2(m))-(sjy-sjx)/6.e0

				uwr = Xmedian(dwr,vr(m),c2(m))
				uwl = Xmedian(dwl,vl(m),c1(m))

				dmr = Xmedian(u3(m),uwr,u4(m))
				dml = Xmedian(u3(m),uwl,u2(m))

				uvr = Xmedian(u3(m),3.e0*u3(m)-2.e0*dml,dmr)         
				uvl = Xmedian(u3(m),3.e0*u3(m)-2.e0*dmr,dml)

				c1(m) = Xmedian(uvl,vl(m),c1(m))         
				c2(m) = Xmedian(uvr,vr(m),c2(m))

            Endif

        Enddo

        CALL AMPLTD(QR,c2,b2)
	 
      Do m=1,7  
		qxr(i,m)=b2(m)
      enddo

      CALL AMPLTD(QR,c1,b1)

      Do m=1,7  
		qxl(i,m)=b1(m)
      enddo

      Endif 

         if(b1(1).lt.1.e-5) then 

            Do m=1,7 

	          uwl=q4(m)-q3(m)
	          uwr=q3(m)-q2(m)

	          u1(m)=0.5e0*(SIGN(1.e0,uwl)+SIGN(1.e0,uwr))*
     1             min(abs(uwl),abs(uwr))

            Enddo
	 
            Do m=1,7  

		      qxl(i,m)=q3(m)-u1(m)/2.e0

	        Enddo

            if(abs(qxl(i,1)-q3(1)).ge.0.8e0*q3(1))then

                Do m=1,7  

		          qxl(i,m)=q3(m)

	            Enddo

            Endif

         Endif

         If(b2(1).lt.1.0e-5) then 

             Do m=1,7  

	          uwl=q4(m)-q3(m)
	          uwr=q3(m)-q2(m)

	          u1(m)=0.5e0*(SIGN(1.e0,uwl)+SIGN(1.e0,uwr))*
     1             min(abs(uwl),abs(uwr))

             Enddo
	 
             Do m=1,7  

 		      qxr(i,m)=q3(m)+u1(m)/2.e0

	         Enddo

             If(abs(qxr(i,1)-q3(1)).ge.0.8e0*q3(1))then

               Do m=1,7  

	            qxr(i,m)=q3(m)

	           Enddo

             Endif

         Endif

         if(abs(q4(1)-q3(1))/max(q4(1),q3(1)).gt.0.8e0.
     1     and.dv2(i).lt.0.)then

          do m=1,7
	        qxr(i,m)=q3(m)
          enddo

         endif

         if(abs(q2(1)-q3(1))/max(q2(1),q3(1)).gt.0.8e0.
     1	     and.dv1(i).lt.0.)then

          do m=1,7
	        qxl(i,m)=q3(m)
          enddo

         endif
		 
		 
        Enddo
		 
      Return
      End
