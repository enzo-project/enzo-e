/* #include "error.def" */

      Subroutine PPML(dn,vx,vy,vz,bx,by,bz,
     1                 dnrx,vxrx,vyrx,vzrx,bxrx,byrx,bzrx,	 
     1                 dnry,vxry,vyry,vzry,bxry,byry,bzry,	 
     1                 dnrz,vxrz,vyrz,vzrz,bxrz,byrz,bzrz,
     1                 dt, CycleNumber, dx,dy,dz,
     1                 nx, ny, nz,
     1                 GridStartIndex, GridEndIndex,
     1                 nsubgrids, lface, rface,
     1                 fistart, fiend, fjstart, fjend,
     1                 array, dnindex,
     1                 vxindex, vyindex, vzindex,
     1                 bxindex, byindex, bzindex)	 

	 Implicit NONE

!	 i1,j1,k1 - initial point of grid
!	 i2,j2,k2 - last point of grid
!    ighost,jghost,kghost - number of ghost zones	 
!    dn,vx,vy,vz,bx,by,bz - central primitive variables
!    dnrx,vxrx,vyrx,vzrx,bxrx,byrx,bzrx - x-right primitive variables
!    dnry,vxry,vyry,vzry,bxry,byry,bzry - y-right primitive variables
!    dnrz,vxrz,vyrz,vzrz,bxrz,byrz,bzrz - z-right primitive variables 
!    f1,f2,f3,f4,f5,f6,f7 - x-fluxes
!	 g1,g2,g3,g4,g5,g6,g7 - y-fluxes
!    h1,h2,h3,h4,h5,h6,h7 - z-fluxes 

      Integer nx,ny,nz,ighost,jghost,kghost,nsubgrids
	  Integer nxb,nxe,nyb,nye,nzb,nze 
	  Integer i1,j1,k1,i2,j2,k2,CycleNumber
	  Integer nxp,nyp,nzp,i,j,k,m,n,idim,offset
      Integer GridStartIndex(3), GridEndIndex(3)
	  Integer m1,m2,m3
	  
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     1        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     1        lface(nsubgrids*3)  , rface(nsubgrids*3)
      integer dnindex(nsubgrids*6), vxindex(nsubgrids*6),
     1        vyindex(nsubgrids*6), vzindex(nsubgrids*6),
     1        bxindex(nsubgrids*6), byindex(nsubgrids*6),
     1        bzindex(nsubgrids*6)
	 
      real    array(1)	 
	  
	  
      real dn(nx,ny,nz)  	 
      real vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)  	 
      real bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)  	 

      real dnrx(nx,ny,nz)  	 
      real vxrx(nx,ny,nz),vyrx(nx,ny,nz),vzrx(nx,ny,nz)  	 
      real bxrx(nx,ny,nz),byrx(nx,ny,nz),bzrx(nx,ny,nz)  	 

      real dnry(nx,ny,nz)  	 
      real vxry(nx,ny,nz),vyry(nx,ny,nz),vzry(nx,ny,nz)  	 
      real bxry(nx,ny,nz),byry(nx,ny,nz),bzry(nx,ny,nz)  	 

      real dnrz(nx,ny,nz)  	 
      real vxrz(nx,ny,nz),vyrz(nx,ny,nz),vzrz(nx,ny,nz)  	 
      real bxrz(nx,ny,nz),byrz(nx,ny,nz),bzrz(nx,ny,nz)  	 

      real f1(nx,ny,nz)  	 
      real f2(nx,ny,nz),f3(nx,ny,nz),f4(nx,ny,nz)  	 
      real f5(nx,ny,nz),f6(nx,ny,nz),f7(nx,ny,nz)  	 

      real g1(nx,ny,nz)  	 
      real g2(nx,ny,nz),g3(nx,ny,nz),g4(nx,ny,nz)  	 
      real g5(nx,ny,nz),g6(nx,ny,nz),g7(nx,ny,nz)  	 

      real h1(nx,ny,nz)  	 
      real h2(nx,ny,nz),h3(nx,ny,nz),h4(nx,ny,nz)  	 
      real h5(nx,ny,nz),h6(nx,ny,nz),h7(nx,ny,nz)  	 

      real ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)  	 
	 
      real qu1(nx,ny,nz)  	 
      real qu2(nx,ny,nz),qu3(nx,ny,nz),qu4(nx,ny,nz)  	 
      real qu5(nx,ny,nz),qu6(nx,ny,nz),qu7(nx,ny,nz)  
	
!bug      real bxr(nx,ny,nz),byr(nx,ny,nz),bzr(nx,ny,nz)  
		 

	  real Qp1(nx,7),Qr1(nx,7)
	  real Qdy1(nx,7),Qdz1(nx,7)
	  real Qlx1(nx,7),Qrx1(nx,7),Qvl1(nx,7),Qvr1(nx,7)
	  real Qrp1(nx,7),Fm1(nx,7)
	  real dx1(nx),dx2(nx)
	  
	  real Qp2(ny,7),Qr2(ny,7)
	  real Qdx2(ny,7),Qdz2(ny,7)
	  real Qlx2(ny,7),Qrx2(ny,7),Qvl2(ny,7),Qvr2(ny,7)
	  real Qrp2(ny,7),Fm2(ny,7)
	  real dy1(ny),dy2(ny)

	  real Qp3(nz,7),Qr3(nz,7)
	  real Qdx3(nz,7),Qdy3(nz,7)
	  real Qlx3(nz,7),Qrx3(nz,7),Qvl3(nz,7),Qvr3(nz,7)
	  real Qrp3(nz,7),Fm3(nz,7)
	  real dz1(nz),dz2(nz)

	real df(nx,ny),vf(nx,ny),uf(nx,ny),wf(nx,ny)
	real bf(nx,ny),tf(nx,ny),pf(nx,ny)
	real dg(ny,nz),vg(ny,nz),ug(ny,nz),wg(ny,nz)
	real bg(ny,nz),tg(ny,nz),pg(ny,nz)	  
	real dh(nz,nx),vh(nz,nx),uh(nz,nx),wh(nz,nx)
	real bh(nz,nx),th(nz,nx),ph(nz,nx)	  
	  
	  real vax,vay,vaz,cga,cgs
	  real dn1,dn2,vx1,vx2,vy1,vy2,vz1,vz2
	  real bx1,bx2,by1,by2,bz1,bz2

      real epsk,v2k,exc,dexy,dexy1,dexy2,exi1,exi2,exk1,exk2
	  real Vleer,dsc,dx,dy,dz,dt 
	  
	  real qh(7),qh5,qh6,qh7

c
c check-it
c
        if (nx.le.0) then
C           ERROR_MESSAGE
        endif
	  
      i1=GridStartIndex(1)+1
      j1=GridStartIndex(2)+1
      k1=GridStartIndex(3)+1

      i2=GridEndIndex(1)+1
      j2=GridEndIndex(2)+1
      k2=GridEndIndex(3)+1

      ighost=GridStartIndex(1)
      jghost=GridStartIndex(2)
      kghost=GridStartIndex(3)	  
!SU      kghost=GridEndIndex(3)	  

        if (i1.eq.0) then
C           ERROR_MESSAGE
        endif

	  
!      ighost=3
!      jghost=3
!      kghost=3 	  
	  
!	  i1=4
!	  j1=4
!	  k1=4
	 
!	  i2=Nx-3
!	  j2=Ny-3
!	  k2=Nz-3

	 
	  nxb=i1-ighost 
	  nxe=i2+ighost 

	  nyb=j1-jghost 
	  nye=j2+jghost 

	  nzb=k1-kghost 
	  nze=k2+kghost 

      nxp=nxe-nxb+1 	  
      nyp=nye-nyb+1 	  
      nzp=nze-nzb+1 	  

c-------------------------------------------------------	  
c
c     Calculation of flux in x - direction
c	  
c--------------------------------------------------------       

      Do k=nzb+1,nze-1
      Do j=nyb+1,nye-1

		Do i=nxb,nxe

			Qp1(i,1)=dn(i,j,k)	 
			Qp1(i,2)=vx(i,j,k)	 
			Qp1(i,3)=vy(i,j,k)	 
			Qp1(i,4)=vz(i,j,k)	 
			Qp1(i,5)=bx(i,j,k)	 
			Qp1(i,6)=by(i,j,k)	 
			Qp1(i,7)=bz(i,j,k)	 

			Qr1(i,1)=dnrx(i,j,k)	 
			Qr1(i,2)=vxrx(i,j,k)	 
			Qr1(i,3)=vyrx(i,j,k)	 
			Qr1(i,4)=vzrx(i,j,k)	 
			Qr1(i,5)=bxrx(i,j,k)	 
			Qr1(i,6)=byrx(i,j,k)	 
			Qr1(i,7)=bzrx(i,j,k)	
           			
            Qdy1(i,1)=VLeer(dn(i,j+1,k)-dn(i,j,k),dn(i,j,k)-dn(i,j-1,k))
            Qdy1(i,2)=VLeer(vx(i,j+1,k)-vx(i,j,k),vx(i,j,k)-vx(i,j-1,k))
            Qdy1(i,3)=VLeer(vy(i,j+1,k)-vy(i,j,k),vy(i,j,k)-vy(i,j-1,k))
            Qdy1(i,4)=VLeer(vz(i,j+1,k)-vz(i,j,k),vz(i,j,k)-vz(i,j-1,k))
            Qdy1(i,5)=VLeer(bx(i,j+1,k)-bx(i,j,k),bx(i,j,k)-bx(i,j-1,k))
            Qdy1(i,6)=VLeer(by(i,j+1,k)-by(i,j,k),by(i,j,k)-by(i,j-1,k))
            Qdy1(i,7)=VLeer(bz(i,j+1,k)-bz(i,j,k),bz(i,j,k)-bz(i,j-1,k))

            Qdz1(i,1)=VLeer(dn(i,j,k+1)-dn(i,j,k),dn(i,j,k)-dn(i,j,k-1))
            Qdz1(i,2)=VLeer(vx(i,j,k+1)-vx(i,j,k),vx(i,j,k)-vx(i,j,k-1))
            Qdz1(i,3)=VLeer(vy(i,j,k+1)-vy(i,j,k),vy(i,j,k)-vy(i,j,k-1))
            Qdz1(i,4)=VLeer(vz(i,j,k+1)-vz(i,j,k),vz(i,j,k)-vz(i,j,k-1))
            Qdz1(i,5)=VLeer(bx(i,j,k+1)-bx(i,j,k),bx(i,j,k)-bx(i,j,k-1))
            Qdz1(i,6)=VLeer(by(i,j,k+1)-by(i,j,k),by(i,j,k)-by(i,j,k-1))
            Qdz1(i,7)=VLeer(bz(i,j,k+1)-bz(i,j,k),bz(i,j,k)-bz(i,j,k-1))

        Enddo 

		Do i=nxb,nxe

			dx1(i)=0.
			dx2(i)=0.

		Enddo	
			
		Do i=nxb+2,nxe-2
 
            vax=0.5*(vx(i+2,j,k)+vx(i+1,j,k))
            vay=0.5*(vx(i,j,k)+vx(i-1,j,k))

            vaz=0.5*(vy(i,j+1,k)+vy(i+1,j+1,k))
            cga=0.5*(vy(i,j-1,k)+vy(i+1,j-1,k))

            dsc=0.5*(vz(i,j,k+1)+vz(i+1,j,k+1))
            cgs=0.5*(vz(i,j,k-1)+vz(i+1,j,k-1))

            dx2(i)=vax-vay+vaz-cga+dsc-cgs

            vax=0.5*(vx(i+1,j,k)+vx(i,j,k))
            vay=0.5*(vx(i-1,j,k)+vx(i-2,j,k))

            vaz=0.5*(vy(i-1,j+1,k)+vy(i,j+1,k))
            cga=0.5*(vy(i-1,j-1,k)+vy(i,j-1,k))

            dsc=0.5*(vz(i-1,j,k+1)+vz(i,j,k+1))
            cgs=0.5*(vz(i-1,j,k-1)+vz(i,j,k-1))

            dx1(i)=vax-vay+vaz-cga+dsc-cgs			
		
		Enddo

        Call Monot(1,nxp,Dx1,Dx2,qp1,qr1,qlx1,qrx1)
		
        Do i=nxb+1,nxe
			Qlx1(i,5)=Qr1(i-1,5)
			Qrx1(i,5)=Qr1(i,5)
		Enddo

        Do i=nxb+1,nxe-1

       		 dn1=Qp1(i,1)
             dn2=dn1

       		 vx1=Qp1(i,2)
             vx2=vx1
       		 vy1=Qp1(i,3)
             vy2=vy1
       		 vz1=Qp1(i,4)
             vz2=vz1
			 
       		 bx1=Qp1(i,5)
             bx2=bx1
       		 by1=Qp1(i,6)
             by2=by1
       		 bz1=Qp1(i,7)
             bz2=bz1

			 do m3=k-1,k+1 
	         do m2=j-1,j+1 
	         do m1=i-1,i+1 

                dn1=min(dn1,dn(m1,m2,m3))   
                dn2=max(dn2,dn(m1,m2,m3))   

                vx1=min(vx1,vx(m1,m2,m3))   
                vx2=max(vx2,vx(m1,m2,m3))   

                vy1=min(vy1,vy(m1,m2,m3))   
                vy2=max(vy2,vy(m1,m2,m3))   

                vz1=min(vz1,vz(m1,m2,m3))   
                vz2=max(vz2,vz(m1,m2,m3))   

                bx1=min(bx1,bx(m1,m2,m3))   
                bx2=max(bx2,bx(m1,m2,m3))   

                by1=min(by1,by(m1,m2,m3))   
                by2=max(by2,by(m1,m2,m3))   

                bz1=min(bz1,bz(m1,m2,m3))   
                bz2=max(bz2,bz(m1,m2,m3))   

			Enddo 
            Enddo 
            Enddo 

			qlx1(i,1)=max(dn1,min(dn2,qlx1(i,1)))
			qrx1(i,1)=max(dn1,min(dn2,qrx1(i,1)))
		
			qlx1(i,2)=max(vx1,min(vx2,qlx1(i,2)))
			qrx1(i,2)=max(vx1,min(vx2,qrx1(i,2)))
		 
			qlx1(i,3)=max(vy1,min(vy2,qlx1(i,3)))
			qrx1(i,3)=max(vy1,min(vy2,qrx1(i,3)))

			qlx1(i,4)=max(vz1,min(vz2,qlx1(i,4)))
			qrx1(i,4)=max(vz1,min(vz2,qrx1(i,4)))
			
			qlx1(i,6)=max(by1,min(by2,qlx1(i,6)))
			qrx1(i,6)=max(by1,min(by2,qrx1(i,6)))

			qlx1(i,7)=max(bz1,min(bz2,qlx1(i,7)))
			qrx1(i,7)=max(bz1,min(bz2,qrx1(i,7)))

		Enddo
		
			CALL QDD6(nxp,Qrx1,Qlx1,Qvr1,Qvl1,Qp1)

			Do i=nxb,nxe
				Qvl1(i,5)=Qlx1(i,5)
				Qvr1(i,5)=Qrx1(i,5)
			Enddo
			
		    CALL POTOKx(nxp,Dx,Dt,Qvl1,Qvr1,Qp1,Qdy1,Qdz1,Qrp1,Fm1)

!bug         Do i=i1,i2
          Do i=i1-1,i2

            dnrx(i,j,k)=Qrp1(i,1)		
            vxrx(i,j,k)=Qrp1(i,2)		
            vyrx(i,j,k)=Qrp1(i,3)		
            vzrx(i,j,k)=Qrp1(i,4)		
            byrx(i,j,k)=Qrp1(i,6)		
            bzrx(i,j,k)=Qrp1(i,7)

            F1(i,j,k)=fm1(i,1)			
            F2(i,j,k)=fm1(i,2)			
            F3(i,j,k)=fm1(i,3)			
            F4(i,j,k)=fm1(i,4)			
            F5(i,j,k)=fm1(i,5)			
            F6(i,j,k)=fm1(i,6)			
            F7(i,j,k)=fm1(i,7)

         Enddo			

	 Enddo
	 Enddo
	 
c-------------------------------------------------------	  
c
c     Calculation of flux in y - direction
c	  
c--------------------------------------------------------       

      Do k=nzb+1,nze-1
      Do i=nxb+1,nxe-1

		Do j=nyb,nye

			Qp2(j,1)=dn(i,j,k)	 
			Qp2(j,2)=vx(i,j,k)	 
			Qp2(j,3)=vy(i,j,k)	 
			Qp2(j,4)=vz(i,j,k)	 
			Qp2(j,5)=bx(i,j,k)	 
			Qp2(j,6)=by(i,j,k)	 
			Qp2(j,7)=bz(i,j,k)	 

			Qr2(j,1)=dnry(i,j,k)	 
			Qr2(j,2)=vxry(i,j,k)	 
			Qr2(j,3)=vyry(i,j,k)	 
			Qr2(j,4)=vzry(i,j,k)	 
			Qr2(j,5)=bxry(i,j,k)	 
			Qr2(j,6)=byry(i,j,k)	 
			Qr2(j,7)=bzry(i,j,k)	
           			
            Qdx2(j,1)=VLeer(dn(i+1,j,k)-dn(i,j,k),dn(i,j,k)-dn(i-1,j,k))
            Qdx2(j,2)=VLeer(vx(i+1,j,k)-vx(i,j,k),vx(i,j,k)-vx(i-1,j,k))
            Qdx2(j,3)=VLeer(vy(i+1,j,k)-vy(i,j,k),vy(i,j,k)-vy(i-1,j,k))
            Qdx2(j,4)=VLeer(vz(i+1,j,k)-vz(i,j,k),vz(i,j,k)-vz(i-1,j,k))
            Qdx2(j,5)=VLeer(bx(i+1,j,k)-bx(i,j,k),bx(i,j,k)-bx(i-1,j,k))
            Qdx2(j,6)=VLeer(by(i+1,j,k)-by(i,j,k),by(i,j,k)-by(i-1,j,k))
            Qdx2(j,7)=VLeer(bz(i+1,j,k)-bz(i,j,k),bz(i,j,k)-bz(i-1,j,k))

            Qdz2(j,1)=VLeer(dn(i,j,k+1)-dn(i,j,k),dn(i,j,k)-dn(i,j,k-1))
            Qdz2(j,2)=VLeer(vx(i,j,k+1)-vx(i,j,k),vx(i,j,k)-vx(i,j,k-1))
            Qdz2(j,3)=VLeer(vy(i,j,k+1)-vy(i,j,k),vy(i,j,k)-vy(i,j,k-1))
            Qdz2(j,4)=VLeer(vz(i,j,k+1)-vz(i,j,k),vz(i,j,k)-vz(i,j,k-1))
            Qdz2(j,5)=VLeer(bx(i,j,k+1)-bx(i,j,k),bx(i,j,k)-bx(i,j,k-1))
            Qdz2(j,6)=VLeer(by(i,j,k+1)-by(i,j,k),by(i,j,k)-by(i,j,k-1))
            Qdz2(j,7)=VLeer(bz(i,j,k+1)-bz(i,j,k),bz(i,j,k)-bz(i,j,k-1))

        Enddo 

		Do j=nyb,nye
    		
			dy1(j)=0.
			dy2(j)=0.

		Enddo	
		
		Do j=nyb+2,nye-2

         vax=0.5*(vy(i,j+2,k)+vy(i,j+1,k))
         vay=0.5*(vy(i,j,k)+vy(i,j-1,k))

         vaz=0.5*(vx(i+1,j,k)+vx(i+1,j+1,k))
         cga=0.5*(vx(i-1,j,k)+vx(i-1,j+1,k))

         dsc=0.5*(vz(i,j,k+1)+vz(i,j+1,k+1))
         cgs=0.5*(vz(i,j,k-1)+vz(i,j+1,k-1))

         dy2(j)=vax-vay+vaz-cga+dsc-cgs

         vax=0.5*(vy(i,j+1,k)+vy(i,j,k))
         vay=0.5*(vy(i,j-1,k)+vy(i,j-2,k))

         vaz=0.5*(vx(i+1,j-1,k)+vx(i+1,j,k))
         cga=0.5*(vx(i-1,j-1,k)+vx(i-1,j,k))

         dsc=0.5*(vz(i,j-1,k+1)+vz(i,j,k+1))
         cgs=0.5*(vz(i,j-1,k-1)+vz(i,j,k-1))

         dy1(j)=vax-vay+vaz-cga+dsc-cgs


		Enddo

        Call Monot(2,nyp,Dy1,Dy2,qp2,qr2,qlx2,qrx2)
		
        Do j=nyb+1,nye
			Qlx2(j,6)=Qr2(j-1,6)
			Qrx2(j,6)=Qr2(j,6)
		Enddo

        Do j=nyb+1,nye-1

       		 dn1=Qp2(j,1)
             dn2=dn1

       		 vx1=Qp2(j,2)
             vx2=vx1
       		 vy1=Qp2(j,3)
             vy2=vy1
       		 vz1=Qp2(j,4)
             vz2=vz1
			 
       		 bx1=Qp2(j,5)
             bx2=bx1
       		 by1=Qp2(j,6)
             by2=by1
       		 bz1=Qp2(j,7)
             bz2=bz1

			 do m3=k-1,k+1 
	         do m2=j-1,j+1 
	         do m1=i-1,i+1 

                dn1=min(dn1,dn(m1,m2,m3))   
                dn2=max(dn2,dn(m1,m2,m3))   

                vx1=min(vx1,vx(m1,m2,m3))   
                vx2=max(vx2,vx(m1,m2,m3))   

                vy1=min(vy1,vy(m1,m2,m3))   
                vy2=max(vy2,vy(m1,m2,m3))   

                vz1=min(vz1,vz(m1,m2,m3))   
                vz2=max(vz2,vz(m1,m2,m3))   

                bx1=min(bx1,bx(m1,m2,m3))   
                bx2=max(bx2,bx(m1,m2,m3))   

                by1=min(by1,by(m1,m2,m3))   
                by2=max(by2,by(m1,m2,m3))   

                bz1=min(bz1,bz(m1,m2,m3))   
                bz2=max(bz2,bz(m1,m2,m3))   

			Enddo 
            Enddo 
            Enddo 

			qlx2(j,1)=max(dn1,min(dn2,qlx2(j,1)))
			qrx2(j,1)=max(dn1,min(dn2,qrx2(j,1)))
		
			qlx2(j,2)=max(vx1,min(vx2,qlx2(j,2)))
			qrx2(j,2)=max(vx1,min(vx2,qrx2(j,2)))
		 
			qlx2(j,3)=max(vy1,min(vy2,qlx2(j,3)))
			qrx2(j,3)=max(vy1,min(vy2,qrx2(j,3)))

			qlx2(j,4)=max(vz1,min(vz2,qlx2(j,4)))
			qrx2(j,4)=max(vz1,min(vz2,qrx2(j,4)))
			
!bug			qlx2(j,6)=max(by1,min(by2,qlx2(j,6)))
!bug			qrx2(j,6)=max(by1,min(by2,qrx2(j,6)))

			qlx2(j,5)=max(bx1,min(bx2,qlx2(j,5)))
			qrx2(j,5)=max(bx1,min(bx2,qrx2(j,5)))

			qlx2(j,7)=max(bz1,min(bz2,qlx2(j,7)))
			qrx2(j,7)=max(bz1,min(bz2,qrx2(j,7)))

		Enddo
		
			CALL QDD6(nyp,Qrx2,Qlx2,Qvr2,Qvl2,Qp2)

			Do j=nyb,nye
				Qvl2(j,6)=Qlx2(j,6)
				Qvr2(j,6)=Qrx2(j,6)
			Enddo
			
		    CALL POTOKy(nyp,Dx,Dt,Qvl2,Qvr2,Qp2,Qdx2,Qdz2,Qrp2,Fm2) 
			

!bug         Do j=j1,j2
         Do j=j1-1,j2

            dnry(i,j,k)=Qrp2(j,1)		
            vxry(i,j,k)=Qrp2(j,2)		
            vyry(i,j,k)=Qrp2(j,3)		
            vzry(i,j,k)=Qrp2(j,4)		
            bxry(i,j,k)=Qrp2(j,5)		
            bzry(i,j,k)=Qrp2(j,7)

            G1(i,j,k)=Fm2(j,1)			
            G2(i,j,k)=Fm2(j,2)			
            G3(i,j,k)=Fm2(j,3)			
            G4(i,j,k)=Fm2(j,4)			
            G5(i,j,k)=Fm2(j,5)			
            G6(i,j,k)=Fm2(j,6)			
            G7(i,j,k)=Fm2(j,7)

         Enddo				

	 Enddo
	 Enddo

c-------------------------------------------------------	  
c
c     Calculation of flux in z - direction
c	  
c--------------------------------------------------------    

      Do j=nyb+1,nye-1
      Do i=nxb+1,nxe-1

		Do k=nzb,nze

			Qp3(k,1)=dn(i,j,k)	 
			Qp3(k,2)=vx(i,j,k)	 
			Qp3(k,3)=vy(i,j,k)	 
			Qp3(k,4)=vz(i,j,k)	 
			Qp3(k,5)=bx(i,j,k)	 
			Qp3(k,6)=by(i,j,k)	 
			Qp3(k,7)=bz(i,j,k)	 

			Qr3(k,1)=dnrz(i,j,k)	 
			Qr3(k,2)=vxrz(i,j,k)	 
			Qr3(k,3)=vyrz(i,j,k)	 
			Qr3(k,4)=vzrz(i,j,k)	 
			Qr3(k,5)=bxrz(i,j,k)	 
			Qr3(k,6)=byrz(i,j,k)	 
			Qr3(k,7)=bzrz(i,j,k)	
           			
            Qdx3(k,1)=VLeer(dn(i+1,j,k)-dn(i,j,k),dn(i,j,k)-dn(i-1,j,k))
            Qdx3(k,2)=VLeer(vx(i+1,j,k)-vx(i,j,k),vx(i,j,k)-vx(i-1,j,k))
            Qdx3(k,3)=VLeer(vy(i+1,j,k)-vy(i,j,k),vy(i,j,k)-vy(i-1,j,k))
            Qdx3(k,4)=VLeer(vz(i+1,j,k)-vz(i,j,k),vz(i,j,k)-vz(i-1,j,k))
            Qdx3(k,5)=VLeer(bx(i+1,j,k)-bx(i,j,k),bx(i,j,k)-bx(i-1,j,k))
            Qdx3(k,6)=VLeer(by(i+1,j,k)-by(i,j,k),by(i,j,k)-by(i-1,j,k))
            Qdx3(k,7)=VLeer(bz(i+1,j,k)-bz(i,j,k),bz(i,j,k)-bz(i-1,j,k))

            Qdy3(k,1)=VLeer(dn(i,j+1,k)-dn(i,j,k),dn(i,j,k)-dn(i,j-1,k))
            Qdy3(k,2)=VLeer(vx(i,j+1,k)-vx(i,j,k),vx(i,j,k)-vx(i,j-1,k))
            Qdy3(k,3)=VLeer(vy(i,j+1,k)-vy(i,j,k),vy(i,j,k)-vy(i,j-1,k))
            Qdy3(k,4)=VLeer(vz(i,j+1,k)-vz(i,j,k),vz(i,j,k)-vz(i,j-1,k))
            Qdy3(k,5)=VLeer(bx(i,j+1,k)-bx(i,j,k),bx(i,j,k)-bx(i,j-1,k))
            Qdy3(k,6)=VLeer(by(i,j+1,k)-by(i,j,k),by(i,j,k)-by(i,j-1,k))
            Qdy3(k,7)=VLeer(bz(i,j+1,k)-bz(i,j,k),bz(i,j,k)-bz(i,j-1,k))

        Enddo 

		Do k=nzb,nze
    		
			dz1(k)=0.
			dz2(k)=0.

		Enddo	
		
		Do k=nzb+2,nze-2

         vax=0.5*(vz(i,j,k+2)+vz(i,j,k+1))
         vay=0.5*(vz(i,j,k)+vz(i,j,k-1))

         vaz=0.5*(vx(i+1,j,k)+vx(i+1,j,k+1))
         cga=0.5*(vx(i-1,j,k)+vx(i-1,j,k+1))

         dsc=0.5*(vy(i,j+1,k)+vy(i,j+1,k+1))
         cgs=0.5*(vy(i,j-1,k)+vy(i,j-1,k+1))

         dz2(k)=vax-vay+vaz-cga+dsc-cgs

         vax=0.5*(vz(i,j,k+1)+vz(i,j,k))
         vay=0.5*(vz(i,j,k-1)+vz(i,j,k-2))

         vaz=0.5*(vx(i+1,j,k-1)+vx(i+1,j,k))
         cga=0.5*(vx(i-1,j,k-1)+vx(i-1,j,k))

         dsc=0.5*(vy(i,j+1,k-1)+vy(i,j+1,k))
         cgs=0.5*(vy(i,j-1,k-1)+vy(i,j-1,k))

         dz1(k)=vax-vay+vaz-cga+dsc-cgs

		 
		Enddo

        Call Monot(3,nzp,Dz1,Dz2,Qp3,Qr3,qlx3,qrx3)
		
        Do k=nzb+1,nze
			Qlx3(k,7)=Qr3(k-1,7)
			Qrx3(k,7)=Qr3(k,7)
		Enddo

        Do k=nzb+1,nze-1

       		 dn1=Qp3(k,1)
             dn2=dn1

       		 vx1=Qp3(k,2)
             vx2=vx1
       		 vy1=Qp3(k,3)
             vy2=vy1
       		 vz1=Qp3(k,4)
             vz2=vz1
			 
       		 bx1=Qp3(k,5)
             bx2=bx1
       		 by1=Qp3(k,6)
             by2=by1
       		 bz1=Qp3(k,7)
             bz2=bz1

			 do m3=k-1,k+1 
	         do m2=j-1,j+1 
	         do m1=i-1,i+1 

                dn1=min(dn1,dn(m1,m2,m3))   
                dn2=max(dn2,dn(m1,m2,m3))   

                vx1=min(vx1,vx(m1,m2,m3))   
                vx2=max(vx2,vx(m1,m2,m3))   

                vy1=min(vy1,vy(m1,m2,m3))   
                vy2=max(vy2,vy(m1,m2,m3))   

                vz1=min(vz1,vz(m1,m2,m3))   
                vz2=max(vz2,vz(m1,m2,m3))   

                bx1=min(bx1,bx(m1,m2,m3))   
                bx2=max(bx2,bx(m1,m2,m3))   

                by1=min(by1,by(m1,m2,m3))   
                by2=max(by2,by(m1,m2,m3))   

                bz1=min(bz1,bz(m1,m2,m3))   
                bz2=max(bz2,bz(m1,m2,m3))   

			Enddo 
            Enddo 
            Enddo 

			qlx3(k,1)=max(dn1,min(dn2,qlx3(k,1)))
			qrx3(k,1)=max(dn1,min(dn2,qrx3(k,1)))
		
			qlx3(k,2)=max(vx1,min(vx2,qlx3(k,2)))
			qrx3(k,2)=max(vx1,min(vx2,qrx3(k,2)))
		 
			qlx3(k,3)=max(vy1,min(vy2,qlx3(k,3)))
			qrx3(k,3)=max(vy1,min(vy2,qrx3(k,3)))

			qlx3(k,4)=max(vz1,min(vz2,qlx3(k,4)))
			qrx3(k,4)=max(vz1,min(vz2,qrx3(k,4)))
			
			qlx3(k,5)=max(bx1,min(bx2,qlx3(k,5)))
			qrx3(k,5)=max(bx1,min(bx2,qrx3(k,5)))

			qlx3(k,6)=max(by1,min(by2,qlx3(k,6)))
			qrx3(k,6)=max(by1,min(by2,qrx3(k,6)))

!bug			qlx3(k,7)=max(bz1,min(bz2,qlx3(k,7)))
!bug			qrx3(k,7)=max(bz1,min(bz2,qrx3(k,7)))

		Enddo
		
			CALL QDD6(nzp,Qrx3,Qlx3,Qvr3,Qvl3,Qp3)

			Do k=nzb,nze
				Qvl3(k,7)=Qlx3(k,7)
				Qvr3(k,7)=Qrx3(k,7)
			Enddo
			
		    CALL POTOKz(nzp,Dx,Dt,Qvl3,Qvr3,Qp3,Qdx3,Qdy3,Qrp3,Fm3) 

!bug         Do k=k1,k2
         Do k=k1-1,k2

            dnrz(i,j,k)=Qrp3(k,1)		
            vxrz(i,j,k)=Qrp3(k,2)		
            vyrz(i,j,k)=Qrp3(k,3)		
            vzrz(i,j,k)=Qrp3(k,4)		
            bxrz(i,j,k)=Qrp3(k,5)		
            byrz(i,j,k)=Qrp3(k,6)		

            H1(i,j,k)=Fm3(k,1)			
            H2(i,j,k)=Fm3(k,2)			
            H3(i,j,k)=Fm3(k,3)			
            H4(i,j,k)=Fm3(k,4)			
            H5(i,j,k)=Fm3(k,5)			
            H6(i,j,k)=Fm3(k,6)			
            H7(i,j,k)=Fm3(k,7)

         Enddo				
			
	 Enddo
	 Enddo   

!---------------------------------------------------------------	 
!
!     Gardiner&Stone CT 
!
!-----------------------------------------------------------------
      epsk=1.0e-12

      DO K=k1-1,k2
      DO J=j1-1,j2
      DO I=i1,i2

          v2k=0.5*(vz(i,j,k)+vz(i,j,k+1))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vy(i,j,k)*bz(i,j,k)-vz(i,j,k)*by(i,j,k))

		  dexy=-g7(i,j,k)-exc

           else 

	      exc=-(vy(i,j,k+1)*bz(i,j,k+1)-vz(i,j,k+1)*by(i,j,k+1))

		  dexy=-g7(i,j,k+1)-exc

           endif

          else

	      exc=-(vy(i,j,k)*bz(i,j,k)-vz(i,j,k)*by(i,j,k))

		  dexy1=-g7(i,j,k)-exc

	      exc=-(vy(i,j,k+1)*bz(i,j,k+1)-vz(i,j,k+1)*by(i,j,k+1))

		  dexy2=-g7(i,j,k+1)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exi1=h6(i,j,k)+dexy  



          v2k=0.5*(vz(i,j+1,k)+vz(i,j+1,k+1))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vy(i,j+1,k)*bz(i,j+1,k)-vz(i,j+1,k)*by(i,j+1,k))
		  dexy=-g7(i,j,k)-exc

           else   
		

	      exc=-(vy(i,j+1,k+1)*bz(i,j+1,k+1)
     $             -vz(i,j+1,k+1)*by(i,j+1,k+1))
		  dexy=-g7(i,j,k+1)-exc

           endif
		else

	      exc=-(vy(i,j+1,k)*bz(i,j+1,k)-vz(i,j+1,k)*by(i,j+1,k))         
		  dexy1=-g7(i,j,k)-exc

	      exc=-(vy(i,j+1,k+1)*bz(i,j+1,k+1)
     $                 -vz(i,j+1,k+1)*by(i,j+1,k+1))
		  dexy2=-g7(i,j,k+1)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exi2=h6(i,j+1,k)+dexy  




          v2k=0.5*(vy(i,j,k)+vy(i,j+1,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vy(i,j,k)*bz(i,j,k)-vz(i,j,k)*by(i,j,k))
		  dexy=h6(i,j,k)-exc

           else 

	      exc=-(vy(i,j+1,k)*bz(i,j+1,k)-vz(i,j+1,k)*by(i,j+1,k))
		  dexy=h6(i,j+1,k)-exc

           endif

		else

	      exc=-(vy(i,j,k)*bz(i,j,k)-vz(i,j,k)*by(i,j,k))
		  dexy1=h6(i,j,k)-exc

	      exc=-(vy(i,j+1,k)*bz(i,j+1,k)-vz(i,j+1,k)*by(i,j+1,k))
		  dexy2=h6(i,j+1,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exk1=-g7(i,j,k)+dexy



		v2k=0.5*(vy(i,j,k+1)+vy(i,j+1,k+1))
	  
	  if(abs(v2k).gt.epsk)then
	   
	    if(v2k.gt.0.)then

	      exc=-(vy(i,j,k+1)*bz(i,j,k+1)-vz(i,j,k+1)*by(i,j,k+1))

		  dexy=h6(i,j,k)-exc
      
	    else

	      exc=-(vy(i,j+1,k+1)*bz(i,j+1,k+1)
     $              -vz(i,j+1,k+1)*by(i,j+1,k+1))

		  dexy=h6(i,j+1,k)-exc

            endif


		else

	      exc=-(vy(i,j,k+1)*bz(i,j,k+1)-vz(i,j,k+1)*by(i,j,k+1))

		  dexy1=h6(i,j,k)-exc

	      exc=-(vy(i,j+1,k+1)*bz(i,j+1,k+1)
     $                 -vz(i,j+1,k+1)*by(i,j+1,k+1))

		  dexy2=h6(i,j+1,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)


	  endif	 

	    exk2=-g7(i,j,k+1)+dexy  

         ex(I,J,K)=0.25*(exi1+exi2+exk1+exk2)

      ENDDO
      ENDDO
      ENDDO

      DO K=k1-1,k2
      DO J=j1,j2
      DO I=i1-1,i2


        v2k=0.5*(vx(i,j,k)+vx(i+1,j,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vz(i,j,k)*bx(i,j,k)-vx(i,j,k)*bz(i,j,k))

		  dexy=-h5(i,j,k)-exc

           else
	    
	      exc=-(vz(i+1,j,k)*bx(i+1,j,k)-vx(i+1,j,k)*bz(i+1,j,k))

		  dexy=-h5(i+1,j,k)-exc

           endif


	  else

	      exc=-(vz(i,j,k)*bx(i,j,k)-vx(i,j,k)*bz(i,j,k))

		  dexy1=-h5(i,j,k)-exc

	      exc=-(vz(i+1,j,k)*bx(i+1,j,k)-vx(i+1,j,k)*bz(i+1,j,k))

		  dexy2=-h5(i+1,j,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	   endif

	    exi1=f7(i,j,k)+dexy  




        v2k=0.5*(vx(i,j,k+1)+vx(i+1,j,k+1))
	  
	  if(abs(v2k).gt.epsk)then
	   
	  if(v2k.gt.0.)then

	      exc=-(vz(i,j,k+1)*bx(i,j,k+1)-vx(i,j,k+1)*bz(i,j,k+1))

		  dexy=-h5(i,j,k)-exc
 
           else 

	      exc=-(vz(i+1,j,k+1)*bx(i+1,j,k+1)
     $             -vx(i+1,j,k+1)*bz(i+1,j,k+1))

		  dexy=-h5(i+1,j,k)-exc

           endif

	  else

	      exc=-(vz(i,j,k+1)*bx(i,j,k+1)-vx(i,j,k+1)*bz(i,j,k+1))

		  dexy1=-h5(i,j,k)-exc

	      exc=-(vz(i+1,j,k+1)*bx(i+1,j,k+1)
     $                 -vx(i+1,j,k+1)*bz(i+1,j,k+1))

		  dexy2=-h5(i+1,j,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exi2=f7(i,j,k+1)+dexy  





        v2k=0.5*(vz(i,j,k)+vz(i,j,k+1))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then
	  
	      exc=-(vz(i,j,k)*bx(i,j,k)-vx(i,j,k)*bz(i,j,k))

		  dexy=f7(i,j,k)-exc

            else

	      exc=-(vz(i,j,k+1)*bx(i,j,k+1)
     $              -vx(i,j,k+1)*bz(i,j,k+1))

		  dexy=f7(i,j,k+1)-exc

            endif

          else

	      exc=-(vz(i,j,k)*bx(i,j,k)-vx(i,j,k)*bz(i,j,k))

		  dexy1=f7(i,j,k)-exc

	      exc=-(vz(i,j,k+1)*bx(i,j,k+1)
     $                 -vx(i,j,k+1)*bz(i,j,k+1))

		  dexy2=f7(i,j,k+1)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exk1=-h5(i,j,k)+dexy  




        v2k=0.5*(vz(i+1,j,k)+vz(i+1,j,k+1))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vz(i+1,j,k)*bx(i+1,j,k)-vx(i+1,j,k)*bz(i+1,j,k))

		  dexy=f7(i,j,k)-exc

           else

	      exc=-(vz(i+1,j,k+1)*bx(i+1,j,k+1)
     $             -vx(i+1,j,k+1)*bz(i+1,j,k+1))

		  dexy=f7(i,j,k+1)-exc

           endif

          else

	      exc=-(vz(i+1,j,k)*bx(i+1,j,k)-vx(i+1,j,k)*bz(i+1,j,k))

		  dexy1=f7(i,j,k)-exc

	      exc=-(vz(i+1,j,k+1)*bx(i+1,j,k+1)
     $                 -vx(i+1,j,k+1)*bz(i+1,j,k+1))

		  dexy2=f7(i,j,k+1)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exk2=-h5(i+1,j,k)+dexy  


        ey(I,J,K)=0.25*(exi1+exi2+exk1+exk2)

      ENDDO
      ENDDO
      ENDDO

	  
	  
      DO K=k1,k2
      DO J=j1-1,j2
      DO I=i1-1,i2

          v2k=0.5*(vy(i,j,k)+vy(i,j+1,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vx(i,j,k)*by(i,j,k)-vy(i,j,k)*bx(i,j,k))

		  dexy=-f6(i,j,k)-exc

           else

	      exc=-(vx(i,j+1,k)*by(i,j+1,k)-vy(i,j+1,k)*bx(i,j+1,k))

		  dexy=-f6(i,j+1,k)-exc

           endif

          else

	      exc=-(vx(i,j,k)*by(i,j,k)-vy(i,j,k)*bx(i,j,k))

		  dexy1=-f6(i,j,k)-exc

	      exc=-(vx(i,j+1,k)*by(i,j+1,k)-vy(i,j+1,k)*bx(i,j+1,k))

		  dexy2=-f6(i,j+1,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	     exi1=g5(i,j,k)+dexy  

		 


        v2k=0.5*(vy(i+1,j,k)+vy(i+1,j+1,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vx(i+1,j,k)*by(i+1,j,k)-vy(i+1,j,k)*bx(i+1,j,k))

		  dexy=-f6(i,j,k)-exc

            else   

	      exc=-(vx(i+1,j+1,k)*by(i+1,j+1,k)
     $              -vy(i+1,j+1,k)*bx(i+1,j+1,k))

		  dexy=-f6(i,j+1,k)-exc

            endif

	  else

	      exc=-(vx(i+1,j,k)*by(i+1,j,k)-vy(i+1,j,k)*bx(i+1,j,k))

		  dexy1=-f6(i,j,k)-exc

	      exc=-(vx(i+1,j+1,k)*by(i+1,j+1,k)
     $                 -vy(i+1,j+1,k)*bx(i+1,j+1,k))

		  dexy2=-f6(i,j+1,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	     exi2=g5(i+1,j,k)+dexy 


		 
		 
		 v2k=0.5*(vx(i,j,k)+vx(i+1,j,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vx(i,j,k)*by(i,j,k)-vy(i,j,k)*bx(i,j,k))

		  dexy=g5(i,j,k)-exc

           else

	      exc=-(vx(i+1,j,k)*by(i+1,j,k)-vy(i+1,j,k)*bx(i+1,j,k))

		  dexy=g5(i+1,j,k)-exc

	   endif	

	  else

	      exc=-(vx(i,j,k)*by(i,j,k)-vy(i,j,k)*bx(i,j,k))

		  dexy1=g5(i,j,k)-exc

	      exc=-(vx(i+1,j,k)*by(i+1,j,k)-vy(i+1,j,k)*bx(i+1,j,k))

		  dexy2=g5(i+1,j,k)-exc

	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exk1=-f6(i,j,k)+dexy 

		
        v2k=0.5*(vx(i,j+1,k)+vx(i+1,j+1,k))
	  
	  if(abs(v2k).gt.epsk)then

	   if(v2k.gt.0.)then

	      exc=-(vx(i,j+1,k)*by(i,j+1,k)-vy(i,j+1,k)*bx(i,j+1,k))

		  dexy=g5(i,j,k)-exc

           else

	      exc=-(vx(i+1,j+1,k)*by(i+1,j+1,k)
     $             -vy(i+1,j+1,k)*bx(i+1,j+1,k))

		  dexy=g5(i+1,j,k)-exc

	   endif	  
	 
	  else

	      exc=-(vx(i,j+1,k)*by(i,j+1,k)-vy(i,j+1,k)*bx(i,j+1,k))

		  dexy1=g5(i,j,k)-exc

                  exc=-(vx(i+1,j+1,k)*by(i+1,j+1,k)-
     $                 vy(i+1,j+1,k)*bx(i+1,j+1,k))

		  dexy2=g5(i+1,j,k)-exc
		  
	  	  dexy=0.5*(dexy1+dexy2)

	  endif	 

	    exk2=-f6(i,j+1,k)+dexy  

        ez(I,J,K)=0.25*(exi1+exi2+exk1+exk2)

      ENDDO
      ENDDO
      ENDDO

c-------------------------------------------------------------------
c	  
c     Evolution of all variables
c
c-------------------------------------------------------------------

	  Do k=k1,k2
       Do j=j1,j2
        Do i=i1-1,i2

!bug            bxr(i,j,k)=bxr(I,J,K)-Dt/Dx*
            bxrx(i,j,k)=bxrx(I,J,K)-Dt/Dx*
     1		  (ez(I,J,K)-ez(I,J-1,K)+
     1                              ey(I,J,K-1)-ey(I,J,K))

	    Enddo
       Enddo
      Enddo

	  Do k=k1,k2
       Do j=j1-1,j2
        Do i=i1,i2

!bug            byr(i,j,k)=byr(I,J,K)-Dt/Dx*
            byry(i,j,k)=byry(I,J,K)-Dt/Dy*
     1		  (ex(I,J,K)-ex(I,J,K-1)+
     1                              ez(I-1,J,K)-ez(I,J,K))

	    Enddo
       Enddo
      Enddo

	  Do k=k1-1,k2
       Do j=j1,j2
        Do i=i1,i2

!bug            bzr(i,j,k)=bzr(I,J,K)-Dt/Dx*
            bzrz(i,j,k)=bzrz(I,J,K)-Dt/Dz*
     1		  (ex(I,J-1,K)-ex(I,J,K)+
     1                              ey(I,J,K)-ey(I-1,J,K))

	    Enddo
       Enddo
      Enddo
	  


c	    open(10,file='res.dat',access='append')
c        write(10,*) i1,i2,j1,j2,k1,k2
c        do i=i1,i2
c          j=(j1+j2)/2 
c          k=(k1+k2)/2 
c          print*,i,f1(i,j,k)
c          write(10,*) i,f1(i,j,k)
c        enddo
c        close(10)

      Call Conserv(nxp,nyp,nzp,dn,vx,vy,vz,bx,by,bz,
     1             qu1,qu2,qu3,qu4,qu5,qu6,qu7) 	  

	  Do k=k1,k2
       Do j=j1,j2
        Do i=i1,i2
		
		  QH(1)=QU1(I,J,K)-Dt/Dx*(F1(I,J,K)-F1(I-1,J,K)
     1	  +G1(I,J,K)-G1(I,J-1,K)+H1(I,J,K)-H1(I,J,K-1))

	      QH(2)=QU2(I,J,K)-Dt/Dx*(F2(I,J,K)-F2(I-1,J,K)
     1	  +G2(I,J,K)-G2(I,J-1,K)+H2(I,J,K)-H2(I,J,K-1))

	      QH(3)=QU3(I,J,K)-Dt/Dx*(F3(I,J,K)-F3(I-1,J,K)
     1	  +G3(I,J,K)-G3(I,J-1,K)+H3(I,J,K)-H3(I,J,K-1))

	      QH(4)=QU4(I,J,K)-Dt/Dx*(F4(I,J,K)-F4(I-1,J,K)
     1	  +G4(I,J,K)-G4(I,J-1,K)+H4(I,J,K)-H4(I,J,K-1))

  	      QH(5)=QU5(I,J,K)-Dt/Dx*(F5(I,J,K)-F5(I-1,J,K)
     1	  +G5(I,J,K)-G5(I,J-1,K)+H5(I,J,K)-H5(I,J,K-1))

	      QH(6)=QU6(I,J,K)-Dt/Dx*(F6(I,J,K)-F6(I-1,J,K)
     1	  +G6(I,J,K)-G6(I,J-1,K)+H6(I,J,K)-H6(I,J,K-1))

	      QH(7)=QU7(I,J,K)-Dt/Dx*(F7(I,J,K)-F7(I-1,J,K)
     1	  +G7(I,J,K)-G7(I,J-1,K)+H7(I,J,K)-H7(I,J,K-1))

!bug          qh(5)=0.5*(bxr(i,j,k)+bxr(i-1,j,k))  
!bug          qh(6)=0.5*(byr(i,j,k)+byr(i,j-1,k))  
!bug          qh(7)=0.5*(bzr(i,j,k)+bzr(i,j,k-1))  

          qh(5)=0.5*(bxrx(i,j,k)+bxrx(i-1,j,k))  
          qh(6)=0.5*(byry(i,j,k)+byry(i,j-1,k))  
          qh(7)=0.5*(bzrz(i,j,k)+bzrz(i,j,k-1))  

          QU1(I,J,K)=QH(1)
          QU2(I,J,K)=QH(2)
          QU3(I,J,K)=QH(3)
          QU4(I,J,K)=QH(4)
          QU5(I,J,K)=QH(5)
          QU6(I,J,K)=QH(6)
          QU7(I,J,K)=QH(7)

	    Enddo
       Enddo
      Enddo
  
          Call Primitiv(nxp,nyp,nzp,qu1,qu2,qu3,qu4,qu5,qu6,qu7,
     1                 dn,vx,vy,vz,bx,by,bz) 	  

	 
c
c  Check this slice against the list of subgrids 
c     (all subgrid quantities are zero based)
c
	  Do k=k1,k2
      do n=0, nsubgrids-1
	   Do j=j1,j2
	    Do i=i1,i2
		  df(i,j)=f1(i,j,k)
		  vf(i,j)=f2(i,j,k)
		  uf(i,j)=f3(i,j,k)
		  wf(i,j)=f4(i,j,k)
		  bf(i,j)=f5(i,j,k)
		  tf(i,j)=f6(i,j,k)
		  pf(i,j)=f7(i,j,k)
	    enddo
	   enddo	
	  if (k .ge. fjstart(n*3+1)+1 .and. k .le. fjend(n*3+1)+1) then
          idim = fiend(n*3+1) - fistart(n*3+1) + 1
          do j=fistart(n*3+1)+1, fiend(n*3+1)+1
             offset = j-fistart(n*3+1) + (k-fjstart(n*3+1)-1)*idim
             array(dnindex(n*6+1)+offset) = df(lface(n*3+1)+1, j)
             array(dnindex(n*6+2)+offset) = df(rface(n*3+1)+2, j)
             array(vxindex(n*6+1)+offset) = vf(lface(n*3+1)+1, j)
             array(vxindex(n*6+2)+offset) = vf(rface(n*3+1)+2, j)
             array(vyindex(n*6+1)+offset) = uf(lface(n*3+1)+1, j)
             array(vyindex(n*6+2)+offset) = uf(rface(n*3+1)+2, j)
             array(vzindex(n*6+1)+offset) = wf(lface(n*3+1)+1, j)
             array(vzindex(n*6+2)+offset) = wf(rface(n*3+1)+2, j)
             array(bxindex(n*6+1)+offset) = bf(lface(n*3+1)+1, j)
             array(bxindex(n*6+2)+offset) = bf(rface(n*3+1)+2, j)
             array(byindex(n*6+1)+offset) = tf(lface(n*3+1)+1, j)
             array(byindex(n*6+2)+offset) =  tf(rface(n*3+1)+2, j)
             array(bzindex(n*6+1)+offset) = pf(lface(n*3+1)+1, j)
             array(bzindex(n*6+2)+offset) = pf(rface(n*3+1)+2, j)
          enddo
        endif
      enddo
      enddo

      Do i=i1,i2
      do n=0, nsubgrids-1
      Do j=j1,j2
	   Do k=k1,k2
		  dg(j,k)=g1(i,j,k)
		  vg(j,k)=g2(i,j,k)
		  ug(j,k)=g3(i,j,k)
		  wg(j,k)=g4(i,j,k)
		  bg(j,k)=g5(i,j,k)
		  tg(j,k)=g6(i,j,k)
		  pg(j,k)=g7(i,j,k)
	    enddo
	   enddo	
        if (i .ge. fistart(n*3+2)+1 .and. i .le. fiend(n*3+2)+1) then
          idim = fiend(n*3+2) - fistart(n*3+2) + 1
          do k=fjstart(n*3+2)+1, fjend(n*3+2)+1
             offset = i-fistart(n*3+2) + (k-fjstart(n*3+2)-1)*idim
             array(dnindex(n*6+3)+offset) = dg(lface(n*3+2)+1, k)
             array(dnindex(n*6+4)+offset) = dg(rface(n*3+2)+2, k)
             array(vxindex(n*6+3)+offset) = vg(lface(n*3+2)+1, k)
             array(vxindex(n*6+4)+offset) = vg(rface(n*3+2)+2, k)
             array(vyindex(n*6+3)+offset) = ug(lface(n*3+2)+1, k)
             array(vyindex(n*6+4)+offset) = ug(rface(n*3+2)+2, k)
             array(vzindex(n*6+3)+offset) = wg(lface(n*3+2)+1, k)
             array(vzindex(n*6+4)+offset) = wg(rface(n*3+2)+2, k)
             array(bxindex(n*6+3)+offset) = bg(lface(n*3+2)+1, k)
             array(bxindex(n*6+4)+offset) = bg(rface(n*3+2)+2, k)
             array(byindex(n*6+3)+offset) = tg(lface(n*3+2)+1, k)
             array(byindex(n*6+4)+offset) = tg(rface(n*3+2)+2, k)
             array(bzindex(n*6+3)+offset) = pg(lface(n*3+2)+1, k)
             array(bzindex(n*6+4)+offset) = pg(rface(n*3+2)+2, k)
          enddo
        endif
      enddo
	  Enddo
	  
      Do j=j1,j2
      do n=0, nsubgrids-1
      Do i=i1,i2
	   Do k=k1,k2
		  dh(k,i)=h1(i,j,k)
		  vh(k,i)=h2(i,j,k)
		  uh(k,i)=h3(i,j,k)
		  wh(k,i)=h4(i,j,k)
		  bh(k,i)=h5(i,j,k)
		  th(k,i)=h6(i,j,k)
		  ph(k,i)=h7(i,j,k)
	    enddo
	   enddo	        
	    if (j .ge. fjstart(n*3+3)+1 .and. j .le. fjend(n*3+3)+1) then
          idim = fiend(n*3+3) - fistart(n*3+3) + 1
          do i=fistart(n*3+3)+1, fiend(n*3+3)+1
             offset = i-fistart(n*3+3) + (j-fjstart(n*3+3)-1)*idim
             array(dnindex(n*6+5)+offset) = dh(lface(n*3+3)+1, i)
             array(dnindex(n*6+6)+offset) = dh(rface(n*3+3)+2, i)
             array(vxindex(n*6+5)+offset) = vh(lface(n*3+3)+1, i)
             array(vxindex(n*6+6)+offset) = vh(rface(n*3+3)+2, i)
             array(vyindex(n*6+5)+offset) = uh(lface(n*3+3)+1, i)
             array(vyindex(n*6+6)+offset) = uh(rface(n*3+3)+2, i)
             array(vzindex(n*6+5)+offset) = wh(lface(n*3+3)+1, i)
             array(vzindex(n*6+6)+offset) = wh(rface(n*3+3)+2, i)
             array(bxindex(n*6+5)+offset) = bh(lface(n*3+3)+1, i)
             array(bxindex(n*6+6)+offset) = bh(rface(n*3+3)+2, i)
             array(byindex(n*6+5)+offset) = th(lface(n*3+3)+1, i)
             array(byindex(n*6+6)+offset) = th(rface(n*3+3)+2, i)
             array(bzindex(n*6+5)+offset) = ph(lface(n*3+3)+1, i)
             array(bzindex(n*6+6)+offset) = ph(rface(n*3+3)+2, i)
          enddo
        endif
      enddo	  
      enddo	  
	  
      Return
	End 
