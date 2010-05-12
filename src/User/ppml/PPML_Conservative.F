c     $Id$
c     See LICENSE_ENZO file for license and copyright information

      Subroutine Conserv(nx,ny,nz,dn,vx,vy,vz,bx,by,bz,
     &                   qu1,qu2,qu3,qu4,qu5,qu6,qu7)
         Implicit NONE

      Integer nx,ny,nz,i,j,k
          REAL dn(nx,ny,nz)           
      REAL vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)           
      REAL bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)           

      REAL qu1(nx,ny,nz)           
      REAL qu2(nx,ny,nz),qu3(nx,ny,nz),qu4(nx,ny,nz)           
      REAL qu5(nx,ny,nz),qu6(nx,ny,nz),qu7(nx,ny,nz)           
          
      Do k=1,nz
       Do j=1,ny
        Do i=1,nx
             QU1(i,j,k)=dn(i,j,k)
             QU2(i,j,k)=dn(i,j,k)*vx(i,j,k)
             QU3(i,j,k)=dn(i,j,k)*vy(i,j,k)
             QU4(i,j,k)=dn(i,j,k)*vz(i,j,k)
             QU5(i,j,k)=bx(i,j,k)
             QU6(i,j,k)=by(i,j,k)
             QU7(i,j,k)=bz(i,j,k)
                Enddo
           Enddo        
      Enddo
           
         RETURN
         END
