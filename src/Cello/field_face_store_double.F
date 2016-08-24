      subroutine field_face_store_double (field_ghost,array,nd3,n3,im3)
      implicit none
      integer nd3(3),n3(3),im3(3)
      real*8 field_ghost(nd3(1),nd3(2),nd3(3))
      real*8 array(n3(1),n3(2),n3(3))
      integer ix,iy,iz
      
      do iz=1,n3(3)
         do iy=1,n3(2)
            do ix=1,n3(1)
               field_ghost(ix,iy,iz) = array(ix,iy,iz)
            end do
         end do
      end do

      return
      end


