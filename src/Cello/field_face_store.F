      subroutine field_face_store_4
     +     (field_ghost,array,nd3,n3,accumulate)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer nd3(3),n3(3),accumulate
      real(kind=real32) field_ghost(nd3(1),nd3(2),nd3(3))
      real(kind=real32) array(n3(1),n3(2),n3(3))
      integer ix,iy,iz
      if (accumulate.eq.0) then
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) = array(ix,iy,iz)
               end do
            end do
         end do
      else
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) =
     +                 field_ghost(ix,iy,iz) + array(ix,iy,iz)
               end do
            end do
         end do
      end if

      return
      end

      subroutine field_face_store_8
     +     (field_ghost,array,nd3,n3,accumulate)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer nd3(3),n3(3),accumulate
      real(kind=real64) field_ghost(nd3(1),nd3(2),nd3(3))
      real(kind=real64) array(n3(1),n3(2),n3(3))
      integer ix,iy,iz
      
      if (accumulate.eq.0) then
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) = array(ix,iy,iz)
               end do
            end do
         end do
      else
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) =
     +                 field_ghost(ix,iy,iz) + array(ix,iy,iz)
               end do
            end do
         end do
      end if

      return
      end

      subroutine field_face_store_16
     +     (field_ghost,array,nd3,n3,accumulate)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer nd3(3),n3(3),accumulate
      real(kind=real128) :: field_ghost(nd3(1),nd3(2),nd3(3))
      real(kind=real128) :: array(n3(1),n3(2),n3(3))
      integer ix,iy,iz
      
      if (accumulate.eq.0) then
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) = array(ix,iy,iz)
               end do
            end do
         end do
      else
         do iz=1,n3(3)
            do iy=1,n3(2)
               do ix=1,n3(1)
                  field_ghost(ix,iy,iz) =
     +                 field_ghost(ix,iy,iz) + array(ix,iy,iz)
               end do
            end do
         end do
      end if

      return
      end


