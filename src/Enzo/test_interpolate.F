      program test_interpolate
      
      implicit none
      integer, parameter :: ndim=2, imethod=1, iposflag=0, ierror=0
      integer, parameter :: mxc=4, myc=8, mzc=1
      integer, parameter :: mxf=4, myf=12, mzf=1
      integer, parameter :: nxc=4, nyc=8, nzc=1
      integer, parameter :: nxf=4, nyf=12, nzf=1
      integer pdims(ndim),pstart(ndim),pend(ndim),refine(ndim)
      integer gdims(ndim),gstart(ndim)
      real parent(mxc,myc), grid(mxf,myf),work(1000)
      integer size,rank
      integer ix,iy,iz
      
      size = 0

      pdims(1) = mxc
      pdims(2) = myc

      pstart(1) = 2
      pstart(2) = 2

      pend(1) = 5
      pend(2) = 13

      refine(1) = 2
      refine(2) = 2

      do ix=1,mxc
         do iy=1,myc
            parent(ix,iy) = ((ix-1)+mxc*(iy-1))
         end do
      end do

      do ix=1,nxc
         print 10, (parent(ix,iy), iy=1,nyc)
      end do
 10   format (1x, 10E10.2)

      gdims(1) = mxf
      gdims(2) = myf

      gstart(1) = 0
      gstart(2) = 0

      do ix=1,mxf
         do iy=1,myf
            do iz=1,mzf
               grid(ix,iy) = (ix-1)+mxf*(iy-1)
            end do
         end do
      end do

      size = size + gdims(1)/2+1
      size = size + gdims(2)/2+1

      rank = 2
      print*, 'pdims',pdims(1),pdims(2)
      print*, 'pstart',pstart(1),pstart(2)
      print*, 'pend',pend(1),pend(2)
      print*, 'gdims',gdims(1),gdims(2)
      print*, 'gstart',gstart(1),gstart(2)
      call interpolate
     +     (rank,parent,pdims,pstart,pend,refine,grid,
     +     gdims,gstart,work,imethod,iposflag,ierror)

      do ix=1,nxf
         print 10, (grid(ix,iy), iy=1,nyf)
      end do
      return
      end
C     2.00  3.00  4.00  5.00  6.00  7.00
C     3.00  4.00  5.00  6.00  7.00  8.00
C     4.00  5.00  6.00  7.00  8.00  9.00
C     5.00  6.00  7.00  8.00  9.00 10.00
C     6.00  7.00  8.00  9.00 10.00 11.00
C     7.00  8.00  9.00 10.00 11.00 12.00
C     
C     3.50  4.00  4.50  5.00  5.50  6.00
C     4.00  4.50  5.00  5.50  6.00  6.50
C     4.50  5.00  5.50  6.00  6.50  7.00
C     5.00  5.50  6.00  6.50  7.00  7.50
C     5.50  6.00  6.50  7.00  7.50  8.00
C     6.00  6.50  7.00  7.50  8.00  8.50
