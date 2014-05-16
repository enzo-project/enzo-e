#include "fortran.def"

!=======================================================================
!////////////////////  SUBROUTINE COOL1D_CLOUDY_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_cloudy_g(d, de, rhoH, metallicity,
     &                in, jn, kn, is, ie, j, k,
     &                logtem, edot, comp2, ispecies, dom, zr,
     &                icmbTfloor, iClHeat, iZscale,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clCooling, clHeating, 
     &                itmask)

!
!  SOLVE CLOUDY METAL COOLING
!
!  written by: Britton Smith
!  date: September, 2009
!
!  PURPOSE:
!    Solve cloudy cooling by interpolating from the data.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!
!    rhoH     - total H mass density
!    metallicity - metallicity
!
!    is,ie    - start and end indices of active region (zero based)
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    logtem   - natural log of temperature values
!
!    dom      - unit conversion to proper number density in code units
!    zr       - current redshift
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    iZscale    - flag to scale cooling by metallicity
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!
!    itmask     - iteration mask
!
!  OUTPUTS:
!    update edot with heating/cooling contributions from metals
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE

!  General Arguments

      integer in, jn, kn, is, ie, j, k, ispecies

      real comp2, dom, zr
      real d(in,jn,kn), de(in,jn,kn),
     &     rhoH(in), metallicity(in), logtem(in)
      double precision edot(in)

!  Cloudy parameters and data

      integer icmbTfloor, iClHeat, iZscale, 
     &     clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3))
      real clCooling(clDataSize), clHeating(clDataSize)

!  Iteration mask

      logical itmask(in)

!  Parameters

!  Locals

      integer i, q, zindex, zmidpt, zhighpt, get_heat
      real dclPar(clGridRank), inv_log10, log10_tCMB
      logical end_int

!  Slice locals

      real e_frac(in), log_e_frac(in), 
     &     cl_e_frac(in), fh(in), log_n_h(in),
     &     log_cool(in), log_cool_cmb(in), log_heat(in),
     &     edot_met(in), log10tem(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      end_int = .false.
      get_heat = iClHeat

      inv_log10 = 1.d0 / log(10.d0)
      log10_tCMB = log10(comp2)

!     Calculate parameter value slopes

      dclPar(1) = (clPar1(clGridDim(1)) - clPar1(1)) / 
     &     real(clGridDim(1) - 1)
      if (clGridRank .gt. 1) then
         dclPar(2) = (clPar2(clGridDim(2)) - clPar2(1)) / 
     &        real(clGridDim(2) - 1)
      endif
      if (clGridRank .gt. 2) then
         dclPar(3) = (clPar3(clGridDim(3)) - clPar3(1)) / 
     &        real(clGridDim(3) - 1)
      endif

      do i=is+1, ie+1
         if ( itmask(i) ) then

            log10tem(i) = logtem(i) * inv_log10

!           Calcualte H mass fraction

            fh(i) = rhoH(i) / d(i,j,k)

!           Calculate proper log(n_H)

            log_n_h(i) = log10(rhoH(i) * dom)

!           Calculate index for redshift dimension

            if (clGridRank .gt. 2) then

!           Get index for redshift dimension via bisection

               if (zr .le. clPar2(1)) then
                  zindex = 1
               else if (zr .ge. clPar2(clGridDim(2)-1)) then
                  zindex = clGridDim(2)
                  end_int = .true.
                  get_heat = 0
               else if (zr .ge. clPar2(clGridDim(2)-2)) then
                  zindex = clGridDim(2) - 2
               else
                  zindex = 1
                  zhighpt = clGridDim(2) - 2
                  do while ((zhighpt - zindex) .gt. 1)
                     zmidpt = int((zhighpt + zindex) / 2)
                     if (zr .ge. clPar2(zmidpt)) then
                        zindex = zmidpt
                     else
                        zhighpt = zmidpt
                     endif
                  enddo
               endif

            endif

!           Call interpolation functions to get heating/cooling

!           Interpolate over temperature.
            if (clGridRank .eq. 1) then
               call interpolate_1D_g(log10tem(i), clGridDim, clPar1,
     &              dclPar(1), clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_1D_g(log10_tCMB, clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clCooling, 
     &                 log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
                  call interpolate_1D_g(log10tem(i), clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clHeating, 
     &                 log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density and temperature.
            else if (clGridRank .eq. 2) then
               call interpolate_2D_g(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.0)) then
                  call interpolate_2D_g(log_n_h(i), log10_tCMB, 
     &                 clGridDim, clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
               call interpolate_2D_g(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density, redshift, and temperature.
            else if (clGridRank .eq. 3) then
               call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &              clGridDim,
     &              clPar1, dclPar(1), 
     &              clPar2, zindex,
     &              clPar3, dclPar(3),
     &              clDataSize, clCooling, 
     &              end_int, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_3Dz_g(log_n_h(i), zr, log10_tCMB,
     &                 clGridDim,
     &                 clPar1, dclPar(1),
     &                 clPar2, zindex,
     &                 clPar3, dclPar(3),
     &                 clDataSize, clCooling, 
     &                 end_int, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
                  call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &                 clGridDim,
     &                 clPar1, dclPar(1),
     &                 clPar2, zindex,
     &                 clPar3, dclPar(3),
     &                 clDataSize, clHeating, 
     &                 end_int, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

            else
               write(*,*) "Maximum cooling data grid rank is 3!"
               return
            endif

!           Scale cooling by metallicity.

            if (iZscale .eq. 1) then
               edot_met(i) = edot_met(i) * metallicity(i)
            endif

            edot(i) = edot(i) + 
     &           (edot_met(i) * rhoH(i) * d(i,j,k))

         end if
      enddo

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_1D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_1D_g(input1, gridDim, gridPar1, dgridPar1, 
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(1)
      real input1, value
      real gridPar1(gridDim(1)), dgridPar1
      real dataField(dataSize)

!  Locals

      integer index1
      real slope

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation index

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))

!     Interpolate over parameter 1

      slope = (dataField(index1+1) - dataField(index1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + dataField(index1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_2D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_2D_g(input1, input2, gridDim, 
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(2)
      real input1, input2, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2
      real dataField(dataSize)

!  Locals

      integer index1, index2, int_index, q
      real slope, value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))

      do q=1, 2

!     interpolate over parameter 2

         int_index = (q+index1-2) * gridDim(2) + index2

         slope = (dataField(int_index+1) - dataField(int_index)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + 
     &        dataField(int_index)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) / 
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_3D_g(input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(3)
      real input1, input2, input3, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3
      real dataField(dataSize)

!  Locals

      integer index1, index2, index3, int_index, q, w
      real slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))

      do q=1, 2

         do w=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) / 
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end
!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_3Dz  \\\\\\\\\\\\\\\\\\\\\
!
!     Similar to interpolate_3D except index2 is calculated
!     ahead of time because it is the redshift and will not 
!     change for the entire grid.
!
      subroutine interpolate_3Dz_g(input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, index2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, 
     &     end_int, value)

!  General Arguments

      integer dataSize, index2
      integer gridDim(3)
      real input1, input2, input3, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)),
     &     gridPar3(gridDim(3)), dgridPar3
      real dataField(dataSize)

!  Locals

      logical end_int
      integer index1, index3, int_index, q, w
      real slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      if (end_int) then
         call interpolate_2Df3D_g(input1, input2, 
     &        input3, gridDim,
     &        gridPar1, dgridPar1,
     &        gridPar2, index2,
     &        gridPar3, dgridPar3,
     &        dataSize, dataField, 
     &        value)
         return
      endif

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))

      do q=1, 2

         do w=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) / 
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!///////////////////  SUBROUTINE INTERPOLATE_2Df3D  \\\\\\\\\\\\\\\\\\\\
!
!     Interpolation in 2 dimensions but with a 3D grid.
!     This is used for interpolating from just the last 
!     slice in the datacube before the redshift where 
!     the UV background turns on.
!
      subroutine interpolate_2Df3D_g(input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, index2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, 
     &     value)

!  General Arguments

      integer dataSize, index2
      integer gridDim(3)
      real input1, input2, input3, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)),
     &     gridPar3(gridDim(3)), dgridPar3
      real dataField(dataSize)

!  Locals

      logical end_int
      integer index1, index3, int_index, q
      real slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))

      do q=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (index2-1)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(q) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

      enddo

!     interpolate over parameter 1

      slope = (value3(2) - value3(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value3(1)

      return
      end
