!
!     C-callable interface to retrieve UB31 station stress/force data
!     and Beam Code Check results.
!
      subroutine getub31stx(nelem, out_array)
      use ub31_module
      implicit none
      integer nelem
      real*8  out_array(6,11)
      integer j, k

      if (allocated(ub31_stx)) then
         do j = 1, 11
            do k = 1, 6
               out_array(k,j) = ub31_stx(k, j, nelem)
            enddo
         enddo
      else
         do j = 1, 11
            do k = 1, 6
               out_array(k,j) = 0.d0
            enddo
         enddo
      endif
      return
      end

      subroutine getub31for(nelem, out_array)
      use ub31_module
      implicit none
      integer nelem
      real*8  out_array(6,11)
      integer j, k

      if (allocated(ub31_for)) then
         do j = 1, 11
            do k = 1, 6
               out_array(k,j) = ub31_for(k, j, nelem)
            enddo
         enddo
      else
         do j = 1, 11
            do k = 1, 6
               out_array(k,j) = 0.d0
            enddo
         enddo
      endif
      return
      end

      subroutine iscodecheckactive(iactive)
      use codecheck_module
      implicit none
      integer iactive

      if (codecheck_active) then
         iactive = 1
      else
         iactive = 0
      endif
      return
      end

      subroutine getcodecheckuc(nelem, out_array)
      use codecheck_module
      implicit none
      integer nelem
      real*8 out_array(6)

      if (allocated(uc_tot) .and. nelem .le. size(uc_tot)) then
         out_array(1) = uc_tot(nelem)
         out_array(2) = uc_ax(nelem)
         out_array(3) = uc_sh(nelem)
         out_array(4) = uc_bnd(nelem)
         out_array(5) = uc_stab(nelem)
         out_array(6) = uc_ltb(nelem)
      else
         out_array = 0.d0
      endif
      return
      end

      subroutine getcodecheckstationuc(nelem, istation, out_array)
      use codecheck_module
      implicit none
      integer nelem, istation
      real*8 out_array(6)
      integer k

      if (allocated(station_uc)) then
         if (nelem .le. size(station_uc, 3) .and.
     &       istation .ge. 1 .and. istation .le. 11) then
            do k = 1, 6
               out_array(k) = station_uc(k, istation, nelem)
            enddo
         else
            out_array = 0.d0
         endif
      else
         out_array = 0.d0
      endif
      return
      end
