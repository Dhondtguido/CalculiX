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

