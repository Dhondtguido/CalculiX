!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine relaxval_al(r,gmatrix,nacti)
!     
      implicit none
!     
      integer i,j,nacti,kcol(3)
!     
      real*8 dssd2,r(*),vssd1(3),gmatrix(nacti,nacti)
!
      logical isdiagdomi
!
      isdiagdomi=.true.
!
!     check for diagonal dominance
!
      do i=1,nacti
        dssd2=sum(abs(gmatrix(:,i))) ! TODO CMT check if rowwise or columnwise is faster (G is symmetric)
        if(2.0d0*abs(gmatrix(i,i))<dssd2) then
          isdiagdomi=.false.
          exit
        endif
      enddo
!
      do i=1,(nacti/3)
!
!       correct column
!
        do j=1,3
          kcol(j)=3*(i-1)+j
        enddo
!        
        if(isdiagdomi) then
          do j=1,3
            vssd1(j)=gmatrix(kcol(j),kcol(j))
          enddo
        else
          do j=1,3
            vssd1(j)=sum(abs(gmatrix(:,kcol(j))))
          enddo
        endif
!        
        dssd2=1.0d0/maxval(vssd1)
!        
        do j=1,3
          r(kcol(j))=dssd2
        enddo
      enddo
!
      return
      end
