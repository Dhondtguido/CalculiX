!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine spooles_write(ad,au,adb,aub,sigma,b,icol,irow,neq,
     &     nzs,symmetryflag,inputformat,nzs3)
!
      implicit none
!
      integer icol(*),irow(*),neq,nzs,symmetryflag,inputformat,nzs3,i
!
      real*8 ad(*),au(*),adb(*),aub(*),sigma,b(*)
!
      open(18,file='spooles_matrix',status='unknown')
!
      write(18,101) neq,nzs,symmetryflag,inputformat,nzs3
      write(18,100) sigma
      do i=1,neq
        write(18,100) ad(i)
      enddo
      do i=1,nzs
        write(18,100) au(i)
      enddo
      do i=1,neq
        write(18,100) b(i)
      enddo
      do i=1,neq
        write(18,101) icol(i)
      enddo
      do i=1,nzs
        write(18,101) irow(i)
      enddo
 100  format(e20.13)
 101  format(5i10)
      close(18)
      return
      end

