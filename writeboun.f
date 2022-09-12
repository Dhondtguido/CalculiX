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
      subroutine writeboun(nodeboun,ndirboun,xboun,typeboun,nboun)
!
!     writes an MPC to standard output (for debugging purposes)
!
      implicit none
!
      character*1 typeboun(*)
      integer nodeboun(*),ndirboun(*),nboun,i
      real*8 xboun(*)
!
      write(*,*)
      write(*,'(''SPC '')') 
      do i=1,nboun
         write(*,'(i5,1x,i10,1x,i5,1x,e11.4,1x,a1)') i,nodeboun(i),
     &             ndirboun(i),xboun(i),typeboun(i)
      enddo
!
      return
      end

