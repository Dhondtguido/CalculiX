!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine init_submodel(nktet,inodfa,ipofa,netet_);
!
      implicit none
!
      integer nktet,i,netet_,inodfa(4,*),ipofa(*)
!
!     initialization of ipofa and inodfa
!
      do i=1,nktet
         ipofa(i)=0
      enddo
      do i=1,4*netet_
         inodfa(4,i)=i+1
      enddo
      inodfa(4,4*netet_)=0
!
      return
      end
