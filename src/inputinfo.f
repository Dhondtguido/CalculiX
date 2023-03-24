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
      subroutine inputinfo(inpc,ipoinpc,iline,text)
!
!     input error message subroutine
!
      implicit none
!
      character*1 inpc(*)
      character*(*) text
!
      integer ipoinpc(0:*),iline,i
!
      write(*,*) '*INFO reading ',text(1:index(text,'%')-1-1),
     &      '. Card image:'
      write(*,'(7x,1320a1)') 
     &     (inpc(i),i=ipoinpc(iline-1)+1,ipoinpc(iline))
      write(*,*)
!
      return
      end
