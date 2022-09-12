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
      subroutine reorderampl(amname,namta,nam)
!
!     reorders amname in alphabetical order
!
      implicit none
!
      character*80 amname(*),amnamecp
!
      integer namta(3,*),nam,id,namtacp(3),i,j
!
      call cident80(amname,amname(nam),nam-1,id)
!
      amnamecp=amname(nam)
      do i=1,3
         namtacp(i)=namta(i,nam)
      enddo
!
      do j=nam,id+2,-1
         amname(j)=amname(j-1)
         do i=1,3
            namta(i,j)=namta(i,j-1)
         enddo
      enddo
!
      amname(id+1)=amnamecp
      do i=1,3
         namta(i,id+1)=namtacp(i)
      enddo
!
      return
      end
