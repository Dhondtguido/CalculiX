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
      subroutine elemperorien(ipoorel,iorel,ielorien,ne,mi)
!
      implicit none
!
      integer ipoorel(*),iorel(2,*),i,ne,mi(*),ielorien(mi(3),*),
     &  iorelfree,iorien
!
!
!
!     determining the elements belonging to the nodes of
!     the elements
!
      iorelfree=1
      do i=1,ne
         iorien=max(0,ielorien(1,i))
         if(iorien.eq.0) cycle
         iorel(1,iorelfree)=i
         iorel(2,iorelfree)=ipoorel(iorien)
         ipoorel(iorien)=iorelfree
         iorelfree=iorelfree+1
      enddo
!
      return
      end
