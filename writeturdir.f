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
      subroutine writeturdir(xn,turdir,nev)
!
!     writes the eigenvalues in the .dat file
!
      implicit none
!
      character*1 turdir(*)
      integer j,nev
      real*8 xn(3)
!
      write(5,*)
      write(5,*)
     &  '    E I G E N M O D E   T U R N I N G   D I R E C T I O N'
      write(5,*)
      write(5,100) (xn(j),j=1,3)
 100  format('    Axis reference direction:',3(1x,e11.4))
      write(5,*)
      write(5,*) 'MODE NO     TURNING DIRECTION (F=FORWARD,B=BACKWARD)'
      write(5,*)
      do j=1,nev
         write(5,'(i7,10x,a1)') j,turdir(j)
      enddo
!
      return
      end

