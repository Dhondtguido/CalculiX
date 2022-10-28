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
      subroutine writedesi(norien,orname)
!
!     writes the orientation design variables in the .dat file
!
      implicit none
!
      character*80 orname(*)
!
      integer j,norien
!
!     
      write(5,*)
      write(5,*) '    D E S I G N   V A R I A B L E S'
      write(5,*)
      write(5,'(a8,1x,a11,62x,a15)') 'VARIABLE','ORIENTATION',
     &    'ROTATION VECTOR'
      write(5,*)
!
      do j=1,norien
         write(5,'(i5,4x,a80,1x,a5)') (j-1)*3+1,orname(j)(1:80),'Rx   '
         write(5,'(i5,4x,a80,1x,a5)') (j-1)*3+2,orname(j)(1:80),'Ry   '
         write(5,'(i5,4x,a80,1x,a5)') (j-1)*3+3,orname(j)(1:80),'Rz   '
      enddo
!
      return
      end

