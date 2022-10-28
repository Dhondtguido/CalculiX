!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine writevector(ad,neq,number)
!
!      writes an vector to file (for debugging purposes)
!
      implicit none
!      
      character*12 name
      character*14 name2
!
      integer neq,i,number
!
      real*8 ad(*)
!
      name='vector_'//char(number+96)//'.out'
      name2='vector_'//char(number+96)//'_t.out'
      open(10,file=name,status='unknown')
      write(10,*) 'vector  number ',number
!      
      do i=1,neq
c         if(ad(i).gt.1.0e-10.or.ad(i).lt.-1.0e-10)then
            write(10,*) 'row ',i,' value ',ad(i)
c         endif
      enddo
!
      close(10)
      return
      end

