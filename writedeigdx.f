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
      subroutine writedeigdx(iev,d,ndesi,orname,dgdx)
!
!     writes the derivative of the eigenfrequencies w.r.t.
!     the orientations in the .dat file
!
      implicit none
!
      character*5 angle
      character*80 orname(*)
!
      integer idesvar,ndesi,iorien,iangle,iev
!
      real*8 dgdx(ndesi,*),d(*)
!
!     
      write(5,*)
      write(5,*) '    E I G E N V A L U E   S E N S I T I V I T Y'
      write(5,*)
      write(5,'(a10,2x,i5,2x,e11.4)') 'EIGENVALUE',iev+1,d(iev+1)
      write(5,*)
!
      do idesvar=1,ndesi
         iorien=(idesvar-1)/3+1
         iangle=idesvar-((idesvar-1)/3)*3
         if(iangle.eq.1) then
            angle='   Rx'
         elseif(iangle.eq.2) then
            angle='   Ry'
         else
            angle='   Rz'
         endif
         write(5,'(a80,1x,a5,1x,e11.4)') orname(iorien),angle,
     &       dgdx(idesvar,1)
      enddo
!
      return
      end

