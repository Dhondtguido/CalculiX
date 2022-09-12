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
      subroutine negativepressure(ne0,ne,mi,stx,pressureratio)
!
!     calculating the ratio of the smallest pressure to the
!     largest pressure for face-to-face contact
!     if the pressure is somewhere negative, this ratio will
!     be negative
!
      implicit none
!
      integer ne0,ne,mi(*),i
!
      real*8 stx(6,mi(1),*),presmin,presmax,pressureratio
!
      presmax=0.d0
      presmin=0.d0
!
      do i=ne0+1,ne
         if(stx(4,1,i).gt.presmax) then
            presmax=stx(4,1,i)
         elseif(stx(4,1,i).lt.presmin) then
            presmin=stx(4,1,i)
         endif
      enddo
      pressureratio=presmin/presmax
!
      return
      end
