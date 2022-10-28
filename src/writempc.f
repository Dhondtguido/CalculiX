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
      subroutine writempc(ipompc,nodempc,coefmpc,labmpc,mpc)
!
!     writes an MPC to standard output (for debugging purposes)
!
      implicit none
!
      character*20 labmpc(*)
      integer ipompc(*),nodempc(3,*),mpc,index,node,idir
      real*8 coefmpc(*),coef
!
      write(*,*)
      write(*,'(''MPC '',i10,1x,a20)') mpc,labmpc(mpc)
      index=ipompc(mpc)
      do
         node=nodempc(1,index)
         idir=nodempc(2,index)
         coef=coefmpc(index)
         write(*,'(i10,1x,i5,1x,e11.4)') node,idir,coef
         index=nodempc(3,index)
         if(index.eq.0) exit
      enddo
!
      return
      end

