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
      subroutine intersectionpoint(pa,pb,xcp,t,xinters)
!     
      implicit none
!
      integer k
!     
      real*8 pa(*),pb(*),xcp(*),t,xinters(*),diff,pab(3),
     &     eplane,tnull
!
      do k=1,3
         pab(k)=pb(k)-pa(k)
       enddo
!       
      diff=0.d0
      tnull=0.d0
!      
      if(abs(eplane(pab,xcp,tnull)).lt.1.d-13)then
         write(*,*) 'SH: IP no intersection point can be found'
         write(*,*) 'SH: IP pab paralell to plane! '
         call exit(201)
      else 
         diff=-eplane(pa,xcp,t)/eplane(pab,xcp,tnull)
      endif  
      do k=1,3
         xinters(k)=pa(k)+diff*pab(k)
      enddo
      return
      end       
      
