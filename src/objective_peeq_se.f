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
      subroutine objective_peeq_se(nk,iobject,mi,depn,objectset,
     &  ialnneigh,naneigh,nbneigh,epn,dksper)
!
!     calculates the 
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer nk,idir,iobject,mi(*),j,k,ialnneigh(*),naneigh,nbneigh
!
      real*8 epn(*),p,rho,xpeeq,dksper,depn(*),mises
!
!     reading rho and the mean peeq for the Kreisselmeier-Steinhauser
!     function
!
      read(objectset(2,iobject)(41:60),'(f20.0)') rho
      read(objectset(2,iobject)(61:80),'(f20.0)') xpeeq
!
      dksper=0.d0
      do j=naneigh,nbneigh        
         k=ialnneigh(j)
!
         dksper=dksper+dexp(rho*depn(k)/xpeeq)*(depn(k)-epn(k))
      enddo
!
      dksper=dksper/xpeeq
!
      return
      end

