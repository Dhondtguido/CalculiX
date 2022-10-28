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
      subroutine distattachline(xig,etg,pneigh,pnode,dist,
     &     nterms,xn,p)
!
!     calculates the distance between a straight line through the node 
!     with coordinates in "pnode" and direction vector "xn" and 
!     the node with local coordinates xig and etg
!     in a face described by "nterms" nodes with coordinates
!     in "pneigh"
!
      implicit none
!
      integer nterms,i,j
!
      real*8 ratio(8),pneigh(3,*),pnode(3),dist,xi,et,xig,etg,p(3),
     &  xn(3),coeff,etm2,xim2,etm,xim,etp,xip,a2,et2,xi2,a
!
!
!
      if(nterms.eq.3) then
         if(xig+etg.le.0.d0) then
            ratio(2)=(xig+1.d0)/2.d0
            ratio(3)=(etg+1.d0)/2.d0
         else
            ratio(2)=(1.d0-etg)/2.d0
            ratio(3)=(1.d0-xig)/2.d0
         endif
         ratio(1)=1.d0-ratio(2)-ratio(3)
      elseif(nterms.eq.4) then
         xip=(1.d0+xig)/4.d0
         xim=(1.d0-xig)/4.d0
         etp=1.d0+etg
         etm=1.d0-etg
         ratio(1)=xim*etm
         ratio(2)=xip*etm
         ratio(3)=xip*etp
         ratio(4)=xim*etp
      elseif(nterms.eq.6) then
         if(xig+etg.le.0.d0) then
            xi=(xig+1.d0)/2.d0
            et=(etg+1.d0)/2.d0
         else
            xi=(1.d0-etg)/2.d0
            et=(1.d0-xig)/2.d0
         endif
         a=1.d0-xi-et
         a2=2.d0*a
         xi2=2.d0*xi
         et2=2.d0*et
         ratio(1)=a*(a2-1.d0)
         ratio(2)=xi*(xi2-1.d0)
         ratio(3)=et*(et2-1.d0)
         ratio(4)=xi2*a2
         ratio(5)=xi2*et2
         ratio(6)=et2*a2
      elseif(nterms.eq.8) then
         xip=1.d0+xig
         xim=1.d0-xig
         xim2=xip*xim/2.d0
         etp=1.d0+etg
         etm=1.d0-etg
         etm2=etp*etm/2.d0
         ratio(5)=xim2*etm
         ratio(6)=xip*etm2
         ratio(7)=xim2*etp
         ratio(8)=xim*etm2
         xim=xim/4.d0
         xip=xip/4.d0
         ratio(1)=xim*etm*(-xig-etp)
         ratio(2)=xip*etm*(xig-etp)
         ratio(3)=xip*etp*(xig-etm)
         ratio(4)=xim*etp*(-xig-etm)
      else
         write(*,*) '*ERROR in distattach: case with ',nterms
         write(*,*) '       terms is not covered'
         call exit(201)
      endif
!
!     calculating the position in the face
!      
      do i=1,3
         p(i)=0.d0
         do j=1,nterms
            p(i)=p(i)+ratio(j)*pneigh(i,j)
         enddo
      enddo
!
!     calculating the distance
!
      coeff=0.d0
      do i=1,3
         coeff=coeff+xn(i)*(p(i)-pnode(i))
      enddo
      dist=(p(1)-pnode(1)-coeff*xn(1))**2+(p(2)-pnode(2)-
     &     coeff*xn(2))**2+(p(3)-pnode(3)-coeff*xn(3))**2
!     
      return
      end
      
