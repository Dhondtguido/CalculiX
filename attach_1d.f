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
      subroutine attach_1d(pneigh,pnode,nterms,ratio,dist,xil)
!
!     ataches node with coordinates in "pnode" to the line containing 
!     "nterms" nodes with coordinates in field "pneigh" (nterms < 9).
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      implicit none
!
      integer nterms,i,k,imin,im
!
      real*8 ratio(3),pneigh(3,3),pnode(3),a,xi(-1:1),
     &  p(3),distmin,d1,dist,xil
!
!
!
      d1=1.d0
!
      xi(0)=0.d0
      call distattach_1d(xi(0),pneigh,pnode,a,p,
     &     ratio,nterms)
      distmin=a
      imin=0
!
      do k=1,8
!
!     initialisation
!
         d1=d1/10.d0
!     
         do i=-1,1
            if(i.eq.0) cycle
!     
            xi(i)=xi(0)+i*d1
!
!              check whether inside the (-1,1) domain
!
            if((xi(i).le.1.d0).and.
     &           (xi(i).ge.-1.d0)) then
               call distattach_1d(xi(i),pneigh,pnode,a,p,
     &              ratio,nterms)
!     
!                 checking for smallest initial distance
!     
               if(a.lt.distmin) then
                  distmin=a
                  imin=i
               endif
            endif
!
         enddo
!     
!     minimizing the distance from the face to the node
!     
         do
!     
!     exit if minimum found
!     
            if(imin.eq.0) exit
!
!           new center of 3 vector
!
            xi(0)=xi(imin)
!
            im=imin
!
            imin=0
!     
            do i=-1,1
               if((i+im.lt.-1).or.(i+im.gt.1)) then
!     
                  xi(i)=xi(0)+i*d1
!
!              check whether inside the (-1,1) domain
!
                  if((xi(i).le.1.d0).and.
     &                 (xi(i).ge.-1.d0)) then
                     call distattach_1d(xi(i),pneigh,
     &                    pnode,a,p,ratio,nterms)
!
!                       check for new minimum
!
                     if(a.lt.distmin) then
                        distmin=a
                        imin=i
                     endif
                  endif
!     
               endif
            enddo
         enddo
      enddo
!
      call distattach_1d(xi(0),pneigh,pnode,a,p,
     &     ratio,nterms)
!
      do i=1,3
        pnode(i)=p(i)
      enddo
!
      dist=dsqrt(a)
      xil=xi(0)
!     
      return
      end
      
