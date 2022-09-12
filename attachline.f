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
      subroutine attachline(pneigh,pnode,nterms,xil,etl,xn,p,dist)
!
!     returns the local coordinates in a face described by
!     the nodal coordinates in pneigh(1..3,*) of the intersection
!     with a straight line through point pnode(1..3) and
!     direction xn(1..3)
!
      implicit none
!
      integer nterms,i,j,k,imin,jmin,im,jm
!
      real*8 pneigh(3,9),pnode(3),xi(-1:1,-1:1),p(3),
     &  et(-1:1,-1:1),distmin,d1,dist,xil,etl,xn(3)
!
!
!
      d1=1.d0
!
      xi(0,0)=0.d0
      et(0,0)=0.d0
      call distattachline(xi(0,0),et(0,0),pneigh,pnode,dist,
     &     nterms,xn,p)
      distmin=dist
      imin=0
      jmin=0
!
      do k=1,6
!
!     initialisation
!
         d1=d1/10.d0
!     
         do i=-1,1
            do j=-1,1
               if((i.eq.0).and.(j.eq.0)) cycle
!
               xi(i,j)=xi(0,0)+i*d1
               et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
               if((xi(i,j).le.1.d0).and.
     &              (xi(i,j).ge.-1.d0).and.
     &              (et(i,j).le.1.d0).and.
     &              (et(i,j).ge.-1.d0)) then
                  call distattachline(xi(i,j),et(i,j),pneigh,pnode,
     &                 dist,nterms,xn,p)
!     
!                 checking for smallest initial distance
!     
                  if(dist.lt.distmin) then
                     distmin=dist
                     imin=i
                     jmin=j
                  endif
               endif
!
            enddo
         enddo
!     
!     minimizing the distance from the face to the node
!     
         do
!     
!     exit if minimum found
!     
            if((imin.eq.0).and.(jmin.eq.0)) exit
!
!           new center of 3x3 matrix
!
            xi(0,0)=xi(imin,jmin)
            et(0,0)=et(imin,jmin)
!
            im=imin
            jm=jmin
!
            imin=0
            jmin=0
!     
            do i=-1,1
               do j=-1,1
                  if((i+im.lt.-1).or.(i+im.gt.1).or.
     &                 (j+jm.lt.-1).or.(j+jm.gt.1)) then
!
                     xi(i,j)=xi(0,0)+i*d1
                     et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
                     if((xi(i,j).le.1.d0).and.
     &                  (xi(i,j).ge.-1.d0).and.
     &                  (et(i,j).le.1.d0).and.
     &                  (et(i,j).ge.-1.d0)) then
                        call distattachline(xi(i,j),et(i,j),pneigh,
     &                       pnode,dist,nterms,xn,p)
!
!                       check for new minimum
!
                        if(dist.lt.distmin) then
                           distmin=dist
                           imin=i
                           jmin=j
                        endif
                     endif
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
      if(nterms.eq.3) then
         if(xi(0,0)+et(0,0).le.0.d0) then
            xil=(xi(0,0)+1.d0)/2.d0
            etl=(et(0,0)+1.d0)/2.d0
         else
            xil=(1.d0-et(0,0))/2.d0
            etl=(1.d0-xi(0,0))/2.d0
         endif
      elseif(nterms.eq.4) then
         xil=xi(0,0)
         etl=et(0,0)
      elseif(nterms.eq.6) then
         if(xi(0,0)+et(0,0).le.0.d0) then
            xil=(xi(0,0)+1.d0)/2.d0
            etl=(et(0,0)+1.d0)/2.d0
         else
            xil=(1.d0-et(0,0))/2.d0
            etl=(1.d0-xi(0,0))/2.d0
         endif
      elseif(nterms.eq.8) then
         xil=xi(0,0)
         etl=et(0,0)
      endif
!     
      return
      end
      
