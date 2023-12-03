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
      subroutine attach_1d_cracklength(pneigh,ratio,dist,xil,
     &     xo,yo,zo,x,y,z,nx,ny,nz,n,pnode)
!
!     searches a position on the line containing 
!     2 nodes with coordinates in field "pneigh" which
!     maximizes the minimum distance from all front nodes
!     of a crack
!
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      implicit none
!
      integer i,k,imin,im,nx(*),ny(*),nz(*),n
!
      real*8 ratio(3),pneigh(3,3),a,xi(-1:1),p(3),distmin,d1,dist,xil,
     &     xo(*),yo(*),zo(*),x(*),y(*),z(*),pnode(*)
!
      d1=1.d0
!
      xi(0)=0.d0
      call distattach_1d_cracklength(xi(0),pneigh,a,p,ratio,
     &     xo,yo,zo,x,y,z,nx,ny,nz,n)
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
               call distattach_1d_cracklength(xi(i),pneigh,a,p,ratio,
     &     xo,yo,zo,x,y,z,nx,ny,nz,n)
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
                     call distattach_1d_cracklength(xi(i),pneigh,
     &                    a,p,ratio,xo,yo,zo,x,y,z,nx,ny,nz,n)
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
      call distattach_1d_cracklength(xi(0),pneigh,a,p,ratio,
     &     xo,yo,zo,x,y,z,nx,ny,nz,n)
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
      
