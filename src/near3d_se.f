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
      subroutine near3d_se(xo,yo,zo,x,y,z,nx,ny,nz,xp,yp,zp,n,
     &     ir,r,nr,radius)     
!
!     determines the nodes out of n within a radius r of
!     the point with coordinates (xp,yp,zp);
!
!
!     INPUT:
!
!     xo                 x-coordinates of cloud of nodes
!     yo                 y-coordinates of cloud of nodes
!     zo                 z-coordinates of cloud of nodes
!     x                  xo ordered in increasing order
!                        (can be done in the calling program 
!                         with dsort)
!     y                  yo ordered in increasing order
!     z                  zo ordered in increasing order
!     nx                 permutations of x-ordering
!     ny                 permutations of y-ordering
!     nz                 permutations of z-ordering
!     xp                 x-coordinate of point of interest
!     yp                 y-coordinate of point of interest
!     zp                 z-coordinate of point of interest
!     n                  number of nodes in cloud
!     radius             radius
!
!     OUTPUT:
!
!     ir                 numbers of the nodes within the given radius
!     r                  distance square of the nodes within the given
!                        radius
!     nr                 number of nodes within the given radius
!
      implicit none
!
      integer n,nx(n),ny(n),nz(n),ir(n+6),nr,nrprev,irnew,
     &  i,j,k,m,id,idx,idy,idz,node
!
      real*8 x(n),y(n),z(n),xo(n),yo(n),zo(n),xp,yp,zp,r(n+6),
     &  xr,yr,zr,c(8),dd,xw,xe,ys,yn,zb,zt,radius,
     &  radius2
!
      radius2=radius*radius
      nrprev=0
!
!     identify position of xp, yp and zp
!      
      call ident(x,xp,n,idx)
      call ident(y,yp,n,idy)
      call ident(z,zp,n,idz)
!
!     initialization of the maximal distance in each direction
!
      xw=0.d0
      xe=0.d0
      ys=0.d0
      yn=0.d0
      zb=0.d0
      zt=0.d0
!
      i=1
!
      do
!
         nr=nrprev
!
!        westp 
!
         id=idx+1-i
         if(id.gt.0) then
            node=nx(id)
            xw=xo(node)-xp
            yr=yo(node)-yp
            zr=zo(node)-zp
            dd=xw*xw+yr*yr+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            xw=1.d30
         endif
!
!        east
!
         id=idx+i
         if(id.le.n) then
            node=nx(id)
            xe=xo(node)-xp
            yr=yo(node)-yp
            zr=zo(node)-zp
            dd=xe*xe+yr*yr+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            xe=1.d30
         endif
!
!        south
!
         id=idy+1-i
         if(id.gt.0) then
            node=ny(id)
            xr=xo(node)-xp
            ys=yo(node)-yp
            zr=zo(node)-zp
            dd=xr*xr+ys*ys+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            ys=1.d30
         endif
!
!        north
!
         id=idy+i
         if(id.le.n) then
            node=ny(id)
            xr=xo(node)-xp
            yn=yo(node)-yp
            zr=zo(node)-zp
            dd=xr*xr+yn*yn+zr*zr
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            yn=1.d30
         endif
!
!        bottom
!
         id=idz+1-i
         if(id.gt.0) then
            node=nz(id)
            xr=xo(node)-xp
            yr=yo(node)-yp
            zb=zo(node)-zp
            dd=xr*xr+yr*yr+zb*zb
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            zb=1.d30
         endif
!
!        top
!
         id=idz+i
         if(id.le.n) then
            node=nz(id)
            xr=xo(node)-xp
            yr=yo(node)-yp
            zt=zo(node)-zp
            dd=xr*xr+yr*yr+zt*zt
            if(dd.lt.radius2) then
               nr=nr+1
               ir(nr)=node
            endif
         else
            zt=1.d30
         endif
!
!        check for new entries
!
         if(nr.gt.nrprev) then
            m=nrprev
            do j=nrprev+1,nr
               irnew=ir(j)
               call nident(ir,irnew,m,id)
               if(id.eq.0) then
                  m=m+1
                  do k=m,2,-1
                     ir(k)=ir(k-1)
                  enddo
                  ir(1)=irnew
               elseif(ir(id).ne.irnew) then
                  m=m+1
                  do k=m,id+2,-1
                     ir(k)=ir(k-1)
                  enddo
                  ir(id+1)=irnew
               endif
            enddo
            nrprev=m
         endif
!
         i=i+1
!
!        check the corners of the box
!
         c(1)=xe*xe+yn*yn+zb*zb
         if(c(1).lt.radius2) cycle
         c(2)=xw*xw+yn*yn+zb*zb
         if(c(2).lt.radius2) cycle
         c(3)=xw*xw+ys*ys+zb*zb
         if(c(3).lt.radius2) cycle
         c(4)=xe*xe+ys*ys+zb*zb
         if(c(4).lt.radius2) cycle
         c(5)=xe*xe+yn*yn+zt*zt
         if(c(5).lt.radius2) cycle
         c(6)=xw*xw+yn*yn+zt*zt
         if(c(6).lt.radius2) cycle
         c(7)=xw*xw+ys*ys+zt*zt
         if(c(7).lt.radius2) cycle
         c(8)=xe*xe+ys*ys+zt*zt
         if(c(8).lt.radius2) cycle
!
!        no new entries possible: finished
!
         nr=nrprev
         do j=1,nr
            node=ir(j)
            xr=xo(node)-xp
            yr=yo(node)-yp
            zr=zo(node)-zp
            r(j)=xr*xr+yr*yr+zr*zr
         enddo
         exit
!
      enddo
!
      return
      end
