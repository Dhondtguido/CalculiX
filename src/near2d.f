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
      subroutine near2d(xo,yo,x,y,nx,ny,xp,yp,n,neighbor,k)
!
!     determines the k closest nodes out of n with coordinates in
!     (xo,yo) to the point with coordinates (xp,yp);
!
      implicit none
!
      integer n,nx(n),ny(n),ir(k+4),nr,neighbor(k),kflag,iflag,
     &  i,j,k,l,m,id,idx,idy,four,node
!
      real*8 x(n),y(n),xo(n),yo(n),xp,yp,r(k+4),xr,yr,c(4),dd,
     &  xw,xe,ys,yn,sqrt_rmaxini
!
      data iflag /1/
      data kflag /2/
      data four /4/
!
      if(k.gt.n) then
         k=n
      endif
!
!     identify position of xp and yp
!      
      call ident(x,xp,n,idx)
      call ident(y,yp,n,idy)
!
!     initialization of r and ir
!
      do i=1,k
         xr=xo(i)-xp
         yr=yo(i)-yp
         r(i)=xr*xr+yr*yr
         ir(i)=i
      enddo
      call dsort(r,ir,k,kflag)
      sqrt_rmaxini=1.d30
!
!     initialization of the maximal distance in each direction
!
      xw=0.d0
      xe=0.d0
      ys=0.d0
      yn=0.d0
!
      i=1
!
      do
!
         nr=k
!
!        west
!
         id=idx+1-i
         if(id.gt.0) then
            node=nx(id)
            xw=xo(node)-xp
            yr=yo(node)-yp
            dd=xw*xw+yr*yr
            if(dd.lt.r(k)) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.0) then
            xw=sqrt_rmaxini
         endif
!
!        east
!
         id=idx+i
         if(id.le.n) then
            node=nx(id)
            xe=xo(node)-xp
            yr=yo(node)-yp
            dd=xe*xe+yr*yr
            if(dd.lt.r(k)) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.n+1) then
            xe=sqrt_rmaxini
         endif
!
!        south
!
         id=idy+1-i
         if(id.gt.0) then
            node=ny(id)
            xr=xo(node)-xp
            ys=yo(node)-yp
            dd=xr*xr+ys*ys
            if(dd.lt.r(k)) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.0) then
            ys=sqrt_rmaxini
         endif
!
!        north
!
         id=idy+i
         if(id.le.n) then
            node=ny(id)
            xr=xo(node)-xp
            yn=yo(node)-yp
            dd=xr*xr+yn*yn
            if(dd.lt.r(k)) then
               nr=nr+1
               ir(nr)=node
               r(nr)=dd
            endif
         elseif(id.eq.n+1) then
            yn=sqrt_rmaxini
         endif
!
!        check the corners of the box
!
         c(1)=xe*xe+yn*yn
         c(2)=xw*xw+yn*yn
         c(3)=xw*xw+ys*ys
         c(4)=xe*xe+ys*ys
         call insertsortd(c,four)
c         call dsort(c,idummy,four,iflag)
!
!        check for new entries
!
         if(nr.gt.k) then
            call dsort(r,ir,nr,kflag)
!     
!     reject equal entries
!     
            m=1
            if(m.lt.k) then
               loop: do j=2,nr
                  do l=m,1,-1
                     if(ir(j).eq.ir(l)) cycle loop
                  enddo
                  m=m+1
                  r(m)=r(j)
                  ir(m)=ir(j)
                  if(m.eq.k) exit
               enddo loop
            endif
         endif
         if(c(1).ge.r(k)) exit
!     
         i=i+1
!
      enddo
!
      do i=1,k
         neighbor(i)=ir(i)
      enddo
!
      return
      end
