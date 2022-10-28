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
!     S.W. Sloan, Adv.Eng.Software,1987,9(1),34-55.
!     Permission for use with the GPL license granted by Prof. Scott
!     Sloan on 17. Nov. 2013
!
      subroutine deltri(numpts,n,x,y,list,bin,v,e,numtri)
!     
      implicit none
!     
      integer n,i,list(*),v(3,*),e(3,*),numtri,bin(*),p,numpts
!     
      real*8 xmin,xmax,ymin,ymax,dmax,c00001,fact,x(*),y(*)
!     
      parameter(c00001=1.d0)
!     
      xmin=x(list(1))
      xmax=xmin
      ymin=y(list(1))
      ymax=ymin
      do 5 i=2,n
         p=list(i)
         xmin=min(xmin,x(p))
         xmax=max(xmax,x(p))
         ymin=min(ymin,y(p))
         ymax=max(ymax,y(p))
 5    continue
      dmax=max(xmax-xmin,ymax-ymin)
      fact=c00001/dmax
      do 10 i=1,n
         p=list(i)
         x(p)=(x(p)-xmin)*fact
         y(p)=(y(p)-ymin)*fact
 10   continue
      call bsort(n,x,y,xmin,xmax,ymin,ymax,dmax,bin,list)
      call delaun(numpts,n,x,y,list,bin,v,e,numtri)
      do 30 i=1,n
         p=list(i)
         x(p)=x(p)*dmax+xmin
         y(p)=y(p)*dmax+ymin
 30   continue
      end
      
