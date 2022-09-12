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
      subroutine bsort(n,x,y,xmin,xmax,ymin,ymax,dmax,bin,list)
!     
      implicit none
!     
      integer list(*),bin(*),n,i,j,k,p,ndiv
!     
      real*8 x(*),y(*),factx,facty,xmin,xmax,ymin,ymax,dmax
!     
      ndiv=nint(real(n)**0.25d0)
      factx=real(ndiv)/((xmax-xmin)*1.01d0/dmax)
      facty=real(ndiv)/((ymax-ymin)*1.01d0/dmax)
      do 10 k=1,n
         p=list(k)
         i=int(y(p)*facty)
         j=int(x(p)*factx)
         if(mod(i,2).eq.0)then
            bin(p)=i*ndiv+j+1
         else
            bin(p)=(i+1)*ndiv-j
         end if
 10   continue
      call qsorti(n,list,bin)
      end
