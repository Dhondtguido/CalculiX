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
      function triloc(xp,yp,x,y,v,e,numtri)
!     
      implicit none
!     
      integer v(3,*),e(3,*),numtri,v1,v2,i,t,triloc
!     
      real*8 x(*),y(*),xp,yp
!     
      t=numtri
 10   continue
      do 20 i=1,3
         v1=v(i,t)
         v2=v(mod(i,3)+1,t)
         if((y(v1)-yp)*(x(v2)-xp).gt.(x(v1)-xp)*(y(v2)-yp)) then
            t=e(i,t)
            goto 10
         end if
 20   continue
      triloc=t
      end
      
