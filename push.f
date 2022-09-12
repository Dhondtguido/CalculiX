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
      subroutine push(item,maxstk,topstk,stack)
!     
      implicit none
!     
      integer topstk,maxstk,stack(*),item
!     
      topstk=topstk+1
      if(topstk.gt.maxstk) then
         write(6,'("0***error in subroutine push***")')
         write(6,'("***stack overflow***")')
         call exit(201)
      else
         stack(topstk)=item
      endif
      end
