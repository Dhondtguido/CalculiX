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
      function pop(topstk,stack)
!     
      implicit none
!     
      integer pop,topstk,stack(*)
!     
      if(topstk.gt.0) then
         pop=stack(topstk)
         topstk=topstk-1
      else
         write(6,'("0***error in function pop***")')
         write(6,'("***stack underflow***")')
         call exit(201)
      end if
      end
      
