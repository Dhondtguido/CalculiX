!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine gennactdofinv3d(nactdof,nactdofinv,nk,mi)
!
!     inverting field nactdof, i.e. creating field nactdofinv
!     listing the node for each independent dof.
!     
      implicit none
!
      integer mi(*),nactdof(0:mi(2),*),nactdofinv(*),nk,i,j,mt
!
!     storing the nodes (in C convention, i.e. starting with 0)
!     in field nactdofinv
!
      mt=mi(2)+1
      do i=1,nk
         do j=0,mi(2)
            if(nactdof(j,i).le.0) cycle
            nactdofinv(nactdof(j,i))=(i-1)*mt+j
         enddo
      enddo
!
      return
      end
