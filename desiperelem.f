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
      subroutine desiperelem(ndesi,istartdesi,ialdesi,ipoeldi,ieldi,
     &      ne,istartelem,ialelem)
!
      implicit none
!
      integer ndesi,istartdesi(*),ialdesi(*),ipoeldi(*),ieldi(2,*),
     &  ieldifree,i,j,nelem,ne,ifree,index,istartelem(*),ialelem(*)
!
!
!
!     storing the design variables per element
!
      ieldifree=1
      do i=1,ndesi
         do j=istartdesi(i),istartdesi(i+1)-1
            nelem=ialdesi(j)
            ieldi(1,ieldifree)=i
            ieldi(2,ieldifree)=ipoeldi(nelem)
            ipoeldi(nelem)=ieldifree
            ieldifree=ieldifree+1
         enddo
      enddo
!
!     adding the zero design variable to all elements with
!     a nonzero ipoeldi value
!
      do i=1,ne
         if(ipoeldi(i).eq.0) cycle
         ieldi(1,ieldifree)=0
         ieldi(2,ieldifree)=ipoeldi(i)
         ipoeldi(i)=ieldifree
         ieldifree=ieldifree+1
      enddo
!
!     determining the design variables belonging to a given 
!     element i. They are stored in ialelem(istartelem(i))..
!     ...up to..... ialdesi(istartelem(i+1)-1)
!
      ifree=1
      do i=1,ne
         istartelem(i)=ifree
         index=ipoeldi(i)
         do
            if(index.eq.0) exit
            ialelem(ifree)=ieldi(1,index)
c            write(*,*) 'desiperelem ',i,ialelem(ifree)
            ifree=ifree+1
            index=ieldi(2,index)
         enddo
      enddo
      istartelem(ne+1)=ifree
!
      return
      end
