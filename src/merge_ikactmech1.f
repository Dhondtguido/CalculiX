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
      subroutine merge_ikactmech1(ikactmech1,nactmech1,neq,ikactmech,
     &     nactmech,num_cpu)
!
!     merging ikactmech1 obtained from the different threads into
!     ikactmech
!
      implicit none
!
      integer ikactmech1(neq,*),nactmech1(*),neq,ikactmech(*),
     &     nactmech,num_cpu,i,j,k,jdof,id
!
!     copying the first thread
!      
      nactmech=nactmech1(1)
      do j=1,nactmech
        ikactmech(j)=ikactmech1(j,1)
      enddo
!
!     adding the other threads
!
      do i=2,num_cpu
        do j=1,nactmech1(i)
          jdof=ikactmech1(j,i)
          call nident(ikactmech,jdof,nactmech,id)
          if(id.gt.0) then
            if(ikactmech(id).eq.jdof) cycle
          endif
          nactmech=nactmech+1
          do k=nactmech,id+2,-1
            ikactmech(k)=ikactmech(k-1)
          enddo
          ikactmech(id+1)=jdof
        enddo
      enddo
!     
      return
      end
