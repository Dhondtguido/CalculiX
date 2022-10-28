!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine characteristiclength(co,istartcrackfro,iendcrackfro,
     &     ncrack,ifront,charlen,datarget)
!     
!     determine the mesh characteristic length for each front
!     
      implicit none
!     
      integer i,k,m,n1,n2,istartcrackfro(*),iendcrackfro(*),
     &     ncrack,ifront(*)
!     
      real*8 co(3,*),dist,sum,charlen(*),datarget
!     
!     first increment: determine for each front a characteristic length
!     
      do i=1,ncrack
!     
!     loop over all nodes belonging to the crack front(s)
!     
        k=0
        sum=0     
        do m=istartcrackfro(i),iendcrackfro(i)-1
!     
!     distance between two adjacent front nodes
!     
          n1=ifront(m)
          n2=ifront(m+1)
          dist=dsqrt((co(1,n2)-co(1,n1))**2+
     &         (co(2,n2)-co(2,n1))**2+
     &         (co(3,n2)-co(3,n2))**2)
          k=k+1
          sum=(sum*(k-1)+dist)/k
        enddo
!     
!     charlen is the mean distance between front nodes for each front
!     however, charlen should not be smaller than the target crack
!     propagation increment (favorizes smooth surfaces)
!     
        charlen(i)=max(datarget,sum)
      enddo
!     
      return
      end
      
      
