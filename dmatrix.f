!     
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
      subroutine dmatrix(ad,au,jqs,irows,icols,ndesi,nodedesi,
     &     add,aud,jqd,irowd,ndesibou,nodedesibou)         
!     
!     calculates the values of the D-matrix
!     
      implicit none
!     
      integer jqs(*),irows(*),icols(*),ndesi,nodedesi(*),idof,j,jdof,
     &     inode1,inode2,jqd(*),irowd(*),ndesibou,nodedesibou(*),ipos1,
     &     ipos2,ipos3
!     
      real*8 ad(*),au(*),add(*),aud(*)
!     
      do idof=1,ndesibou
        inode1=nodedesibou(idof)
        call nident(nodedesi,inode1,ndesi,ipos1)
        add(idof)=ad(ipos1)
        do j=jqd(idof),jqd(idof+1)-1
          jdof=irowd(j)
          inode2=nodedesibou(jdof)
          call nident(nodedesi,inode2,ndesi,ipos2)
          call nident(irows(jqs(ipos1)),ipos2,icols(ipos1),ipos3)
!     
!     assign the value to the D-matrix
!     
          aud(j)=au(jqs(ipos1)-1+ipos3)     
        enddo
      enddo
!     
      return        
      end
