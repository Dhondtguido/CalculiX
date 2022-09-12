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
      subroutine cmatrix(ad,au,jqs,irows,icols,ndesi,nodedesi,
     &     auc,jqc,irowc,nodedesibou)         
!     
!     calculates the values of the C-matrix
!     
      implicit none
!     
      integer jqs(*),irows(*),icols(*),ndesi,nodedesi(*),
     &     inode1,inode2,jqc(*),irowc(*),nodedesibou(*),
     &     kk,jj,irow,ipos1,ipos2
!     
      real*8 ad(*),au(*),auc(*)
!     
      do kk=1,ndesi
        inode1=nodedesi(kk)
        do jj=jqc(kk),jqc(kk+1)-1
          irow=irowc(jj)
          inode2=nodedesibou(irow)
          if(inode1.eq.inode2) then
            auc(jj)=ad(kk)
          elseif(inode1.lt.inode2) then
            call nident(nodedesi,inode2,ndesi,ipos1)
            call nident(irows(jqs(kk)),ipos1,icols(kk),ipos2)
            auc(jj)=au(jqs(kk)-1+ipos2)
          else
            call nident(nodedesi,inode2,ndesi,ipos1)
            call nident(irows(jqs(ipos1)),kk,icols(ipos1),ipos2)
            auc(jj)=au(jqs(ipos1)-1+ipos2)
          endif
        enddo
      enddo
!     
      return        
      end
