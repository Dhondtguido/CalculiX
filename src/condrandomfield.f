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
      subroutine condrandomfield(ad,au,jqs,irows,ndesi,
     &  rhs,vector,idesvar,jqc,auc,irowc)         
!
!     computation of the conditional randomfield entries
!
      implicit none
!
      integer jqs(*),irows(*),ndesi,i,j,irowc(*),idesvar,jqc(*),irow
!
      real*8 ad(*),au(*),auc(*),rhs(*),vector(*)
!
      do i=1,ndesi
         do j=jqc(i),jqc(i+1)-1
            irow=irowc(j)
            vector(i)=vector(i)+rhs(irow)*auc(j)
         enddo
      enddo
!
!     subtraction of diagonal entry  
!
      ad(idesvar)=ad(idesvar)-vector(idesvar)
!
!     subtraction of subdiagonal entries
!
      do i=jqs(idesvar),jqs(idesvar+1)-1
         irow=irows(i)
         au(i)=au(i)-vector(irow)
      enddo      
!
      return        
      end
