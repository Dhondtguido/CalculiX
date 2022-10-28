!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine mulmatvec_asym(au,jq,irow,ncol,x,y,itranspose)
!     
!     asymmetric sparse matrix vector multiplication in
!     Compressed Sparse Column (CSC) format: y = a*x
!      
      implicit none
!     
      integer i,j,itranspose,jq(*),irow(*),ncol
!     
      real*8 au(*),x(*),y(*)
!     
!     itranspose=0: non transposed
!     itranspose=1: transposed
!     
      if(itranspose.eq.0) then
!     
!     NONtransposed multiplication
!     
        do i=1,ncol
          do j=jq(i),jq(i+1)-1
            y(irow(j))=y(irow(j))+au(j)*x(i)
          enddo
        enddo
!     
      else
!     
!     transposed multiplication 
!     
        do i=1,ncol
          do j=jq(i),jq(i+1)-1
            y(i)=y(i)+au(j)*x(irow(j))
          enddo
        enddo
      endif
!     
      return
      end
