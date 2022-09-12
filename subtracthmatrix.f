!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine subtracthmatrix(neqp,aubh,adbh,aux,dp,jqp,irowp,b,
     &     theta1,dtimef)
!     
!     transfer effect of H-matrix for compressible fluids to the
!     rhs
!     
      implicit none
!
      integer i,j,neqp,ir,ic,jqp(*),irowp(*)
!
      real*8 aux(*),aubh(*),adbh(*),dp(*),value,constant,theta1,dtimef,
     &     b(*)
!
      do i=1,neqp
        aux(i)=adbh(i)*dp(i)
      enddo
!
      do ic=1,neqp
        do j=jqp(ic),jqp(ic+1)-1
          ir=irowp(j)
          value=aubh(j)
          aux(ir)=aux(ir)+value*dp(ic)
          aux(ic)=aux(ic)+value*dp(ir)
        enddo
      enddo
!
      constant=dtimef*dtimef*theta1
!
      do i=1,neqp
        b(i)=b(i)-constant*aux(i)
      enddo
!     
      return
      end
      
