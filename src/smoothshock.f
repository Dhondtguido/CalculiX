!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine smoothshock(aub,adl,sol,aux,irow,jq,
     &  neqa,neqb,sa)
!
!     smoothing the finite element solution
!
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer irow(*),jq(*),neqa,neqb,i,j
!
      real*8 aub(*),adl(*),sol(*),aux(*),sa(*)
!
      do i=neqa,neqb
        do j=jq(i),jq(i+1)-1
          aux(i)=aux(i)+aub(j)*sol(irow(j))
        enddo
         sol(i)=sol(i)+sa(i)*aux(i)*adl(i)
      enddo
!
      return
      end
