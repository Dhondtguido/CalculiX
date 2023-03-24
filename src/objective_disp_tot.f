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
      subroutine objective_disp_tot(dgdx,df,ndesi,iobject,jqs,
     &   irows,dgdu)
!
      implicit none
!
      integer ndesi,iobject,idesvar,j,jqs(*),irows(*),idof
!      
      real*8 dgdx(ndesi,*),df(*),dgdu(*)
!
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     non-linear:  dgdx = dgdx + dgdu * ( df )
!     ----------------------------------------------------------------
!
!     Calculation of the total differential:    
!
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject) 
     &           +dgdu(idof)*df(j)
         enddo
      enddo
!      
      return
      end
