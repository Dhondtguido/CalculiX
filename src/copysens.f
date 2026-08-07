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
      subroutine copysens(rhs,dgdxglob,iobject,icopy,
     &   nk,ndesi,nodedesi)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      integer iobject,icopy,nk,ndesi,nodedesi(*),idesvar,inode,istart
!     
      real*8 rhs(*),dgdxglob(2,nk,*)
!  
!     icopy=0: copy unfiltered sensitivities to dgdxglob
!     icopy=1: copy filtered sensitivities to dgdxglob
!
!     FORTRAN convention for iobject
!
      istart=iobject+1
!
!     copy unfiltered sensitivities
!
      if(icopy.eq.0) then
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
            dgdxglob(1,inode,istart)=rhs(idesvar)
         enddo
!   
!     copy filtered sensitivities
!
      elseif(icopy.eq.1) then
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
            dgdxglob(2,inode,istart)=rhs(idesvar)
         enddo
      endif
!      
      return        
      end
