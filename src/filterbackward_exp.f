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
      subroutine filterbackward_exp(adf,auf,jqf,irowf,ndesi,
     &   nodedesi,dgdxglob,dgdx,nobject,nk,nobjectstart,
     &   weighting)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      integer jqf(*),irowf(*),ndesi,nodedesi(*),nobject,
     &   nk,nobjectstart,idesvar,jdesvar,j,iobject,inode,jnode
!     
      real*8 auf(*),dgdxglob(2,nk,*),dgdx(ndesi,*),weighting(*),adf(*)
!
!     Explicit filter matrix A consists of the following matrices:
!     A = V^(-1) * W * M
!     V = weighting matrix to normalize the filter functions stored in 
!         the field weighting(*)
!     W = filter matrix 
!     M = mass matrix
!     In backward filtering s = A^T * x --> A^T = M * W * V^(-1)
!   
!     Step 1: multiplication of V^(-1) with sensitivities in dgdx
!     loop over all design responses
!   
      do iobject=1+nobjectstart,nobject
         do idesvar=1,ndesi
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject)
     &         /weighting(idesvar) 
         enddo
      enddo
!
!     Step 2: multiplication of W with sensitivities in dgdx
!     loop over all design responses
!   
      do iobject=1+nobjectstart,nobject
!   
!        loop over all rows for main- and sub-diagonal
!   
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
            dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
     &         +dgdx(idesvar,iobject)*adf(idesvar)
!
            do j=jqf(idesvar),jqf(idesvar+1)-1
               jdesvar=irowf(j)
               jnode=nodedesi(jdesvar)
!
!              entries for all rows in sub-diagonal
!   
               dgdxglob(2,jnode,iobject)=dgdxglob(2,jnode,iobject)
     &            +dgdx(idesvar,iobject)*auf(j)
!
!              entries for all columns in top-diagonal
!   
               dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
     &            +dgdx(jdesvar,iobject)*auf(j)
            enddo
         enddo
      enddo
!
      return        
      end
