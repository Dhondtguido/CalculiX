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
      subroutine filterforward_exp(adf,auf,jqf,irowf,ndesi,
     &   nodedesi,gradproj,feasdir,weighting,temparray,adb,aub,
     &   jq,irow)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      integer jqf(*),irowf(*),ndesi,nodedesi(*),idesvar,jdesvar,j,
     &   inode,jnode,jq(*),irow(*)
!     
      real*8 auf(*),gradproj(3,*),feasdir(2,*),weighting(*),adf(*),
     &   temparray(*),adb(*),aub(*)
!
!     copy unfiltered sensitivities from gradproj to feasdir
!  
      do idesvar=1,ndesi
         inode=nodedesi(idesvar)
         feasdir(1,inode)=gradproj(3,inode)
      enddo
!
!     Explicit filter matrix A consitst for the following matrices:
!     A = V^(-1) * W * M
!     V = weighting matrix to normalize the filter functions stored in 
!         the field weighting(*)
!     W = filter matrix 
!     M = mass matrix
!     In forward filtering s = A * x --> A = V^(-1) * W * M
!   
!     Step 1: multiplication of M with sensitivities in temparray
!     loop over all rows for main- and sub-diagonal
!  
      do idesvar=1,ndesi
         inode=nodedesi(idesvar)
         temparray(idesvar)=temparray(idesvar)
     &      +gradproj(3,inode)*adb(idesvar)
!
         do j=jq(idesvar),jq(idesvar+1)-1
            jdesvar=irow(j)
            jnode=nodedesi(jdesvar)
!
!           loop over all rows in sub-diagonal
!  
            temparray(jdesvar)=temparray(jdesvar)
     &         +gradproj(3,inode)*aub(j)
!
!           loop over all columns in top-diagonal
!  
            temparray(idesvar)=temparray(idesvar)
     &         +gradproj(3,jnode)*aub(j)
         enddo
      enddo
!
!     Step 2: multiplication of W with sensitivities in dgdx
!     loop over all rows for main- and sub-diagonal
!  
      do idesvar=1,ndesi
         inode=nodedesi(idesvar)
         feasdir(2,inode)=feasdir(2,inode)
     &      +temparray(idesvar)*adf(idesvar)
!
         do j=jqf(idesvar),jqf(idesvar+1)-1
            jdesvar=irowf(j)
            jnode=nodedesi(jdesvar)
!
!           loop over all rows in sub-diagonal
!  
            feasdir(2,jnode)=feasdir(2,jnode)
     &         +temparray(idesvar)*auf(j)
!
!           loop over all columns in top-diagonal
!  
            feasdir(2,inode)=feasdir(2,inode)
     &         +temparray(jdesvar)*auf(j)
         enddo
      enddo
!
!     Step 3: multiplication of V^(-1) with sensitivities in dgdx
!  
      do idesvar=1,ndesi
         inode=nodedesi(idesvar)
         feasdir(2,inode)=feasdir(2,inode)/weighting(idesvar) 
      enddo
!
      return        
      end
