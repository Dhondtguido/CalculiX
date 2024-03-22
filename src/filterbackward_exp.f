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
      subroutine filterbackward_exp(adf,auf,jqf,irowf,icolf,ndesi,
     &   nodedesi,dgdxglob,dgdx,nobject,nk,nobjectstart,objectset,
     &   weighting,temparray,adb,aub,jq,irow,icol)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      character*81 objectset(5,*)
!
      integer jqf(*),irowf(*),icolf(*),ndesi,nodedesi(*),nobject,ipos,
     &   inode1,inode2,kk,jj,ii,nk,nobjectstart,idesvar,jdesvar,j,
     &   iobject,inode,jnode,jq(*),irow(*),icol(*),iobjloc
!     
      real*8 auf(*),dgdxglob(2,nk,*),dgdx(ndesi,*),weighting(*),adf(*),
     &   temparray(ndesi,*),adb(*),aub(*)
!  
!
!     Explicit fitler matrix A consitst fo the following matrices:
!     A = V^(-1) * W * M
!     V = weighting matrix to normailze the filter functions stored in 
!         the field weighting(*)
!     W = filter matrix 
!     M = mass matrix
!     In backward filtering s = A^T * x --> A^T = M * W * V^(-1)
!   
!     Step1: multiplication of V^(-1) with sensitivities in dgdx
!     loop over all design responses
      do iobject=1+nobjectstart,nobject
         iobjloc=iobject-nobjectstart
         do idesvar=1,ndesi
            dgdx(idesvar,iobjloc)=dgdx(idesvar,iobjloc)
     &         /weighting(idesvar) 
         enddo
      enddo
!
!     Step2: multiplication of W with sensitivities in dgdx
!     loop over all design responses
      do iobject=1+nobjectstart,nobject
         iobjloc=iobject-nobjectstart
!   
!        loop over all rows for main- and sub-diagonal
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
            dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
     &         +dgdx(idesvar,iobjloc)*adf(idesvar)
!
            do j=jqf(idesvar),jqf(idesvar+1)-1
               jdesvar=irowf(j)
               jnode=nodedesi(jdesvar)
!
!              entries for all rows in sub-diagonal
               dgdxglob(2,jnode,iobject)=dgdxglob(2,jnode,iobject)
     &            +dgdx(idesvar,iobjloc)*auf(j)
!
!              entries for all columns in top-diagonal
               dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
     &            +dgdx(jdesvar,iobjloc)*auf(j)
            enddo
         enddo
      enddo
!
!     Step3: multiplication of M with sensitivities in temparray
!     loop over all design responses
C      do iobject=1+nobjectstart,nobject
C         iobjloc=iobject-nobjectstart
!
!        loop over all rows for main- and sub-diagonal
C         do idesvar=1,ndesi
C            inode=nodedesi(idesvar)
C      dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
C     &         +temparray(idesvar,iobjloc)*adb(idesvar)
!
C            do j=jq(idesvar),jq(idesvar+1)-1
C         jdesvar=irow(j)
C         jnode=nodedesi(jdesvar)
!
!              entries for all rows in sub-diagonal
C               dgdxglob(2,jnode,iobject)=dgdxglob(2,jnode,iobject)
C     &            +temparray(idesvar,iobjloc)*aub(j)
!
!              entries for all columns in top-diagonal
C               dgdxglob(2,inode,iobject)=dgdxglob(2,inode,iobject)
C     &            +temparray(jdesvar,iobjloc)*aub(j)
C      enddo
C   enddo
C      enddo
! 
C      do iobject=1+nobjectstart,nobject
C   do idesvar=1,ndesi
C      inode=nodedesi(idesvar)
C     dgdxglob(2,inode,iobject)=temparray(idesvar,iobject)
C  enddo
C      enddo
!
      return        
      end
