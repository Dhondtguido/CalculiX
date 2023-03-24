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
      subroutine objective_freq_cs(dgdx,
     &  df,vold,ndesi,iobject,mi,nactdofinv,
     &  jqs,irows,nk,nzss)
!
      implicit none
!
      integer ndesi,iobject,mi(*),idesvar,j,idir,nk,nzss,
     &  jqs(*),irows(*),nactdofinv(*),node,idof,inode,mt
!      
      real*8 dgdx(ndesi,*),df(*),vold(0:mi(2),*)
!
!
!
!     ----------------------------------------------------------------
!     Calculation of the frequency sensitivity in the cyclic symmetric case
!
!     dgdx = -vold_R^T*df_R+vold_I^T*df_I
!
!     notice that df was calculated in mafillsmcsse as 
!     - (K-lambda*M)*eigenvector
!     ----------------------------------------------------------------
!     
      mt=mi(2)+1
!
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            inode=nactdofinv(idof)
            node=inode/mt+1
            idir=inode-mt*(inode/mt)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject)
     &          -vold(idir,node)*df(j)-vold(idir,nk+node)*df(nzss+j)
         enddo
      enddo
!      
      return
      end
