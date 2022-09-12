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
      subroutine calcmach(vold,v,nk,ntmat_,shcon,nshcon,physcon,
     &  inomat,mi)
!
      implicit none
!
      integer ntmat_,mi(*),nshcon(*),nk,i,imat,inomat(*)
!
      real*8 v(nk,0:mi(2)),vold(0:mi(2),*),cp,r,temp,
     &  shcon(0:3,ntmat_,*),physcon(*)
!     
!     calculate kappa (cp/cv) and store it in v(*,0)
!     calculate the Mach number and store it in v(*,1)
!     
      do i=1,nk
        imat=inomat(i)
        if(imat.eq.0) cycle
         temp=vold(0,i)
         call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &        nshcon,cp,physcon)
         r=shcon(3,1,imat)
         v(i,0)=cp/(cp-r)
         v(i,1)=dsqrt((vold(1,i)**2+vold(2,i)**2+vold(3,i)**2)
     &        /(v(i,0)*r*(temp-physcon(1))))
      enddo
!     
      return
      end
      
