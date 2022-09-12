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
!     diagonal preconditioning of a matrix
!     
      subroutine inviscidpre(nk,inomat,ntmat_,shcon,nshcon,physcon,
     &     xmach2,nactdoh,b,vold,v,mi)
!     
      implicit none
!     
      integer i,nk,imat,inomat(*),ntmat_,nshcon(*),nactdoh(0:4,*),
     &     idof,mi(*)
!     
      real*8 temp,shcon(0:3,ntmat_,*),cp,physcon(*),r,xkappa,vel2,
     &     xmach,xmach2(*),constant,b(*),vold(0:mi(2),*),
     &     v(0:mi(2),*)
!     
      do i=1,nk
        imat=inomat(i)
        if(imat.eq.0) cycle
!     
        idof=nactdoh(0,i)
        temp=vold(0,i)
!     
!     material properties r, cp and kappa
!     
        call materialdata_cp_sec(imat,ntmat_,temp,shcon,
     &       nshcon,cp,physcon)
        r=shcon(3,1,imat)
        xkappa=cp/(cp-r)
!     
!     size of the velocity
!     
        vel2=vold(1,i)**2+vold(2,i)**2+vold(3,i)**2
!     
        xmach=dsqrt(vel2)/(xkappa*r*(temp-physcon(1)))
!     
!     determine the relevant Mach number
!     
        if(xmach.lt.1.d-5) then
          xmach2(idof)=1.d-10
        elseif(xmach.lt.1.d0) then
          xmach2(idof)=xmach**2
        else
          xmach2(idof)=1.d0
        endif
!     
        constant=1.d0-1.d0/xmach2(idof)
!     
!     modifying the right hand side
!     
        b(idof)=b(idof)+vel2*constant/3.d0*v(0,i)-constant*
     &       -(vold(1,i)*v(1,i)+vold(2,i)*v(2,i)+vold(3,i)*v(3,i))
!     
      enddo
!     
      return
      end
