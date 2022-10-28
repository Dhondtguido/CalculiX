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
      subroutine radresult(ntr,xloadact,bcr,nloadtr,tarea,
     &     tenv,physcon,erad,auview,fenv,irowrad,jqrad,
     &     nzsrad,q)
!     
      implicit none
!     
      integer i,j,k,ntr,nloadtr(*),irowrad(*),jqrad(*),nzsrad
!     
      real*8 xloadact(2,*), tarea(*),tenv(*),auview(*),
     &     erad(*),fenv(*),physcon(*),bcr(ntr),q(*)
!     
!     calculating the flux and transforming the flux into an
!     equivalent temperature
!     
      write(*,*) ''
!      
      do i=1,ntr
         q(i)=bcr(i)
      enddo
!
!        lower triangle
!
      do i=1,ntr
         do j=jqrad(i),jqrad(i+1)-1
            k=irowrad(j)
            q(k)=q(k)-auview(j)*bcr(i)
!     
!        upper triangle
!
            q(i)=q(i)-auview(nzsrad+j)*bcr(k)
         enddo
      enddo
!
      do i=1,ntr
         j=nloadtr(i)
         q(i)=q(i)-fenv(i)*physcon(2)*tenv(i)**4
         xloadact(2,j)=
     &        max(tarea(i)**4-q(i)/(erad(i)*physcon(2)),0.d0)
         xloadact(2,j)=(xloadact(2,j))**0.25d0+physcon(1)
      enddo
!     
      return
      end
