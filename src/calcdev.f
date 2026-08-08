!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine calcdev(vold,vcon,v,nk,iturbulent,mi,vconmax,vmax,
     &     iexplicit,num_cpus)
!
!     calculates the change in solution
!     
      implicit none
!
      integer iturbulent,mi(*),nk,i,j,k,l,iexplicit,num_cpus
!
      real*8 v(nk,0:mi(2)),vold(0:mi(2),*),vcon(nk,0:mi(2)),vmax(0:6),
     &     vconmax(0:6),sum1,sum2
!
      vmax=0.d0
      vconmax=0.d0
      sum1=0.d0
      sum2=0.d0
!
!$omp parallel private(j,k,l), num_threads(num_cpus)
!
      do j=0,3
!$omp do reduction(+:sum1)
         do i=1,nk
            sum1=sum1+vcon(i,j)**2
         end do
!$omp end do
         vconmax(j)=vconmax(j)+sum1
      end do
!
      sum1=0.d0
      if(iexplicit.eq.1) then
!$omp do reduction(+:sum1)
         do i=1,nk
!           for incompressible fluids the density is stored
!           in vcon(4,*), the change in density in v(*,4)
            sum1=sum1+vcon(i,4)**2
         end do
!$omp end do
      else
!$omp do reduction(+:sum1)
         do i=1,nk
!           for incompressible fluids the pressure is stored
!           in vold(4,*), the change in pressure in v(*,4)
            sum1=sum1+vold(4,i)**2
         end do
!$omp end do
      end if
      vconmax(4)=vconmax(4)+sum1
!
      sum1=0.d0
      do k=0,4
!$omp do reduction(+:sum1)
         do i=1,nk
            sum1=sum1+v(i,k)**2
         end do
!$omp end do
         vmax(k)=vmax(k)+sum1
      end do
!
      if(iturbulent.ne.0) then
         do l=5,6
            sum1=0.d0
            sum2=0.d0
!$omp do reduction(+:sum1,sum2)
            do i=1,nk
               sum1=sum1+vcon(i,l)**2
               sum2=sum2+v(i,l)**2
            enddo
!$omp end do
            vconmax(l)=vconmax(l)+sum1
            vmax(l)=vmax(l)+sum2
         enddo
      endif
!
!$omp end parallel
!
      return
      end
