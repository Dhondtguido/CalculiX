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
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!
!     here, an ABAQUS umat routine can be inserted
!
!     note that reals should be double precision (REAL*8)
!
      implicit none
!
      character*80 cmname
!
      integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,
     &  kstep,kinc
!
      real*8 stress(ntens),statev(nstatv),
     &  ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     &  stran(ntens),dstran(ntens),time(2),celent,
     &  props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     &  sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,predef,dpred,
     &  pnewdt
!
!     START EXAMPLE LINEAR ELASTIC MATERIAL
!
      integer i,j
      real*8 e,un,al,um,am1,am2
!
c      write(*,*) 'noel,npt ',noel,npt
c      write(*,*) 'stress ',(stress(i),i=1,6)
c      write(*,*) 'stran ',(stran(i),i=1,6)
c      write(*,*) 'dstran ',(dstran(i),i=1,6)
c      write(*,*) 'drot ',((drot(i,j),i=1,3),j=1,3)
      e=props(1)
      un=props(2)
      al=un*e/(1.d0+un)/(1.d0-2.d0*un)
      um=e/2.d0/(1.d0+un)
      am1=al+2.d0*um
      am2=um
!
!     stress
!      
      stress(1)=stress(1)+am1*dstran(1)+al*(dstran(2)+dstran(3))
      stress(2)=stress(2)+am1*dstran(2)+al*(dstran(1)+dstran(3))
      stress(3)=stress(3)+am1*dstran(3)+al*(dstran(1)+dstran(2))
      stress(4)=stress(4)+am2*dstran(4)
      stress(5)=stress(5)+am2*dstran(5)
      stress(6)=stress(6)+am2*dstran(6)
!
!     stiffness
!
      do i=1,6
         do j=1,6
            ddsdde(i,j)=0.d0
         enddo
      enddo
      ddsdde(1,1)=al+2.d0*um
      ddsdde(1,2)=al
      ddsdde(2,1)=al
      ddsdde(2,2)=al+2.d0*um
      ddsdde(1,3)=al
      ddsdde(3,1)=al
      ddsdde(2,3)=al
      ddsdde(3,2)=al
      ddsdde(3,3)=al+2.d0*um
      ddsdde(4,4)=um
      ddsdde(5,5)=um
      ddsdde(6,6)=um
!
!     END EXAMPLE LINEAR ELASTIC MATERIAL
!
      return
      end
