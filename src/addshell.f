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
      subroutine addshell(nactdof,node,b,mi,iperturb,nmethod,cam,
     &     v)
!
!     updates the translational and rotational dofs
!     for true shells
!
!     the translational dofs are simply added;
!      
!     the rotational dofs
!     are at the start of this routine available as rotational
!     vectors. These vectors are transformed into a rotational matrix,
!     both matrices are multiplied and reverted into a vector;
!
      implicit none
!
      integer mi(*),nactdof(0:mi(2),*),node,iperturb(*),nmethod,
     &     i,j,k,l
!
      real*8 b(*),bnac(6),cam(*),v(0:mi(2),*),xn(3),ww,c1,c2,c3,
     &     r0(3,3),dr(3,3),r(3,3),th,d(3,3),e(3,3,3)
!
!     d(3,3): delta Dirac function
!     e(3,3,3): alternating symbol
!
      data d /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data e /0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,
     &     0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,0.d0,
     &     0.d0,-1.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
!
!     storing the change in solution in bnac
!
      do j=1,6
        if(nactdof(j,node).gt.0) then
          bnac(j)=b(nactdof(j,node))
        else
          bnac(j)=0.d0
        endif
        if((iperturb(1).ne.0).and.(abs(nmethod).eq.1)) then
          if(dabs(bnac(j)).gt.cam(1)) then
            cam(1)=dabs(bnac(j))
            cam(4)=nactdof(j,node)-0.5d0
          endif
        endif
      enddo
!
!     updata translational dofs
!
      do j=1,3
        v(j,node)=v(j,node)+bnac(j)
      enddo
!
!     update rotational dofs
!      
!     previous solution      
!
      xn(1)=v(4,node)
      xn(2)=v(5,node)
      xn(3)=v(6,node)
!
      ww=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
!     
      c1=dcos(ww)
      if(ww.gt.1.d-10) then
        c2=dsin(ww)/ww
      else
        c2=1.d0
      endif
      if(ww.gt.1.d-5) then
        c3=(1.d0-c1)/ww**2
      else
        c3=0.5d0
      endif
!     
!     rotation matrix r0 (Buch by Dhondt, Wiley(2004), Eq. (3.76))
!     
      do k=1,3
        do l=1,3
          r0(k,l)=c1*d(k,l)+
     &         c2*(e(k,1,l)*xn(1)+e(k,2,l)*xn(2)+
     &         e(k,3,l)*xn(3))+c3*xn(k)*xn(l)
        enddo
      enddo
!
!     change in solution
!
      xn(1)=bnac(4)
      xn(2)=bnac(5)
      xn(3)=bnac(6)
!
      ww=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
!     
      c1=dcos(ww)
      if(ww.gt.1.d-10) then
        c2=dsin(ww)/ww
      else
        c2=1.d0
      endif
      if(ww.gt.1.d-5) then
        c3=(1.d0-c1)/ww**2
      else
        c3=0.5d0
      endif
!     
!     rotation matrix dr (Buch by Dhondt, Wiley(2004), Eq. (3.76))
!     
      do k=1,3
        do l=1,3
          dr(k,l)=c1*d(k,l)+
     &         c2*(e(k,1,l)*xn(1)+e(k,2,l)*xn(2)+
     &         e(k,3,l)*xn(3))+c3*xn(k)*xn(l)
        enddo
      enddo
!
!     multiplying the matrices
!
      do i=1,3
        do j=1,3
          r(i,j)=0.d0
          do k=1,3
            r(i,j)=r(i,j)+dr(i,k)*r0(k,j)
          enddo
        enddo
      enddo
!
!     convert into a rotational vector (inverse formulas obtained by:
!     C_ii=1+2*cos(theta)
!     (C_ij-C_ji)/2=sin(theta)*e_ikj*n_k
!
      th=dacos((r(1,1)+r(2,2)+r(3,3)-1.d0)/2.d0)
!
      if(dabs(th).lt.1.d-10) then
!
!       default: rotation of 1.d-10 about the x-axis
!
        v(4,node)=1.d-10
        v(5,node)=0.d0
        v(6,node)=0.d0
      else
        th=th/(2.d0*dsin(th))
        v(4,node)=th*(r(3,2)-r(2,3))
        v(5,node)=th*(r(1,3)-r(3,1))
        v(6,node)=th*(r(2,1)-r(1,2))
      endif
!     
      return
      end
      
