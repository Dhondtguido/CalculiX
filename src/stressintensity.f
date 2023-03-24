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
      subroutine stressintensity(nfront,ifrontrel,stress,xt,xn,xa,
     &     dk1,dk2,dk3,xkeq,phi,psi,acrack,shape,nstep)
!     
!     calculate the stress intensity factors along the crack fronts
!     
      implicit none
!     
      integer nfront,ifrontrel(*),i,noderel,nstep,m
!     
      real*8 s(3,3),stress(6,nstep,*),xt(3,*),xn(3,*),xa(3,*),
     &     dk1(nstep,*),dk2(nstep,*),dk3(nstep,*),xkeq(nstep,*),
     &     phi(nstep,*),psi(nstep,*),pi,c2,c3,c4,term,
     &     acrack(*),ratio,constant,t(3),shape(3,*),c1
!     
      pi=4.d0*datan(1.d0)
      c2=70.d0*pi/180.d0
      c3=78.d0*pi/180.d0
      c4=33.d0*pi/180.d0
!     
      do i=1,nfront
!     
!     loop over all nodes belonging to the crack front(s)
!     
        noderel=ifrontrel(i)
!     
        do m=1,nstep
          s(1,1)=stress(1,m,noderel)
          s(1,2)=stress(4,m,noderel)
          s(1,3)=stress(6,m,noderel)
          s(2,1)=s(1,2)
          s(2,2)=stress(2,m,noderel)
          s(2,3)=stress(5,m,noderel)
          s(3,1)=s(1,3)
          s(3,2)=s(2,3)
          s(3,3)=stress(3,m,noderel)
!     
!     calculating the stress vector on the crack plane
!     
          t(1)=s(1,1)*xn(1,i)+s(1,2)*xn(2,i)+s(1,3)*xn(3,i)
          t(2)=s(2,1)*xn(1,i)+s(2,2)*xn(2,i)+s(2,3)*xn(3,i)
          t(3)=s(3,1)*xn(1,i)+s(3,2)*xn(2,i)+s(3,3)*xn(3,i)
!     
!     calculating the stress intensity factors
!     
          dk1(m,i)=t(1)*xn(1,i)+t(2)*xn(2,i)+t(3)*xn(3,i)
          dk2(m,i)=t(1)*xa(1,i)+t(2)*xa(2,i)+t(3)*xa(3,i)
          dk3(m,i)=t(1)*xt(1,i)+t(2)*xt(2,i)+t(3)*xt(3,i)
!     
!     taking the formula for the subsurface circular crack
!     
c          write(*,*) 'sh1 ',shape(1,i),shape(2,i),shape(3,i)
c          write(*,*)
          constant=dsqrt(pi*acrack(i))
          dk1(m,i)=dk1(m,i)*shape(1,i)*constant
          dk2(m,i)=dk2(m,i)*shape(2,i)*constant
          dk3(m,i)=dk3(m,i)*shape(3,i)*constant
        enddo
      enddo
!     
!     calculating the equivalent K-factor, the deflection angle
!     and twist angle (formulas by Hans Richard, University of Paderborn;
!     slightly modified to accomodate negative dk1 as well)
!     
      do i=1,nfront
        do m=1,nstep
          if(dk1(m,i).ge.0.d0) then
            xkeq(m,i)=(dk1(m,i)+dsqrt(dk1(m,i)*dk1(m,i)+
     &           5.3361*dk2(m,i)*dk2(m,i)+
     &           4.d0*dk3(m,i)*dk3(m,i)))/2.d0
            if(xkeq(m,i).gt.1.d-20) then
              term=dk1(m,i)+dabs(dk2(m,i))+dabs(dk3(m,i))
              ratio=dabs(dk2(m,i))/term
              phi(m,i)=-c2*ratio*(2.d0-ratio)*dk2(m,i)/dabs(dk2(m,i))
              ratio=dabs(dk3(m,i))/term
              psi(m,i)=-ratio*(c3-c4*ratio)*dk3(m,i)/dabs(dk3(m,i))
            else
              phi(m,i)=0.d0
              psi(m,i)=0.d0
            endif
          else
            xkeq(m,i)=-(-dk1(m,i)+dsqrt(dk1(m,i)*dk1(m,i)+
     &           5.3361*dk2(m,i)*dk2(m,i)+
     &           4.d0*dk3(m,i)*dk3(m,i)))/2.d0
!
!     for negative equivalent K-factor the formulas of Richard
!     probably do not apply: to be checked!
!
            phi(m,i)=0.d0
            psi(m,i)=0.d0
          endif
        enddo
      enddo
!     
      write(*,*) 'stressintensity k1 k2 k3'
      write(*,*)
      c1=1.d0/dsqrt(1000.d0)
      do m=1,nstep
        do i=1,nfront
          write(*,100) m,i,c1*dk1(m,i),c1*dk2(m,i),c1*dk3(m,i)
        enddo
        write(*,*)
      enddo
!     
      write(*,*) 'stressintensity keq phi psi'
      write(*,*)
      do m=1,nstep
        do i=1,nfront
          write(*,100) m,i,c1*xkeq(m,i),phi(m,i)*180.d0/pi,
     &         psi(m,i)*180.d0/pi
        enddo
        write(*,*)
      enddo
 100  format(2i10,3(1x,e11.4))
!     
      return
      end
