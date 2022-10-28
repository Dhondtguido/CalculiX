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
      subroutine objective_stress_dx_dy(nodeset,istartset,iendset,
     &  ialset,nk,idesvar1c,idesvar2c,iobject,dgdx,ndesi,nobject,stn,
     &  dstn1,dstn2,dstn12,objectset,g0,dgdxdy)
!
!     calculates derivative of the sum of the square of the 
!     von Mises stress of a node set w.r.t. the coordinates of the mesh
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer nk,istartset(*),iendset(*),ialset(*),nodeset,idesvar2c,
     &  idesvar1c,iobject,j,k,ndesi,nobject,idesvar1,idesvar2
!
      real*8 dgdx(ndesi,nobject),stn(6,*),dstn1(6,*),g0(nobject),
     &  rho,xstress,svm,dsvmd1,p,dstn2(6,*),dstn12(6,*),dsvmd2,
     &  dsvmd12,dgdxdy(ndesi,ndesi,nobject),argument
!
      idesvar1=idesvar1c+1
      idesvar2=idesvar2c+1
!
!     reading rho and the mean stress for the Kreisselmeier-Steinhauser
!     function
!
      read(objectset(2,iobject)(41:60),'(f20.0)') rho
      read(objectset(2,iobject)(61:80),'(f20.0)') xstress
!
!     check for the existence of a set, else take the complete mesh
!
      dgdx(idesvar1,iobject)=0.d0
      dgdx(idesvar2,iobject)=0.d0
!
      if(nodeset.eq.0) then
         do j=1,nk
            p=-(stn(1,j)+stn(2,j)+stn(3,j))/3.d0
            svm=dsqrt(1.5d0*
     &        ((stn(1,j)+p)**2+(stn(2,j)+p)**2+(stn(3,j)+p)**2+
     &         2.d0*(stn(4,j)**2+stn(5,j)**2+stn(6,j)**2)))
            dsvmd1=
     &          ((2.d0*stn(1,j)-stn(2,j)-stn(3,j))*dstn1(1,j)+
     &           (2.d0*stn(2,j)-stn(1,j)-stn(3,j))*dstn1(2,j)+
     &           (2.d0*stn(3,j)-stn(1,j)-stn(2,j))*dstn1(3,j)+
     &           6.d0*(stn(4,j)*dstn1(4,j)+stn(5,j)*dstn1(5,j)+
     &                 stn(6,j)*dstn1(6,j)))/(2.d0*svm)
            argument=rho*svm/xstress
            if(argument.gt.600.d0) then
               write(*,*) '*ERROR in objective_stress: argument'
               write(*,*) '       ',argument,
     &               ' of exponential function is too big;'
               write(*,*) 
     &             '       choose other Kreisselmeier-Steinhauser'
               write(*,*) '       coefficients'
               call exit(201)
            endif
            dgdx(idesvar1,iobject)=dgdx(idesvar1,iobject)+
     &           dexp(argument)*dsvmd1
            dsvmd2=
     &          ((2.d0*stn(1,j)-stn(2,j)-stn(3,j))*dstn2(1,j)+
     &           (2.d0*stn(2,j)-stn(1,j)-stn(3,j))*dstn2(2,j)+
     &           (2.d0*stn(3,j)-stn(1,j)-stn(2,j))*dstn2(3,j)+
     &           6.d0*(stn(4,j)*dstn2(4,j)+stn(5,j)*dstn2(5,j)+
     &                 stn(6,j)*dstn2(6,j)))/(2.d0*svm)
            dgdx(idesvar2,iobject)=dgdx(idesvar2,iobject)+
     &           dexp(argument)*dsvmd2
            dsvmd12=-dsvmd1*dsvmd2/svm+
     &         ((2.d0*dstn2(1,j)-dstn2(2,j)-dstn2(3,j))*dstn1(1,j)+
     &          (2.d0*dstn2(2,j)-dstn2(3,j)-dstn2(1,j))*dstn1(2,j)+
     &          (2.d0*dstn2(3,j)-dstn2(1,j)-dstn2(2,j))*dstn1(3,j)+
     &          (2.d0*stn(1,j)-stn(2,j)-stn(3,j))*dstn12(1,j)+
     &          (2.d0*stn(2,j)-stn(3,j)-stn(1,j))*dstn12(1,j)+
     &          (2.d0*stn(3,j)-stn(1,j)-stn(2,j))*dstn12(1,j)+
     &          6.d0*(dstn1(4,j)*dstn2(4,j)+stn(4,j)*dstn12(4,j)+
     &                dstn1(5,j)*dstn2(5,j)+stn(5,j)*dstn12(5,j)+
     &                dstn1(6,j)+dstn2(6,j)+stn(6,j)*dstn12(6,j)))/
     &            (2.d0*svm)
            dgdxdy(idesvar1,idesvar2,iobject)=
     &        dgdxdy(idesvar1,idesvar2,iobject)+dexp(argument)*
     &        (dsvmd12+rho*dsvmd1*dsvmd2/xstress)
         enddo
         dgdx(idesvar1,iobject)=dgdx(idesvar1,iobject)*
     &         dexp(-rho*g0(iobject))/xstress
         dgdx(idesvar2,iobject)=dgdx(idesvar2,iobject)*
     &         dexp(-rho*g0(iobject))/xstress
         dgdxdy(idesvar1,idesvar2,iobject)=
     &      dgdxdy(idesvar1,idesvar2,iobject)
     &      *dexp(-rho*g0(iobject))/xstress
     &      -rho*dgdx(idesvar1,iobject)*dgdx(idesvar2,iobject)
      else
         do j=istartset(nodeset),iendset(nodeset)
            if(ialset(j).gt.0) then
               k=ialset(j)
               p=-(stn(1,k)+stn(2,k)+stn(3,k))/3.d0
               svm=dsqrt(1.5d0*
     &              ((stn(1,k)+p)**2+(stn(2,k)+p)**2+(stn(3,k)+p)**2+
     &              2.d0*(stn(4,k)**2+stn(5,k)**2+stn(6,k)**2)))
               dsvmd1=
     &             ((2.d0*stn(1,k)-stn(2,k)-stn(3,k))*dstn1(1,k)+
     &              (2.d0*stn(2,k)-stn(1,k)-stn(3,k))*dstn1(2,k)+
     &              (2.d0*stn(3,k)-stn(1,k)-stn(2,k))*dstn1(3,k)+
     &              6.d0*(stn(4,k)*dstn1(4,k)+stn(5,k)*dstn1(5,k)+
     &              stn(6,k)*dstn1(6,j)))/(2.d0*svm)
               argument=rho*svm/xstress
               if(argument.gt.600.d0) then
                  write(*,*) '*ERROR in objective_stress: argument'
                  write(*,*) '       ',argument,
     &               ' of exponential function is too big;'
                  write(*,*) 
     &             '       choose other Kreisselmeier-Steinhauser'
                  write(*,*) '       coefficients'
                  call exit(201)
               endif
               dgdx(idesvar1,iobject)=dgdx(idesvar1,iobject)+
     &              dexp(argument)*dsvmd1
               dsvmd2=
     &             ((2.d0*stn(1,k)-stn(2,k)-stn(3,k))*dstn2(1,k)+
     &              (2.d0*stn(2,k)-stn(1,k)-stn(3,k))*dstn2(2,k)+
     &              (2.d0*stn(3,k)-stn(1,k)-stn(2,k))*dstn2(3,k)+
     &              6.d0*(stn(4,k)*dstn2(4,k)+stn(5,k)*dstn2(5,k)+
     &              stn(6,k)*dstn2(6,j)))/(2.d0*svm)
               dgdx(idesvar2,iobject)=dgdx(idesvar2,iobject)+
     &              dexp(argument)*dsvmd2
               dsvmd12=-dsvmd1*dsvmd2/svm+
     &         ((2.d0*dstn2(1,j)-dstn2(2,j)-dstn2(3,j))*dstn1(1,j)+
     &          (2.d0*dstn2(2,j)-dstn2(3,j)-dstn2(1,j))*dstn1(2,j)+
     &          (2.d0*dstn2(3,j)-dstn2(1,j)-dstn2(2,j))*dstn1(3,j)+
     &          (2.d0*stn(1,j)-stn(2,j)-stn(3,j))*dstn12(1,j)+
     &          (2.d0*stn(2,j)-stn(3,j)-stn(1,j))*dstn12(1,j)+
     &          (2.d0*stn(3,j)-stn(1,j)-stn(2,j))*dstn12(1,j)+
     &          6.d0*(dstn1(4,j)*dstn2(4,j)+stn(4,j)*dstn12(4,j)+
     &                dstn1(5,j)*dstn2(5,j)+stn(5,j)*dstn12(5,j)+
     &                dstn1(6,j)+dstn2(6,j)+stn(6,j)*dstn12(6,j)))/
     &            (2.d0*svm)
               dgdxdy(idesvar1,idesvar2,iobject)=
     &          dgdxdy(idesvar1,idesvar2,iobject)+dexp(argument)*
     &          (dsvmd12+rho*dsvmd1*dsvmd2/xstress)
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  p=-(stn(1,k)+stn(2,k)+stn(3,k))/3.d0
                  svm=dsqrt(1.5d0*
     &                 ((stn(1,k)+p)**2+(stn(2,k)+p)**2+(stn(3,k)+p)**2+
     &                 2.d0*(stn(4,k)**2+stn(5,k)**2+stn(6,k)**2)))
                  dsvmd1=
     &                ((2.d0*stn(1,k)-stn(2,k)-stn(3,k))*dstn1(1,k)+
     &                 (2.d0*stn(2,k)-stn(1,k)-stn(3,k))*dstn1(2,k)+
     &                 (2.d0*stn(3,k)-stn(1,k)-stn(2,k))*dstn1(3,k)+
     &                 6.d0*(stn(4,k)*dstn1(4,k)+stn(5,k)*dstn1(5,k)+
     &                 stn(6,k)*dstn1(6,j)))/(2.d0*svm)
                  argument=rho*svm/xstress
                  if(argument.gt.600.d0) then
                     write(*,*) '*ERROR in objective_stress: argument'
                     write(*,*) '       ',argument,
     &               ' of exponential function is too big;'
                     write(*,*) 
     &             '       choose other Kreisselmeier-Steinhauser'
                     write(*,*) '       coefficients'
                     call exit(201)
                  endif
                  dgdx(idesvar1,iobject)=dgdx(idesvar1,iobject)+
     &                 dexp(argument)*dsvmd1
                  dsvmd2=
     &                ((2.d0*stn(1,k)-stn(2,k)-stn(3,k))*dstn2(1,k)+
     &                 (2.d0*stn(2,k)-stn(1,k)-stn(3,k))*dstn2(2,k)+
     &                 (2.d0*stn(3,k)-stn(1,k)-stn(2,k))*dstn2(3,k)+
     &                 6.d0*(stn(4,k)*dstn2(4,k)+stn(5,k)*dstn2(5,k)+
     &                 stn(6,k)*dstn2(6,j)))/(2.d0*svm)
                  dgdx(idesvar2,iobject)=dgdx(idesvar2,iobject)+
     &                 dexp(argument)*dsvmd2
                  dsvmd12=-dsvmd1*dsvmd2/svm+
     &            ((2.d0*dstn2(1,j)-dstn2(2,j)-dstn2(3,j))*dstn1(1,j)+
     &             (2.d0*dstn2(2,j)-dstn2(3,j)-dstn2(1,j))*dstn1(2,j)+
     &             (2.d0*dstn2(3,j)-dstn2(1,j)-dstn2(2,j))*dstn1(3,j)+
     &             (2.d0*stn(1,j)-stn(2,j)-stn(3,j))*dstn12(1,j)+
     &             (2.d0*stn(2,j)-stn(3,j)-stn(1,j))*dstn12(1,j)+
     &             (2.d0*stn(3,j)-stn(1,j)-stn(2,j))*dstn12(1,j)+
     &              6.d0*(dstn1(4,j)*dstn2(4,j)+stn(4,j)*dstn12(4,j)+
     &                dstn1(5,j)*dstn2(5,j)+stn(5,j)*dstn12(5,j)+
     &                dstn1(6,j)+dstn2(6,j)+stn(6,j)*dstn12(6,j)))/
     &             (2.d0*svm)
                  dgdxdy(idesvar1,idesvar2,iobject)=
     &                dgdxdy(idesvar1,idesvar2,iobject)+
     &                dexp(argument)*
     &                (dsvmd12+rho*dsvmd1*dsvmd2/xstress)
               enddo
            endif
         enddo
         dgdx(idesvar1,iobject)=dgdx(idesvar1,iobject)*
     &         dexp(-rho*g0(iobject))/xstress
         dgdx(idesvar2,iobject)=dgdx(idesvar2,iobject)*
     &         dexp(-rho*g0(iobject))/xstress
         dgdxdy(idesvar1,idesvar2,iobject)=
     &      dgdxdy(idesvar1,idesvar2,iobject)
     &      *dexp(-rho*g0(iobject))/xstress
     &      -rho*dgdx(idesvar1,iobject)*dgdx(idesvar2,iobject)
      endif
!     
      return
      end
      
