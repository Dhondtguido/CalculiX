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
      subroutine rectcyl(co,v,fn,stn,qfn,een,cs,n,icntrl,t,filab,
     &  imag,mi,emn)
!
!     icntrl=1:  rectangular to cylindrical coordinates for nodal
!                coordinates in field co
!     icntrl=-1: cylindrical to rectangular coordinates for nodal
!                coordinates in field co
!     icntrl=2:  rectangular to cylindrical coordinates for fields
!                v,fn,stn,een and emn
!     icntrl=-2: cylindrical to rectangular coordinates for fields
!                v,fn,stn, een and emn
!
!     the axis of the cylindrical coordinates is defined by points
!     a with coordinates csab(1..3) and b with coordinates csab(4..6).
!     Theta=0 (2nd cylindrical coordinate) is defined by the vector t,
!     which is perpendicular to the axis. The subroutine should be called
!     with icntrl=1 before calling it with icntrl=-1.
!
!     for icntrl=2 the imaginary part is extra taken into account if
!     imag=1
!
      implicit none
!
      character*87 filab(*)
      integer i,j,n,icntrl,imag,mi(*)
      real*8 co(3,*),v(0:mi(2),*),fn(0:mi(2),*),stn(6,*),een(6,*),
     &  a(3,3),emn(6,*),
     &  xr,xt,xz,b(3,3),cs(17,*),t(3),u(3),qfn(3,*),csab(7),
     &  xn(3),r(3),z,theta,rr,c(3,3),ctm,ct,st,ddx,ddy,dd
!
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
!
      if(icntrl.eq.1) then
!
!        normal along the cylindrical axis
!
         xn(1)=csab(4)-csab(1)
         xn(2)=csab(5)-csab(2)
         xn(3)=csab(6)-csab(3)
         dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
         do i=1,3
            xn(i)=xn(i)/dd
         enddo
!
!        normal to the cylindrical axis (vector t)
!
         if(dabs(xn(1)).gt.1.d-10) then
            t(2)=1.d0
            t(3)=0.d0
            t(1)=-xn(2)/xn(1)
         elseif(dabs(xn(2)).gt.1.d-10) then
            t(3)=1.d0
            t(1)=0.d0
            t(2)=-xn(3)/xn(2)
         else
            t(1)=1.d0
            t(2)=0.d0
            t(3)=-xn(1)/xn(3)
         endif
         dd=dsqrt(t(1)*t(1)+t(2)*t(2)+t(3)*t(3))
         do i=1,3
            t(i)=t(i)/dd
         enddo
!
!        normal to xn and t
!
         u(1)=xn(2)*t(3)-xn(3)*t(2)
         u(2)=-xn(1)*t(3)+xn(3)*t(1)
         u(3)=xn(1)*t(2)-xn(2)*t(1)
!
!        loop over all nodes to convert
!
         do i=1,n
            do j=1,3
               r(j)=co(j,i)-csab(j)
            enddo
            z=r(1)*xn(1)+r(2)*xn(2)+r(3)*xn(3)
            do j=1,3
               r(j)=r(j)-z*xn(j)
            enddo
            rr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            if(dabs(rr).lt.1.d-10) then
               theta=0.d0
            else
               do j=1,3
                  r(j)=r(j)/rr
               enddo
               ddx=t(1)*r(1)+t(2)*r(2)+t(3)*r(3)
               ddy=u(1)*r(1)+u(2)*r(2)+u(3)*r(3)
               theta=datan2(ddy,ddx)
            endif
            co(1,i)=rr
            co(2,i)=theta
            co(3,i)=z
         enddo
      elseif(icntrl.eq.-1) then
!
!        normal along the cylindrical axis
!
         xn(1)=csab(4)-csab(1)
         xn(2)=csab(5)-csab(2)
         xn(3)=csab(6)-csab(3)
         dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
         do i=1,3
            xn(i)=xn(i)/dd
         enddo
!
!        loop over all nodes to convert
!         
         do i=1,n
            rr=co(1,i)
            theta=co(2,i)
c            write(*,*) 'rectcyl',i,co(2,i)
            z=co(3,i)
            ct=dcos(theta)
            st=dsin(theta)
            ctm=1.d0-ct
!
!           rotation matrix
!
            c(1,1)=ct+ctm*xn(1)*xn(1)
            c(1,2)=-st*xn(3)+ctm*xn(1)*xn(2)
            c(1,3)=st*xn(2)+ctm*xn(1)*xn(3)
            c(2,1)=st*xn(3)+ctm*xn(2)*xn(1)
            c(2,2)=ct+ctm*xn(2)*xn(2)
            c(2,3)=-st*xn(1)+ctm*xn(2)*xn(3)
            c(3,1)=-st*xn(2)+ctm*xn(3)*xn(1)
            c(3,2)=st*xn(1)+ctm*xn(3)*xn(2)
            c(3,3)=ct+ctm*xn(3)*xn(3)
!
            co(1,i)=csab(1)+z*xn(1)+
     &              rr*(c(1,1)*t(1)+c(1,2)*t(2)+c(1,3)*t(3))
            co(2,i)=csab(2)+z*xn(2)+
     &              rr*(c(2,1)*t(1)+c(2,2)*t(2)+c(2,3)*t(3))
            co(3,i)=csab(3)+z*xn(3)+
     &              rr*(c(3,1)*t(1)+c(3,2)*t(2)+c(3,3)*t(3))
         enddo
      elseif(icntrl.eq.2) then
         do i=1,n
            j=i
            call transformatrix(csab,co(1,i),a)
!
            if((filab(1)(1:3).eq.'U  ').or.
     &         (filab(11)(1:4).eq.'PU'))  then 
               xr=v(1,j)*a(1,1)+v(2,j)*a(2,1)+v(3,j)*a(3,1)
               xt=v(1,j)*a(1,2)+v(2,j)*a(2,2)+v(3,j)*a(3,2)
               xz=v(1,j)*a(1,3)+v(2,j)*a(2,3)+v(3,j)*a(3,3)
               v(1,j)=xr
               v(2,j)=xt
               v(3,j)=xz
            endif
!
            if((filab(3)(1:4).eq.'S   ').or.
     &         (filab(18)(1:4).eq.'PHS ')) then
               b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(2,1)+stn(5,j)*a(3,1)
               b(1,2)=stn(1,j)*a(1,2)+stn(4,j)*a(2,2)+stn(5,j)*a(3,2)
               b(1,3)=stn(1,j)*a(1,3)+stn(4,j)*a(2,3)+stn(5,j)*a(3,3)
               b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(2,1)+stn(6,j)*a(3,1)
               b(2,2)=stn(4,j)*a(1,2)+stn(2,j)*a(2,2)+stn(6,j)*a(3,2)
               b(2,3)=stn(4,j)*a(1,3)+stn(2,j)*a(2,3)+stn(6,j)*a(3,3)
               b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(2,1)+stn(3,j)*a(3,1)
               b(3,2)=stn(5,j)*a(1,2)+stn(6,j)*a(2,2)+stn(3,j)*a(3,2)
               b(3,3)=stn(5,j)*a(1,3)+stn(6,j)*a(2,3)+stn(3,j)*a(3,3)
!
               stn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               stn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               stn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               stn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               stn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               stn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
            if(filab(4)(1:4).eq.'E   ') then
               b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(2,1)+een(5,j)*a(3,1)
               b(1,2)=een(1,j)*a(1,2)+een(4,j)*a(2,2)+een(5,j)*a(3,2)
               b(1,3)=een(1,j)*a(1,3)+een(4,j)*a(2,3)+een(5,j)*a(3,3)
               b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(2,1)+een(6,j)*a(3,1)
               b(2,2)=een(4,j)*a(1,2)+een(2,j)*a(2,2)+een(6,j)*a(3,2)
               b(2,3)=een(4,j)*a(1,3)+een(2,j)*a(2,3)+een(6,j)*a(3,3)
               b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(2,1)+een(3,j)*a(3,1)
               b(3,2)=een(5,j)*a(1,2)+een(6,j)*a(2,2)+een(3,j)*a(3,2)
               b(3,3)=een(5,j)*a(1,3)+een(6,j)*a(2,3)+een(3,j)*a(3,3)
!
               een(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               een(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               een(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               een(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               een(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               een(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
            if(filab(5)(1:4).eq.'RF  ') then
               xr=fn(1,j)*a(1,1)+fn(2,j)*a(2,1)+fn(3,j)*a(3,1)
               xt=fn(1,j)*a(1,2)+fn(2,j)*a(2,2)+fn(3,j)*a(3,2)
               xz=fn(1,j)*a(1,3)+fn(2,j)*a(2,3)+fn(3,j)*a(3,3)
               fn(1,j)=xr
               fn(2,j)=xt
               fn(3,j)=xz
            endif
!
            if(filab(9)(1:4).eq.'HFL ') then
               xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(2,1)+qfn(3,j)*a(3,1)
               xt=qfn(1,j)*a(1,2)+qfn(2,j)*a(2,2)+qfn(3,j)*a(3,2)
               xz=qfn(1,j)*a(1,3)+qfn(2,j)*a(2,3)+qfn(3,j)*a(3,3)
               qfn(1,j)=xr
               qfn(2,j)=xt
               qfn(3,j)=xz
            endif
!
            if(filab(32)(1:4).eq.'ME  ') then
               b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(2,1)+emn(5,j)*a(3,1)
               b(1,2)=emn(1,j)*a(1,2)+emn(4,j)*a(2,2)+emn(5,j)*a(3,2)
               b(1,3)=emn(1,j)*a(1,3)+emn(4,j)*a(2,3)+emn(5,j)*a(3,3)
               b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(2,1)+emn(6,j)*a(3,1)
               b(2,2)=emn(4,j)*a(1,2)+emn(2,j)*a(2,2)+emn(6,j)*a(3,2)
               b(2,3)=emn(4,j)*a(1,3)+emn(2,j)*a(2,3)+emn(6,j)*a(3,3)
               b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(2,1)+emn(3,j)*a(3,1)
               b(3,2)=emn(5,j)*a(1,2)+emn(6,j)*a(2,2)+emn(3,j)*a(3,2)
               b(3,3)=emn(5,j)*a(1,3)+emn(6,j)*a(2,3)+emn(3,j)*a(3,3)
!
               emn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
               emn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
               emn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
               emn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
               emn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
               emn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
            endif
!
!           imaginary part for cyclic symmetry frequency calculations
!
            if(imag.eq.1) then
!
               j=i+n
!
               if((filab(1)(1:3).eq.'U  ').or.
     &            (filab(11)(1:4).eq.'PU'))  then 
                  xr=v(1,j)*a(1,1)+v(2,j)*a(2,1)+v(3,j)*a(3,1)
                  xt=v(1,j)*a(1,2)+v(2,j)*a(2,2)+v(3,j)*a(3,2)
                  xz=v(1,j)*a(1,3)+v(2,j)*a(2,3)+v(3,j)*a(3,3)
                  v(1,j)=xr
                  v(2,j)=xt
                  v(3,j)=xz
               endif
!
               if((filab(3)(1:4).eq.'S   ').or.
     &            (filab(18)(1:4).eq.'PHS ')) then
                  b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(2,1)+stn(5,j)*a(3,1)
                  b(1,2)=stn(1,j)*a(1,2)+stn(4,j)*a(2,2)+stn(5,j)*a(3,2)
                  b(1,3)=stn(1,j)*a(1,3)+stn(4,j)*a(2,3)+stn(5,j)*a(3,3)
                  b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(2,1)+stn(6,j)*a(3,1)
                  b(2,2)=stn(4,j)*a(1,2)+stn(2,j)*a(2,2)+stn(6,j)*a(3,2)
                  b(2,3)=stn(4,j)*a(1,3)+stn(2,j)*a(2,3)+stn(6,j)*a(3,3)
                  b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(2,1)+stn(3,j)*a(3,1)
                  b(3,2)=stn(5,j)*a(1,2)+stn(6,j)*a(2,2)+stn(3,j)*a(3,2)
                  b(3,3)=stn(5,j)*a(1,3)+stn(6,j)*a(2,3)+stn(3,j)*a(3,3)
!
                  stn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  stn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  stn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  stn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  stn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  stn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
!
               if(filab(4)(1:4).eq.'E   ') then
                  b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(2,1)+een(5,j)*a(3,1)
                  b(1,2)=een(1,j)*a(1,2)+een(4,j)*a(2,2)+een(5,j)*a(3,2)
                  b(1,3)=een(1,j)*a(1,3)+een(4,j)*a(2,3)+een(5,j)*a(3,3)
                  b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(2,1)+een(6,j)*a(3,1)
                  b(2,2)=een(4,j)*a(1,2)+een(2,j)*a(2,2)+een(6,j)*a(3,2)
                  b(2,3)=een(4,j)*a(1,3)+een(2,j)*a(2,3)+een(6,j)*a(3,3)
                  b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(2,1)+een(3,j)*a(3,1)
                  b(3,2)=een(5,j)*a(1,2)+een(6,j)*a(2,2)+een(3,j)*a(3,2)
                  b(3,3)=een(5,j)*a(1,3)+een(6,j)*a(2,3)+een(3,j)*a(3,3)
!
                  een(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  een(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  een(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  een(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  een(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  een(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
!
               if(filab(5)(1:4).eq.'RF  ') then
                  xr=fn(1,j)*a(1,1)+fn(2,j)*a(2,1)+fn(3,j)*a(3,1)
                  xt=fn(1,j)*a(1,2)+fn(2,j)*a(2,2)+fn(3,j)*a(3,2)
                  xz=fn(1,j)*a(1,3)+fn(2,j)*a(2,3)+fn(3,j)*a(3,3)
                  fn(1,j)=xr
                  fn(2,j)=xt
                  fn(3,j)=xz
               endif
!
               if(filab(9)(1:4).eq.'HFL ') then
                  xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(2,1)+qfn(3,j)*a(3,1)
                  xt=qfn(1,j)*a(1,2)+qfn(2,j)*a(2,2)+qfn(3,j)*a(3,2)
                  xz=qfn(1,j)*a(1,3)+qfn(2,j)*a(2,3)+qfn(3,j)*a(3,3)
                  qfn(1,j)=xr
                  qfn(2,j)=xt
                  qfn(3,j)=xz
               endif
!
               if(filab(32)(1:4).eq.'ME  ') then
                  b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(2,1)+emn(5,j)*a(3,1)
                  b(1,2)=emn(1,j)*a(1,2)+emn(4,j)*a(2,2)+emn(5,j)*a(3,2)
                  b(1,3)=emn(1,j)*a(1,3)+emn(4,j)*a(2,3)+emn(5,j)*a(3,3)
                  b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(2,1)+emn(6,j)*a(3,1)
                  b(2,2)=emn(4,j)*a(1,2)+emn(2,j)*a(2,2)+emn(6,j)*a(3,2)
                  b(2,3)=emn(4,j)*a(1,3)+emn(2,j)*a(2,3)+emn(6,j)*a(3,3)
                  b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(2,1)+emn(3,j)*a(3,1)
                  b(3,2)=emn(5,j)*a(1,2)+emn(6,j)*a(2,2)+emn(3,j)*a(3,2)
                  b(3,3)=emn(5,j)*a(1,3)+emn(6,j)*a(2,3)+emn(3,j)*a(3,3)
!
                  emn(1,j)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
                  emn(2,j)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
                  emn(3,j)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
                  emn(4,j)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
                  emn(5,j)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
                  emn(6,j)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
               endif
            endif
         enddo
      elseif(icntrl.eq.-2) then
         do i=1,n
            j=i
            call transformatrix(csab,co(1,i),a)
!
            if((filab(1)(1:3).eq.'U  ').or.
     &         (filab(11)(1:4).eq.'PU'))  then 
               xr=v(1,j)*a(1,1)+v(2,j)*a(1,2)+v(3,j)*a(1,3)
               xt=v(1,j)*a(2,1)+v(2,j)*a(2,2)+v(3,j)*a(2,3)
               xz=v(1,j)*a(3,1)+v(2,j)*a(3,2)+v(3,j)*a(3,3)
               v(1,j)=xr
               v(2,j)=xt
               v(3,j)=xz
            endif
!
            if((filab(3)(1:4).eq.'S   ').or.
     &         (filab(18)(1:4).eq.'PHS ')) then
               b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(1,2)+stn(5,j)*a(1,3)
               b(1,2)=stn(1,j)*a(2,1)+stn(4,j)*a(2,2)+stn(5,j)*a(2,3)
               b(1,3)=stn(1,j)*a(3,1)+stn(4,j)*a(3,2)+stn(5,j)*a(3,3)
               b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(1,2)+stn(6,j)*a(1,3)
               b(2,2)=stn(4,j)*a(2,1)+stn(2,j)*a(2,2)+stn(6,j)*a(2,3)
               b(2,3)=stn(4,j)*a(3,1)+stn(2,j)*a(3,2)+stn(6,j)*a(3,3)
               b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(1,2)+stn(3,j)*a(1,3)
               b(3,2)=stn(5,j)*a(2,1)+stn(6,j)*a(2,2)+stn(3,j)*a(2,3)
               b(3,3)=stn(5,j)*a(3,1)+stn(6,j)*a(3,2)+stn(3,j)*a(3,3)
!
               stn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
               stn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
               stn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
               stn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
               stn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
               stn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
            endif
!
            if(filab(4)(1:4).eq.'E   ') then
               b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(1,2)+een(5,j)*a(1,3)
               b(1,2)=een(1,j)*a(2,1)+een(4,j)*a(2,2)+een(5,j)*a(2,3)
               b(1,3)=een(1,j)*a(3,1)+een(4,j)*a(3,2)+een(5,j)*a(3,3)
               b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(1,2)+een(6,j)*a(1,3)
               b(2,2)=een(4,j)*a(2,1)+een(2,j)*a(2,2)+een(6,j)*a(2,3)
               b(2,3)=een(4,j)*a(3,1)+een(2,j)*a(3,2)+een(6,j)*a(3,3)
               b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(1,2)+een(3,j)*a(1,3)
               b(3,2)=een(5,j)*a(2,1)+een(6,j)*a(2,2)+een(3,j)*a(2,3)
               b(3,3)=een(5,j)*a(3,1)+een(6,j)*a(3,2)+een(3,j)*a(3,3)
!
               een(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
               een(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
               een(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
               een(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
               een(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
               een(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
            endif
!
            if(filab(5)(1:4).eq.'RF  ') then
               xr=fn(1,j)*a(1,1)+fn(2,j)*a(1,2)+fn(3,j)*a(1,3)
               xt=fn(1,j)*a(2,1)+fn(2,j)*a(2,2)+fn(3,j)*a(2,3)
               xz=fn(1,j)*a(3,1)+fn(2,j)*a(3,2)+fn(3,j)*a(3,3)
               fn(1,j)=xr
               fn(2,j)=xt
               fn(3,j)=xz
            endif
!
            if(filab(9)(1:4).eq.'HFL ') then
               xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(1,2)+qfn(3,j)*a(1,3)
               xt=qfn(1,j)*a(2,1)+qfn(2,j)*a(2,2)+qfn(3,j)*a(2,3)
               xz=qfn(1,j)*a(3,1)+qfn(2,j)*a(3,2)+qfn(3,j)*a(3,3)
               qfn(1,j)=xr
               qfn(2,j)=xt
               qfn(3,j)=xz
            endif
!
            if(filab(32)(1:4).eq.'ME  ') then
               b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(1,2)+emn(5,j)*a(1,3)
               b(1,2)=emn(1,j)*a(2,1)+emn(4,j)*a(2,2)+emn(5,j)*a(2,3)
               b(1,3)=emn(1,j)*a(3,1)+emn(4,j)*a(3,2)+emn(5,j)*a(3,3)
               b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(1,2)+emn(6,j)*a(1,3)
               b(2,2)=emn(4,j)*a(2,1)+emn(2,j)*a(2,2)+emn(6,j)*a(2,3)
               b(2,3)=emn(4,j)*a(3,1)+emn(2,j)*a(3,2)+emn(6,j)*a(3,3)
               b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(1,2)+emn(3,j)*a(1,3)
               b(3,2)=emn(5,j)*a(2,1)+emn(6,j)*a(2,2)+emn(3,j)*a(2,3)
               b(3,3)=emn(5,j)*a(3,1)+emn(6,j)*a(3,2)+emn(3,j)*a(3,3)
!
               emn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
               emn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
               emn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
               emn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
               emn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
               emn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
            endif
!
!           imaginary part for cyclic symmetry frequency calculations
!
            if(imag.eq.1) then
!
               j=i+n
!
               if((filab(1)(1:3).eq.'U  ').or.
     &            (filab(11)(1:4).eq.'PU'))  then 
                  xr=v(1,j)*a(1,1)+v(2,j)*a(1,2)+v(3,j)*a(1,3)
                  xt=v(1,j)*a(2,1)+v(2,j)*a(2,2)+v(3,j)*a(2,3)
                  xz=v(1,j)*a(3,1)+v(2,j)*a(3,2)+v(3,j)*a(3,3)
                  v(1,j)=xr
                  v(2,j)=xt
                  v(3,j)=xz
               endif
!     
               if((filab(3)(1:4).eq.'S   ').or.
     &            (filab(18)(1:4).eq.'PHS ')) then
                  b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(1,2)+stn(5,j)*a(1,3)
                  b(1,2)=stn(1,j)*a(2,1)+stn(4,j)*a(2,2)+stn(5,j)*a(2,3)
                  b(1,3)=stn(1,j)*a(3,1)+stn(4,j)*a(3,2)+stn(5,j)*a(3,3)
                  b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(1,2)+stn(6,j)*a(1,3)
                  b(2,2)=stn(4,j)*a(2,1)+stn(2,j)*a(2,2)+stn(6,j)*a(2,3)
                  b(2,3)=stn(4,j)*a(3,1)+stn(2,j)*a(3,2)+stn(6,j)*a(3,3)
                  b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(1,2)+stn(3,j)*a(1,3)
                  b(3,2)=stn(5,j)*a(2,1)+stn(6,j)*a(2,2)+stn(3,j)*a(2,3)
                  b(3,3)=stn(5,j)*a(3,1)+stn(6,j)*a(3,2)+stn(3,j)*a(3,3)
!     
                  stn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  stn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  stn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  stn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  stn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  stn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
!     
               if(filab(4)(1:4).eq.'E   ') then
                  b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(1,2)+een(5,j)*a(1,3)
                  b(1,2)=een(1,j)*a(2,1)+een(4,j)*a(2,2)+een(5,j)*a(2,3)
                  b(1,3)=een(1,j)*a(3,1)+een(4,j)*a(3,2)+een(5,j)*a(3,3)
                  b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(1,2)+een(6,j)*a(1,3)
                  b(2,2)=een(4,j)*a(2,1)+een(2,j)*a(2,2)+een(6,j)*a(2,3)
                  b(2,3)=een(4,j)*a(3,1)+een(2,j)*a(3,2)+een(6,j)*a(3,3)
                  b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(1,2)+een(3,j)*a(1,3)
                  b(3,2)=een(5,j)*a(2,1)+een(6,j)*a(2,2)+een(3,j)*a(2,3)
                  b(3,3)=een(5,j)*a(3,1)+een(6,j)*a(3,2)+een(3,j)*a(3,3)
!     
                  een(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  een(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  een(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  een(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  een(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  een(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
!     
               if(filab(5)(1:4).eq.'RF  ') then
                  xr=fn(1,j)*a(1,1)+fn(2,j)*a(1,2)+fn(3,j)*a(1,3)
                  xt=fn(1,j)*a(2,1)+fn(2,j)*a(2,2)+fn(3,j)*a(2,3)
                  xz=fn(1,j)*a(3,1)+fn(2,j)*a(3,2)+fn(3,j)*a(3,3)
                  fn(1,j)=xr
                  fn(2,j)=xt
                  fn(3,j)=xz
               endif
!     
               if(filab(9)(1:4).eq.'HFL ') then
                  xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(1,2)+qfn(3,j)*a(1,3)
                  xt=qfn(1,j)*a(2,1)+qfn(2,j)*a(2,2)+qfn(3,j)*a(2,3)
                  xz=qfn(1,j)*a(3,1)+qfn(2,j)*a(3,2)+qfn(3,j)*a(3,3)
                  qfn(1,j)=xr
                  qfn(2,j)=xt
                  qfn(3,j)=xz
               endif
!     
               if(filab(32)(1:4).eq.'ME  ') then
                  b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(1,2)+emn(5,j)*a(1,3)
                  b(1,2)=emn(1,j)*a(2,1)+emn(4,j)*a(2,2)+emn(5,j)*a(2,3)
                  b(1,3)=emn(1,j)*a(3,1)+emn(4,j)*a(3,2)+emn(5,j)*a(3,3)
                  b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(1,2)+emn(6,j)*a(1,3)
                  b(2,2)=emn(4,j)*a(2,1)+emn(2,j)*a(2,2)+emn(6,j)*a(2,3)
                  b(2,3)=emn(4,j)*a(3,1)+emn(2,j)*a(3,2)+emn(6,j)*a(3,3)
                  b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(1,2)+emn(3,j)*a(1,3)
                  b(3,2)=emn(5,j)*a(2,1)+emn(6,j)*a(2,2)+emn(3,j)*a(2,3)
                  b(3,3)=emn(5,j)*a(3,1)+emn(6,j)*a(3,2)+emn(3,j)*a(3,3)
!     
                  emn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
                  emn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
                  emn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
                  emn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
                  emn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
                  emn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
               endif
            endif
!
         enddo
      endif
!
      return
      end













