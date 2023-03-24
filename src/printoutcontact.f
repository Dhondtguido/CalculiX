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
      subroutine printoutcontact(co,vold,lakon,ne0,ne,pslavsurf,stx,
     &  prset,ttime,nprint,prlab,mi,ipkon,kon,springarea,
     &  time,tieset,itiefac,ntie,pmastsurf)
!
!     calculation and printout of the contact forces
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*81 prset(*),tieset(3,*)
!
      integer ii,nprint,i,j,k,i1,nope,nopes,iposslave,ne0,ne,nopem,
     &  iflag,indexe,ipos,itie,ipkon(*),kon(*),mi(*),igauss,jfaces,
     &  node,istart,iend,ntie,itiefac(2,*)
!
      real*8 co(3,*),f(0:3),time,pslavsurf(3,*),vold(0:mi(2),*),xi,et,
     &  weight,ttime,coords(3),xm(3),df(3),cg(3),
     &  area,xn(3),xnormforc,shearforc,xntot(3),t1(3),t2(3),pl(3,16),
     &  xsj2s(3),xs2s(3,7),darea,dt1,springarea(2,*),shp2s(7,9),
     &  pmastsurf(6,*),stx(6,mi(1),*)
!
!
      include "gauss.f"
!
!     flag to check whether sof and/or som output was already
!     printed
!
      do ii=1,nprint
!
!        check whether there are contact print requests
!
!        CF: forces and moments on a contact slave surface
!
         if(prlab(ii)(1:2).eq.'CF')      
     &      then
!     
!     printing the header
!     
            write(5,*)
!
!              initialize total force and total moment
!
            do i=1,3
               f(i)=0.d0
               xm(i)=0.d0
               cg(i)=0.d0
               xntot(i)=0.d0
            enddo
            area=0.d0
!
!           identifying the appropriate tie
!
            read(prset(ii)(1:10),'(i10)') itie
            istart=itiefac(1,itie)
            iend=itiefac(2,itie)
!
!           loop over all contact elements
!
            do i=ne0+1,ne
!
!              check whether contact element
!               
               if((lakon(i)(1:1).ne.'E').or.
     &            (lakon(i)(7:7).ne.'C')) cycle
!
               indexe=ipkon(i)
               nope=kon(indexe)
!
!              check whether contact element belongs to correct slave surface
!
               jfaces=kon(indexe+nope+2)
               if((jfaces.lt.istart).or.(jfaces.gt.iend)) cycle
!
!              location of integration point information
!               
               igauss=kon(indexe+nope+1)
!
!              number of master and slave nodes
!               
               nopem=ichar(lakon(i)(8:8))-48
               nopes=nope-nopem
!
!              coordinates of slave nodes
!
               do k=nopem+1,nope
                  node=kon(indexe+k)
                  do j=1,3
                     pl(j,k)=co(j,node)+vold(j,node)
!                            +clearini(j,i-nopem,jfaces)*reltime ???
                  enddo
               enddo
!
!              location and weight of slave location
!               
               xi=pslavsurf(1,igauss)
               et=pslavsurf(2,igauss)
               weight=pslavsurf(3,igauss)
!
!              shape functions
!               
               iflag=1
               if(nopes.eq.8) then
                  call shape8q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &             iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &             iflag)
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &             iflag)
               else
                  call shape3tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &             iflag)
               endif
!               
!              global coordinates of the slave location
!
               do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                     coords(k)=coords(k)+shp2s(4,j)*pl(k,nopem+j)
                  enddo
               enddo
!
!              normal on the master face
!
               xn(1)=pmastsurf(4,igauss)
               xn(2)=pmastsurf(5,igauss)
               xn(3)=pmastsurf(6,igauss)
!     
!              calculating the local directions on master surface
!     
               if(1.d0 - dabs(xn(1)).lt.1.5231d-6) then       
                  t1(1)=-xn(3)*xn(1)
                  t1(2)=-xn(3)*xn(2)
                  t1(3)=1.d0-xn(3)*xn(3)
               else
                  t1(1)=1.d0-xn(1)*xn(1)
                  t1(2)=-xn(1)*xn(2)
                  t1(3)=-xn(1)*xn(3)
               endif
               dt1=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
               do i1=1,3
                  t1(i1)=t1(i1)/dt1
               enddo
               t2(1)=xn(2)*t1(3)-xn(3)*t1(2)
               t2(2)=xn(3)*t1(1)-xn(1)*t1(3)
               t2(3)=xn(1)*t1(2)-xn(2)*t1(1)
!
               do i1=1,3
                  xn(i1)=-xn(i1)
                  t1(i1)=-t1(i1)
                  t2(i1)=-t2(i1)
               enddo
!               
!              area corresponding to the slave location
!
               darea=springarea(1,igauss)
!
!              force
!
               if(prlab(ii)(1:3).eq.'CFN') then
                  do i1=1,3
                     df(i1)=(-stx(4,1,i)*xn(i1))*darea
                  enddo
               elseif(prlab(ii)(1:3).eq.'CFS') then
                  do i1=1,3
                     df(i1)=(stx(5,1,i)*t1(i1)
     &                      +stx(6,1,i)*t2(i1))*darea
                  enddo
               else
                  do i1=1,3
                     df(i1)=(-stx(4,1,i)*xn(i1)
     &                       +stx(5,1,i)*t1(i1)
     &                       +stx(6,1,i)*t2(i1))*darea
                  enddo
               endif
!
!              update total force and total moment
!
               do i1=1,3
                  f(i1)=f(i1)+df(i1)
               enddo
               xm(1)=xm(1)+coords(2)*df(3)-coords(3)*df(2)
               xm(2)=xm(2)+coords(3)*df(1)-coords(1)*df(3)
               xm(3)=xm(3)+coords(1)*df(2)-coords(2)*df(1)
!     
               area=area+darea
!               
               do i1=1,3
                  cg(i1)=cg(i1)+coords(i1)*darea
                  xntot(i1)=xntot(i1)+xn(i1)*darea
               enddo
!
            enddo
!
!              surface set statistics
!
            write(5,*)
            write(5,130) tieset(2,itie)(1:index(tieset(2,itie),' ')-2),
     &                   tieset(3,itie)(1:index(tieset(3,itie),' ')-2),
     &                   ttime+time
 130        format(' statistics for slave set ',A,', master set ',A,
     &           ' and time ',e14.7)
!
!           total force and moment about the origin
!
            write(5,*)
            if(prlab(ii)(1:2).eq.'CF') then     
              write(5,126)
            elseif(prlab(ii)(1:3).eq.'CFN') then     
              write(5,131)
            elseif(prlab(ii)(1:3).eq.'CFS') then
              write(5,132)
            endif
!     
 126        format('   total surface force (fx,fy,fz) ',
     &           'and moment about the origin (mx,my,mz)')
 131        format('   total normal surface force (fx,fy,fz) ',
     &           'and its moment about the origin (mx,my,mz)')
 132        format('   total shear surface force (fx,fy,fz) ',
     &           'and its moment about the origin (mx,my,mz)')
            write(5,*)
            write(5,'(2x,1p,6(1x,e13.6))') (f(j),j=1,3),(xm(j),j=1,3)
!
!           center of gravity and mean normal
!
            do i1=1,3
               cg(i1)=cg(i1)/area
               xntot(i1)=xntot(i1)/area
            enddo
            write(5,*)
            write(5,127)
 127        format('   center of gravity and mean normal')
            write(5,*)
            write(5,'(2x,1p,6(1x,e13.6))')(cg(j),j=1,3),(xntot(j),j=1,3)
!
!           moment about the center of gravity
!
            write(5,*)
            write(5,129)
 129        format(
     &           '   moment about the center of gravity(mx,my,mz)')
            write(5,*)
            write(5,'(2x,1p,6(1x,e13.6))') 
     &           xm(1)-cg(2)*f(3)+cg(3)*f(2),
     &           xm(2)-cg(3)*f(1)+cg(1)*f(3),
     &           xm(3)-cg(1)*f(2)+cg(2)*f(1)
!
!           area, normal and shear force
!
            xnormforc=f(1)*xntot(1)+f(2)*xntot(2)+f(3)*xntot(3)
            shearforc=sqrt((f(1)-xnormforc*xntot(1))**2+
     &           (f(2)-xnormforc*xntot(2))**2+
     &           (f(3)-xnormforc*xntot(3))**2)
            write(5,*)
            write(5,128)
 128        format('   area, ',
     &           ' normal force (+ = tension) and shear force (size)')
            write(5,*)
            write(5,'(2x,1p,3(1x,e13.6))') area,xnormforc,shearforc
!     
         endif
      enddo
!     
      return
      end
