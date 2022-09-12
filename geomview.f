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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine geomview(vold,co,pmid,e1,e2,e3,kontri,area,
     &     cs,mcs,inocs,ntrit,nk,mi,sidemean)
!     
      implicit none
!
!     change following line if nlabel is increased
!     
      character*87 label(55)
!     
      integer i,j,l,mi(*),kontri(4,*),i1,mcs,inocs(*),i2,i3,
     &     ntrit,jj,is,m,nkt,icntrl,imag,nk,nlabel
!     
      real*8 vold(0:mi(2),*),co(3,*),
     &     pmid(3,*),e3(4,*),e1(3,*),e2(3,*),p1(3),p2(3),p3(3),
     &     cs(17,*),area(*),dd,p21(3),p32(3),
     &     fn,stn,qfn,een,t(3),sidemean,emn
!
      write(*,*) 'Calculating the viewfactors'
      write(*,*)
!     
!     change following line if nlabel is increased and the dimension
!     of field label above!
!     
      nlabel=55
!     
!     updating the displacements for cyclic symmetric structures
!     
      if(mcs.gt.0) then
         nkt=0
         do i=1,mcs
            if(int(cs(1,i)).gt.nkt) nkt=int(cs(1,i))
         enddo
         nkt=nk*nkt
         do i=1,nlabel
            do l=1,87
               label(i)(l:l)=' '
            enddo
         enddo
         label(1)(1:1)='U'
         imag=0
         icntrl=2
         call rectcyl(co,vold,fn,stn,qfn,een,cs,nk,icntrl,t,
     &        label,imag,mi,emn)
         
         do jj=0,mcs-1
           is=int(cs(1,jj+1))
c           write(*,*) 'geomview ',cs(1,jj+1)
c            is=cs(1,jj+1)
!     
            do i=1,is-1
               do l=1,nk
                  if(inocs(l).ne.jj) cycle
                  do m=1,mi(2)
                     vold(m,l+nk*i)=vold(m,l)
                  enddo
               enddo
            enddo
         enddo
         icntrl=-2
         call rectcyl(co,vold,fn,stn,qfn,een,cs,nkt,icntrl,t,
     &        label,imag,mi,emn)
      endif
!     
!     calculating the momentaneous center of the triangles,
!     area of the triangles and normal to the triangles
!     
      sidemean=0.d0
      do i=1,ntrit
         i1=kontri(1,i)
         if(i1.eq.0) cycle
         i2=kontri(2,i)
         i3=kontri(3,i)
         do j=1,3
            p1(j)=co(j,i1)+vold(j,i1)
            p2(j)=co(j,i2)+vold(j,i2)
            p3(j)=co(j,i3)+vold(j,i3)
            pmid(j,i)=(p1(j)+p2(j)+p3(j))/3.d0
            p21(j)=p2(j)-p1(j)
            p32(j)=p3(j)-p2(j)
         enddo
!     
!     normal to the triangle
!     
         e3(1,i)=p21(2)*p32(3)-p32(2)*p21(3)
         e3(2,i)=p21(3)*p32(1)-p32(3)*p21(1)
         e3(3,i)=p21(1)*p32(2)-p32(1)*p21(2)
!     
         dd=dsqrt(e3(1,i)*e3(1,i)+e3(2,i)*e3(2,i)+
     &        e3(3,i)*e3(3,i))
!     
!     check for degenerated triangles
!     
         if(dd.lt.1.d-20) then
            area(i)=0.d0
            cycle
         endif
!     
         do j=1,3
            e3(j,i)=e3(j,i)/dd
         enddo
!     
!     area of the triangle
!     
         area(i)=dd/2.d0
!     
!     unit vector parallel to side 1-2
!     
         dd=dsqrt(p21(1)*p21(1)+p21(2)*p21(2)+p21(3)*p21(3))
         sidemean=sidemean+dd
         do j=1,3
            e1(j,i)=p21(j)/dd
         enddo
!     
!     unit vector orthogonal to e1 and e3
!     
         e2(1,i)=e3(2,i)*e1(3,i)-e3(3,i)*e1(2,i)
         e2(2,i)=e3(3,i)*e1(1,i)-e3(1,i)*e1(3,i)
         e2(3,i)=e3(1,i)*e1(2,i)-e3(2,i)*e1(1,i)
!     
!     the fourth component in e3 is the constant term in the
!     equation of the triangle plane in the form
!     e3(1)*x+e3(2)*y+e3(3)*z+e3(4)=0
!     
         e3(4,i)=-(e3(1,i)*p1(1)+e3(2,i)*p1(2)
     &        +e3(3,i)*p1(3))
      enddo
      sidemean=sidemean/ntrit
!     
      return
      end
      
      
