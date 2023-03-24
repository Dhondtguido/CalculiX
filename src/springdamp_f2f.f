!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine springdamp_f2f(xl,elas,voldl,s,imat,elcon,
     &  ncmat_,ntmat_,nope,lakonl,iperturb,springarea,nmethod,mi,
     &  reltime,nasym,jfaces,igauss,pslavsurf,pmastsurf,clearini)
!
!     calculates the contact damping matrix (face-to-face penalty)
!
      implicit none
!
      character*8 lakonl
!
      integer i,j,imat,ncmat_,ntmat_,k,l,nope,iflag,iperturb(*),
     &  nmethod,mi(*),nasym,jfaces,igauss,nopem,nopes,nopep
!
      real*8 xl(3,19),elas(21),pproj(3),shp2m(7,9),al(3),s(60,60),
     &  voldl(0:mi(2),19),pl(3,19),xn(3),c3,elcon(0:ncmat_,ntmat_,*),
     &  xm(3),dval(3,19),fpu(3,3,19),xi,et,dal(3,3,19),xs2(3,7),xk,
     &  stickslope,springarea(2),tu(3,3,19),clear,reltime,
     &  xsj2s(3),xs2s(3,7),shp2s(7,9),weight,pslavsurf(3,*),
     &  pmastsurf(6,*),clearini(3,9,*)
!
      data iflag /1/
!     
!     # of master nodes
!
      read(lakonl(8:8),'(i1)') nopem
!
!     # of slave nodes
!
      nopes=nope-nopem
!
!     actual positions of the nodes belonging to the contact spring
!     (otherwise no contact force)
!
!     master nodes
!
      do i=1,nopem
         do j=1,3
            pl(j,i)=xl(j,i)+voldl(j,i)
         enddo
      enddo
!
!     slave nodes
!
      do i=nopem+1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+voldl(j,i)+clearini(j,i-nopem,jfaces)
     &             *reltime
         enddo
      enddo
!
!     contact springs
!
      read(lakonl(8:8),'(i1)') nopem
      nopes = nope - nopem
!
      xi=pslavsurf(1,igauss)
      et=pslavsurf(2,igauss)
      weight=pslavsurf(3,igauss)
!
      if(nopes.eq.8) then
          call shape8q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      elseif(nopes.eq.4) then
          call shape4q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      elseif(nopes.eq.6) then
          call shape6tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      else
          call shape3tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
      endif
!
      nopep=nope+1
!
!     position and displacements of the integration point in the
!     slave face
!
      do k=1,3
          pl(k,nopep)=0.d0
          voldl(k,nopep)=0.d0
          do j=1,nopes
              pl(k,nopep)=pl(k,nopep)+shp2s(4,j)*pl(k,nopem+j)
              voldl(k,nopep)=voldl(k,nopep)+shp2s(4,j)*voldl(k,nopem+j)
          enddo
      enddo   
!
      xi=pmastsurf(1,igauss)
      et=pmastsurf(2,igauss)
!
!     determining the jacobian vector on the surface 
!
      if(nopem.eq.8) then
         call shape8q(xi,et,pl,xm,xs2,shp2m,iflag)
      elseif(nopem.eq.4) then
         call shape4q(xi,et,pl,xm,xs2,shp2m,iflag)
      elseif(nopem.eq.6) then
         call shape6tri(xi,et,pl,xm,xs2,shp2m,iflag)
      else
         call shape3tri(xi,et,pl,xm,xs2,shp2m,iflag)
      endif
!
!     position of the projection of the slave integration point
!     on the master faces (done at the start of the increment)
!
      do i=1,3
         pproj(i)=0.d0
         do j=1,nopem
            pproj(i)=pproj(i)+shp2m(4,j)*pl(i,j)
         enddo
      enddo
!
!     vector connecting the integration point with its projection
!     on the master face
!
      do i=1,3
         al(i)=pl(i,nopep)-pproj(i)
      enddo
!
!     dal(i,j,k) is the derivative of al(i) w.r.t pl(j,k)
!     ( d al / d u_k)
!
      do k=1,nopem
         do i=1,3
            do j=1,3
               dal(i,j,k)=0.d0
            enddo
         enddo
      enddo
      do i=1,3
         do j=1,3
            dal(i,j,nopep)=0.d0
         enddo
      enddo
!
      do i=1,nopem
         do j=1,3
            dal(j,j,i)=-shp2m(4,i)
         enddo
      enddo
      do j=1,3
         dal(j,j,nopep)=1.d0
      enddo
!
!     normal vector on master face
!
      xn(1)=pmastsurf(4,igauss)
      xn(2)=pmastsurf(5,igauss)
      xn(3)=pmastsurf(6,igauss)
!
!     distance from surface along normal (= clearance)
!
      clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
      if(nmethod.eq.1) then
         clear=clear-springarea(2)*(1.d0-reltime)
      endif
!
      elas(2)=-springarea(1)*elcon(5,1,imat)
      elas(1)=elas(2)*clear
!     
!     derivatives of the size of the jacobian vector w.r.t. the
!     displacement vectors (d ||m||/d u_k)
!
      do k=1,nopem
!
!        auxiliary variable: (d clear d u_k)*||m||
!
         do i=1,3
            dval(i,k)=xn(1)*dal(1,i,k)+xn(2)*dal(2,i,k)+xn(3)*dal(3,i,k)
         enddo
!
      enddo
      do i=1,3
         dval(i,nopep)=xn(1)*dal(1,i,nopep)+xn(2)*dal(2,i,nopep)+
     &        xn(3)*dal(3,i,nopep)
      enddo
!
      c3=elas(2)
!
!     derivatives of the forces w.r.t. the displacement vectors
!
      do k=1,nopem
         do j=1,3
            do i=1,3
               fpu(i,j,k)=-c3*xn(i)*dval(j,k)
            enddo
         enddo
      enddo
      do j=1,3
         do i=1,3
            fpu(i,j,nopep)=-c3*xn(i)*dval(j,nopep)
         enddo
      enddo
!
!     tangential damping
!    
      if(elcon(8,1,imat).gt.0.d0) then
!     
         stickslope=elcon(8,1,imat)*elcon(5,1,imat)
!     
!        stiffness of shear stress versus slip curve
!     
         xk=stickslope*springarea(1)
!     
!       (d al/d u_k).||m|| -> (d al^*/d u_k).||m||
!     
         do k=1,nopem
            do i=1,3
               dval(i,k)=xn(1)*dal(1,i,k)+xn(2)*dal(2,i,k)
     &              +xn(3)*dal(3,i,k)
            enddo
         enddo
         do i=1,3
            dval(i,nopep)=xn(1)*dal(1,i,nopep)+xn(2)*dal(2,i,nopep)
     &           +xn(3)*dal(3,i,nopep)
         enddo
!     
!        d t/d u_k
!     
         do k=1,nopem
            do j=1,3
               do i=1,3
                  tu(i,j,k)=dal(i,j,k)-xn(i)*dval(j,k)
               enddo
            enddo
         enddo
         do j=1,3
            do i=1,3
               tu(i,j,nopep)=dal(i,j,nopep)-xn(i)*dval(j,nopep)
            enddo
         enddo
!     
!        damping matrix
!     
         do k=1,nopem
            do j=1,3
               do i=1,3
                  fpu(i,j,k)=fpu(i,j,k)+xk*tu(i,j,k)
               enddo
            enddo
         enddo  
         do j=1,3
            do i=1,3
               fpu(i,j,nopep)=fpu(i,j,nopep)+xk*tu(i,j,nopep)
            enddo
         enddo
      endif
!     
!     determining the stiffness matrix contributions
!     
!     dFkm/dUlm
!     
      do k=1,nopem
         do i=1,3
            do l=1,nopem
               do j=1,3
                  s(i+(k-1)*3,j+(l-1)*3)=-shp2m(4,k)*fpu(i,j,l)
               enddo
            enddo
         enddo
      enddo
!     
!     dFks/dUls
!     
      do k=nopem+1,nopem+nopes
         do i=1,3
            do l=nopem+1,nopem+nopes
               do j=1,3
                  s(i+(k-1)*3,j+(l-1)*3)=shp2s(4,k-nopem)*
     &                 shp2s(4,l-nopem)*fpu(i,j,nopep)
               enddo
            enddo
         enddo
      enddo
!     
!     dFkm/dUls
!     
      do k=1,nopem
         do i=1,3
            do l=nopem+1,nopem+nopes
               do j=1,3
                  s(i+(k-1)*3,j+(l-1)*3)=-shp2s(4,l-nopem)*
     &                shp2m(4,k)*fpu(i,j,nopep)
               enddo
            enddo
         enddo
      enddo
!
!     dFks/dUlm
!
      do k=nopem+1,nopem+nopes
         do i=1,3
            do l=1,nopem
               do j=1,3
                  s(i+(k-1)*3,j+(l-1)*3)=shp2s(4,k-nopem)*fpu(i,j,l)
               enddo
            enddo
         enddo
      enddo
!
!     symmetrizing the matrix
!     this is done in the absence of friction or for modal dynamic
!     calculations
!
      if((nasym.eq.0).or.((nmethod.eq.4).and.(iperturb(1).le.1))) then
         do j=1,3*nope
            do i=1,j-1
               s(i,j)=(s(i,j)+s(j,i))/2.d0
            enddo
         enddo
      endif
!
      return
      end

