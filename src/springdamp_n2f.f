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
      subroutine springdamp_n2f(xl,elas,voldl,s,imat,elcon,
     &  ncmat_,ntmat_,nope,iperturb,springarea,nmethod,mi,
     &  reltime,nasym)
!
!     calculates the contact damping matrix (node-to-face penalty)
!     (User's manual -> theory -> boundary conditions -> 
!      node-to-face penalty contact)
!
      implicit none
!
      integer i,j,imat,ncmat_,ntmat_,k,l,nope,nterms,iflag,i1,
     &  iperturb(*),nmethod,mi(*),nasym
!
      real*8 xl(3,10),elas(21),ratio(9),pproj(3),val,shp2(7,9),
     &  al(3),s(60,60),voldl(0:mi(2),10),pl(3,10),xn(3),dm,
     &  c1,c2,c3,c4,elcon(0:ncmat_,ntmat_,*),xm(3),xmu(3,3,10),
     &  dxmu(3,10),dval(3,10),fpu(3,3,10),xi,et,xs2(3,7),xk,
     &  a11,a12,a22,b1(3,10),b2(3,10),dal(3,3,10),qxxy(3),fnl(3),
     &  qxyy(3),dxi(3,10),det(3,10),determinant,c11,c12,c22,
     &  qxyx(3),qyxy(3),springarea(2),dist,tu(3,3,10),
     &  clear,reltime,alnew(3)
!
      data iflag /4/
!
!     actual positions of the nodes belonging to the contact spring
!     (otherwise no contact force)
!
      do i=1,nope
         do j=1,3
            pl(j,i)=xl(j,i)+voldl(j,i)
         enddo
      enddo
!
!     contact springs
!
      nterms=nope-1
!
!     vector al connects the actual position of the slave node 
!     with its projection on the master face = vec_r (User's
!     manual -> theory -> boundary conditions -> node-to-face
!     penalty contact)
!
      do i=1,3
         pproj(i)=pl(i,nope)
      enddo
      call attach_2d(pl,pproj,nterms,ratio,dist,xi,et)
      do i=1,3
         al(i)=pl(i,nope)-pproj(i)
      enddo
!
!     determining the jacobian vector on the surface 
!
      if(nterms.eq.8) then
         call shape8q(xi,et,pl,xm,xs2,shp2,iflag)
      elseif(nterms.eq.4) then
         call shape4q(xi,et,pl,xm,xs2,shp2,iflag)
      elseif(nterms.eq.6) then
         call shape6tri(xi,et,pl,xm,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,pl,xm,xs2,shp2,iflag)
      endif
!
!     d xi / d vec_u_j
!     d eta / d vec_u_j
!
!     dxi and det are determined from the orthogonality
!     condition
!
      a11=-(xs2(1,1)*xs2(1,1)+xs2(2,1)*xs2(2,1)+xs2(3,1)*xs2(3,1))
     &    +al(1)*xs2(1,5)+al(2)*xs2(2,5)+al(3)*xs2(3,5)
      a12=-(xs2(1,1)*xs2(1,2)+xs2(2,1)*xs2(2,2)+xs2(3,1)*xs2(3,2))
     &    +al(1)*xs2(1,6)+al(2)*xs2(2,6)+al(3)*xs2(3,6)
      a22=-(xs2(1,2)*xs2(1,2)+xs2(2,2)*xs2(2,2)+xs2(3,2)*xs2(3,2))
     &    +al(1)*xs2(1,7)+al(2)*xs2(2,7)+al(3)*xs2(3,7)
!
      do i=1,3
         do j=1,nterms
            b1(i,j)=shp2(4,j)*xs2(i,1)-shp2(1,j)*al(i)
            b2(i,j)=shp2(4,j)*xs2(i,2)-shp2(2,j)*al(i)
         enddo
         b1(i,nope)=-xs2(i,1)
         b2(i,nope)=-xs2(i,2)
      enddo
!
      determinant=a11*a22-a12*a12
      c11=a22/determinant
      c12=-a12/determinant
      c22=a11/determinant
!
      do i=1,3
         do j=1,nope
            dxi(i,j)=c11*b1(i,j)+c12*b2(i,j)
            det(i,j)=c12*b1(i,j)+c22*b2(i,j)
         enddo
      enddo
!
!     d vec_r / d vec_u_k
!
      do i=1,nope
         do j=1,3
            do k=1,3
               dal(j,k,i)=-xs2(j,1)*dxi(k,i)-xs2(j,2)*det(k,i)
            enddo
         enddo
      enddo
      do i=1,nterms
         do j=1,3
            dal(j,j,i)=dal(j,j,i)-shp2(4,i)
         enddo
      enddo
      do j=1,3
         dal(j,j,nope)=dal(j,j,nope)+1.d0
      enddo
!
!     d2q/dxx x dq/dy
!
      qxxy(1)=xs2(2,5)*xs2(3,2)-xs2(3,5)*xs2(2,2)
      qxxy(2)=xs2(3,5)*xs2(1,2)-xs2(1,5)*xs2(3,2)
      qxxy(3)=xs2(1,5)*xs2(2,2)-xs2(2,5)*xs2(1,2)
!
!     dq/dx x d2q/dyy
!
      qxyy(1)=xs2(2,1)*xs2(3,7)-xs2(3,1)*xs2(2,7)
      qxyy(2)=xs2(3,1)*xs2(1,7)-xs2(1,1)*xs2(3,7)
      qxyy(3)=xs2(1,1)*xs2(2,7)-xs2(2,1)*xs2(1,7)
!
!     Modified by Stefan Sicklinger
!
!     dq/dx x d2q/dxy
!
      qxyx(1)=xs2(2,1)*xs2(3,6)-xs2(3,1)*xs2(2,6)
      qxyx(2)=xs2(3,1)*xs2(1,6)-xs2(1,1)*xs2(3,6)
      qxyx(3)=xs2(1,1)*xs2(2,6)-xs2(2,1)*xs2(1,6)
!
!     d2q/dxy x dq/dy
!
      qyxy(1)=xs2(2,6)*xs2(3,2)-xs2(3,6)*xs2(2,2)
      qyxy(2)=xs2(3,6)*xs2(1,2)-xs2(1,6)*xs2(3,2)
      qyxy(3)=xs2(1,6)*xs2(2,2)-xs2(2,6)*xs2(1,2)
!
!
!     End modifications
!
!     normal on the surface
!
      dm=dsqrt(xm(1)*xm(1)+xm(2)*xm(2)+xm(3)*xm(3))
      do i=1,3
         xn(i)=xm(i)/dm
      enddo
!
!     distance from surface along normal (= clearance)
!
      clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
      if(nmethod.eq.1) then
         clear=clear-springarea(2)*(1.d0-reltime)
      endif
!
!     representative area: usually the slave surface stored in
!     springarea; however, if no area was assigned because the
!     node does not belong to any element, the master surface
!     is used
!
      if(springarea(1).le.0.d0) then
         if(nterms.eq.3) then
            springarea(1)=dm/2.d0
         else
            springarea(1)=dm*4.d0
         endif
      endif
!
      elas(2)=-springarea(1)*elcon(5,1,imat)
      elas(1)=elas(2)*clear
!
!     contact force
!
      do i=1,3
         fnl(i)=-elas(1)*xn(i)
      enddo
!
!     derivatives of the jacobian vector w.r.t. the displacement
!     vectors (d m / d u_k)
!
      do k=1,nterms
         xmu(1,1,k)=0.d0
         xmu(2,2,k)=0.d0
         xmu(3,3,k)=0.d0
         xmu(1,2,k)=shp2(1,k)*xs2(3,2)-shp2(2,k)*xs2(3,1)
         xmu(2,3,k)=shp2(1,k)*xs2(1,2)-shp2(2,k)*xs2(1,1)
         xmu(3,1,k)=shp2(1,k)*xs2(2,2)-shp2(2,k)*xs2(2,1)
         xmu(1,3,k)=-xmu(3,1,k)
         xmu(2,1,k)=-xmu(1,2,k)
         xmu(3,2,k)=-xmu(2,3,k)
      enddo
      do i=1,3
         do j=1,3
            xmu(i,j,nope)=0.d0
         enddo
      enddo
!
!     correction due to change of xi and eta
!
      do k=1,nope
         do j=1,3
            do i=1,3
!
!     modified by Stefan Sicklinger
!
               xmu(i,j,k)=xmu(i,j,k)+(qxxy(i)+qxyx(i))*dxi(j,k)
     &                                +(qxyy(i)+qyxy(i))*det(j,k)
            enddo
         enddo
      enddo
!     
!     derivatives of the size of the jacobian vector w.r.t. the
!     displacement vectors (d ||m||/d u_k)
!
      do k=1,nope
         do i=1,3
            dxmu(i,k)=xn(1)*xmu(1,i,k)+xn(2)*xmu(2,i,k)+
     &           xn(3)*xmu(3,i,k)
         enddo
!
!        auxiliary variable: (d clear d u_k)*||m||
!
         do i=1,3
            dval(i,k)=al(1)*xmu(1,i,k)+al(2)*xmu(2,i,k)+
     &               al(3)*xmu(3,i,k)-clear*dxmu(i,k)+
     &               xm(1)*dal(1,i,k)+xm(2)*dal(2,i,k)+xm(3)*dal(3,i,k)
         enddo
!
      enddo
!
      c1=1.d0/dm
      c2=c1*c1
      c3=elas(2)*c2
      c4=elas(1)*c1
!
!     derivatives of the forces w.r.t. the displacement vectors
!
      do k=1,nope
         do j=1,3
            do i=1,3
               fpu(i,j,k)=-c3*xm(i)*dval(j,k)
     &                    +c4*(xn(i)*dxmu(j,k)-xmu(i,j,k))
            enddo
         enddo
      enddo    
!
!     tangential damping
!    
      if(elcon(8,1,imat).gt.0.d0) then
!     
!     stiffness of shear stress versus slip curve
!     
         xk=elcon(8,1,imat)*elcon(5,1,imat)*springarea(1)
!     
!     calculating the relative displacement between the slave node
!     and its projection on the master surface
!     
         do i=1,3
            alnew(i)=voldl(i,nope)
            do j=1,nterms
               alnew(i)=alnew(i)-ratio(j)*voldl(i,j)
            enddo
         enddo
!     
!     calculating the difference in relative displacement since
!     the start of the increment = vec_s
!     
         do i=1,3
            al(i)=alnew(i)
         enddo
!     
!     s=||vec_s||
!     
         val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!     
!     d vec_s/ d vec_u_k
!     notice: xi & et are const.
!     
         do k=1,nope
            do i=1,3
               do j=1,3
                  dal(i,j,k)=0.d0
               enddo
            enddo
         enddo
!     
         do i=1,nterms
            do j=1,3
               dal(j,j,i)=-shp2(4,i)
            enddo
         enddo
!     
         do j=1,3
            dal(j,j,nope)=1.d0
         enddo
!     
!     d s/ dvec_u_k.||m||
!     
         do k=1,nope
            do i=1,3
               dval(i,k)=al(1)*xmu(1,i,k)+al(2)*xmu(2,i,k)
     &              +al(3)*xmu(3,i,k)-val*dxmu(i,k)
     &              +xm(1)*dal(1,i,k)+xm(2)*dal(2,i,k)
     &              +xm(3)*dal(3,i,k)
            enddo
         enddo
!     
!     d vec_t/d vec_u_k
!     
         do k=1,nope
            do j=1,3
               do i=1,3
                  tu(i,j,k)=dal(i,j,k)
     &                 -c1*(xn(i)*(dval(j,k)-val*dxmu(j,k))
     &                 +val*xmu(i,j,k))
               enddo
            enddo
         enddo
!     
!     damping matrix
!     
         do k=1,nope
            do j=1,3
               do i=1,3
                  fpu(i,j,k)=fpu(i,j,k)+xk*tu(i,j,k)
               enddo
            enddo
         enddo  
      endif
!     
!     determining the damping matrix contributions
!     
!     complete field shp2 
!
      shp2(1,nope)=0.d0
      shp2(2,nope)=0.d0
      shp2(4,nope)=-1.d0
!
      do k=1,nope
         do l=1,nope
            do i=1,3
               i1=i+(k-1)*3
               do j=1,3
                  s(i1,j+(l-1)*3)=-shp2(4,k)*fpu(i,j,l)
     &                            -shp2(1,k)*fnl(i)*dxi(j,l)
     &                            -shp2(2,k)*fnl(i)*det(j,l)
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

