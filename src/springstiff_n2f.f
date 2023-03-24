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
      subroutine springstiff_n2f(xl,elas,konl,voldl,s,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,plicon,
     &  nplicon,npmat_,iperturb,springarea,nmethod,mi,ne0,
     &  nstate_,xstateini,xstate,reltime,nasym,ielorien,orab,
     &  norien,nelem)
!
!     calculates the stiffness of a spring (node-to-face penalty)
!     (User's manual -> theory -> boundary conditions -> 
!      node-to-face penalty contact)
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),i,j,imat,ncmat_,ntmat_,k,l,nope,nterms,iflag,
     &  i1,kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),
     &  iperturb(*),nmethod,mi(*),ne0,nstate_,nasym,ielorien(mi(3),*),
     &  norien,nelem,idof,idof1,idof2,iorien
!
      real*8 xl(3,10),elas(21),ratio(9),pproj(3),val,shp2(7,9),
     &  al(3),s(60,60),voldl(0:mi(2),10),pl(3,10),xn(3),dm,
     &  c1,c2,c3,c4,alpha,beta,elcon(0:ncmat_,ntmat_,*),xm(3),
     &  xmu(3,3,10),dxmu(3,10),dval(3,10),fpu(3,3,10),xi,et,
     &  xs2(3,7),t1l,elconloc(*),plconloc(802),xk,fk,
     &  xiso(200),yiso(200),dd0,plicon(0:2*npmat_,ntmat_,*),
     &  a11,a12,a22,b1(3,10),b2(3,10),dal(3,3,10),qxxy(3),fnl(3),
     &  qxyy(3),dxi(3,10),det(3,10),determinant,c11,c12,c22,
     &  qxyx(3),qyxy(3),springarea(2),dd,dist,t(3),tu(3,3,10),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  um,eps,pi,dftdt(3,3),tp(3),te(3),ftrial(3),clear,
     &  dftrial,dfnl,dfshear,dg,dte,alnew(3),dfn(3,10),reltime,
     &  overlap,pres,dpresdoverlap,overclosure,orab(7,*),
     &  xn1(3),xn2(3),a(3,3)
!
      iflag=4
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
      if(lakonl(7:7).eq.'A') then
         if(lakonl(4:6).eq.'RNG') then
!
!           SPRINGA-element
!
            dd0=dsqrt((xl(1,2)-xl(1,1))**2
     &           +(xl(2,2)-xl(2,1))**2
     &           +(xl(3,2)-xl(3,1))**2)
            dd=dsqrt((pl(1,2)-pl(1,1))**2
     &           +(pl(2,2)-pl(2,1))**2
     &           +(pl(3,2)-pl(3,1))**2)
            if(dd.le.0.d0) then
               write(*,*) '*ERROR in springstiff_n2f: spring connects'
               write(*,*) '       two nodes at the same location:'
               write(*,*) '       spring length is zero'
               call exit(201)
            endif
            do i=1,3
               xn(i)=(pl(i,2)-pl(i,1))/dd
            enddo
            val=dd-dd0
!     
!     interpolating the material data
!     
            call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &           elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!     
!     calculating the spring force and the spring constant
!     
            if(kode.eq.2)then
               xk=elconloc(1)
               fk=xk*val
            else
               niso=int(plconloc(801))
               do i=1,niso
                  xiso(i)=plconloc(2*i-1)
                  yiso(i)=plconloc(2*i)
               enddo
               call ident(xiso,val,niso,id)
               if(id.eq.0) then
                  xk=0.d0
                  fk=yiso(1)
               elseif(id.eq.niso) then
                  xk=0.d0
                  fk=yiso(niso)
               else
                  xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
                  fk=yiso(id)+xk*(val-xiso(id))
               endif
            endif
!     
            c1=fk/dd
            c2=xk-c1
            do i=1,3
               do j=1,3
                  s(i,j)=c2*xn(i)*xn(j)
               enddo
               s(i,i)=s(i,i)+c1
            enddo
         else
!
!           GAP-element
!           interpolating the material data
!     
            call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &           elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!
            dd=elconloc(1)
            xn(1)=elconloc(2)
            xn(2)=elconloc(3)
            xn(3)=elconloc(4)
            xk=elconloc(5)
            pi=4.d0*datan(1.d0)
            eps=pi*elconloc(6)/xk
            overclosure=-dd-xn(1)*(voldl(1,2)-voldl(1,1))
     &                     -xn(2)*(voldl(2,2)-voldl(2,1))
     &                     -xn(3)*(voldl(3,2)-voldl(3,1))
            fk=-xk*overclosure*(0.5d0+datan(overclosure/eps)/pi)
            c2=xk*((0.5d0+datan(overclosure/eps)/pi)
     &        +overclosure/(pi*eps*(1.d0+(overclosure/eps)**2)))
            do i=1,3
               do j=1,3
                  s(i,j)=c2*xn(i)*xn(j)
               enddo
            enddo
         endif
!
         do i=1,3
            do j=1,3
               s(i+3,j)=-s(i,j)
               s(i,j+3)=-s(i,j)
               s(i+3,j+3)=s(i,j)
            enddo
         enddo
         return
      elseif(lakonl(7:7).eq.'1') then
!
!        spring1-element
!        determine the direction of action of the spring
!
         idof=nint(elcon(3,1,imat))
!
         if(norien.gt.0) then
            iorien=max(0,ielorien(1,nelem))
         else
            iorien=0
         endif
!
         if(iorien.eq.0) then
            do i=1,3
               xn(i)=0.d0
            enddo
            xn(idof)=1.d0
         else
            call transformatrix(orab(1,iorien),xl(1,1),a)
            do i=1,3
               xn(i)=a(i,idof)
            enddo
         endif
!     
!     interpolating the material data
!     
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &           elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!     
!     calculating the spring constant
!     
         if(kode.eq.2)then
            xk=elconloc(1)
         else
            niso=int(plconloc(801))
            do i=1,niso
               xiso(i)=plconloc(2*i-1)
               yiso(i)=plconloc(2*i)
            enddo
            call ident(xiso,val,niso,id)
            if(id.eq.0) then
               xk=0.d0
            elseif(id.eq.niso) then
               xk=0.d0
            else
               xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            endif
         endif
!
!        assembling the stiffness matrix
!
         do i=1,3
            do j=1,3
               s(i,j)=xk*xn(i)*xn(j)
            enddo
         enddo
         return
      elseif(lakonl(7:7).eq.'2') then
!
!        spring2-element
!        determine the direction of action of the spring
!
         idof1=nint(elcon(3,1,imat))
         idof2=nint(elcon(4,1,imat))
!
         if(norien.gt.0) then
            iorien=max(0,ielorien(1,nelem))
         else
            iorien=0
         endif
!
         if(iorien.eq.0) then
            do i=1,3
               xn1(i)=0.d0
               xn2(i)=0.d0
            enddo
            xn1(idof1)=1.d0
            xn2(idof2)=1.d0
         else
            call transformatrix(orab(1,iorien),xl(1,1),a)
            do i=1,3
               xn1(i)=a(i,idof1)
            enddo
            call transformatrix(orab(1,iorien),xl(1,2),a)
            do i=1,3
               xn2(i)=a(i,idof2)
            enddo
         endif
!     
!     interpolating the material data
!     
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &           elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
!     
!        calculating the spring constant
!     
         if(kode.eq.2)then
            xk=elconloc(1)
         else
            niso=int(plconloc(801))
            do i=1,niso
               xiso(i)=plconloc(2*i-1)
               yiso(i)=plconloc(2*i)
            enddo
            call ident(xiso,val,niso,id)
            if(id.eq.0) then
               xk=0.d0
            elseif(id.eq.niso) then
               xk=0.d0
            else
               xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            endif
         endif
!
!        assembling the stiffness matrix
!
         do i=1,3
            do j=1,3
               s(i,j)=xk*xn1(i)*xn1(j)
               s(i,j+3)=-xk*xn1(i)*xn2(j)
               s(i+3,j)=-xk*xn2(i)*xn1(j)
               s(i+3,j+3)=xk*xn2(i)*xn2(j)
            enddo
         enddo
         return
!
      endif
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
!     alpha and beta, taking the representative area into account
!     (conversion of pressure into force)
!
      if(int(elcon(3,1,imat)).eq.1) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            elas(1)=0.d0
            elas(2)=0.d0
         else
            alpha=elcon(2,1,imat)*springarea(1)
            beta=elcon(1,1,imat)
            if(-beta*clear.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/clear
            endif
            elas(1)=dexp(-beta*clear+dlog(alpha))
            elas(2)=-beta*elas(1)
         endif
      elseif(int(elcon(3,1,imat)).eq.2) then
!     
!        linear overclosure
!
         pi=4.d0*datan(1.d0)
         eps=elcon(1,1,imat)*pi/elcon(2,1,imat)
         elas(1)=(-springarea(1)*elcon(2,1,imat)*clear*
     &            (0.5d0+datan(-clear/eps)/pi))
         elas(2)=-springarea(1)*elcon(2,1,imat)*
     &            ((0.5d0+datan(-clear/eps)/pi)-
     &             clear/(pi*eps*(1.d0+(clear/eps)**2)))
      elseif(int(elcon(3,1,imat)).eq.3) then
!     
!        tabular overclosure
!
!        interpolating the material data
!
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
         overlap=-clear
         niso=int(plconloc(801))
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,overlap,niso,id)
         if(id.eq.0) then
            dpresdoverlap=0.d0
            pres=yiso(1)
         elseif(id.eq.niso) then
            dpresdoverlap=0.d0
            pres=yiso(niso)
         else
            dpresdoverlap=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            pres=yiso(id)+dpresdoverlap*(overlap-xiso(id))
         endif
         elas(1)=springarea(1)*pres
         elas(2)=-springarea(1)*dpresdoverlap
      endif
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
!     Coulomb friction for static calculations
!    
      if(ncmat_.ge.7) then
         um=elcon(6,1,imat)
         if(um.gt.0.d0) then
!     
!     stiffness of shear stress versus slip curve
!     
            xk=elcon(7,1,imat)*springarea(1)
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
               al(i)=alnew(i)-xstateini(3+i,1,ne0+konl(nope+1))
            enddo
!     
!     s=||vec_s||
!     
            val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!     
!     update the relative tangential displacement vec_t
!     
            do i=1,3
               t(i)=xstateini(6+i,1,ne0+konl(nope+1))+al(i)-val*xn(i)
            enddo
!     
!     store the actual relative displacement and
!     the actual relative tangential displacement
!     
            do i=1,3
               xstate(3+i,1,ne0+konl(nope+1))=alnew(i)
               xstate(6+i,1,ne0+konl(nope+1))=t(i)
            enddo
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
     &                 +al(3)*xmu(3,i,k)-val*dxmu(i,k)
     &                 +xm(1)*dal(1,i,k)+xm(2)*dal(2,i,k)
     &                 +xm(3)*dal(3,i,k)
               enddo
            enddo
!     
!     d vec_t/d vec_u_k
!     
            do k=1,nope
               do j=1,3
                  do i=1,3
                     tu(i,j,k)=dal(i,j,k)
     &                    -c1*(xn(i)*(dval(j,k)-val*dxmu(j,k))
     &                    +val*xmu(i,j,k))
                  enddo
               enddo
            enddo
!     
!     size of normal force
!     
            dfnl=dsqrt(fnl(1)**2+fnl(2)**2+fnl(3)**2)
!     
!     maximum size of shear force
!     
            dfshear=um*dfnl       
!     
!     plastic and elastic slip
!     
            do i=1,3
               tp(i)=xstateini(i,1,ne0+konl(nope+1))
               te(i)=t(i)-tp(i)
            enddo
!     
!     the force due to normal contact is in -xn
!     direction (internal force) -> minus signs
!     
            do k=1,nope
               do i=1,3
                  dfn(i,k)=-xn(1)*fpu(1,i,k)-xn(2)*fpu(2,i,k)-
     &                 xn(3)*fpu(3,i,k)  
               enddo
            enddo
!     
            dte=dsqrt(te(1)*te(1)+te(2)*te(2)+te(3)*te(3))
!     
!     trial force
!     
            do i=1,3
               ftrial(i)=xk*te(i)
            enddo
            dftrial=dsqrt(ftrial(1)**2+ftrial(2)**2+ftrial(3)**2)
!     
!     check whether stick or slip
!     
            if((dftrial.lt.dfshear) .or. (dftrial.le.0.d0)) then
!     
!     stick force
!     
               do i=1,3
                  fnl(i)=fnl(i)+ftrial(i)
                  xstate(i,1,ne0+konl(nope+1))=tp(i)
               enddo
!     
!     stick stiffness
!     
               do k=1,nope
                  do j=1,3
                     do i=1,3
                        fpu(i,j,k)=fpu(i,j,k)+xk*tu(i,j,k)
                     enddo
                  enddo
               enddo  
            else
!     
!     slip force
!     
               dg=(dftrial-dfshear)/xk
               do i=1,3
                  ftrial(i)=te(i)/dte
                  fnl(i)=fnl(i)+dfshear*ftrial(i)
                  xstate(i,1,ne0+konl(nope+1))=tp(i)+dg*ftrial(i)
               enddo
!     
!     slip stiffness
!     
               c1=xk*dfshear/dftrial
               do i=1,3
                  do j=1,3
                     dftdt(i,j)=-c1*ftrial(i)*ftrial(j)
                  enddo
                  dftdt(i,i)=dftdt(i,i)+c1
               enddo
!     
               do k=1,nope
                  do j=1,3
                     do i=1,3
                        do l=1,3
                           fpu(i,j,k)=fpu(i,j,k)+dftdt(i,l)*tu(l,j,k)
                        enddo
                        if((nmethod.ne.4).or.(iperturb(1).gt.1)) then
                           fpu(i,j,k)=fpu(i,j,k)+um*ftrial(i)*dfn(j,k)
                        endif
                     enddo
                  enddo
               enddo
            endif
         endif
      endif
!     
!     determining the stiffness matrix contributions
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

