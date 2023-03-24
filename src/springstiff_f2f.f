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
      subroutine springstiff_f2f(xl,elas,voldl,s,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,plicon,
     &  nplicon,npmat_,iperturb,springarea,nmethod,mi,ne0,
     &  nstate_,xstateini,xstate,reltime,nasym,
     &  jfaces,igauss,pslavsurf,pmastsurf,clearini,kscale)
!
!     calculates the stiffness of a spring (face-to-face penalty)
!
      implicit none
!
      character*8 lakonl
!
      integer i,j,imat,ncmat_,ntmat_,k,l,nope,iflag,
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),
     &  iperturb(*),nmethod,mi(*),ne0,nstate_,nasym,
     &  jfaces,igauss,nopem,nopes,nopep,kscale
!
      real*8 xl(3,19),elas(21),pproj(3),val,shp2m(7,9),
     &  al(3),s(60,60),voldl(0:mi(2),19),pl(3,19),xn(3),
     &  c1,c3,alpha,beta,elcon(0:ncmat_,ntmat_,*),xm(3),
     &  fpu(3,3),xi,et,fnl(3),
     &  xs2(3,7),t1l,elconloc(*),plconloc(82),xk,stickslope,
     &  xiso(20),yiso(20),plicon(0:2*npmat_,ntmat_,*),
     &  springarea(2),t(3),tu(3,3),overlap,pres,dpresdoverlap,
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  um,dftdt(3,3),tp(3),te(3),ftrial(3),clear,
     &  dftrial,dfnl,dfshear,dg,dte,alnew(3),dfn(3),reltime,
     &  xsj2s(3),xs2s(3,7),shp2s(7,9),weight,pslavsurf(3,*),
     &  pmastsurf(6,*),clearini(3,9,*)
!
      iflag=1
!     
!     # of master nodes
!
      nopem=ichar(lakonl(8:8))-48
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
      elseif((int(elcon(3,1,imat)).eq.2).or.
     &       (int(elcon(3,1,imat)).eq.4)) then
!     
!        linear overclosure/tied overclosure
!
c        write(*,*) 'springstiff_f2f',springarea(1),elcon(2,1,imat),
c     &       kscale
         elas(2)=-springarea(1)*elcon(2,1,imat)/kscale
         elas(1)=elas(2)*clear
      elseif(int(elcon(3,1,imat)).eq.3) then
!     
!        tabular overclosure
!
!        interpolating the material data
!
         call materialdata_sp(elcon,nelcon,imat,ntmat_,i,t1l,
     &     elconloc,kode,plicon,nplicon,npmat_,plconloc,ncmat_)
         overlap=-clear
         niso=int(plconloc(81))
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
      c3=elas(2)
!
!     derivatives of the forces w.r.t. the displacement vectors
!
      do j=1,3
         do i=1,3
            fpu(i,j)=-c3*xn(i)*xn(j)
         enddo
      enddo
!
!     Coulomb friction for static calculations
!    
      if((ncmat_.ge.7).or.(int(elcon(3,1,imat)).eq.4)) then
!
!        tied contact
!
         if(int(elcon(3,1,imat)).eq.4) then
            um=1.d30
         else
            um=elcon(6,1,imat)
         endif
         stickslope=elcon(7,1,imat)/kscale
!
         if(um.gt.0.d0) then
!     
!     stiffness of shear stress versus slip curve
!     
            xk=stickslope*springarea(1)
!     
!     calculating the relative displacement between the slave node
!     and its projection on the master surface
!     
            do i=1,3
               alnew(i)=voldl(i,nopep)
               do j=1,nopem
                  alnew(i)=alnew(i)-shp2m(4,j)*voldl(i,j)
               enddo
            enddo
!     
!     calculating the difference in relative displacement since
!     the start of the increment = lamda^*
!           
            do i=1,3
               al(i)=alnew(i)-xstateini(3+i,1,ne0+igauss)
            enddo
!     
!     ||lambda^*||
!     
            val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!     
!     update the relative tangential displacement
!     
            do i=1,3
               t(i)=xstateini(6+i,1,ne0+igauss)+al(i)-val*xn(i)
            enddo
!     
!     store the actual relative displacement and
!     the actual relative tangential displacement
!     
            do i=1,3
               xstate(3+i,1,ne0+igauss)=alnew(i)
               xstate(6+i,1,ne0+igauss)=t(i)
            enddo
!     
!     d t/d u_k
!     
            do j=1,3
               do i=1,3
                  tu(i,j)=-xn(i)*xn(j)
               enddo
               tu(j,j)=tu(j,j)+1.d0
            enddo
!     
!     size of normal force
!     
            dfnl=dsqrt(fnl(1)**2+fnl(2)**2+fnl(3)**2)
!     
!     maximum size of shear force
!     
            if(int(elcon(3,1,imat)).eq.4) then
               dfshear=1.d30
            else
               dfshear=um*dfnl 
            endif
c            dfshear=um*dfnl       
!     
!     plastic and elastic slip
!     
            do i=1,3
               tp(i)=xstateini(i,1,ne0+igauss)
               te(i)=t(i)-tp(i)
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
                  xstate(i,1,ne0+igauss)=tp(i)
               enddo
!     
!     stick stiffness
!     
               do j=1,3
                  do i=1,3
                     fpu(i,j)=fpu(i,j)+xk*tu(i,j)
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
                  xstate(i,1,ne0+igauss)=tp(i)+dg*ftrial(i)
               enddo
!     
!     slip stiffness
!     
               do i=1,3
                  dfn(i)=-xn(1)*fpu(1,i)-xn(2)*fpu(2,i)-
     &                 xn(3)*fpu(3,i)  
               enddo
!
               c1=xk*dfshear/dftrial
               do i=1,3
                  do j=1,3
                     dftdt(i,j)=-c1*ftrial(i)*ftrial(j)
                  enddo
                  dftdt(i,i)=dftdt(i,i)+c1
               enddo
!     
               do j=1,3
                  do i=1,3
                     do l=1,3
                        fpu(i,j)=fpu(i,j)+
     &                                 dftdt(i,l)*tu(l,j)
                     enddo
                     if((nmethod.ne.4).or.(iperturb(1).gt.1)) then
                        fpu(i,j)=fpu(i,j)+
     &                                 um*ftrial(i)*dfn(j)
                     endif
                  enddo
               enddo
            endif
         endif
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
                  s(i+(k-1)*3,j+(l-1)*3)=
     &                 shp2m(4,k)*shp2m(4,l)*fpu(i,j)
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
     &                                shp2s(4,l-nopem)*fpu(i,j)
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
     &                shp2m(4,k)*fpu(i,j)
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
                  s(i+(k-1)*3,j+(l-1)*3)=
     &                -shp2s(4,k-nopem)*shp2m(4,l)*fpu(i,j)
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

