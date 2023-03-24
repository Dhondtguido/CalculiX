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
      subroutine springforc_f2f(xl,vl,imat,elcon,nelcon,
     &  elas,fnl,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &  plicon,nplicon,npmat_,senergy,nener,cstr,mi,
     &  springarea,nmethod,ne0,nstate_,xstateini,
     &  xstate,reltime,ielas,jfaces,igauss,
     &  pslavsurf,pmastsurf,clearini,venergy,kscale,
     &  konl,iout,nelem,smscale,mscalmethod)
!
!     calculates the force of the spring (face-to-face penalty)
!
      implicit none
!
      character*8 lakonl
!
      integer i,j,k,imat,ncmat_,ntmat_,nope,iflag,mi(*),
     &  kode,niso,id,nplicon(0:ntmat_,*),npmat_,nelcon(2,*),nener,
     &  nmethod,ne0,nstate_,ielas,jfaces,kscale,konl(26),
     &  igauss,nopes,nopem,nopep,iout,nelem,mscalmethod
!
      real*8 xl(3,10),elas(21),t1l,al(3),vl(0:mi(2),19),stickslope,
     &  pl(3,19),xn(3),alpha,beta,fnl(3,19),tp(3),te(3),ftrial(3),
     &  t(3),dftrial,elcon(0:ncmat_,ntmat_,*),pproj(3),clear,
     &  xi,et,elconloc(*),plconloc(82),xk,val,xiso(20),yiso(20),
     &  plicon(0:2*npmat_,ntmat_,*),um,senergy,cstr(6),dg,venergy,
     &  dfshear,dfnl,springarea(2),overlap,pres,clearini(3,9,*),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),t1(3),t2(3),
     &  dt1,dte,alnew(3),reltime,weight,xsj2m(3),xs2m(3,7),shp2m(7,9),
     &  xsj2s(3),xs2s(3,7),shp2s(7,9),pslavsurf(3,*),pmastsurf(6,*),
     &  smscale(*)
!
      include "gauss.f"
!
      iflag=2
!     
!     # of master nodes
!
      nopem=ichar(lakonl(8:8))-48
!
!     # of slave nodes
!
      nopes=nope-nopem
!
!     actual positions of the master nodes belonging to the contact spring
!     
      do i=1,nopem
         do j=1,3
            pl(j,i)=xl(j,i)+vl(j,i)
         enddo
      enddo
!
!     actual positions of the slave nodes belonging to the contact spring
!     
      if(nmethod.ne.2) then
         do i=nopem+1,nope
            do j=1,3
               pl(j,i)=xl(j,i)+clearini(j,i-nopem,jfaces)*reltime
     &                        +vl(j,i)
            enddo
         enddo
      else
!
!        for frequency calculations the eigenmodes are freely
!        scalable, leading to problems with contact finding 
!
         do i=nopem+1,nope
            do j=1,3
               pl(j,i)=xl(j,i)+clearini(j,i-nopem,jfaces)*reltime
            enddo
         enddo
      endif
!
!     location of integration point in slave face
!
      xi=pslavsurf(1,igauss)
      et=pslavsurf(2,igauss)
      weight=pslavsurf(3,igauss)
!
      iflag=1
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
      do k=1,3
          pl(k,nopep)=0.d0
          vl(k,nopep)=0.d0
          do j=1,nopes
              pl(k,nopep)=pl(k,nopep)+shp2s(4,j)*pl(k,nopem+j)
              vl(k,nopep)=vl(k,nopep)+shp2s(4,j)*vl(k,nopem+j)
          enddo
      enddo    
!
!     corresponding location in the master face
!
      xi=pmastsurf(1,igauss)
      et=pmastsurf(2,igauss)
!
!     determining the jacobian vector on the surface 
!
      iflag=2
      if(nopem.eq.8) then
         call shape8q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      elseif(nopem.eq.4) then
         call shape4q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      elseif(nopem.eq.6) then
         call shape6tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      else
         call shape3tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
      endif
!
      do i=1,3
         pproj(i)=0.d0
         do j=1,nopem
            pproj(i)=pproj(i)+shp2m(4,j)*pl(i,j)
         enddo
      enddo
!
!     distance vector between both
!
      do i=1,3
         al(i)=pl(i,nopep)-pproj(i)
      enddo
!
!     normal on the master face
!
      xn(1)=pmastsurf(4,igauss)
      xn(2)=pmastsurf(5,igauss)
      xn(3)=pmastsurf(6,igauss)
!
!     distance from surface along normal (= clearance)
!
      clear=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!
!     check for a reduction of the initial penetration, if any
!
      if(nmethod.eq.1) then
         clear=clear-springarea(2)*(1.d0-reltime)
      endif
      cstr(1)=clear
!
!
      if(int(elcon(3,1,imat)).eq.1) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            elas(1)=0.d0
            beta=1.d0
         else
!     
            alpha=elcon(2,1,imat)*springarea(1)
            beta=elcon(1,1,imat)
            if(-beta*clear.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/clear
            endif
            elas(1)=dexp(-beta*clear+dlog(alpha))
         endif
         if(nener.eq.1) then
            senergy=elas(1)/beta
         endif
      elseif((int(elcon(3,1,imat)).eq.2).or.
     &       (int(elcon(3,1,imat)).eq.4)) then
!     
!        linear overclosure/tied overclosure
!     
         elas(1)=-springarea(1)*elcon(2,1,imat)*clear/kscale
!
!     spring scaling for explicit dynamics
!     
         if((mscalmethod.eq.2).or.(mscalmethod.eq.3)) then
            elas(1)=elas(1)*smscale(nelem)
         endif
!         
         if(nener.eq.1) then
            senergy=-elas(1)*clear/2.d0
         endif
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
            pres=yiso(1)
            if(nener.eq.1) then
               senergy=yiso(1)*overlap
            endif
         elseif(id.eq.niso) then
            pres=yiso(niso)
            if(nener.eq.1) then
               senergy=yiso(1)*xiso(1)
               do i=2,niso
                  senergy=senergy+(xiso(i)-xiso(i-1))*
     &                            (yiso(i)+yiso(i-1))/2.d0
               enddo
               senergy=senergy+(pres-xiso(niso))*yiso(niso)
            endif
         else
            xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            pres=yiso(id)+xk*(overlap-xiso(id))
            if(nener.eq.1) then
               senergy=yiso(1)*xiso(1)
               do i=2, id
                  senergy=senergy+(xiso(i)-xiso(i-1))*
     &                 (yiso(i)+yiso(i-1))/2.d0
               enddo
               senergy=senergy+(overlap-xiso(id))*(pres+yiso(id))/2.d0
            endif
         endif
         elas(1)=springarea(1)*pres
         if(nener.eq.1) senergy=springarea(1)*senergy
      endif
!
!     forces in the nodes of the contact element
!
      do i=1,3
         fnl(i,nopep)=-elas(1)*xn(i)
      enddo
      cstr(4)=elas(1)/springarea(1)
!
!     Coulomb friction for static calculations
!
      if((int(elcon(3,1,imat)).eq.4).and.(ncmat_.lt.7)) then
         write(*,*) '*ERROR in springforc_f2f: for contact'
         write(*,*) '       with pressure-overclosure=tied'
         write(*,*) '       a stick slope using the *FRICTION'
         write(*,*) '       card is mandatory'
         call exit(201)
      endif
!
      if(ncmat_.ge.7) then
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
!     spring scaling for explicit dynamics
!     
         if((mscalmethod.eq.2).or.(mscalmethod.eq.3)) then
            stickslope=stickslope*smscale(nelem)
         endif
!
         if(um.gt.0.d0) then
            if(1.d0 - dabs(xn(1)).lt.1.5231d-6) then       
!     
!     calculating the local directions on master surface
!     
               t1(1)=-xn(3)*xn(1)
               t1(2)=-xn(3)*xn(2)
               t1(3)=1.d0-xn(3)*xn(3)
            else
               t1(1)=1.d0-xn(1)*xn(1)
               t1(2)=-xn(1)*xn(2)
               t1(3)=-xn(1)*xn(3)
            endif
            dt1=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
            do i=1,3
               t1(i)=t1(i)/dt1
            enddo
            t2(1)=xn(2)*t1(3)-xn(3)*t1(2)
            t2(2)=xn(3)*t1(1)-xn(1)*t1(3)
            t2(3)=xn(1)*t1(2)-xn(2)*t1(1)           
!     
!     linear stiffness of the shear stress versus
!     slip curve
!     
            xk=stickslope*springarea(1)
!     
!     calculating the relative displacement between the slave node
!     and its projection on the master surface
!     
            do i=1,3
               alnew(i)=vl(i,nopep)
               do j=1,nopem
                  alnew(i)=alnew(i)-shp2m(4,j)*vl(i,j)
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
!     size of normal force
!     
            dfnl=dsqrt(fnl(1,nopep)**2+fnl(2,nopep)**2+
     &           fnl(3,nopep)**2)
!     
!     maximum size of shear force
! 
            if(int(elcon(3,1,imat)).eq.4) then
               dfshear=1.d30
            else
               dfshear=um*dfnl 
            endif
!     
!     plastic and elastic slip
!     
            do i=1,3
               tp(i)=xstateini(i,1,ne0+igauss)
               te(i)=t(i)-tp(i)
            enddo
            dte=dsqrt(te(1)*te(1)+te(2)*te(2)+te(3)*te(3))
!     
!     trial force
!     
            do i=1,3
               ftrial(i)=xk*te(i)
            enddo
            dftrial=xk*dte
c            dftrial=dsqrt(ftrial(1)**2+ftrial(2)**2+ftrial(3)**2)
!     
!     check whether stick or slip
!     
            if((dftrial.lt.dfshear).or.(dftrial.le.0.d0).or.
     &           (ielas.eq.1)) then
!     
!     stick
!     
               do i=1,3
                  fnl(i,nopep)=fnl(i,nopep)+ftrial(i)
                  xstate(i,1,ne0+igauss)=tp(i)
               enddo
               cstr(5)=(ftrial(1)*t1(1)+ftrial(2)*t1(2)+
     &              ftrial(3)*t1(3))/springarea(1)
               cstr(6)=(ftrial(1)*t2(1)+ftrial(2)*t2(2)+
     &              ftrial(3)*t2(3))/springarea(1)
!
!              shear elastic energy
!
               if(nener.eq.1) senergy=senergy+dftrial*dte/2.d0
            else
!     
!     slip
!     
               dg=(dftrial-dfshear)/xk
               do i=1,3
                  ftrial(i)=te(i)/dte
                  fnl(i,nopep)=fnl(i,nopep)+dfshear*ftrial(i)
                  xstate(i,1,ne0+igauss)=tp(i)+dg*ftrial(i)
               enddo
               cstr(5)=(dfshear*ftrial(1)*t1(1)+
     &              dfshear*ftrial(2)*t1(2)+
     &              dfshear*ftrial(3)*t1(3))/springarea(1)
               cstr(6)=(dfshear*ftrial(1)*t2(1)+
     &              dfshear*ftrial(2)*t2(2)+
     &              dfshear*ftrial(3)*t2(3))/springarea(1)
!
!              shear elastic and viscous energy
!
               if(nener.eq.1) then
                  senergy=senergy+dfshear*dfshear/(xk*2.d0)
                  venergy=venergy+dg*dfshear
               endif
!
            endif
         endif
!     
!     storing the tangential displacements
!     
         cstr(2)=t(1)*t1(1)+t(2)*t1(2)+t(3)*t1(3)
         cstr(3)=t(1)*t2(1)+t(2)*t2(2)+t(3)*t2(3)
      endif
!
!     force in the master nodes
!
      do i=1,3
         do j=1,nopem
            fnl(i,j)=-shp2m(4,j)*fnl(i,nopep)
         enddo
      enddo
!
!     forces in the nodes of the slave face
!
      do i=1,3
          do j=1,nopes
              fnl(i,nopem+j)=shp2s(4,j)*fnl(i,nopep)
          enddo
      enddo
!
!     write statements for Malte Krack
!
c      if(iout.gt.0) then
c         write(*,*) 'contact element: ',nelem
c         write(*,*) 'undeformed location of the integration point'
c         write(*,*) ((pl(j,nopep)-vl(j,nopep)),j=1,3)
c         write(*,*) 'deformed location of the integration point'
c         write(*,*) (pl(j,nopep),j=1,3)
c         write(*,*) 'nodes and shape values'
c         do j=1,nopem
c            write(*,*) konl(j),-shp2m(4,j)
c         enddo
c         do j=1,nopes
c            write(*,*) konl(nopem+j),shp2s(4,j)
c         enddo
c         write(*,*) 'contact force'
c         write(*,*) fnl(1,nopep),fnl(2,nopep),fnl(3,nopep)
c         write(*,*) 'slave area'
c         write(*,*) springarea(1)
c      endif
!
      return
      end

