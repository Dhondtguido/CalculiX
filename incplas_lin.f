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
      subroutine incplas_lin(elconloc,plconloc,xstate,xstateini,
     &  elas,emec,ithermal,icmd,beta,stre,vj,kode,
     &  ielas,amat,t1l,dtime,time,ttime,iel,iint,nstate_,mi,
     &  eloc,pgauss,nmethod,pnewdt,depvisc)
!
!     calculates stiffness and stresses for the incremental plasticity
!     material law
!
!     icmd=3: calculates stress at mechanical strain
!     else: calculates stress at mechanical strain and the stiffness
!           matrix
!
!     This routine is meant for small strains. Reference: 
!     Dhondt, Guido; The Finite Element Method for Three-dimensional
!     Thermomechanical Applications (Wiley, 2004), Section 5.3
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,i,k,l,m,n,kode,ivisco,ielastic,kel(4,21),
     &  niso,nkin,ielas,iel,iint,nstate_,mi(*),id,leximp,lend,layer,
     &  kspt,kstep,kinc,iloop,nmethod,user_hardening,user_creep
!
      real*8 elconloc(*),elas(21),emec(6),beta(6),stre(6),
     &  vj,plconloc(802),stbl(6),stril(6),xitril(6),
     &  ee,un,um,al,cop,dxitril,xn(3,3),epl(6),c1,c2,c3,c4,c7,
     &  c8,ftrial,xiso(200),yiso(200),xkin(200),ykin(200),
     &  fiso,dfiso,fkin,dfkin,fiso0,fkin0,ep,t1l,dtime,
     &  epini,a1,dsvm,xxa,xxn,dkl(3,3),el(6),tracee,traces,
     &  dcop,time,ttime,eloc(6),xstate(nstate_,mi(1),*),
     &  xstateini(nstate_,mi(1),*),decra(5),deswa(5),serd,
     &  esw(2),ec(2),p,qtild,predef(1),dpred(1),timeabq(2),pgauss(3),
     &  dtemp,pnewdt,um2,depvisc
!
!
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &          1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &          3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &          1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
      dkl=reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),
     &      (/3,3/))
!
      leximp=1
      lend=2
      user_creep=0
!
!     localizing the plastic fields
!
      do i=1,6
         epl(i)=xstateini(1+i,iint,iel)
         stbl(i)=xstateini(7+i,iint,iel)
      enddo
      epini=xstateini(1,iint,iel)
!
      ee=elconloc(1)
      un=elconloc(2)
      um2=ee/(1.d0+un)
      al=um2*un/(1.d0-2.d0*un)
      um=um2/2.d0
!
      ep=epini
!
!     check whether the user activated a viscous calculation
!     (*VISCO instead of *STATIC)
!
      if((nmethod.ne.1).or.(ithermal(1).eq.3)) then
         ivisco=1
      else
         ivisco=0
      endif
!
!     check for user subroutines
!
      if((plconloc(801).lt.0.8d0).and.(plconloc(802).lt.0.8d0)) then
         user_hardening=1
      else
         user_hardening=0
      endif
      if((kode.eq.-52).and.(ivisco.eq.1)) then
         if(elconloc(3).lt.0.d0) then
            user_creep=1
         else
            user_creep=0
            xxa=elconloc(3)*(ttime+time)**elconloc(5)
            if(xxa.lt.1.d-20) xxa=1.d-20
            xxn=elconloc(4)
            a1=xxa*dtime
         endif
      endif
!
!     trial elastic strain
!     
      do i=1,6
         el(i)=emec(i)-epl(i)
      enddo
!
!     trial stress
!
      tracee=el(1)+el(2)+el(3)
      do i=1,6
         stre(i)=um2*el(i)-beta(i)
      enddo
      do i=1,3
         stre(i)=stre(i)+al*tracee
      enddo
!
!     trial deviatoric stress
!
      traces=(stre(1)+stre(2)+stre(3))/3.d0
      do i=1,3
         stril(i)=stre(i)-traces
      enddo
      do i=4,6
         stril(i)=stre(i)
      enddo
!
!     subtracting the back stress -> trial radius vector
!
      do i=1,6
         xitril(i)=stril(i)-stbl(i)
      enddo
!
!     size of trial radius vector
!
      dxitril=0.d0
      do i=1,3
         dxitril=dxitril+xitril(i)*xitril(i)
      enddo
      do i=4,6
         dxitril=dxitril+2.d0*xitril(i)*xitril(i)
      enddo
      dxitril=dsqrt(dxitril)
!
!        restoring the hardening curves for the actual temperature
!        plconloc contains the true stresses. By multiplying by
!        the Jacobian, yiso and ykin are Kirchhoff stresses, as
!        required by the hyperelastic theory (cf. Simo, 1988).
!
      niso=int(plconloc(801))
      nkin=int(plconloc(802))
      if(niso.ne.0) then
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
      endif
      if(nkin.ne.0) then
         do i=1,nkin
            xkin(i)=plconloc(399+2*i)
            ykin(i)=plconloc(400+2*i)
         enddo
      endif
!
!     if no viscous calculation is performed a pure creep calculatino
!     (without plasticity) is reduced to an elastic calculation
!
      ielastic=0
      if(ivisco.eq.0) then
         if(niso.eq.2) then
            if((dabs(yiso(1)).lt.1.d-10).and.
     &         (dabs(yiso(2)).lt.1.d-10)) then
               ielastic=1
            endif
         endif
      endif
!
!     check for yielding
!
      if(user_hardening.eq.1) then
         call uhardening(amat,iel,iint,t1l,epini,ep,dtime,
     &        fiso,dfiso,fkin,dfkin)
      else
         if(niso.ne.0) then
            call ident(xiso,ep,niso,id)
            if(id.eq.0) then
               fiso=yiso(1)
            elseif(id.eq.niso) then
               fiso=yiso(niso)
            else
               dfiso=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
               fiso=yiso(id)+dfiso*(ep-xiso(id))
            endif
         elseif(nkin.ne.0) then
            fiso=ykin(1)
         else
            fiso=0.d0
         endif
      endif
!
      c1=2.d0/3.d0
      c2=dsqrt(c1)
!
      ftrial=dxitril-c2*fiso
      if((ftrial.le.1.d-10).or.(ielas.eq.1).or.(ielastic.eq.1)
     &    .or.(dtime.lt.1.d-30)) then
!
!        updating the plastic fields
!
         do i=1,6
            xstate(1+i,iint,iel)=epl(i)
            xstate(7+i,iint,iel)=stbl(i)
         enddo
         xstate(1,iint,iel)=ep
!
         if(icmd.ne.3) then
            elas(1)=al+um2
            elas(2)=al
            elas(3)=al+um2
            elas(4)=al
            elas(5)=al
            elas(6)=al+um2
            elas(7)=0.d0
            elas(8)=0.d0
            elas(9)=0.d0
            elas(10)=um
            elas(11)=0.d0
            elas(12)=0.d0
            elas(13)=0.d0
            elas(14)=0.d0
            elas(15)=um
            elas(16)=0.d0
            elas(17)=0.d0
            elas(18)=0.d0
            elas(19)=0.d0
            elas(20)=0.d0
            elas(21)=um
         endif
!
         return
      endif
!
!        plastic deformation
!
!        calculating the consistency parameter
!
      iloop=0
      cop=0.d0
!
      loop: do
         iloop=iloop+1
         ep=epini+c2*cop
!
         if(user_hardening.eq.1) then
            call uhardening(amat,iel,iint,t1l,epini,ep,dtime,
     &           fiso,dfiso,fkin,dfkin)
         else
            if(niso.ne.0) then
               call ident(xiso,ep,niso,id)
               if(id.eq.0) then
                  fiso=yiso(1)
                  dfiso=0.d0
               elseif(id.eq.niso) then
                  fiso=yiso(niso)
                  dfiso=0.d0
               else
                  dfiso=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
                  fiso=yiso(id)+dfiso*(ep-xiso(id))
               endif
            elseif(nkin.ne.0) then
               fiso=ykin(1)
               dfiso=0.d0
            else
               fiso=0.d0
               dfiso=0.d0
            endif
!
            if(nkin.ne.0) then
               call ident(xkin,ep,nkin,id)
               if(id.eq.0) then
                  fkin=ykin(1)
                  dfkin=0.d0
               elseif(id.eq.nkin) then
                  fkin=ykin(nkin)
                  dfkin=0.d0
               else
                  dfkin=(ykin(id+1)-ykin(id))/(xkin(id+1)-xkin(id))
                  fkin=ykin(id)+dfkin*(ep-xkin(id))
               endif
            elseif(niso.ne.0) then
               fkin=yiso(1)
               dfkin=0.d0
            else
               fkin=0.d0
               dfkin=0.d0
            endif
         endif
!
         if(dabs(cop).lt.1.d-10) then
            fiso0=fiso
            fkin0=fkin
         endif
!
         if((kode.eq.-51).or.(ivisco.eq.0)) then
            dcop=(dxitril-c2*(fiso+fkin-fkin0)-um2*cop)/
     &           (um2+c1*(dfiso+dfkin))
         else
            if(user_creep.eq.1) then
               if(ithermal(1).eq.0) then
                  write(*,*) '*ERROR in incplas: no temperature defined'
                  call exit(201)
               endif
               timeabq(1)=time
               timeabq(2)=ttime+time
!
!              Eqn. (5.139)
!     
               qtild=(dxitril-um2*cop)/c2-fiso-(fkin-fkin0)
!
!              the Von Mises stress must be positive
!
               if(qtild.lt.1.d-10) qtild=1.d-10
               ec(1)=epini
               call creep(decra,deswa,xstateini(1,iint,iel),serd,ec,
     &             esw,p,qtild,t1l,dtemp,predef,dpred,timeabq,dtime,
     &             amat,leximp,lend,pgauss,nstate_,iel,iint,layer,kspt,
     &             kstep,kinc)
               dsvm=1.d0/decra(5)
               dcop=(decra(1)-c2*cop)/
     &                 (c2*(decra(5)*(dfiso+um*(3.d0+dfkin/um))+1.d0))
            else
               qtild=(dxitril-um2*cop)/c2-fiso-(fkin-fkin0)
!
!              the Von Mises stress must be positive
!
               if(qtild.lt.1.d-10) qtild=1.d-10
               decra(1)=a1*qtild**xxn
               decra(5)=xxn*decra(1)/qtild
               dsvm=1.d0/decra(5)
               dcop=(decra(1)-c2*cop)/
     &                 (c2*(decra(5)*(dfiso+um*(3.d0+dfkin/um))+1.d0))
            endif
         endif
         cop=cop+dcop
!
         if((dabs(dcop).lt.cop*1.d-4).or.
     &      (dabs(dcop).lt.1.d-10)) exit
!
!        check for endless loops or a negative consistency
!        parameter
!
         if((iloop.gt.15).or.(cop.le.0.d0)) then
            pnewdt=0.25d0
            return
         endif
!
      enddo loop
!
!        updating the equivalent plastic strain
!
      ep=epini+c2*cop
!
!     updating the back stress
!
      c7=c2*(fkin-fkin0)/dxitril
      do i=1,6
         stbl(i)=stbl(i)-c7*xitril(i)
      enddo
!
!     updating the plastic strain tensor
!
      c7=cop/dxitril
      do i=1,6
         epl(i)=epl(i)+c7*xitril(i)
      enddo
!
!     trial elastic strain
!     
      do i=1,6
         el(i)=emec(i)-epl(i)
      enddo
!
!     trial stress
!
      tracee=el(1)+el(2)+el(3)
      do i=1,6
         stre(i)=um2*el(i)-beta(i)
      enddo
      do i=1,3
         stre(i)=stre(i)+al*tracee
      enddo
!
!     calculating the local stiffness matrix
!
      if(icmd.ne.3) then
         c7=um2**2*cop/dxitril
!
         if((kode.eq.-51).or.(ivisco.eq.0)) then
            c8=um2**2/(um2+c1*(dfiso+dfkin))
         else
            c8=um2**2/(um2+c1*(dfiso+dfkin+1.d0/decra(5)))
         endif
!
         c1=um-c7/2.d0
         c2=al+c7/3.d0
         c3=(c7-c8)/(dxitril**2)
!
         xn(1,1)=xitril(1)
         xn(2,2)=xitril(2)
         xn(3,3)=xitril(3)
         xn(1,2)=xitril(4)
         xn(1,3)=xitril(5)
         xn(2,3)=xitril(6)
         xn(2,1)=xn(1,2)
         xn(3,1)=xn(1,3)
         xn(3,2)=xn(2,3)
!
         do i=1,21
            k=kel(1,i)
            l=kel(2,i)
            m=kel(3,i)
            n=kel(4,i)
            elas(i)=c1*(dkl(k,m)*dkl(l,n)+dkl(k,n)*dkl(l,m))
     &             +c2*dkl(k,l)*dkl(m,n)
     &             +c3*xn(k,l)*xn(m,n)
         enddo
      endif
!
!        updating the plastic fields
!
      do i=1,6
         xstate(1+i,iint,iel)=epl(i)
         xstate(7+i,iint,iel)=stbl(i)
      enddo
      xstate(1,iint,iel)=ep
!
!     maximum difference in the equivalent viscoplastic strain 
!     in this increment based on the viscoplastic strain rate at
!     the start and the end of the increment
!
      if(ivisco.eq.1) then
         depvisc=max(depvisc,abs(decra(1)-dtime*xstateini(14,iint,iel)))
         xstate(14,iint,iel)=decra(1)/dtime
      endif
!
!     maximum difference in the equivalent viscoplastic strain 
!     in this increment based on the viscoplastic strain rate at
!     the start and the end of the increment
!
      if(ivisco.eq.1) then
         depvisc=max(depvisc,abs(decra(1)-dtime*xstateini(14,iint,iel)))
         xstate(14,iint,iel)=decra(1)/dtime
      endif
c      depvisc=max(depvisc,ep-epini)
!
      return
      end
