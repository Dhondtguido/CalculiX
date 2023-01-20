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
      subroutine incplas(elconloc,plconloc,xstate,xstateini,
     &  elas,emec,ithermal,icmd,beta,stre,vj,kode,
     &  ielas,amat,t1l,dtime,time,ttime,iel,iint,nstate_,mi,
     &  eloc,pgauss,nmethod,pnewdt,depvisc)
!
!     calculates stiffness and stresses for the incremental plasticity
!     material law (Ref: J.C. Simo, A framework for finite strain
!     elastoplasticity, Comp. Meth. Appl. Mech. Engng., 66(1988)199-219
!     and 68(1988)1-31)
!
!     icmd=3: calculates stress at mechanical strain
!     else: calculates stress at mechanical strain and the stiffness
!           matrix
!
!     the stresses in the routine proposed by Simo are Kirchhoff 
!     stresses. Since the stress in the hardening laws are Chauchy
!     stresses, they are converted into Kirchhoff stress by 
!     multiplication with the Jacobian determinant
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,i,j,kode,ivisco,ielastic,
     &  niso,nkin,ielas,iel,iint,nstate_,mi(*),id,leximp,lend,layer,
     &  kspt,kstep,kinc,iloop,nmethod,user_hardening,user_creep
!
      real*8 elconloc(*),elas(21),emec(6),beta(6),stre(6),
     &  vj,plconloc(802),stbl(6),stril(6),xitril(6),
     &  ee,un,um,al,xk,cop,umb,umbb,dxitril,f0,d0,f1,d1,d2,xg(3,3),
     &  xs(3,3),xx(3,3),xn(3,3),xd(3,3),cpl(6),c(6),ci(6),
     &  c1,c2,c3,c4,c5,c6,c7,c8,c9,cplb(6),stblb(6),
     &  ftrial,xiso(200),yiso(200),xkin(200),ykin(200),
     &  fiso,dfiso,fkin,dfkin,fiso0,fkin0,ep,t1l,dtime,
     &  epini,a1,dsvm,xxa,xxn,vj2,vj23,sc(3,3),depvisc,
     &  cop1,cop2,fu1,fu2,fu,dcop,time,ttime,eloc(6),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,
     &  g17,g18,g28,g29,g30,g31,g32,g33,decra(5),deswa(5),serd,
     &  esw(2),ec(2),p,qtild,predef(1),dpred(1),timeabq(2),pgauss(3),
     &  dtemp,pnewdt
!
      leximp=1
      lend=2
      user_creep=0
!
!     localizing the plastic fields
!
      do i=1,6
         cpl(i)=-2.d0*xstateini(1+i,iint,iel)
         stbl(i)=xstateini(7+i,iint,iel)
      enddo
      do i=1,3
         cpl(i)=cpl(i)+1.d0
      enddo
      epini=xstateini(1,iint,iel)
!
      ee=elconloc(1)
      un=elconloc(2)
      um=ee/(1.d0+un)
      al=um*un/(1.d0-2.d0*un)
      xk=al+um/3.d0
      um=um/2.d0
!
      ep=epini
!
!     right Cauchy-Green tensor (eloc contains the Lagrange strain,
!     including thermal strain)
!
      c(1)=2.d0*emec(1)+1.d0
      c(2)=2.d0*emec(2)+1.d0
      c(3)=2.d0*emec(3)+1.d0
      c(4)=2.d0*emec(4)
      c(5)=2.d0*emec(5)
      c(6)=2.d0*emec(6)
!
!     calculating the Jacobian
!
      vj=c(1)*(c(2)*c(3)-c(6)*c(6))
     &  -c(4)*(c(4)*c(3)-c(6)*c(5))
     &  +c(5)*(c(4)*c(6)-c(2)*c(5))
      if(vj.gt.1.d-30) then
         vj=dsqrt(vj)
      else
ccc         write(*,*) '*WARNING in incplas: deformation inside-out'
!
!        deformation is reset to zero in order to continue the
!        calculation. Alternatively, a flag could be set forcing
!        a reiteration of the increment with a smaller size (to
!        be done)
!
         c(1)=1.d0
         c(2)=1.d0
         c(3)=1.d0
         c(4)=0.d0
         c(5)=0.d0
         c(6)=0.d0
         vj=1.d0
      endif
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
c      if(kode.eq.-52) then
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
!        inversion of the right Cauchy-Green tensor
!
      vj2=vj*vj
      ci(1)=(c(2)*c(3)-c(6)*c(6))/vj2
      ci(2)=(c(1)*c(3)-c(5)*c(5))/vj2
      ci(3)=(c(1)*c(2)-c(4)*c(4))/vj2
      ci(4)=(c(5)*c(6)-c(4)*c(3))/vj2
      ci(5)=(c(4)*c(6)-c(2)*c(5))/vj2
      ci(6)=(c(4)*c(5)-c(1)*c(6))/vj2
!
!        reducing the plastic right Cauchy-Green tensor and
!        the back stress to "isochoric" quantities (b stands
!        for bar)
!
      vj23=vj**(2.d0/3.d0)
      do i=1,6
         cplb(i)=cpl(i)/vj23
         stblb(i)=stbl(i)/vj23
      enddo
!
!        calculating the (n+1) trace and the (n+1) deviation of 
!        the (n) "isochoric" plastic right Cauchy-Green tensor
!
      umb=(c(1)*cplb(1)+c(2)*cplb(2)+c(3)*cplb(3)+
     &     2.d0*(c(4)*cplb(4)+c(5)*cplb(5)+c(6)*cplb(6)))/3.d0
      do i=1,6
         cplb(i)=cplb(i)-umb*ci(i)
      enddo
!
!        calculating the (n+1) trace and the (n+1) deviation of 
!        the (n) "isochoric" back stress tensor
!
      umbb=(c(1)*stblb(1)+c(2)*stblb(2)+c(3)*stblb(3)+
     &     2.d0*(c(4)*stblb(4)+c(5)*stblb(5)+c(6)*stblb(6)))/3.d0
      do i=1,6
         stblb(i)=stblb(i)-umbb*ci(i)
      enddo
!
!        calculating the trial stress
!
      do i=1,6
         stril(i)=um*cplb(i)-beta(i)
      enddo
!
!        calculating the trial radius vector of the yield surface
!
      do i=1,6
         xitril(i)=stril(i)-stblb(i)
      enddo
!
      sc(1,1)=xitril(1)*c(1)+xitril(4)*c(4)+xitril(5)*c(5)
      sc(1,2)=xitril(1)*c(4)+xitril(4)*c(2)+xitril(5)*c(6)
      sc(1,3)=xitril(1)*c(5)+xitril(4)*c(6)+xitril(5)*c(3)
      sc(2,1)=xitril(4)*c(1)+xitril(2)*c(4)+xitril(6)*c(5)
      sc(2,2)=xitril(4)*c(4)+xitril(2)*c(2)+xitril(6)*c(6)
      sc(2,3)=xitril(4)*c(5)+xitril(2)*c(6)+xitril(6)*c(3)
      sc(3,1)=xitril(5)*c(1)+xitril(6)*c(4)+xitril(3)*c(5)
      sc(3,2)=xitril(5)*c(4)+xitril(6)*c(2)+xitril(3)*c(6)
      sc(3,3)=xitril(5)*c(5)+xitril(6)*c(6)+xitril(3)*c(3)
!
      dxitril=sc(1,1)*sc(1,1)+sc(2,2)*sc(2,2)+sc(3,3)*sc(3,3)+
     &  2.d0*(sc(1,2)*sc(2,1)+sc(1,3)*sc(3,1)+sc(2,3)*sc(3,2))
c      g1=c(6)
c      g2=xitril(6)
c      g3=xitril(3)
c      g4=xitril(2)
c      g5=c(5)
c      g6=xitril(5)
c      g7=xitril(4)
c      g8=c(4)
c      g9=c(3)
c      g10=c(2)
c      g11=c(1)
c      g12=xitril(1)
c      g13=g12*g11
c      g14=g10*g4
c      g15=g9*g3
c      g16=g8*g7
c      g17=g6*g5
c      g18=g2*g1
c      g28=4*(g16 + g15)
c      g29=4*g13
c      g30=4*g14
c      g31=4*g6*g1
c      g32=4*g8*g2
c      g33=4*g7*g5
c      dxitril=(g13*g13 + g14*g14 + g15*g15 + g16*(g30 + g29 + 2*
c     &     g16) + g17*(g29 + g28 + 2*g17) + g18*(g30 + g28 + 2*
c     &     g18 + 4*g17) + g11*g7*(g31 + 2*g10*g7) + g9*g6*(g32 + 
c     &     2*g11*g6) + g10*g2*(g33 + 2*g9*g2) + g8*g4*(g31 + 2*
c     &     g12*g8) + g12*g5*(g32 + 2*g5*g3) + g3*g1*(g33 + 2*g4*
c     &     g1))
!
!     in principal dxitril is a norm and cannot be less than zero
!     however, due to round-off it can happen anyway.
!     In order to avoid this, one would have to use the deformation
!     gradient and calculate tau_ij*tau_kl*g_ik*g_jl instead of
!     S_IJ*S_KL*C_IK*C_JL, however, this would increase the
!     calculational effort substantially.
!
      if(dxitril.lt.0.d0) then
ccc         write(*,*) '*WARNING in incplas: dxitril < 0'
         dxitril=0.d0
      else
         dxitril=dsqrt(dxitril)
      endif
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
            yiso(i)=vj*plconloc(2*i)
         enddo
      endif
      if(nkin.ne.0) then
         do i=1,nkin
            xkin(i)=plconloc(399+2*i)
            ykin(i)=vj*plconloc(400+2*i)
         enddo
      endif
!
!     if no viscous calculation is performed a pure creep calculation
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
         fiso=fiso*vj
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
      ftrial=dxitril-dsqrt(2.d0/3.d0)*fiso
      if((ftrial.le.1.d-10).or.(ielas.eq.1).or.(ielastic.eq.1)
     &    .or.(dtime.lt.1.d-30)) then
!
!        no plastic deformation
!        beta contains the Cauchy residual stresses
!
c         write(*,*) 'no plastic deformation'
         c8=xk*(vj2-1.d0)/2.d0
!
!           residual stresses are de facto PK2 stresses
!           (Piola-Kirchhoff of the second kind)
!
         stre(1)=c8*ci(1)+stril(1)-beta(1)
         stre(2)=c8*ci(2)+stril(2)-beta(2)
         stre(3)=c8*ci(3)+stril(3)-beta(3)
         stre(4)=c8*ci(4)+stril(4)-beta(4)
         stre(5)=c8*ci(5)+stril(5)-beta(5)
         stre(6)=c8*ci(6)+stril(6)-beta(6)
!
!        updating the plastic fields
!
         do i=1,3
            cpl(i)=cpl(i)-1.d0
         enddo
         do i=1,6
            xstate(1+i,iint,iel)=-cpl(i)/2.d0
            xstate(7+i,iint,iel)=stbl(i)
         enddo
         xstate(1,iint,iel)=ep
!
         if(icmd.ne.3) then
!
            umb=um*umb
!
!           calculating the local stiffness matrix
!
            xg(1,1)=(c(2)*c(3)-c(6)*c(6))/vj2
            xg(2,2)=(c(1)*c(3)-c(5)*c(5))/vj2
            xg(3,3)=(c(1)*c(2)-c(4)*c(4))/vj2
            xg(1,2)=(c(5)*c(6)-c(4)*c(3))/vj2
            xg(1,3)=(c(4)*c(6)-c(2)*c(5))/vj2
            xg(2,3)=(c(4)*c(5)-c(1)*c(6))/vj2
            xg(2,1)=xg(1,2)
            xg(3,1)=xg(1,3)
            xg(3,2)=xg(2,3)
!
            xs(1,1)=stril(1)
            xs(2,2)=stril(2)
            xs(3,3)=stril(3)
            xs(1,2)=stril(4)
            xs(2,1)=stril(4)
            xs(1,3)=stril(5)
            xs(3,1)=stril(5)
            xs(2,3)=stril(6)
            xs(3,2)=stril(6)
!
            elas(1)=umb*(xg(1,1)*xg(1,1)+xg(1,1)*xg(1,1)-
     &           2.d0*xg(1,1)*xg(1,1)/3.d0)
     &           -2.d0*(xs(1,1)*xg(1,1)+xg(1,1)*xs(1,1))/3.d0
     &           +xk*vj2*xg(1,1)*xg(1,1)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(1,1)
     &           +xg(1,1)*xg(1,1))/2.d0
            elas(2)=umb*(xg(1,2)*xg(1,2)+xg(1,2)*xg(1,2)-
     &           2.d0*xg(1,1)*xg(2,2)/3.d0)
     &           -2.d0*(xs(1,1)*xg(2,2)+xg(1,1)*xs(2,2))/3.d0
     &           +xk*vj2*xg(1,1)*xg(2,2)
     &           -xk*(vj2-1.d0)*(xg(1,2)*xg(1,2)
     &           +xg(1,2)*xg(1,2))/2.d0
            elas(3)=umb*(xg(2,2)*xg(2,2)+xg(2,2)*xg(2,2)-
     &           2.d0*xg(2,2)*xg(2,2)/3.d0)
     &           -2.d0*(xs(2,2)*xg(2,2)+xg(2,2)*xs(2,2))/3.d0
     &           +xk*vj2*xg(2,2)*xg(2,2)
     &           -xk*(vj2-1.d0)*(xg(2,2)*xg(2,2)
     &           +xg(2,2)*xg(2,2))/2.d0
            elas(4)=umb*(xg(1,3)*xg(1,3)+xg(1,3)*xg(1,3)-
     &           2.d0*xg(1,1)*xg(3,3)/3.d0)
     &           -2.d0*(xs(1,1)*xg(3,3)+xg(1,1)*xs(3,3))/3.d0
     &           +xk*vj2*xg(1,1)*xg(3,3)
     &           -xk*(vj2-1.d0)*(xg(1,3)*xg(1,3)
     &           +xg(1,3)*xg(1,3))/2.d0
            elas(5)=umb*(xg(2,3)*xg(2,3)+xg(2,3)*xg(2,3)-
     &           2.d0*xg(2,2)*xg(3,3)/3.d0)
     &           -2.d0*(xs(2,2)*xg(3,3)+xg(2,2)*xs(3,3))/3.d0
     &           +xk*vj2*xg(2,2)*xg(3,3)
     &           -xk*(vj2-1.d0)*(xg(2,3)*xg(2,3)
     &           +xg(2,3)*xg(2,3))/2.d0
            elas(6)=umb*(xg(3,3)*xg(3,3)+xg(3,3)*xg(3,3)-
     &           2.d0*xg(3,3)*xg(3,3)/3.d0)
     &           -2.d0*(xs(3,3)*xg(3,3)+xg(3,3)*xs(3,3))/3.d0
     &           +xk*vj2*xg(3,3)*xg(3,3)
     &           -xk*(vj2-1.d0)*(xg(3,3)*xg(3,3)
     &           +xg(3,3)*xg(3,3))/2.d0
            elas(7)=umb*(xg(1,1)*xg(1,2)+xg(1,2)*xg(1,1)-
     &           2.d0*xg(1,1)*xg(1,2)/3.d0)
     &           -2.d0*(xs(1,1)*xg(1,2)+xg(1,1)*xs(1,2))/3.d0
     &           +xk*vj2*xg(1,1)*xg(1,2)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(1,2)
     &           +xg(1,2)*xg(1,1))/2.d0
            elas(8)=umb*(xg(2,1)*xg(2,2)+xg(2,2)*xg(2,1)-
     &           2.d0*xg(2,2)*xg(1,2)/3.d0)
     &           -2.d0*(xs(2,2)*xg(1,2)+xg(2,2)*xs(1,2))/3.d0
     &           +xk*vj2*xg(2,2)*xg(1,2)
     &           -xk*(vj2-1.d0)*(xg(2,1)*xg(2,2)
     &           +xg(2,2)*xg(2,1))/2.d0
            elas(9)=umb*(xg(3,1)*xg(3,2)+xg(3,2)*xg(3,1)-
     &           2.d0*xg(3,3)*xg(1,2)/3.d0)
     &           -2.d0*(xs(3,3)*xg(1,2)+xg(3,3)*xs(1,2))/3.d0
     &           +xk*vj2*xg(3,3)*xg(1,2)
     &           -xk*(vj2-1.d0)*(xg(3,1)*xg(3,2)
     &           +xg(3,2)*xg(3,1))/2.d0
            elas(10)=umb*(xg(1,1)*xg(2,2)+xg(1,2)*xg(2,1)-
     &           2.d0*xg(1,2)*xg(1,2)/3.d0)
     &           -2.d0*(xs(1,2)*xg(1,2)+xg(1,2)*xs(1,2))/3.d0
     &           +xk*vj2*xg(1,2)*xg(1,2)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(2,2)
     &           +xg(1,2)*xg(2,1))/2.d0
            elas(11)=umb*(xg(1,1)*xg(1,3)+xg(1,3)*xg(1,1)-
     &           2.d0*xg(1,1)*xg(1,3)/3.d0)
     &           -2.d0*(xs(1,1)*xg(1,3)+xg(1,1)*xs(1,3))/3.d0
     &           +xk*vj2*xg(1,1)*xg(1,3)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(1,3)
     &           +xg(1,3)*xg(1,1))/2.d0
            elas(12)=umb*(xg(2,1)*xg(2,3)+xg(2,3)*xg(2,1)-
     &           2.d0*xg(2,2)*xg(1,3)/3.d0)
     &           -2.d0*(xs(2,2)*xg(1,3)+xg(2,2)*xs(1,3))/3.d0
     &           +xk*vj2*xg(2,2)*xg(1,3)
     &           -xk*(vj2-1.d0)*(xg(2,1)*xg(2,3)
     &           +xg(2,3)*xg(2,1))/2.d0
            elas(13)=umb*(xg(3,1)*xg(3,3)+xg(3,3)*xg(3,1)-
     &           2.d0*xg(3,3)*xg(1,3)/3.d0)
     &           -2.d0*(xs(3,3)*xg(1,3)+xg(3,3)*xs(1,3))/3.d0
     &           +xk*vj2*xg(3,3)*xg(1,3)
     &           -xk*(vj2-1.d0)*(xg(3,1)*xg(3,3)
     &           +xg(3,3)*xg(3,1))/2.d0
            elas(14)=umb*(xg(1,1)*xg(2,3)+xg(1,3)*xg(2,1)-
     &           2.d0*xg(1,2)*xg(1,3)/3.d0)
     &           -2.d0*(xs(1,2)*xg(1,3)+xg(1,2)*xs(1,3))/3.d0
     &           +xk*vj2*xg(1,2)*xg(1,3)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(2,3)
     &           +xg(1,3)*xg(2,1))/2.d0
            elas(15)=umb*(xg(1,1)*xg(3,3)+xg(1,3)*xg(3,1)-
     &           2.d0*xg(1,3)*xg(1,3)/3.d0)
     &           -2.d0*(xs(1,3)*xg(1,3)+xg(1,3)*xs(1,3))/3.d0
     &           +xk*vj2*xg(1,3)*xg(1,3)
     &           -xk*(vj2-1.d0)*(xg(1,1)*xg(3,3)
     &           +xg(1,3)*xg(3,1))/2.d0
            elas(16)=umb*(xg(1,2)*xg(1,3)+xg(1,3)*xg(1,2)-
     &           2.d0*xg(1,1)*xg(2,3)/3.d0)
     &           -2.d0*(xs(1,1)*xg(2,3)+xg(1,1)*xs(2,3))/3.d0
     &           +xk*vj2*xg(1,1)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(1,2)*xg(1,3)
     &           +xg(1,3)*xg(1,2))/2.d0
            elas(17)=umb*(xg(2,2)*xg(2,3)+xg(2,3)*xg(2,2)-
     &           2.d0*xg(2,2)*xg(2,3)/3.d0)
     &           -2.d0*(xs(2,2)*xg(2,3)+xg(2,2)*xs(2,3))/3.d0
     &           +xk*vj2*xg(2,2)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(2,2)*xg(2,3)
     &           +xg(2,3)*xg(2,2))/2.d0
            elas(18)=umb*(xg(3,2)*xg(3,3)+xg(3,3)*xg(3,2)-
     &           2.d0*xg(3,3)*xg(2,3)/3.d0)
     &           -2.d0*(xs(3,3)*xg(2,3)+xg(3,3)*xs(2,3))/3.d0
     &           +xk*vj2*xg(3,3)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(3,2)*xg(3,3)
     &           +xg(3,3)*xg(3,2))/2.d0
            elas(19)=umb*(xg(1,2)*xg(2,3)+xg(1,3)*xg(2,2)-
     &           2.d0*xg(1,2)*xg(2,3)/3.d0)
     &           -2.d0*(xs(1,2)*xg(2,3)+xg(1,2)*xs(2,3))/3.d0
     &           +xk*vj2*xg(1,2)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(1,2)*xg(2,3)
     &           +xg(1,3)*xg(2,2))/2.d0
            elas(20)=umb*(xg(1,2)*xg(3,3)+xg(1,3)*xg(3,2)-
     &           2.d0*xg(1,3)*xg(2,3)/3.d0)
     &           -2.d0*(xs(1,3)*xg(2,3)+xg(1,3)*xs(2,3))/3.d0
     &           +xk*vj2*xg(1,3)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(1,2)*xg(3,3)
     &           +xg(1,3)*xg(3,2))/2.d0
            elas(21)=umb*(xg(2,2)*xg(3,3)+xg(2,3)*xg(3,2)-
     &           2.d0*xg(2,3)*xg(2,3)/3.d0)
     &           -2.d0*(xs(2,3)*xg(2,3)+xg(2,3)*xs(2,3))/3.d0
     &           +xk*vj2*xg(2,3)*xg(2,3)
     &           -xk*(vj2-1.d0)*(xg(2,2)*xg(3,3)
     &           +xg(2,3)*xg(3,2))/2.d0
!
         endif
!
         return
      endif
!
!        plastic deformation
!
      umb=um*umb
      umbb=umb-umbb
!
!        calculating the consistency parameter
!
      c1=2.d0/3.d0
      c2=dsqrt(c1)
      c3=c1/um
      c4=c2/um
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
            fiso=fiso*vj
            dfiso=dfiso*vj
            fkin=fkin*vj
            dfkin=dfkin*vj
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
            dcop=(ftrial-c2*(fiso-fiso0)
     &           -umbb*(2.d0*cop+c4*(fkin-fkin0)))/
     &           (-c1*dfiso-umbb*(2.d0+c3*dfkin))
         else
            if(user_creep.eq.1) then
               if(ithermal(1).eq.0) then
ccc                  write(*,*) '*ERROR in incplas: no temperature defined'
                  call exit(201)
               endif
               timeabq(1)=time
               timeabq(2)=ttime+time
               qtild=(ftrial-c2*(fiso-fiso0)
     &              -umbb*(2.d0*cop+c4*(fkin-fkin0)))/(c2*vj)
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
               dcop=-(decra(1)-c2*cop)/
     &                 (c2*(decra(5)*(dfiso+umbb*(3.d0+dfkin/um))+1.d0))
            else
               qtild=(ftrial-c2*(fiso-fiso0)
     &              -umbb*(2.d0*cop+c4*(fkin-fkin0)))/(c2*vj)
!
!              the Von Mises stress must be positive
!
               if(qtild.lt.1.d-10) qtild=1.d-10
               decra(1)=a1*qtild**xxn
               decra(5)=xxn*decra(1)/qtild
               dsvm=1.d0/decra(5)
               dcop=-(decra(1)-c2*cop)/
     &                 (c2*(decra(5)*(dfiso+umbb*(3.d0+dfkin/um))+1.d0))
            endif
         endif
         cop=cop-dcop
!
         if((dabs(dcop).lt.cop*1.d-4).or.
     &      (dabs(dcop).lt.1.d-10)) exit
!
!        check for endless loops or a negative consistency
!        parameter
!
         if((iloop.gt.15).or.(cop.le.0.d0)) then
            iloop=1
            cop=0.d0
            do
               if(iloop.gt.100) then
c                  NOTE: write statements cause problems for
c                        parallellized execution
c                  write(*,*) '*WARNING in incplas: material loop'
c                  write(*,*) '         did not converge in integration'
c                  write(*,*) '         point',iint,'in element',iel,';'
c                  write(*,*) '         the increment size is reduced'
c                  write(*,*)
                  pnewdt=0.25d0
                  return
               endif
               ep=epini+c2*cop
!
               if(user_hardening.eq.1) then
                  call uhardening(amat,iel,iint,t1l,epini,ep,dtime,
     &                 fiso,dfiso,fkin,dfkin)
                  fiso=fiso*vj
                  fkin=fkin*vj
               else
                  if(niso.ne.0) then
                     call ident(xiso,ep,niso,id)
                     if(id.eq.0) then
                        fiso=yiso(1)
                     elseif(id.eq.niso) then
                        fiso=yiso(niso)
                     else
                        dfiso=(yiso(id+1)-yiso(id))/
     &                        (xiso(id+1)-xiso(id))
                        fiso=yiso(id)+dfiso*(ep-xiso(id))
                     endif
                  elseif(nkin.ne.0) then
                     fiso=ykin(1)
                  else
                     fiso=0.d0
                  endif
!
                  if(nkin.ne.0) then
                     call ident(xkin,ep,nkin,id)
                     if(id.eq.0) then
                        fkin=ykin(1)
                     elseif(id.eq.nkin) then
                        fkin=ykin(nkin)
                     else
                        dfkin=(ykin(id+1)-ykin(id))/
     &                        (xkin(id+1)-xkin(id))
                        fkin=ykin(id)+dfkin*(ep-xkin(id))
                     endif
                  elseif(niso.ne.0) then
                     fkin=yiso(1)
                  else
                     fkin=0.d0
                  endif
               endif
!
               if(dabs(cop).lt.1.d-10) then
                  fiso0=fiso
                  fkin0=fkin
               endif
!
               if((kode.eq.-51).or.(ivisco.eq.0)) then
                  fu=(ftrial-c2*(fiso-fiso0)
     &                 -umbb*(2.d0*cop+c4*(fkin-fkin0)))
               else
                  if(user_creep.eq.1) then
                     timeabq(1)=time
                     timeabq(2)=ttime+time
                     qtild=(ftrial-c2*(fiso-fiso0)
     &                    -umbb*(2.d0*cop+c4*(fkin-fkin0)))/(c2*vj)
!
!                    the Von Mises stress must be positive
!
                     if(qtild.lt.1.d-10) qtild=1.d-10
                     ec(1)=epini
                     call creep(decra,deswa,xstateini(1,iint,iel),serd,
     &                    ec,esw,p,qtild,t1l,dtemp,predef,dpred,timeabq,
     &                    dtime,amat,leximp,lend,pgauss,nstate_,iel,
     &                    iint,layer,kspt,kstep,kinc)
                     dsvm=1.d0/decra(5)
                     fu=decra(1)-c2*cop
                  else
                     qtild=(ftrial-c2*(fiso-fiso0)
     &                    -umbb*(2.d0*cop+c4*(fkin-fkin0)))/(c2*vj)
!
!                    the Von Mises stress must be positive
!
                     if(qtild.lt.1.d-10) qtild=1.d-10
                     decra(1)=a1*qtild**xxn
                     decra(5)=xxn*decra(1)/qtild
                     dsvm=1.d0/decra(5)
                     fu=decra(1)-c2*cop
                  endif
               endif
!
               if(iloop.eq.1) then
c                  write(*,*) 'cop,fu ',cop,fu
                  cop1=0.d0
                  fu1=fu
                  iloop=2
                  cop=1.d-10
               elseif(iloop.eq.2) then
                  if(fu*fu1.le.0.d0) then
c                     write(*,*) cop,fu
                     iloop=3
                     fu2=fu
                     cop2=cop
                     cop=(cop1+cop2)/2.d0
                     dcop=(cop2-cop1)/2.d0
                  else
c                     write(*,*) cop,fu
                     cop=cop*10.d0
                     if(cop.gt.100.d0) then
                        write(*,*) '*ERROR: no convergence in incplas'
                        pnewdt=0.25d0
                        return
c                        call exit(201)
                     endif
                  endif
               else
c                     write(*,*) cop,fu
                  if(fu*fu1.ge.0.d0) then
                     cop1=cop
                     fu1=fu
                  else
                     cop2=cop
                     fu2=fu
                  endif
                  cop=(cop1+cop2)/2.d0
                  dcop=(cop2-cop1)/2.d0
                  if((dabs(dcop).lt.cop*1.d-4).or.
     &                 (dabs(dcop).lt.1.d-10)) exit loop
               endif
            enddo
         endif
!
      enddo loop
!
!        updating the equivalent plastic strain
!
      ep=epini+c2*cop
!
!        updating the back stress
!
      c5=2.d0*umbb*cop/dxitril
      c6=c5/(3.d0*um)
      c7=c6*dfkin*vj23
      do i=1,6
         stbl(i)=stbl(i)+c7*xitril(i)
      enddo
!
!        updating the stress
!        vj: Jacobian of the total deformation gradient
!
      c8=xk*(vj2-1.d0)/2.d0
!
      do i=1,6
         stre(i)=c8*ci(i)-beta(i)+stril(i)-c5*xitril(i)
      enddo
!
!        updating the plastic right Cauchy-Green tensor
!
      c9=c6*3.d0*vj23
      do i=1,6
         cpl(i)=cpl(i)-c9*xitril(i)
      enddo
!
      if(icmd.ne.3) then
!
!        calculating the local stiffness matrix
!
         xg(1,1)=(c(2)*c(3)-c(6)*c(6))/vj2
         xg(2,2)=(c(1)*c(3)-c(5)*c(5))/vj2
         xg(3,3)=(c(1)*c(2)-c(4)*c(4))/vj2
         xg(1,2)=(c(5)*c(6)-c(4)*c(3))/vj2
         xg(1,3)=(c(4)*c(6)-c(2)*c(5))/vj2
         xg(2,3)=(c(4)*c(5)-c(1)*c(6))/vj2
         xg(2,1)=xg(1,2)
         xg(3,1)=xg(1,3)
         xg(3,2)=xg(2,3)
!                                
         xs(1,1)=stril(1)
         xs(2,2)=stril(2)
         xs(3,3)=stril(3)
         xs(1,2)=stril(4)
         xs(2,1)=stril(4)
         xs(1,3)=stril(5)
         xs(3,1)=stril(5)
         xs(2,3)=stril(6)
         xs(3,2)=stril(6)
!
         f0=2.d0*umbb*cop/dxitril
         d0=1.d0+(dfkin/um+dfiso/umbb)/3.d0
!
!        creep contribution
!
         if((kode.eq.-52).and.(ivisco.eq.1)) then
            d0=d0+dsvm/(3.d0*umbb)
         endif
!
         f1=1.d0/d0-f0
         d1=2.d0*f1*umbb-((1.d0+dfkin/(3.d0*um))/d0-1.d0)*
     &      4.d0*cop*dxitril/3.d0
         d2=2d0*dxitril*f1
!
         xx(1,1)=xitril(1)
         xx(2,2)=xitril(2)
         xx(3,3)=xitril(3)
         xx(1,2)=xitril(4)
         xx(2,1)=xitril(4)
         xx(1,3)=xitril(5)
         xx(3,1)=xitril(5)
         xx(2,3)=xitril(6)
         xx(3,2)=xitril(6)
!
         xn(1,1)=xitril(1)/dxitril
         xn(2,2)=xitril(2)/dxitril
         xn(3,3)=xitril(3)/dxitril
         xn(1,2)=xitril(4)/dxitril
         xn(2,1)=xitril(4)/dxitril
         xn(1,3)=xitril(5)/dxitril
         xn(3,1)=xitril(5)/dxitril
         xn(2,3)=xitril(6)/dxitril
         xn(3,2)=xitril(6)/dxitril
!
         do i=1,3
            do j=i,3
               xd(i,j)=xn(i,1)*xn(1,j)*c(1)+xn(i,1)*xn(2,j)*c(4)+
     &                 xn(i,1)*xn(3,j)*c(5)+xn(i,2)*xn(1,j)*c(4)+
     &                 xn(i,2)*xn(2,j)*c(2)+xn(i,2)*xn(3,j)*c(6)+
     &                 xn(i,3)*xn(1,j)*c(5)+xn(i,3)*xn(2,j)*c(6)+
     &                 xn(i,3)*xn(3,j)*c(3)
            enddo
         enddo
         xd(2,1)=xd(1,2)
         xd(3,1)=xd(1,3)
         xd(3,2)=xd(2,3)
!
!        deviatoric part
!
         c1=(xd(1,1)*c(1)+xd(2,2)*c(2)+xd(3,3)*c(3)+
     &      2.d0*(xd(1,2)*c(4)+xd(1,3)*c(5)+xd(2,3)*c(6)))/3.d0
         do i=1,3
            do j=i,3
               xd(i,j)=xd(i,j)-c1*xg(i,j)
            enddo
         enddo
         xd(2,1)=xd(1,2)
         xd(3,1)=xd(1,3)
         xd(3,2)=xd(2,3)
!
         elas(1)=(umb-f0*umbb)*(xg(1,1)*xg(1,1)+xg(1,1)*xg(1,1)-
     &        2.d0*xg(1,1)*xg(1,1)/3.d0)
     &        -2.d0*(xs(1,1)*xg(1,1)+xg(1,1)*xs(1,1))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(1,1)+xg(1,1)*xx(1,1))/3.d0
     &        -d1*xn(1,1)*xn(1,1)-d2*(xn(1,1)*xd(1,1)+
     &        xd(1,1)*xn(1,1))/2.d0+xk*vj2*xg(1,1)*xg(1,1)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(1,1)+xg(1,1)*xg(1,1))/2.d0
         elas(2)=(umb-f0*umbb)*(xg(1,2)*xg(1,2)+xg(1,2)*xg(1,2)-
     &        2.d0*xg(1,1)*xg(2,2)/3.d0)
     &        -2.d0*(xs(1,1)*xg(2,2)+xg(1,1)*xs(2,2))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(2,2)+xg(1,1)*xx(2,2))/3.d0
     &        -d1*xn(1,1)*xn(2,2)-d2*(xn(1,1)*xd(2,2)+
     &        xd(1,1)*xn(2,2))/2.d0+xk*vj2*xg(1,1)*xg(2,2)
     &        -xk*(vj2-1.d0)*(xg(1,2)*xg(1,2)+xg(1,2)*xg(1,2))/2.d0
         elas(3)=(umb-f0*umbb)*(xg(2,2)*xg(2,2)+xg(2,2)*xg(2,2)-
     &        2.d0*xg(2,2)*xg(2,2)/3.d0)
     &        -2.d0*(xs(2,2)*xg(2,2)+xg(2,2)*xs(2,2))/3.d0
     &        +f0*2.d0*(xx(2,2)*xg(2,2)+xg(2,2)*xx(2,2))/3.d0
     &        -d1*xn(2,2)*xn(2,2)-d2*(xn(2,2)*xd(2,2)+
     &        xd(2,2)*xn(2,2))/2.d0+xk*vj2*xg(2,2)*xg(2,2)
     &        -xk*(vj2-1.d0)*(xg(2,2)*xg(2,2)+xg(2,2)*xg(2,2))/2.d0
         elas(4)=(umb-f0*umbb)*(xg(1,3)*xg(1,3)+xg(1,3)*xg(1,3)-
     &        2.d0*xg(1,1)*xg(3,3)/3.d0)
     &        -2.d0*(xs(1,1)*xg(3,3)+xg(1,1)*xs(3,3))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(3,3)+xg(1,1)*xx(3,3))/3.d0
     &        -d1*xn(1,1)*xn(3,3)-d2*(xn(1,1)*xd(3,3)+
     &        xd(1,1)*xn(3,3))/2.d0+xk*vj2*xg(1,1)*xg(3,3)
     &        -xk*(vj2-1.d0)*(xg(1,3)*xg(1,3)+xg(1,3)*xg(1,3))/2.d0
         elas(5)=(umb-f0*umbb)*(xg(2,3)*xg(2,3)+xg(2,3)*xg(2,3)-
     &        2.d0*xg(2,2)*xg(3,3)/3.d0)
     &        -2.d0*(xs(2,2)*xg(3,3)+xg(2,2)*xs(3,3))/3.d0
     &        +f0*2.d0*(xx(2,2)*xg(3,3)+xg(2,2)*xx(3,3))/3.d0
     &        -d1*xn(2,2)*xn(3,3)-d2*(xn(2,2)*xd(3,3)+
     &        xd(2,2)*xn(3,3))/2.d0+xk*vj2*xg(2,2)*xg(3,3)
     &        -xk*(vj2-1.d0)*(xg(2,3)*xg(2,3)+xg(2,3)*xg(2,3))/2.d0
         elas(6)=(umb-f0*umbb)*(xg(3,3)*xg(3,3)+xg(3,3)*xg(3,3)-
     &        2.d0*xg(3,3)*xg(3,3)/3.d0)
     &        -2.d0*(xs(3,3)*xg(3,3)+xg(3,3)*xs(3,3))/3.d0
     &        +f0*2.d0*(xx(3,3)*xg(3,3)+xg(3,3)*xx(3,3))/3.d0
     &        -d1*xn(3,3)*xn(3,3)-d2*(xn(3,3)*xd(3,3)+
     &        xd(3,3)*xn(3,3))/2.d0+xk*vj2*xg(3,3)*xg(3,3)
     &        -xk*(vj2-1.d0)*(xg(3,3)*xg(3,3)+xg(3,3)*xg(3,3))/2.d0
         elas(7)=(umb-f0*umbb)*(xg(1,1)*xg(1,2)+xg(1,2)*xg(1,1)-
     &        2.d0*xg(1,1)*xg(1,2)/3.d0)
     &        -2.d0*(xs(1,1)*xg(1,2)+xg(1,1)*xs(1,2))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(1,2)+xg(1,1)*xx(1,2))/3.d0
     &        -d1*xn(1,1)*xn(1,2)-d2*(xn(1,1)*xd(1,2)+
     &        xd(1,1)*xn(1,2))/2.d0+xk*vj2*xg(1,1)*xg(1,2)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(1,2)+xg(1,2)*xg(1,1))/2.d0
         elas(8)=(umb-f0*umbb)*(xg(2,1)*xg(2,2)+xg(2,2)*xg(2,1)-
     &        2.d0*xg(2,2)*xg(1,2)/3.d0)
     &        -2.d0*(xs(2,2)*xg(1,2)+xg(2,2)*xs(1,2))/3.d0
     &        +f0*2.d0*(xx(2,2)*xg(1,2)+xg(2,2)*xx(1,2))/3.d0
     &        -d1*xn(2,2)*xn(1,2)-d2*(xn(2,2)*xd(1,2)+
     &        xd(2,2)*xn(1,2))/2.d0+xk*vj2*xg(2,2)*xg(1,2)
     &        -xk*(vj2-1.d0)*(xg(2,1)*xg(2,2)+xg(2,2)*xg(2,1))/2.d0
         elas(9)=(umb-f0*umbb)*(xg(3,1)*xg(3,2)+xg(3,2)*xg(3,1)-
     &        2.d0*xg(3,3)*xg(1,2)/3.d0)
     &        -2.d0*(xs(3,3)*xg(1,2)+xg(3,3)*xs(1,2))/3.d0
     &        +f0*2.d0*(xx(3,3)*xg(1,2)+xg(3,3)*xx(1,2))/3.d0
     &        -d1*xn(3,3)*xn(1,2)-d2*(xn(3,3)*xd(1,2)+
     &        xd(3,3)*xn(1,2))/2.d0+xk*vj2*xg(3,3)*xg(1,2)
     &        -xk*(vj2-1.d0)*(xg(3,1)*xg(3,2)+xg(3,2)*xg(3,1))/2.d0
         elas(10)=(umb-f0*umbb)*(xg(1,1)*xg(2,2)+xg(1,2)*xg(2,1)-
     &        2.d0*xg(1,2)*xg(1,2)/3.d0)
     &        -2.d0*(xs(1,2)*xg(1,2)+xg(1,2)*xs(1,2))/3.d0
     &        +f0*2.d0*(xx(1,2)*xg(1,2)+xg(1,2)*xx(1,2))/3.d0
     &        -d1*xn(1,2)*xn(1,2)-d2*(xn(1,2)*xd(1,2)+
     &        xd(1,2)*xn(1,2))/2.d0+xk*vj2*xg(1,2)*xg(1,2)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(2,2)+xg(1,2)*xg(2,1))/2.d0
         elas(11)=(umb-f0*umbb)*(xg(1,1)*xg(1,3)+xg(1,3)*xg(1,1)-
     &        2.d0*xg(1,1)*xg(1,3)/3.d0)
     &        -2.d0*(xs(1,1)*xg(1,3)+xg(1,1)*xs(1,3))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(1,3)+xg(1,1)*xx(1,3))/3.d0
     &        -d1*xn(1,1)*xn(1,3)-d2*(xn(1,1)*xd(1,3)+
     &        xd(1,1)*xn(1,3))/2.d0+xk*vj2*xg(1,1)*xg(1,3)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(1,3)+xg(1,3)*xg(1,1))/2.d0
         elas(12)=(umb-f0*umbb)*(xg(2,1)*xg(2,3)+xg(2,3)*xg(2,1)-
     &        2.d0*xg(2,2)*xg(1,3)/3.d0)
     &        -2.d0*(xs(2,2)*xg(1,3)+xg(2,2)*xs(1,3))/3.d0
     &        +f0*2.d0*(xx(2,2)*xg(1,3)+xg(2,2)*xx(1,3))/3.d0
     &        -d1*xn(2,2)*xn(1,3)-d2*(xn(2,2)*xd(1,3)+
     &        xd(2,2)*xn(1,3))/2.d0+xk*vj2*xg(2,2)*xg(1,3)
     &        -xk*(vj2-1.d0)*(xg(2,1)*xg(2,3)+xg(2,3)*xg(2,1))/2.d0
         elas(13)=(umb-f0*umbb)*(xg(3,1)*xg(3,3)+xg(3,3)*xg(3,1)-
     &        2.d0*xg(3,3)*xg(1,3)/3.d0)
     &        -2.d0*(xs(3,3)*xg(1,3)+xg(3,3)*xs(1,3))/3.d0
     &        +f0*2.d0*(xx(3,3)*xg(1,3)+xg(3,3)*xx(1,3))/3.d0
     &        -d1*xn(3,3)*xn(1,3)-d2*(xn(3,3)*xd(1,3)+
     &        xd(3,3)*xn(1,3))/2.d0+xk*vj2*xg(3,3)*xg(1,3)
     &        -xk*(vj2-1.d0)*(xg(3,1)*xg(3,3)+xg(3,3)*xg(3,1))/2.d0
         elas(14)=(umb-f0*umbb)*(xg(1,1)*xg(2,3)+xg(1,3)*xg(2,1)-
     &        2.d0*xg(1,2)*xg(1,3)/3.d0)
     &        -2.d0*(xs(1,2)*xg(1,3)+xg(1,2)*xs(1,3))/3.d0
     &        +f0*2.d0*(xx(1,2)*xg(1,3)+xg(1,2)*xx(1,3))/3.d0
     &        -d1*xn(1,2)*xn(1,3)-d2*(xn(1,2)*xd(1,3)+
     &        xd(1,2)*xn(1,3))/2.d0+xk*vj2*xg(1,2)*xg(1,3)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(2,3)+xg(1,3)*xg(2,1))/2.d0
         elas(15)=(umb-f0*umbb)*(xg(1,1)*xg(3,3)+xg(1,3)*xg(3,1)-
     &        2.d0*xg(1,3)*xg(1,3)/3.d0)
     &        -2.d0*(xs(1,3)*xg(1,3)+xg(1,3)*xs(1,3))/3.d0
     &        +f0*2.d0*(xx(1,3)*xg(1,3)+xg(1,3)*xx(1,3))/3.d0
     &        -d1*xn(1,3)*xn(1,3)-d2*(xn(1,3)*xd(1,3)+
     &        xd(1,3)*xn(1,3))/2.d0+xk*vj2*xg(1,3)*xg(1,3)
     &        -xk*(vj2-1.d0)*(xg(1,1)*xg(3,3)+xg(1,3)*xg(3,1))/2.d0
         elas(16)=(umb-f0*umbb)*(xg(1,2)*xg(1,3)+xg(1,3)*xg(1,2)-
     &        2.d0*xg(1,1)*xg(2,3)/3.d0)
     &        -2.d0*(xs(1,1)*xg(2,3)+xg(1,1)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(1,1)*xg(2,3)+xg(1,1)*xx(2,3))/3.d0
     &        -d1*xn(1,1)*xn(2,3)-d2*(xn(1,1)*xd(2,3)+
     &        xd(1,1)*xn(2,3))/2.d0+xk*vj2*xg(1,1)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(1,2)*xg(1,3)+xg(1,3)*xg(1,2))/2.d0
         elas(17)=(umb-f0*umbb)*(xg(2,2)*xg(2,3)+xg(2,3)*xg(2,2)-
     &        2.d0*xg(2,2)*xg(2,3)/3.d0)
     &        -2.d0*(xs(2,2)*xg(2,3)+xg(2,2)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(2,2)*xg(2,3)+xg(2,2)*xx(2,3))/3.d0
     &        -d1*xn(2,2)*xn(2,3)-d2*(xn(2,2)*xd(2,3)+
     &        xd(2,2)*xn(2,3))/2.d0+xk*vj2*xg(2,2)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(2,2)*xg(2,3)+xg(2,3)*xg(2,2))/2.d0
         elas(18)=(umb-f0*umbb)*(xg(3,2)*xg(3,3)+xg(3,3)*xg(3,2)-
     &        2.d0*xg(3,3)*xg(2,3)/3.d0)
     &        -2.d0*(xs(3,3)*xg(2,3)+xg(3,3)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(3,3)*xg(2,3)+xg(3,3)*xx(2,3))/3.d0
     &        -d1*xn(3,3)*xn(2,3)-d2*(xn(3,3)*xd(2,3)+
     &        xd(3,3)*xn(2,3))/2.d0+xk*vj2*xg(3,3)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(3,2)*xg(3,3)+xg(3,3)*xg(3,2))/2.d0
         elas(19)=(umb-f0*umbb)*(xg(1,2)*xg(2,3)+xg(1,3)*xg(2,2)-
     &        2.d0*xg(1,2)*xg(2,3)/3.d0)
     &        -2.d0*(xs(1,2)*xg(2,3)+xg(1,2)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(1,2)*xg(2,3)+xg(1,2)*xx(2,3))/3.d0
     &        -d1*xn(1,2)*xn(2,3)-d2*(xn(1,2)*xd(2,3)+
     &        xd(1,2)*xn(2,3))/2.d0+xk*vj2*xg(1,2)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(1,2)*xg(2,3)+xg(1,3)*xg(2,2))/2.d0
         elas(20)=(umb-f0*umbb)*(xg(1,2)*xg(3,3)+xg(1,3)*xg(3,2)-
     &        2.d0*xg(1,3)*xg(2,3)/3.d0)
     &        -2.d0*(xs(1,3)*xg(2,3)+xg(1,3)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(1,3)*xg(2,3)+xg(1,3)*xx(2,3))/3.d0
     &        -d1*xn(1,3)*xn(2,3)-d2*(xn(1,3)*xd(2,3)+
     &        xd(1,3)*xn(2,3))/2.d0+xk*vj2*xg(1,3)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(1,2)*xg(3,3)+xg(1,3)*xg(3,2))/2.d0
         elas(21)=(umb-f0*umbb)*(xg(2,2)*xg(3,3)+xg(2,3)*xg(3,2)-
     &        2.d0*xg(2,3)*xg(2,3)/3.d0)
     &        -2.d0*(xs(2,3)*xg(2,3)+xg(2,3)*xs(2,3))/3.d0
     &        +f0*2.d0*(xx(2,3)*xg(2,3)+xg(2,3)*xx(2,3))/3.d0
     &        -d1*xn(2,3)*xn(2,3)-d2*(xn(2,3)*xd(2,3)+
     &        xd(2,3)*xn(2,3))/2.d0+xk*vj2*xg(2,3)*xg(2,3)
     &        -xk*(vj2-1.d0)*(xg(2,2)*xg(3,3)+xg(2,3)*xg(3,2))/2.d0
!     
      endif
!
!        updating the plastic fields
!
      do i=1,3
         cpl(i)=cpl(i)-1.d0
      enddo
      do i=1,6
         xstate(1+i,iint,iel)=-cpl(i)/2.d0
         xstate(7+i,iint,iel)=stbl(i)
      enddo
      xstate(1,iint,iel)=ep
!
!     maximum difference in the equivalent viscoplastic strain 
!     in this increment based on the viscoplastic strain rate at
!     the start and the end of the increment
!
      if((kode.eq.-52).and.(ivisco.eq.1)) then
c      if(ivisco.eq.1) then
         depvisc=max(depvisc,abs(decra(1)-dtime*xstateini(14,iint,iel)))
         xstate(14,iint,iel)=decra(1)/dtime
      endif
c      depvisc=max(depvisc,ep-epini)
!
      return
      end
