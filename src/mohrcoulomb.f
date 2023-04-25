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
      subroutine mohrcoulomb(elconloc,plconloc,xstate,xstateini,
     &  elas,emec,ithermal,icmd,beta,stre,vj,kode,
     &  ielas,amat,t1l,dtime,time,ttime,iel,iint,nstate_,mi,
     &  eloc,pgauss,nmethod,pnewdt,depvisc)
!
!     calculates stiffness and stresses for the Mohr-Coulomb
!     material law
!
!     icmd=3: calculates stress at mechanical strain
!     else: calculates stress at mechanical strain and the stiffness
!           matrix
!
!     This routine is meant for small strains. Reference: 
!     Elasto-Plastic Strain Hardening Mohr-Coulomb Model
!     Emil Smed Sorensen, Aalborg University
!     M.Sc. 4th Semester, 8 June 2012
!
      implicit none
!     
      character*80 amat
!     
      integer ithermal(*),icmd,i,k,l,m,n,kode,ivisco,ielastic,kel(4,21),
     &     niso,nkin,ielas,iel,iint,nstate_,mi(*),id,leximp,lend,layer,
     &     kspt,kstep,kinc,iloop,nmethod,user_hardening,user_creep,ier,
     &     iregion,matz
!     
      real*8 elconloc(*),elas(21),emec(6),beta(6),stre(6),
     &     vj,plconloc(802),stril(6),xitril(6),xk,xm,sa,
     &     ee,un,um,al,cop,dxitril,xn(3,3),epl(6),c1,c2,c3,c4,c7,
     &     c8,ftrial,xiso(200),yiso(200),xkin(200),ykin(200),
     &     fiso,dfiso,fkin,dfkin,fiso0,fkin0,ep,t1l,dtime,
     &     epini,a1,dsvm,xxa,xxn,dkl(3,3),el(6),tracee,traces,
     &     dcop,time,ttime,eloc(6),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),decra(5),deswa(5),serd,
     &     esw(2),ec(2),p,qtild,predef(1),dpred(1),timeabq(2),pgauss(3),
     &     dtemp,pnewdt,um2,depvisc,c,fv1(3),fv2(3),pf1l1,pf1l6,
     &     pl1ra,pl6f1,pl6ra,traceb,z(3,3),s(3,3),b1(3),b2(3),b6(3),
     &     rl6(3),s1xrl1(3),sb(3),s1xrl6(3),rl1(3),s1(3),s2(3),s6(3),
     &     s1xs2(3),s6xs1(3),dlambda,ddlambda,dlambda2(2),ddlambda2(2),
     &     h,dh,h2(2),dh2(2,2),dk,dm,det,dlambda6(2),ddlambda6(2),
     &     h6(2),dh6(2,2),dlambdar(3),ddlambdar(3),hr(3),dhr(3,3),
     &     s3(3)
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
      dkl=reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),
     &     (/3,3/))
!     
      leximp=1
      lend=2
      user_creep=0
!     
!     localizing the plastic fields
!     
      do i=1,6
        epl(i)=xstateini(1+i,iint,iel)
      enddo
      epini=xstateini(1,iint,iel)
!     
      ee=elconloc(1)
      un=elconloc(2)
      um2=ee/(1.d0+un)
      al=um2*un/(1.d0-2.d0*un)
      um=um2/2.d0
!     
      xk=elconloc(3)
      xm=elconloc(4)
      sa=elconloc(5)
!     
      ep=epini
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
      s(1,1)=stre(1)
      s(2,2)=stre(2)
      s(3,3)=stre(3)
      s(1,2)=stre(4)
      s(1,3)=stre(5)
      s(2,3)=stre(6)
      s(2,1)=s(1,2)
      s(3,1)=s(1,3)
      s(3,2)=s(2,3)
!     
      n=3
      matz=1
      call rs(n,n,s,sb,matz,z,fv1,fv2,ier)
!     
      if(ier.ne.0) then
        write(*,*) '
     &*ERROR calculating the eigenvalues/vectors in umat_abaqusnl'
        call exit(201)
      endif
!     
      ftrial=k*sb(1)-sb(3)-2.d0*c*dsqrt(xk)
      if((ftrial.le.1.d-10).or.(ielas.eq.1).or.(ielastic.eq.1)
     &     .or.(dtime.lt.1.d-30)) then
!     
!     updating the plastic fields
!     
        do i=1,6
          xstate(1+i,iint,iel)=epl(i)
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
!     plastic deformation
!
!     check for hardening
!
      niso=int(plconloc(801))
      if(niso.ne.0) then
        do i=1,niso
          xiso(i)=plconloc(2*i-1)
          yiso(i)=plconloc(2*i)
        enddo
      endif
!
!     D.b for sector I
!     
      b1(1)=xm
      b1(2)=0.d0
      b1(3)=-1.d0
      traceb=xm-1.d0
      do i=1,3
        s1(i)=um2*b1(i)+al*traceb
      enddo
!     
!     vector along yield line between sector I and II
!     
      rl1(1)=1.d0
      rl1(2)=1.d0
      rl1(3)=xk
!     
!     vector along yield line between sector I and VI
!     
      rl6(1)=1.d0
      rl6(2)=xk
      rl6(3)=xk
!     
      sa=2.d0*c*dsqrt(xk)/(xk-1.d0)
!     
!     s1 x rl1
!     
      s1xrl1(1)=s1(2)*rl1(3)-s1(3)*rl1(2)
      s1xrl1(2)=s1(3)*rl1(1)-s1(1)*rl1(3)
      s1xrl1(3)=s1(1)*rl1(2)-s1(2)*rl1(1)
!     
!     s1 x rl6
!     
      s1xrl6(1)=s1(2)*rl6(3)-s1(3)*rl6(2)
      s1xrl6(2)=s1(3)*rl6(1)-s1(1)*rl6(3)
      s1xrl6(3)=s1(1)*rl6(2)-s1(2)*rl6(1)
!     
!     boundary plane between sector I and II
!     (s1 x rl1).(sb-sa)
!     
      pf1l1=s1xrl1(1)*(sb(1)-sa)+
     &     s1xrl1(2)*(sb(2)-sa)+
     &     s1xrl1(3)*(sb(3)-sa)
!     
!     boundary plane between sector I and VI
!     (s1 x rl6).(sb-sa)
!     
      pf1l6=s1xrl6(1)*(sb(1)-sa)+
     &     s1xrl6(2)*(sb(2)-sa)+
     &     s1xrl6(3)*(sb(3)-sa)
!     
      if((pf1l1.le.0.d0).and.(pl6f1.ge.0.d0)) then
        iregion=1
      else
!     
!     D.b for sector II
!     
        b2(1)=0.d0
        b2(2)=xm
        b2(3)=-1.d0
        do i=1,3
          s2(i)=um2*b2(i)+al*traceb
        enddo
!     
!     s1 x s2
!     
        s1xs2(1)=s1(2)*s2(3)-s1(3)*s2(2)
        s1xs2(2)=s1(3)*s2(1)-s1(1)*s2(3)
        s1xs2(3)=s1(1)*s2(2)-s1(2)*s2(1)
!     
!     boundary plane between sector I and "ra"
!     (s1 x s2).(sb-sa)
!     
        pl1ra=s1xs2(1)*(sb(1)-sa)+
     &       s1xs2(2)*(sb(2)-sa)+
     &       s1xs2(3)*(sb(3)-sa)
        if((pf1l1.ge.0.d0).and.(pl1ra.le.0.d0)) then
          iregion=2
!     
!     D.b for sector VI
!     
          b6(1)=xm
          b6(2)=-1.d0
          b6(3)=0.d0
          do i=1,3
            s6(i)=um2*b6(i)+al*traceb
          enddo
        else
!     
!     s6 x s1
!     
          s6xs1(1)=s6(2)*s1(3)-s6(3)*s1(2)
          s6xs1(2)=s6(3)*s1(1)-s6(1)*s1(3)
          s6xs1(3)=s6(1)*s1(2)-s6(2)*s1(1)
!     
!     boundary plane between sector I and "ra"
!     (s6 x s1).(sb-sa)
!     
          pl6ra=s6xs1(1)*(sb(1)-sa)+
     &         s6xs1(2)*(sb(2)-sa)+
     &         s6xs1(3)*(sb(3)-sa)
          if((pl6f1.le.0.d0).and.(pl6ra.le.0.d0)) then
            iregion=3
          else
            iregion=4
          endif
        endif
      endif
!
!     calculate the change in lambda
!
      dk=2.d0*dsqrt(xk)
      dm=dsqrt(2.d0*(xm*xm+1.d0)/3.d0)
      iloop=0
      dlambda=0.d0
!      
      if(iregion.eq.1) then
        do
          iloop=iloop+1
          ep=epini+dm*dlambda
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
          else
!            
!           needed?            
!            
            fiso=0.d0
            dfiso=0.d0
          endif
!
          h=xk*(sb(1)-dlambda*s1(1))-(sb(3)-dlambda*s1(3))-dk*fiso
          dh=-xk*s1(1)+s1(3)-dk*dm*dfiso
          ddlambda=-h/dh
!
          if((ddlambda.lt.1.d-10).or.(ddlambda.lt.1.d-4*dlambda)) exit
          dlambda=dlambda+ddlambda
          if((iloop.gt.15).or.(dlambda.le.0.d0)) then
            pnewdt=0.25d0
            return
          endif
        enddo
      elseif(iregion.eq.2) then
        iloop=0
        dlambda2(1)=0.d0
        dlambda2(2)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambda2(1)+dlambda2(2))
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
          else
!            
!           needed?            
!            
            fiso=0.d0
            dfiso=0.d0
          endif
!
!         setting up the 2x2 equation system: right hand side
!
          h2(1)=xk*(sb(1)-dlambda2(1)*s1(1)-dlambda2(2)*s2(1))
     &            -(sb(3)-dlambda2(1)*s1(3)-dlambda2(2)*s2(3))
     &            -dk*fiso
          h2(2)=xk*(sb(2)-dlambda2(1)*s1(2)-dlambda2(2)*s2(2))
     &            -(sb(3)-dlambda2(1)*s1(3)-dlambda2(2)*s2(3))
     &            -dk*fiso
!
!         setting up the 2x2 equation system: left hand side
!
          dh2(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dh2(1,2)=-xk*s2(1)+s2(3)-dk*dm*dfiso
          dh2(2,1)=-xk*s1(2)+s1(3)-dk*dm*dfiso
          dh2(2,2)=-xk*s2(2)+s2(3)-dk*dm*dfiso
!
          det=dh2(1,1)*dh2(2,2)-dh2(2,1)*dh2(1,2)
!
!         solving the system
!
          ddlambda2(1)=(h2(1)*dh2(2,2)-h2(2)*dh2(1,2))/det
          ddlambda2(2)=(dh2(1,1)*h2(2)-dh2(2,1)*h2(1))/det
!
          if(((ddlambda2(1).lt.1.d-10).or.
     &        (ddlambda2(1).lt.1.d-4*dlambda2(1))).and.
     &       ((ddlambda2(2).lt.1.d-10).or.
     &        (ddlambda2(2).lt.1.d-4*dlambda2(2)))) exit
!
          dlambda2(1)=dlambda2(1)+ddlambda2(1)
          dlambda2(2)=dlambda2(2)+ddlambda2(2)
!
          if((iloop.gt.15).or.(dlambda2(1).le.0.d0).or.
     &       (dlambda2(2).le.0.d0)) then
            pnewdt=0.25d0
            return
          endif
        enddo
      elseif(iregion.eq.3) then
        iloop=0
        dlambda6(1)=0.d0
        dlambda6(2)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambda6(1)+dlambda6(2))
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
          else
!            
!           needed?            
!            
            fiso=0.d0
            dfiso=0.d0
          endif
!
!         setting up the 2x2 equation system: right hand side
!
          h6(1)=xk*(sb(1)-dlambda6(1)*s1(1)-dlambda6(2)*s6(1))
     &            -(sb(3)-dlambda6(1)*s1(3)-dlambda6(2)*s6(3))
     &            -dk*fiso
          h6(2)=xk*(sb(2)-dlambda6(1)*s1(2)-dlambda6(2)*s6(2))
     &            -(sb(3)-dlambda6(1)*s1(3)-dlambda6(2)*s6(3))
     &            -dk*fiso
!
!         setting up the 2x2 equation system: left hand side
!
          dh6(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dh6(1,2)=-xk*s6(1)+s6(3)-dk*dm*dfiso
          dh6(2,1)=-xk*s1(2)+s1(3)-dk*dm*dfiso
          dh6(2,2)=-xk*s6(2)+s6(3)-dk*dm*dfiso
!
          det=dh6(1,1)*dh6(2,2)-dh6(2,1)*dh6(1,2)
!
!         solving the system
!
          ddlambda6(1)=(h6(1)*dh6(2,2)-h6(2)*dh6(1,2))/det
          ddlambda6(2)=(dh6(1,1)*h6(2)-dh6(2,1)*h6(1))/det
!
          if(((ddlambda6(1).lt.1.d-10).or.
     &        (ddlambda6(1).lt.1.d-4*dlambda6(1))).and.
     &       ((ddlambda6(2).lt.1.d-10).or.
     &        (ddlambda6(2).lt.1.d-4*dlambda6(2)))) exit
!
          dlambda6(1)=dlambda6(1)+ddlambda6(1)
          dlambda6(2)=dlambda6(2)+ddlambda6(2)
!
          if((iloop.gt.15).or.(dlambda6(1).le.0.d0).or.
     &       (dlambda6(2).le.0.d0)) then
            pnewdt=0.25d0
            return
          endif
        enddo
      else
!
!       region IV
!
        iloop=0
        dlambdar(1)=0.d0
        dlambdar(2)=0.d0
        dlambdar(3)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambdar(1)+dlambdar(2)+dlambdar(3))
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
          else
!            
!           needed?            
!            
            fiso=0.d0
            dfiso=0.d0
          endif
!
!         setting up the 3x3 equation system: right hand side
!
          hr(1)=xk*(sb(1)-dlambdar(1)*s1(1)-dlambdar(2)*s2(1)
     &             -dlambdar(3)*s3(1))
     &            -(sb(3)-dlambdar(1)*s1(3)-dlambdar(2)*s2(3)
     &             -dlambdar(3)*s3(3))
     &            -dk*fiso
          hr(1)=xk*(sb(2)-dlambdar(1)*s1(2)-dlambdar(2)*s2(2)
     &             -dlambdar(3)*s3(2))
     &            -(sb(3)-dlambdar(1)*s1(3)-dlambdar(2)*s2(3)
     &             -dlambdar(3)*s3(3))
     &            -dk*fiso
          hr(1)=xk*(sb(1)-dlambdar(1)*s1(1)-dlambdar(2)*s2(1)
     &             -dlambdar(3)*s3(1))
     &            -(sb(2)-dlambdar(1)*s1(2)-dlambdar(2)*s2(2)
     &             -dlambdar(3)*s3(2))
     &            -dk*fiso
!
!         setting up the 2x2 equation system: left hand side
!
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dhr(1,2)=-xk*s2(1)+s2(3)-dk*dm*dfiso
          dhr(1,3)=-xk*s3(1)+s3(3)-dk*dm*dfiso
          dhr(2,1)=-xk*s1(2)+s1(3)-dk*dm*dfiso
!          ******
!          *****
          dhr(2,2)=-xk*s1(1)+s2(2)-dk*dm*dfiso
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
!
          det=dh2(1,1)*dh2(2,2)-dh2(2,1)*dh2(1,2)
!
!         solving the system
!
          ddlambdar(1)=(h2(1)*dh2(2,2)-h2(2)*dh2(1,2))/det
          ddlambdar(2)=(dh2(1,1)*h2(2)-dh2(2,1)*h2(1))/det
!
          if(((ddlambdar(1).lt.1.d-10).or.
     &        (ddlambdar(1).lt.1.d-4*dlambdar(1))).and.
     &       ((ddlambdar(2).lt.1.d-10).or.
     &        (ddlambdar(2).lt.1.d-4*dlambdar(2)))) exit
!
          dlambdar(1)=dlambdar(1)+ddlambdar(1)
          dlambdar(2)=dlambdar(2)+ddlambdar(2)
!
          if((iloop.gt.15).or.(dlambdar(1).le.0.d0).or.
     &       (dlambdar(2).le.0.d0)) then
            pnewdt=0.25d0
            return
          endif
        enddo
      endif
          






      
!     
      return
      end
