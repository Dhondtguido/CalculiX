!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
     &  elas,emec,icmd,beta,stre,
     &  ielas,dtime,time,ttime,iel,iint,nstate_,mi,pnewdt)
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
      integer icmd,i,j,k,n,niso,ielas,iel,iint,nstate_,mi(*),id,
     &     kstep,kinc,iloop,ier,iregion,matz,j1,j2,j3,j4,jj,kal(2,6),
     &     kel(4,21),j5,j6,j7,j8
!     
      real*8 elconloc(*),elas(21),emec(6),beta(6),stre(6),sc(6),
     &     plconloc(802),xk,xm,sa,stiff(6,6),ttime,ee,un,um,al,epl(6),
     &     ftrial,xiso(200),yiso(200),da6(3),fiso,dfiso,ep,dtime,denom,
     &     epini,el(6),tracee,a2(3),a6(3),time,xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),tracea,da1(3),delas(21),
     &     pnewdt,um2,fv1(3),fv2(3),ps1r1,ps1r6,
     &     ps1s2,ps6s1,traceb,z(3,3),s(3,3),b1(3),b2(3),b6(3),
     &     r6(3),s1xr1(3),sb(3),s1xr6(3),r1(3),s1(3),s2(3),s6(3),
     &     s1xs2(3),s6xs1(3),dlambda,ddlambda,dlambda2(2),ddlambda2(2),
     &     h,dh,h2(2),dh2(2,2),dk,dm,det,dlambda6(2),ddlambda6(2),
     &     h6(2),dh6(2,2),dlambdar(3),ddlambdar(3),hr(3),dhr(3,3),
     &     a(3,3),b(3,3),a1(3),da2(3),t(6,6),dum(6,6),dum1,ya(3,3,3,3)
!     
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
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
!     
      ep=epini
!     
!     hardening
!     
      niso=int(plconloc(801))
      do i=1,niso
        xiso(i)=plconloc(2*i-1)
        yiso(i)=plconloc(2*i)
      enddo
!     
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
!     switching eigenvalues and eigenvectors such that sb(1) is the
!     biggest and sb(3) is the smallest eigenvector
!
      dum1=sb(1)
      sb(1)=sb(3)
      sb(3)=dum1
!
      do i=1,3
        dum1=z(i,1)
        z(i,1)=z(i,3)
        z(i,3)=dum1
      enddo
!
      ftrial=xk*sb(1)-sb(3)-2.d0*fiso*dsqrt(xk)
      if((ftrial.le.1.d-10).or.(ielas.eq.1).or.(dtime.lt.1.d-30)) then
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
      r1(1)=1.d0
      r1(2)=1.d0
      r1(3)=xk
!     
!     vector along yield line between sector I and VI
!     
      r6(1)=1.d0
      r6(2)=xk
      r6(3)=xk
!     
      sa=2.d0*fiso*dsqrt(xk)/(xk-1.d0)
!     
!     s1 x r1
!     
      s1xr1(1)=s1(2)*r1(3)-s1(3)*r1(2)
      s1xr1(2)=s1(3)*r1(1)-s1(1)*r1(3)
      s1xr1(3)=s1(1)*r1(2)-s1(2)*r1(1)
!     
!     s1 x r6
!     
      s1xr6(1)=s1(2)*r6(3)-s1(3)*r6(2)
      s1xr6(2)=s1(3)*r6(1)-s1(1)*r6(3)
      s1xr6(3)=s1(1)*r6(2)-s1(2)*r6(1)
!     
!     boundary plane between sector I and II
!     (s1 x r1).(sb-sa)
!     
      ps1r1=s1xr1(1)*(sb(1)-sa)+
     &     s1xr1(2)*(sb(2)-sa)+
     &     s1xr1(3)*(sb(3)-sa)
!     
!     boundary plane between sector I and VI
!     (s1 x r6).(sb-sa)
!     
      ps1r6=s1xr6(1)*(sb(1)-sa)+
     &     s1xr6(2)*(sb(2)-sa)+
     &     s1xr6(3)*(sb(3)-sa)
!     
      if((ps1r1.ge.0.d0).and.(ps1r6.le.0.d0)) then
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
        ps1s2=s1xs2(1)*(sb(1)-sa)+
     &       s1xs2(2)*(sb(2)-sa)+
     &       s1xs2(3)*(sb(3)-sa)
        if((ps1r1.le.0.d0).and.(ps1s2.le.0.d0)) then
          iregion=2
        else
!     
!     D.b for sector VI
!     
          b6(1)=xm
          b6(2)=-1.d0
          b6(3)=0.d0
          do i=1,3
            s6(i)=um2*b6(i)+al*traceb
          enddo
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
          ps6s1=s6xs1(1)*(sb(1)-sa)+
     &         s6xs1(2)*(sb(2)-sa)+
     &         s6xs1(3)*(sb(3)-sa)
          if((ps1r6.ge.0.d0).and.(ps6s1.le.0.d0)) then
            iregion=3
          else
            iregion=4
          endif
        endif
      endif
c      write(*,*) 'mohrcoulomb ',iel,iint,iregion
!     
!     calculate the change in lambda
!     
      dk=2.d0*dsqrt(xk)
      dm=dsqrt(2.d0*(xm*xm+1.d0)/3.d0)
!     
      if(iregion.eq.1) then
        iloop=0
        dlambda=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*dlambda
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
!     
          h=xk*(sb(1)-dlambda*s1(1))-(sb(3)-dlambda*s1(3))-dk*fiso
          dh=-xk*s1(1)+s1(3)-dk*dm*dfiso
          ddlambda=-h/dh
!     
          if((dabs(ddlambda).lt.1.d-10).or.
     &         (dabs(ddlambda).lt.1.d-4*dlambda)) exit
          dlambda=dlambda+ddlambda
          if((iloop.gt.15).or.(dlambda.le.-1.d-10)) then
            pnewdt=0.25d0
            return
          endif
        enddo
        dlambda=max(dlambda,0.d0)
!     
!     calculate the stress at C
!     
        do i=1,3
          sc(i)=sb(i)-dlambda*s1(i)
        enddo
        do i=4,6
          sc(i)=0.d0
        enddo
!     
        if(icmd.ne.3) then
!     
!     calculate the tangent stiffness matrix
!     
          a1(1)=xk
          a1(2)=0.d0
          a1(3)=-1.d0
          tracea=xk-1.d0
          do i=1,3
            da1(i)=um2*a1(i)+al*tracea
          enddo
          denom=a1(1)*s1(1)+a1(2)*s1(2)+a1(3)*s1(3)+dk*dm*dfiso
          do i=1,3
            do j=1,3
              stiff(i,j)=al-s1(i)*da1(j)/denom
            enddo
            stiff(i,i)=stiff(i,i)+um2
          enddo
          do i=1,3
            do j=4,6
              stiff(i,j)=0.d0
              stiff(j,i)=0.d0
            enddo
          enddo
          do i=4,6
            do j=4,6
              stiff(i,j)=0.d0
            enddo
          enddo
        endif
        stiff(4,4)=(sc(1)-sc(2))/(sb(1)-sb(2))*um
        stiff(5,5)=(sc(1)-sc(3))/(sb(1)-sb(3))*um
        stiff(6,6)=(sc(2)-sc(3))/(sb(2)-sb(3))*um
!     
      elseif(iregion.eq.2) then
        iloop=0
        dlambda2(1)=0.d0
        dlambda2(2)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambda2(1)+dlambda2(2))
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
!     
!     setting up the 2x2 equation system: right hand side
!     
          h2(1)=xk*(sb(1)-dlambda2(1)*s1(1)-dlambda2(2)*s2(1))
     &         -(sb(3)-dlambda2(1)*s1(3)-dlambda2(2)*s2(3))
     &         -dk*fiso
          h2(2)=xk*(sb(2)-dlambda2(1)*s1(2)-dlambda2(2)*s2(2))
     &         -(sb(3)-dlambda2(1)*s1(3)-dlambda2(2)*s2(3))
     &         -dk*fiso
!     
!     setting up the 2x2 equation system: left hand side
!     
          dh2(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dh2(1,2)=-xk*s2(1)+s2(3)-dk*dm*dfiso
          dh2(2,1)=-xk*s1(2)+s1(3)-dk*dm*dfiso
          dh2(2,2)=-xk*s2(2)+s2(3)-dk*dm*dfiso
!     
          det=dh2(1,1)*dh2(2,2)-dh2(2,1)*dh2(1,2)
!     
!     solving the system
!     
          ddlambda2(1)=-(h2(1)*dh2(2,2)-h2(2)*dh2(1,2))/det
          ddlambda2(2)=-(dh2(1,1)*h2(2)-dh2(2,1)*h2(1))/det
!     
          if(((dabs(ddlambda2(1)).lt.1.d-10).or.
     &         (dabs(ddlambda2(1)).lt.1.d-4*dlambda2(1))).and.
     &         ((dabs(ddlambda2(2)).lt.1.d-10).or.
     &         (dabs(ddlambda2(2)).lt.1.d-4*dlambda2(2)))) exit
!     
          dlambda2(1)=dlambda2(1)+ddlambda2(1)
          dlambda2(2)=dlambda2(2)+ddlambda2(2)
!     
          if((iloop.gt.15).or.(dlambda2(1).le.-1.d-10).or.
     &         (dlambda2(2).le.-1.d-10)) then
            pnewdt=0.25d0
            return
          endif
        enddo
        dlambda2(1)=max(dlambda2(1),0.d0)
        dlambda2(2)=max(dlambda2(2),0.d0)
!     
!     calculate the stress at C
!     
        do i=1,3
          sc(i)=sb(i)-dlambda2(1)*s1(i)-dlambda2(2)*s2(i)
        enddo
        do i=4,6
          sc(i)=0.d0
        enddo
!     
        if(icmd.ne.3) then
!     
!     calculate the tangent stiffness matrix
!     
          a1(1)=xk
          a1(2)=0.d0
          a1(3)=-1.d0
          a2(1)=0.d0
          a2(2)=xk
          a2(3)=-1.d0
          tracea=xk-1.d0
          do i=1,3
            da1(i)=um2*a1(i)+al*tracea
            da2(i)=um2*a2(i)+al*tracea
          enddo
!     
!     setting up lhs matrix a(*,*)
!     
          a(1,1)=a1(1)*s1(1)+a1(2)*s1(2)+a1(3)*s1(3)+dk*dm*dfiso
          a(1,2)=a1(1)*s2(1)+a1(2)*s2(2)+a1(3)*s2(3)+dk*dm*dfiso
          a(2,1)=a2(1)*s1(1)+a2(2)*s1(2)+a2(3)*s1(3)+dk*dm*dfiso
          a(2,2)=a2(1)*s2(1)+a2(2)*s2(2)+a2(3)*s2(3)+dk*dm*dfiso
!     
!     inverting the matrix -> b(*,*)
!     
          det=a(1,1)*a(2,2)-a(2,1)*a(1,2)
          b(1,1)=a(2,2)/det
          b(1,2)=-a(1,2)/det
          b(2,1)=-a(2,1)/det
          b(2,2)=a(1,1)/det
          do i=1,3
            do j=1,3
              stiff(i,j)=al-b(1,1)*s1(i)*da1(j)
     &             -b(1,2)*s1(i)*da2(j)
     &             -b(2,1)*s2(i)*da1(j)
     &             -b(2,2)*s2(i)*da2(j)
            enddo
            stiff(i,i)=stiff(i,i)+um2
          enddo
          do i=1,3
            do j=4,6
              stiff(i,j)=0.d0
              stiff(j,i)=0.d0
            enddo
          enddo
          do i=4,6
            do j=4,6
              stiff(i,j)=0.d0
            enddo
          enddo
          stiff(4,4)=(sc(1)-sc(2))/(sb(1)-sb(2))*um
          stiff(5,5)=(sc(1)-sc(3))/(sb(1)-sb(3))*um
          stiff(6,6)=(sc(2)-sc(3))/(sb(2)-sb(3))*um
        endif
      elseif(iregion.eq.3) then
        iloop=0
        dlambda6(1)=0.d0
        dlambda6(2)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambda6(1)+dlambda6(2))
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
!     
!     setting up the 2x2 equation system: right hand side
!     
          h6(1)=xk*(sb(1)-dlambda6(1)*s1(1)-dlambda6(2)*s6(1))
     &         -(sb(3)-dlambda6(1)*s1(3)-dlambda6(2)*s6(3))
     &         -dk*fiso
          h6(2)=xk*(sb(1)-dlambda6(1)*s1(1)-dlambda6(2)*s6(1))
     &         -(sb(2)-dlambda6(1)*s1(2)-dlambda6(2)*s6(2))
     &         -dk*fiso
!     
!     setting up the 2x2 equation system: left hand side
!     
          dh6(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dh6(1,2)=-xk*s6(1)+s6(3)-dk*dm*dfiso
          dh6(2,1)=-xk*s1(1)+s1(2)-dk*dm*dfiso
          dh6(2,2)=-xk*s6(1)+s6(2)-dk*dm*dfiso
!     
          det=dh6(1,1)*dh6(2,2)-dh6(2,1)*dh6(1,2)
!     
!     solving the system
!     
          ddlambda6(1)=-(h6(1)*dh6(2,2)-h6(2)*dh6(1,2))/det
          ddlambda6(2)=-(dh6(1,1)*h6(2)-dh6(2,1)*h6(1))/det
!     
          if(((dabs(ddlambda6(1)).lt.1.d-10).or.
     &         (dabs(ddlambda6(1)).lt.1.d-4*dlambda6(1))).and.
     &         ((dabs(ddlambda6(2)).lt.1.d-10).or.
     &         (dabs(ddlambda6(2)).lt.1.d-4*dlambda6(2)))) exit
!     
          dlambda6(1)=dlambda6(1)+ddlambda6(1)
          dlambda6(2)=dlambda6(2)+ddlambda6(2)
!     
          if((iloop.gt.15).or.(dlambda6(1).le.-1.d-10).or.
     &         (dlambda6(2).le.-1.e-10)) then
            pnewdt=0.25d0
            return
          endif
        enddo
        dlambda6(1)=max(dlambda6(1),0.d0)
        dlambda6(2)=max(dlambda6(2),0.d0)
!     
!     calculate the stress at C
!     
        do i=1,3
          sc(i)=sb(i)-dlambda6(1)*s1(i)-dlambda6(2)*s6(i)
        enddo
        do i=4,6
          sc(i)=0.d0
        enddo
!     
        if(icmd.ne.3) then
!     
!     calculate the tangent stiffness matrix
!     
          a1(1)=xk
          a1(2)=0.d0
          a1(3)=-1.d0
          a6(1)=xk
          a6(2)=-1.d0
          a6(3)=0.d0
          tracea=xk-1.d0
          do i=1,3
            da1(i)=um2*a1(i)+al*tracea
            da6(i)=um2*a6(i)+al*tracea
          enddo
!     
!     setting up lhs matrix a(*,*)
!     
          a(1,1)=a1(1)*s1(1)+a1(2)*s1(2)+a1(3)*s1(3)+dk*dm*dfiso
          a(1,2)=a1(1)*s6(1)+a1(2)*s6(2)+a1(3)*s6(3)+dk*dm*dfiso
          a(2,1)=a6(1)*s1(1)+a6(2)*s1(2)+a6(3)*s1(3)+dk*dm*dfiso
          a(2,2)=a6(1)*s6(1)+a6(2)*s6(2)+a6(3)*s6(3)+dk*dm*dfiso
!     
!     inverting the matrix -> b(*,*)
!     
          det=a(1,1)*a(2,2)-a(2,1)*a(1,2)
          b(1,1)=a(2,2)/det
          b(1,2)=-a(1,2)/det
          b(2,1)=-a(2,1)/det
          b(2,2)=a(1,1)/det
          do i=1,3
            do j=1,3
              stiff(i,j)=al-b(1,1)*s1(i)*da1(j)
     &             -b(1,2)*s1(i)*da6(j)
     &             -b(2,1)*s6(i)*da1(j)
     &             -b(2,2)*s6(i)*da6(j)
            enddo
            stiff(i,i)=stiff(i,i)+um2
          enddo
          do i=1,3
            do j=4,6
              stiff(i,j)=0.d0
              stiff(j,i)=0.d0
            enddo
          enddo
          do i=4,6
            do j=4,6
              stiff(i,j)=0.d0
            enddo
          enddo
          stiff(4,4)=(sc(1)-sc(2))/(sb(1)-sb(2))*um
          stiff(5,5)=(sc(1)-sc(3))/(sb(1)-sb(3))*um
          stiff(6,6)=(sc(2)-sc(3))/(sb(2)-sb(3))*um
        endif
      else
!     
!     region IV
!     
        iloop=0
        dlambdar(1)=0.d0
        dlambdar(2)=0.d0
        dlambdar(3)=0.d0
        do
          iloop=iloop+1
          ep=epini+dm*(dlambdar(1)+dlambdar(2)+dlambdar(3))
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
!     
!     setting up the 3x3 equation system: right hand side
!     
          hr(1)=xk*(sb(1)-dlambdar(1)*s1(1)-dlambdar(2)*s2(1)
     &         -dlambdar(3)*s6(1))
     &         -(sb(3)-dlambdar(1)*s1(3)-dlambdar(2)*s2(3)
     &         -dlambdar(3)*s6(3))
     &         -dk*fiso
          hr(2)=xk*(sb(2)-dlambdar(1)*s1(2)-dlambdar(2)*s2(2)
     &         -dlambdar(3)*s6(2))
     &         -(sb(3)-dlambdar(1)*s1(3)-dlambdar(2)*s2(3)
     &         -dlambdar(3)*s6(3))
     &         -dk*fiso
          hr(3)=xk*(sb(1)-dlambdar(1)*s1(1)-dlambdar(2)*s2(1)
     &         -dlambdar(3)*s6(1))
     &         -(sb(2)-dlambdar(1)*s1(2)-dlambdar(2)*s2(2)
     &         -dlambdar(3)*s6(2))
     &         -dk*fiso
!     
!     setting up the 3x3 equation system: left hand side
!     
          dhr(1,1)=-xk*s1(1)+s1(3)-dk*dm*dfiso
          dhr(1,2)=-xk*s2(1)+s2(3)-dk*dm*dfiso
          dhr(1,3)=-xk*s6(1)+s6(3)-dk*dm*dfiso
          dhr(2,1)=-xk*s1(2)+s1(3)-dk*dm*dfiso
          dhr(2,2)=-xk*s2(2)+s2(3)-dk*dm*dfiso
          dhr(2,3)=-xk*s6(2)+s6(3)-dk*dm*dfiso
          dhr(3,1)=-xk*s1(1)+s1(2)-dk*dm*dfiso
          dhr(3,2)=-xk*s2(1)+s2(2)-dk*dm*dfiso
          dhr(3,3)=-xk*s6(1)+s6(2)-dk*dm*dfiso
!     
          det=dhr(1,1)*(dhr(2,2)*dhr(3,3)-dhr(2,3)*dhr(3,2))
     &         -dhr(1,2)*(dhr(2,1)*dhr(3,3)-dhr(2,3)*dhr(3,1))
     &         +dhr(1,3)*(dhr(2,1)*dhr(3,2)-dhr(2,2)*dhr(3,1))
!     
!     solving the system
!     
          ddlambdar(1)=-(hr(1)*(dhr(2,2)*dhr(3,3)-dhr(2,3)*dhr(3,2))
     &         -dhr(1,2)*(hr(2)*dhr(3,3)-dhr(2,3)*hr(3))
     &         +dhr(1,3)*(hr(2)*dhr(3,2)-dhr(2,2)*hr(3)))/det
          ddlambdar(2)=-(dhr(1,1)*(hr(2)*dhr(3,3)-dhr(2,3)*hr(3))
     &         -hr(1)*(dhr(2,1)*dhr(3,3)-dhr(2,3)*dhr(3,1))
     &         +dhr(1,3)*(dhr(2,1)*hr(3)-hr(2)*dhr(3,1)))/det
          ddlambdar(3)=-(dhr(1,1)*(dhr(2,2)*hr(3)-hr(2)*dhr(3,2))
     &         -dhr(1,2)*(dhr(2,1)*hr(3)-hr(2)*dhr(3,1))
     &         +hr(1)*(dhr(2,1)*dhr(3,2)-dhr(2,2)*dhr(3,1)))/det
!     
          if(((dabs(ddlambdar(1)).lt.1.d-10).or.
     &         (dabs(ddlambdar(1)).lt.1.d-4*dlambdar(1))).and.
     &         ((dabs(ddlambdar(2)).lt.1.d-10).or.
     &         (dabs(ddlambdar(2)).lt.1.d-4*dlambdar(2))).and.
     &         ((dabs(ddlambdar(3)).lt.1.d-10).or.
     &         (dabs(ddlambdar(3)).lt.1.d-4*dlambdar(3)))) exit
!     
          dlambdar(1)=dlambdar(1)+ddlambdar(1)
          dlambdar(2)=dlambdar(2)+ddlambdar(2)
          dlambdar(3)=dlambdar(3)+ddlambdar(3)
!     
          if((iloop.gt.15).or.(dlambdar(1).le.-1.d-10).or.
     &       (dlambdar(2).le.-1.d-10).or.(dlambdar(3).le.-1.d-10)) then
            pnewdt=0.25d0
            return
          endif
        enddo
        dlambdar(1)=max(dlambdar(1),0.d0)
        dlambdar(2)=max(dlambdar(2),0.d0)
        dlambdar(3)=max(dlambdar(3),0.d0)
!     
!     calculate the stress at C
!     
        do i=1,3
          sc(i)=sb(i)
     &         -dlambdar(1)*s1(i)-dlambdar(2)*s2(i)-dlambdar(3)*s6(i)
        enddo
        do i=4,6
          sc(i)=0.d0
        enddo
!     
        if(icmd.ne.3) then
!     
!     calculate the tangent stiffness matrix
!     
          a1(1)=xk
          a1(2)=0.d0
          a1(3)=-1.d0
          a2(1)=0.d0
          a2(2)=xk
          a2(3)=-1.d0
          a6(1)=xk
          a6(2)=-1.d0
          a6(3)=0.d0
          tracea=xk-1.d0
          do i=1,3
            da1(i)=um2*a1(i)+al*tracea
            da2(i)=um2*a2(i)+al*tracea
            da6(i)=um2*a6(i)+al*tracea
          enddo
!     
!     setting up lhs matrix a(*,*)
!     
          a(1,1)=a1(1)*s1(1)+a1(2)*s1(2)+a1(3)*s1(3)+dk*dm*dfiso
          a(1,2)=a1(1)*s2(1)+a1(2)*s2(2)+a1(3)*s2(3)+dk*dm*dfiso
          a(1,3)=a1(1)*s6(1)+a1(2)*s6(2)+a1(3)*s6(3)+dk*dm*dfiso
          a(2,1)=a2(1)*s1(1)+a2(2)*s1(2)+a2(3)*s1(3)+dk*dm*dfiso
          a(2,2)=a2(1)*s2(1)+a2(2)*s2(2)+a2(3)*s2(3)+dk*dm*dfiso
          a(2,3)=a2(1)*s6(1)+a2(2)*s6(2)+a2(3)*s6(3)+dk*dm*dfiso
          a(3,1)=a6(1)*s1(1)+a6(2)*s1(2)+a6(3)*s1(3)+dk*dm*dfiso
          a(3,2)=a6(1)*s2(1)+a6(2)*s2(2)+a6(3)*s2(3)+dk*dm*dfiso
          a(3,3)=a6(1)*s6(1)+a6(2)*s6(2)+a6(3)*s6(3)+dk*dm*dfiso
!     
!     inverting the matrix -> b(*,*)
!     
          det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &         -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &         +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
          b(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/det
          b(1,2)=-(a(1,2)*a(3,3)-a(3,2)*a(1,3))/det
          b(1,3)=-(a(1,2)*a(2,3)-a(2,2)*a(1,3))/det
          b(2,1)=-(a(2,1)*a(3,3)-a(3,1)*a(2,3))/det
          b(2,2)=(a(1,1)*a(3,3)-a(3,1)*a(1,3))/det
          b(2,3)=-(a(1,1)*a(2,3)-a(2,1)*a(1,3))/det
          b(3,1)=-(a(2,1)*a(3,2)-a(3,1)*a(2,2))/det
          b(3,2)=-(a(1,1)*a(3,2)-a(3,1)*a(1,2))/det
          b(3,3)=(a(1,1)*a(2,2)-a(2,1)*a(1,2))/det
!     
          do i=1,3
            do j=1,3
              stiff(i,j)=al-b(1,1)*s1(i)*da1(j)
     &             -b(1,2)*s1(i)*da2(j)
     &             -b(1,3)*s1(i)*da6(j)
     &             -b(2,1)*s2(i)*da1(j)
     &             -b(2,2)*s2(i)*da2(j)
     &             -b(2,3)*s2(i)*da6(j)
     &             -b(3,1)*s6(i)*da1(j)
     &             -b(3,2)*s6(i)*da2(j)
     &             -b(3,3)*s6(i)*da6(j)
            enddo
            stiff(i,i)=stiff(i,i)+um2
          enddo
          do i=1,3
            do j=4,6
              stiff(i,j)=0.d0
              stiff(j,i)=0.d0
            enddo
          enddo
          do i=4,6
            do j=4,6
              stiff(i,j)=0.d0
            enddo
          enddo
          stiff(4,4)=(sc(1)-sc(2))/(sb(1)-sb(2))*um
          stiff(5,5)=(sc(1)-sc(3))/(sb(1)-sb(3))*um
          stiff(6,6)=(sc(2)-sc(3))/(sb(2)-sb(3))*um
        endif
      endif
!     
!     BACKTRANSFORMATION FROM PRINCIPAL AXES INTO GLOBAL AXES
!     
      t(1,1)=z(1,1)*z(1,1)
      t(1,2)=z(2,1)*z(2,1)
      t(1,3)=z(3,1)*z(3,1)
      t(1,4)=z(1,1)*z(2,1)
      t(1,5)=z(3,1)*z(1,1)
      t(1,6)=z(2,1)*z(3,1)
      t(2,1)=z(1,2)*z(1,2)
      t(2,2)=z(2,2)*z(2,2)
      t(2,3)=z(3,2)*z(3,2)
      t(2,4)=z(1,2)*z(2,2)
      t(2,5)=z(3,2)*z(1,2)
      t(2,6)=z(2,2)*z(3,2)
      t(3,1)=z(1,3)*z(1,3)
      t(3,2)=z(2,3)*z(2,3)
      t(3,3)=z(3,3)*z(3,3)
      t(3,4)=z(1,3)*z(2,3)
      t(3,5)=z(3,3)*z(1,3)
      t(3,6)=z(2,3)*z(3,3)
      t(4,1)=2.d0*z(1,1)*z(1,2)
      t(4,2)=2.d0*z(2,1)*z(2,2)
      t(4,3)=2.d0*z(3,1)*z(3,2)
      t(4,4)=z(1,1)*z(2,2)+z(1,2)*z(2,1)
      t(4,5)=z(3,1)*z(1,2)+z(3,2)*z(1,1)
      t(4,6)=z(2,1)*z(3,2)+z(2,2)*z(3,1)
      t(5,1)=2.d0*z(1,3)*z(1,1)
      t(5,2)=2.d0*z(2,3)*z(2,1)
      t(5,3)=2.d0*z(3,3)*z(3,1)
      t(5,4)=z(1,3)*z(2,1)+z(1,1)*z(2,3)
      t(5,5)=z(3,3)*z(1,1)+z(3,1)*z(1,3)
      t(5,6)=z(2,3)*z(3,1)+z(2,1)*z(3,3)
      t(6,1)=2.d0*z(1,2)*z(1,3)
      t(6,2)=2.d0*z(2,2)*z(2,3)
      t(6,3)=2.d0*z(3,2)*z(3,3)
      t(6,4)=z(1,2)*z(2,3)+z(1,3)*z(2,2)
      t(6,5)=z(3,2)*z(1,3)+z(3,3)*z(1,2)
      t(6,6)=z(2,2)*z(3,3)+z(2,3)*z(3,2)
!     
!     transforming the stress into the global system
!     
      do i=1,6
        stre(i)=0.d0
        do j=1,6
          stre(i)=stre(i)+t(j,i)*sc(j)
        enddo
      enddo
c      write(*,*) 'mohrcoulomb s'
c      write(*,*) (stre(i),i=1,6)
!
!     
cccc
c      s(1,1)=sc(1)
c      s(2,2)=sc(2)
c      s(3,3)=sc(3)
c      s(1,2)=sc(4)
c      s(1,3)=sc(5)
c      s(2,3)=sc(6)
c      s(2,1)=s(1,2)
c      s(3,1)=s(1,3)
c      s(3,2)=s(2,3)
c      do jj=1,6
c        stre(jj)=0.d0
c        j1=kal(1,jj)
c        j2=kal(2,jj)
c        do j3=1,3
c          do j4=1,3
c            stre(jj)=stre(jj)+
c     &           s(j3,j4)*z(j1,j3)*z(j2,j4)
c          enddo
c        enddo
c      enddo
cccc
c      write(*,*) 'mohrcoulomb s'
c      write(*,*) (stre(i),i=1,6)
!     
      if(icmd.ne.3) then
!     
!     transforming the stiffness matrix into the global system
!
c        stiff(2,1)=(stiff(1,2)+stiff(2,1))/2.d0
c        stiff(3,1)=(stiff(1,3)+stiff(3,1))/2.d0
c        stiff(3,2)=(stiff(2,3)+stiff(3,2))/2.d0
c        stiff(4,1)=(stiff(1,4)+stiff(4,1))/2.d0
c        stiff(4,2)=(stiff(2,4)+stiff(4,2))/2.d0
c        stiff(4,3)=(stiff(3,4)+stiff(4,3))/2.d0
c        stiff(5,1)=(stiff(1,5)+stiff(5,1))/2.d0
c        stiff(5,2)=(stiff(2,5)+stiff(5,2))/2.d0
c        stiff(5,3)=(stiff(3,5)+stiff(5,3))/2.d0
c        stiff(5,4)=(stiff(4,5)+stiff(5,4))/2.d0
c        stiff(6,1)=(stiff(1,6)+stiff(6,1))/2.d0
c        stiff(6,2)=(stiff(2,6)+stiff(6,2))/2.d0
c        stiff(6,3)=(stiff(3,6)+stiff(6,3))/2.d0
c        stiff(6,4)=(stiff(4,6)+stiff(6,4))/2.d0
c        stiff(6,5)=(stiff(5,6)+stiff(6,5))/2.d0
c        do j=2,6
c          do i=1,j-1
c            stiff(i,j)=stiff(j,i)
c          enddo
c        enddo
!        
        do i=1,6
          do j=1,6
            dum(i,j)=0.d0
            do k=1,6
c     ERROR in Appendix A of Ph.D by Johan Clausen, p 1058
c      it should be C'=A^T.C.A               
c              dum(i,j)=dum(i,j)+stiff(i,k)*t(j,k)
              dum(i,j)=dum(i,j)+stiff(i,k)*t(k,j)
            enddo
          enddo
        enddo
!     
        do i=1,6
          do j=1,6
            stiff(i,j)=0.d0
            do k=1,6
              stiff(i,j)=stiff(i,j)+t(k,i)*dum(k,j)
c              stiff(i,j)=stiff(i,j)+t(i,k)*dum(k,j)
            enddo
          enddo
        enddo
!     
!     symmetrizing the matrix
!     
        elas(1)=stiff(1,1)
        elas(2)=(stiff(1,2)+stiff(2,1))/2.d0
        elas(3)=stiff(2,2)
        elas(4)=(stiff(1,3)+stiff(3,1))/2.d0
        elas(5)=(stiff(2,3)+stiff(3,2))/2.d0
        elas(6)=stiff(3,3)
        elas(7)=(stiff(1,4)+stiff(4,1))/2.d0
        elas(8)=(stiff(2,4)+stiff(4,2))/2.d0
        elas(9)=(stiff(3,4)+stiff(4,3))/2.d0
        elas(10)=stiff(4,4)
        elas(11)=(stiff(1,5)+stiff(5,1))/2.d0
        elas(12)=(stiff(2,5)+stiff(5,2))/2.d0
        elas(13)=(stiff(3,5)+stiff(5,3))/2.d0
        elas(14)=(stiff(4,5)+stiff(5,4))/2.d0
        elas(15)=stiff(5,5)
        elas(16)=(stiff(1,6)+stiff(6,1))/2.d0
        elas(17)=(stiff(2,6)+stiff(6,2))/2.d0
        elas(18)=(stiff(3,6)+stiff(6,3))/2.d0
        elas(19)=(stiff(4,6)+stiff(6,4))/2.d0
        elas(20)=(stiff(5,6)+stiff(6,5))/2.d0
        elas(21)=stiff(6,6)
ccccc
c        call anisotropic(elas,ya)
c!     
c        do jj=1,21
c          j1=kel(1,jj)
c          j2=kel(2,jj)
c          j3=kel(3,jj)
c          j4=kel(4,jj)
c          delas(jj)=0.d0
c          do j5=1,3
c            do j6=1,3
c              do j7=1,3
c                do j8=1,3
c                  delas(jj)=delas(jj)+ya(j5,j6,j7,j8)*
c     &                 z(j1,j5)*z(j2,j6)*z(j3,j7)*z(j4,j8)
c                enddo
c              enddo
c            enddo
c          enddo
c        enddo
c        do jj=1,21
c          elas(jj)=delas(jj)
c        enddo
ccccc
c        write(*,*) 'mohrcoulomb stiff'
c          write(*,*) (elas(i),i=1,1)
c          write(*,*) (elas(i),i=2,3)
c          write(*,*) (elas(i),i=4,6)
c          write(*,*) (elas(i),i=7,10)
c          write(*,*) (elas(i),i=11,15)
c          write(*,*) (elas(i),i=16,21)
      endif
!     
!     updating the plastic fields
!     
      do i=1,6
        xstate(1+i,iint,iel)=epl(i)
      enddo
      xstate(1,iint,iel)=ep
!     
      return
      end
