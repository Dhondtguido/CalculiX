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
      subroutine umat_aniso_creep(amat,iel,iint,kode,elconloc,emec,
     &        emec0,beta,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,iorien,
     &        pgauss,orab,nmethod,pnewdt,depvisc)
!
!     calculates stiffness and stresses for orthotropic elasticity
!     in combination with isotropic plasticity
!
!     Ref: Dhondt,G., The Finite Element Method for Three-Dimensional
!                     Thermomechanical Applications, Wiley (2004)
!                     Section 5.5
!
!     icmd=3: calculates stress at mechanical strain
!     else: calculates stress at mechanical strain and the stiffness
!           matrix
!
!     INPUT:
!
!     amat               material name
!     iel                element number
!     iint               integration point number
!
!     kode               material type (-100-#of constants entered
!                        under *USER MATERIAL): can be used for materials
!                        with varying number of constants
!
!     elconloc(*)        user defined constants defined by the keyword
!                        card *USER MATERIAL (max. 21, actual # =
!                        -kode-100), interpolated for the
!                        actual temperature t1l
!
!     emec(6)            Lagrange mechanical strain tensor (component order:
!                        11,22,33,12,13,23) at the end of the increment
!                        (thermal strains are subtracted)
!     emec0(6)           Lagrange mechanical strain tensor at the start of the
!                        increment (thermal strains are subtracted)
!     beta(6)            residual stress tensor (the stress entered under
!                        the keyword *INITIAL CONDITIONS,TYPE=STRESS)
!
!     xkl(3,3)           deformation gradient at the end of the increment
!     vj                 Jacobian at the end of the increment
!
!     ithermal           0: no thermal effects are taken into account
!                        >0: thermal effects are taken into account (triggered
!                        by the keyword *INITIAL CONDITIONS,TYPE=TEMPERATURE)
!     t1l                temperature at the end of the increment
!     dtime              time length of the increment
!     time               step time at the end of the current increment
!     ttime              total time at the start of the current step
!
!     icmd               not equal to 3: calculate stress and stiffness
!                        3: calculate only stress
!     ielas              0: no elastic iteration: irreversible effects
!                        are allowed
!                        1: elastic iteration, i.e. no irreversible
!                           deformation allowed
!
!     mi(1)              max. # of integration points per element in the
!                        model
!     nstate_            max. # of state variables in the model
!
!     xstateini(nstate_,mi(1),# of elements)
!                        state variables at the start of the increment
!     xstate(nstate_,mi(1),# of elements)
!                        state variables at the end of the increment
!
!     stre(6)            Piola-Kirchhoff stress of the second kind
!                        at the start of the increment
!
!     iorien             number of the local coordinate axis system
!                        in the integration point at stake (takes the value
!                        0 if no local system applies)
!     pgauss(3)          global coordinates of the integration point
!     orab(7,*)          description of all local coordinate systems.
!                        If a local coordinate system applies the global 
!                        tensors can be obtained by premultiplying the local
!                        tensors with skl(3,3). skl is  determined by calling
!                        the subroutine transformatrix: 
!                        call transformatrix(orab(1,iorien),pgauss,skl)
!
!
!     OUTPUT:
!
!     xstate(nstate_,mi(1),# of elements)
!                        updated state variables at the end of the increment
!     stre(6)            Piola-Kirchhoff stress of the second kind at the
!                        end of the increment
!     stiff(21):         consistent tangent stiffness matrix in the material
!                        frame of reference at the end of the increment. In
!                        other words: the derivative of the PK2 stress with
!                        respect to the Lagrangian strain tensor. The matrix
!                        is supposed to be symmetric, only the upper half is
!                        to be given in the same order as for a fully
!                        anisotropic elastic material (*ELASTIC,TYPE=ANISO).
!                        Notice that the matrix is an integral part of the 
!                        fourth order material tensor, i.e. the Voigt notation
!                        is not used.
!
      implicit none
!     
      character*80 amat
!     
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &     i,j,ipiv(6),info,neq,lda,ldb,j1,j2,j3,j4,j5,j6,j7,j8,
     &     nrhs,iplas,kel(4,21),nmethod,iloop,leximp,lend,layer,kspt,
     &     kstep,kinc,exitcriterion
!     
      real*8 ep0(6),eeq,ep(6),b,Pn(6),ec(2),timeabq(2),dtemp,
     &     dg,ddg,c(21),x(21),cm1(21),svm,stri(6),htri,sg(6),r(6),
     &     ee(6),dd,gl(6,6),gr(6,6),c0,c1,c2,pnewdt,
     &     skl(3,3),gcreep,gm1,ya(3,3,3,3),dsg,detc,strinv,
     &     elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &     vj,t1l,dtime,xkl(3,3),pgauss(3),orab(7,*),
     &     time,ttime,xstate(nstate_,mi(1),*),depvisc,
     &     xstateini(nstate_,mi(1),*),decra(5),deswa(5),serd,esw(2),
     &     p,predef(1),dpred(1)
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!     
      iloop=0
!     
      c0=dsqrt(2.d0/3.d0)
      c1=2.d0/3.d0
      c2=-1.d0/3.d0
!     
!     elastic constants
!     
      if(iorien.gt.0) then
!     
        call transformatrix(orab(1,iorien),pgauss,skl)
!     
        call orthotropic(elconloc,ya)
!     
        do j=1,21
          j1=kel(1,j)
          j2=kel(2,j)
          j3=kel(3,j)
          j4=kel(4,j)
          c(j)=0.d0
          do j5=1,3
            do j6=1,3
              do j7=1,3
                do j8=1,3
                  c(j)=c(j)+ya(j5,j6,j7,j8)*
     &                 skl(j1,j5)*skl(j2,j6)*skl(j3,j7)*skl(j4,j8)
                enddo
              enddo
            enddo
          enddo
        enddo
!     
      else
        do i=1,9
          c(i)=elconloc(i)
        enddo
      endif
!     
!     state variables
!     
!     equivalent plastic strain
!     
      eeq=xstateini(1,iint,iel)
!     
!     plastic strain
!     
      do i=1,6
        ep0(i)=xstateini(1+i,iint,iel)
      enddo
!     
!     elastic strains
!     
      do i=1,6
        ee(i)=emec(i)-ep0(i)
      enddo
!     
!     global trial stress tensor
!     
      if(iorien.gt.0) then
        stri(1)=c(1)*ee(1)+c(2)*ee(2)+c(4)*ee(3)+
     &       2.d0*(c(7)*ee(4)+c(11)*ee(5)+c(16)*ee(6))
     &       -beta(1)
        stri(2)=c(2)*ee(1)+c(3)*ee(2)+c(5)*ee(3)+
     &       2.d0*(c(8)*ee(4)+c(12)*ee(5)+c(17)*ee(6))
     &       -beta(2)
        stri(3)=c(4)*ee(1)+c(5)*ee(2)+c(6)*ee(3)+
     &       2.d0*(c(9)*ee(4)+c(13)*ee(5)+c(18)*ee(6))
     &       -beta(3)
        stri(4)=c(7)*ee(1)+c(8)*ee(2)+c(9)*ee(3)+
     &       2.d0*(c(10)*ee(4)+c(14)*ee(5)+c(19)*ee(6))
     &       -beta(4)
        stri(5)=c(11)*ee(1)+c(12)*ee(2)+c(13)*ee(3)+
     &       2.d0*(c(14)*ee(4)+c(15)*ee(5)+c(20)*ee(6))
     &       -beta(5)
        stri(6)=c(16)*ee(1)+c(17)*ee(2)+c(18)*ee(3)+
     &       2.d0*(c(19)*ee(4)+c(20)*ee(5)+c(21)*ee(6))
     &       -beta(6)
      else
        stri(1)=c(1)*ee(1)+c(2)*ee(2)+c(4)*ee(3)-beta(1)
        stri(2)=c(2)*ee(1)+c(3)*ee(2)+c(5)*ee(3)-beta(2)
        stri(3)=c(4)*ee(1)+c(5)*ee(2)+c(6)*ee(3)-beta(3)
        stri(4)=2.d0*c(7)*ee(4)-beta(4)
        stri(5)=2.d0*c(8)*ee(5)-beta(5)
        stri(6)=2.d0*c(9)*ee(6)-beta(6)
      endif
!     
!     stress radius (only deviatoric part of stress enters)
!     
      strinv=(stri(1)+stri(2)+stri(3))/3.d0
      do i=1,3
        sg(i)=stri(i)-strinv
      enddo
      do i=4,6
        sg(i)=stri(i)
      enddo
      dsg=dsqrt(sg(1)*sg(1)+sg(2)*sg(2)+sg(3)*sg(3)+
     &     2.d0*(sg(4)*sg(4)+sg(5)*sg(5)+sg(6)*sg(6)))
!     
!     evaluation of the yield surface
!     
      htri=dsg
!     
!     check whether plasticity occurs
!     
      if(htri.gt.0.d0) then
        iplas=1
      else
        iplas=0
      endif
!     
      if((iplas.eq.0).or.(ielas.eq.1).or.(dtime.lt.1.d-30).or.
     &                   ((nmethod.eq.1).and.(ithermal(1).ne.3))) then
!     
!     elastic stress
!     
        do i=1,6
          stre(i)=stri(i)
        enddo
!     
!     updating the state variables
!     
        xstate(1,iint,iel)=eeq
        do i=1,6
          xstate(1+i,iint,iel)=ep0(i)
        enddo
!     
!     elastic stiffness
!     
        if(icmd.ne.3) then
          if(iorien.gt.0) then
            do i=1,21
              stiff(i)=c(i)
            enddo
          else
            stiff(1)=c(1)
            stiff(2)=c(2)
            stiff(3)=c(3)
            stiff(4)=c(4)
            stiff(5)=c(5)
            stiff(6)=c(6)
            stiff(7)=0.d0
            stiff(8)=0.d0
            stiff(9)=0.d0
            stiff(10)=c(7)
            stiff(11)=0.d0
            stiff(12)=0.d0
            stiff(13)=0.d0
            stiff(14)=0.d0
            stiff(15)=c(8)
            stiff(16)=0.d0
            stiff(17)=0.d0
            stiff(18)=0.d0
            stiff(19)=0.d0
            stiff(20)=0.d0
            stiff(21)=c(9)
          endif
        endif
!     
        return
      endif
!     
!     plastic deformation
!     
      neq=6
      nrhs=1
      lda=6
      ldb=6
!     
!     initializing the state variables
!     
      do i=1,6
        ep(i)=ep0(i)
      enddo
      dg=0.d0
!     
!     determining the inverse of C (field cm)
!     
      if(iorien.gt.0) then
!     
!     solve gl:C=gr
!     
        gl(1,1)=c(1) 
        gl(1,2)=c(2) 
        gl(2,2)=c(3) 
        gl(1,3)=c(4) 
        gl(2,3)=c(5) 
        gl(3,3)=c(6) 
        gl(1,4)=c(7) 
        gl(2,4)=c(8) 
        gl(3,4)=c(9) 
        gl(4,4)=c(10)
        gl(1,5)=c(11)
        gl(2,5)=c(12)
        gl(3,5)=c(13)
        gl(4,5)=c(14)
        gl(5,5)=c(15)
        gl(1,6)=c(16)
        gl(2,6)=c(17)
        gl(3,6)=c(18)
        gl(4,6)=c(19)
        gl(5,6)=c(20)
        gl(6,6)=c(21)
        do i=1,6
          do j=1,i-1
            gl(i,j)=gl(j,i)
          enddo
        enddo
        do i=1,6
          do j=4,6
            gl(i,j)=2.d0*gl(i,j)
          enddo
        enddo
        do i=1,6
          do j=1,6
            gr(i,j)=0.d0
          enddo
          if(i.le.3) then
            gr(i,i)=1.d0
          else
            gr(i,i)=0.5d0
          endif
        enddo
        nrhs=6
        call dgesv(neq,nrhs,gl,lda,ipiv,gr,ldb,info)
        if(info.ne.0) then
          write(*,*) '*ERROR in sc.f: linear equation solver'
          write(*,*) '       exited with error: info = ',info
          call exit(201)
        endif
        nrhs=1
        cm1(1)=gr(1,1)
        cm1(2)=gr(1,2)
        cm1(3)=gr(2,2)
        cm1(4)=gr(1,3)
        cm1(5)=gr(2,3)
        cm1(6)=gr(3,3)
        cm1(7)=gr(1,4)
        cm1(8)=gr(2,4)
        cm1(9)=gr(3,4)
        cm1(10)=gr(4,4)
        cm1(11)=gr(1,5)
        cm1(12)=gr(2,5)
        cm1(13)=gr(3,5)
        cm1(14)=gr(4,5)
        cm1(15)=gr(5,5)
        cm1(16)=gr(1,6)
        cm1(17)=gr(2,6)
        cm1(18)=gr(3,6)
        cm1(19)=gr(4,6)
        cm1(20)=gr(5,6)
        cm1(21)=gr(6,6)
      else
        detc=c(1)*(c(3)*c(6)-c(5)*c(5))-
     &       c(2)*(c(2)*c(6)-c(4)*c(5))+
     &       c(4)*(c(2)*c(5)-c(4)*c(3))
        cm1(1)=(c(3)*c(6)-c(5)*c(5))/detc
        cm1(2)=(c(5)*c(4)-c(2)*c(6))/detc
        cm1(3)=(c(1)*c(6)-c(4)*c(4))/detc
        cm1(4)=(c(2)*c(5)-c(3)*c(4))/detc
        cm1(5)=(c(2)*c(4)-c(1)*c(5))/detc
        cm1(6)=(c(1)*c(3)-c(2)*c(2))/detc
        cm1(7)=1.d0/(4.d0*c(7))
        cm1(8)=1.d0/(4.d0*c(8))
        cm1(9)=1.d0/(4.d0*c(9))
      endif
!     
!     loop
!     
      do
!     
        iloop=iloop+1
        if(iloop.gt.100) then
c     NOTE: write statements cause problems for
c     parallellized execution
c     write(*,*) '*WARNING in umat_aniso_plas: material loop'
c     write(*,*) '         did not converge in integration'
c     write(*,*) '         point',iint,'in element',iel,';'
c     write(*,*) '         the increment size is reduced'
c     write(*,*)
          pnewdt=0.25d0
          return
        endif
!
!       for iloop=1 the values were calculated before the loop
!
        if(iloop.gt.1) then
!     
!     elastic strains
!     
          do i=1,6
            ee(i)=emec(i)-ep(i)
          enddo
!     
!     global trial stress tensor
!     
          if(iorien.gt.0) then
            stri(1)=c(1)*ee(1)+c(2)*ee(2)+c(4)*ee(3)+
     &           2.d0*(c(7)*ee(4)+c(11)*ee(5)+c(16)*ee(6))
     &           -beta(1)
            stri(2)=c(2)*ee(1)+c(3)*ee(2)+c(5)*ee(3)+
     &           2.d0*(c(8)*ee(4)+c(12)*ee(5)+c(17)*ee(6))
     &           -beta(2)
            stri(3)=c(4)*ee(1)+c(5)*ee(2)+c(6)*ee(3)+
     &           2.d0*(c(9)*ee(4)+c(13)*ee(5)+c(18)*ee(6))
     &           -beta(3)
            stri(4)=c(7)*ee(1)+c(8)*ee(2)+c(9)*ee(3)+
     &           2.d0*(c(10)*ee(4)+c(14)*ee(5)+c(19)*ee(6))
     &           -beta(4)
            stri(5)=c(11)*ee(1)+c(12)*ee(2)+c(13)*ee(3)+
     &           2.d0*(c(14)*ee(4)+c(15)*ee(5)+c(20)*ee(6))
     &           -beta(5)
            stri(6)=c(16)*ee(1)+c(17)*ee(2)+c(18)*ee(3)+
     &           2.d0*(c(19)*ee(4)+c(20)*ee(5)+c(21)*ee(6))
     &           -beta(6)
          else
            stri(1)=c(1)*ee(1)+c(2)*ee(2)+c(4)*ee(3)-beta(1)
            stri(2)=c(2)*ee(1)+c(3)*ee(2)+c(5)*ee(3)-beta(1)
            stri(3)=c(4)*ee(1)+c(5)*ee(2)+c(6)*ee(3)-beta(1)
            stri(4)=2.d0*c(7)*ee(4)-beta(4)
            stri(5)=2.d0*c(8)*ee(5)-beta(5)
            stri(6)=2.d0*c(9)*ee(6)-beta(6)
          endif
!     
!     stress radius (only deviatoric part of stress enters)
!     
          strinv=(stri(1)+stri(2)+stri(3))/3.d0
          do i=1,3
            sg(i)=stri(i)-strinv
          enddo
          do i=4,6
            sg(i)=stri(i)
          enddo
          dsg=dsqrt(sg(1)*sg(1)+sg(2)*sg(2)+sg(3)*sg(3)+
     &         2.d0*(sg(4)*sg(4)+sg(5)*sg(5)+sg(6)*sg(6)))
        endif
!     
!     evaluation of the yield surface
!     
         ec(1)=eeq
         decra(1)=c0*dg
         timeabq(1)=time
         timeabq(2)=ttime+time
         call creep(decra,deswa,xstateini(1,iint,iel),serd,ec,
     &        esw,p,svm,t1l,dtemp,predef,dpred,timeabq,dtime,
     &        amat,leximp,lend,pgauss,nstate_,iel,iint,layer,kspt,
     &        kstep,kinc)
!
!        if the creep routine returns an increased value of decra(1)
!        it means that there is a lower cut-off for decra(1);
!        if the routine stays in a range lower than this cut-off,
!        it will never leave it and the exit conditions are
!        assumed to be satisfied.
!
         if(decra(1).gt.c0*dg) then
            dg=decra(1)/c0
            if(iloop.gt.1) exitcriterion=1
         endif
!     
         htri=dsg-c0*svm
!
        do i=1,6
          sg(i)=sg(i)/dsg
        enddo
!     
!     determining the residual matrix
!     
        do i=1,6  
          r(i)=ep0(i)-ep(i)+dg*sg(i)
        enddo
!     
!     check convergence (for iloop=1 the remaining parts of the loop
!     have not been executed: => P:n ... needed for the stiffness
!     matrix are not calculated)
!
        if(exitcriterion.eq.1) exit
!
        if(iloop.gt.1) then
          if((dabs(htri).le.1.d-5).or.(dabs(ddg).lt.1.d-3*dabs(dg)))
     &         then
          dd=0.d0
          do i=1,6
            dd=dd+r(i)*r(i)
          enddo
          dd=sqrt(dd)
          if(dd.le.1.d-10) then
            exit
          endif
        endif
      endif
!     
!     determining b.x
!     
        b=dg/dsg
!     
        x(1)=b*(c1-sg(1)*sg(1))
        x(2)=b*(c2-sg(1)*sg(2))
        x(3)=b*(c1-sg(2)*sg(2))
        x(4)=b*(c2-sg(1)*sg(3))
        x(5)=b*(c2-sg(2)*sg(3))
        x(6)=b*(c1-sg(3)*sg(3))
        x(7)=-b*sg(1)*sg(4)
        x(8)=-b*sg(2)*sg(4)
        x(9)=-b*sg(3)*sg(4)
        x(10)=b*(.5d0-sg(4)*sg(4))
        x(11)=-b*sg(1)*sg(5)
        x(12)=-b*sg(2)*sg(5)
        x(13)=-b*sg(3)*sg(5)
        x(14)=-b*sg(4)*sg(5)
        x(15)=b*(.5d0-sg(5)*sg(5))
        x(16)=-b*sg(1)*sg(6)
        x(17)=-b*sg(2)*sg(6)
        x(18)=-b*sg(3)*sg(6)
        x(19)=-b*sg(4)*sg(6)
        x(20)=-b*sg(5)*sg(6)
        x(21)=b*(.5d0-sg(6)*sg(6))
!     
!     filling the LHS (Equation 5.336)
!     
        if(iorien.gt.0) then
          gl(1,1)=cm1(1)+x(1)
          gl(1,2)=cm1(2)+x(2)
          gl(2,2)=cm1(3)+x(3)
          gl(1,3)=cm1(4)+x(4)
          gl(2,3)=cm1(5)+x(5)
          gl(3,3)=cm1(6)+x(6)
          gl(1,4)=cm1(7)+x(7)
          gl(2,4)=cm1(8)+x(8)
          gl(3,4)=cm1(9)+x(9)
          gl(4,4)=cm1(10)+x(10)
          gl(1,5)=cm1(11)+x(11)
          gl(2,5)=cm1(12)+x(12)
          gl(3,5)=cm1(13)+x(13)
          gl(4,5)=cm1(14)+x(14)
          gl(5,5)=cm1(15)+x(15)
          gl(1,6)=cm1(16)+x(16)
          gl(2,6)=cm1(17)+x(17)
          gl(3,6)=cm1(18)+x(18)
          gl(4,6)=cm1(19)+x(19)
          gl(5,6)=cm1(20)+x(20)
          gl(6,6)=cm1(21)+x(21)
          do i=1,6
            do j=1,i-1
              gl(i,j)=gl(j,i)
            enddo
          enddo
          do i=1,6
            do j=4,6
              gl(i,j)=2.d0*gl(i,j)
            enddo
          enddo
        else
          gl(1,1)=cm1(1)+x(1)
          gl(1,2)=cm1(2)+x(2)
          gl(2,2)=cm1(3)+x(3)
          gl(1,3)=cm1(4)+x(4)
          gl(2,3)=cm1(5)+x(5)
          gl(3,3)=cm1(6)+x(6)
          gl(1,4)=x(7)
          gl(2,4)=x(8)
          gl(3,4)=x(9)
          gl(4,4)=cm1(7)+x(10)
          gl(1,5)=x(11)
          gl(2,5)=x(12)
          gl(3,5)=x(13)
          gl(4,5)=x(14)
          gl(5,5)=cm1(8)+x(15)
          gl(1,6)=x(16)
          gl(2,6)=x(17)
          gl(3,6)=x(18)
          gl(4,6)=x(19)
          gl(5,6)=x(20)
          gl(6,6)=cm1(9)+x(21)
          do i=1,6
            do j=1,i-1
              gl(i,j)=gl(j,i)
            enddo
          enddo
          do i=1,6
            do j=4,6
              gl(i,j)=2.d0*gl(i,j)
            enddo
          enddo
        endif
!     
!     filling the RHS
!     
        do i=1,6
          gr(i,1)=sg(i)
        enddo
!     
!     solve gl:(P:n)=gr
!     
        call dgesv(neq,nrhs,gl,lda,ipiv,gr,ldb,info)
        if(info.ne.0) then
          write(*,*) '*ERROR in sc.f: linear equation solver'
          write(*,*) '       exited with error: info = ',info
          call exit(201)
        endif
!     
!       Pn = P:n     
!     
        do i=1,6
          Pn(i)=gr(i,1)
        enddo
!     
!     calculating the creep contribution
!     
        if(dg.le.0.d0) then
          ec(1)=eeq
          decra(1)=c0*1.d-10
          timeabq(1)=time
          timeabq(2)=ttime+time
          call creep(decra,deswa,xstateini(1,iint,iel),serd,ec,
     &         esw,p,svm,t1l,dtemp,predef,dpred,timeabq,dtime,
     &         amat,leximp,lend,pgauss,nstate_,iel,iint,layer,kspt,
     &         kstep,kinc)
        endif
        gcreep=c1/decra(5)
!
!     calculating the correction to the consistency parameter
!     
        gm1=Pn(1)*sg(1)+Pn(2)*sg(2)+Pn(3)*sg(3)+
     &       2.d0*(Pn(4)*sg(4)+Pn(5)*sg(5)+Pn(6)*sg(6))
        gm1=1.d0/(gm1+gcreep)
        ddg=gm1*(htri-(Pn(1)*r(1)+Pn(2)*r(2)+Pn(3)*r(3)+
     &       2.d0*(Pn(4)*r(4)+Pn(5)*r(5)+Pn(6)*r(6))))
!     
!     updating the residual matrix
!     
        do i=1,6
          r(i)=r(i)+ddg*sg(i)
        enddo
!     
!     update the plastic strain
!     
        gr(1,1)=r(1)
        gr(2,1)=r(2)
        gr(3,1)=r(3)
        gr(4,1)=r(4)
        gr(5,1)=r(5)
        gr(6,1)=r(6)
!     
!       Solving Equation (5.349)    
!     
        call dgetrs('No transpose',neq,nrhs,gl,lda,ipiv,gr,ldb,info)
        if(info.ne.0) then
          write(*,*) '*ERROR in sc.f: linear equation solver'
          write(*,*) '       exited with error: info = ',info
          call exit(201)
        endif
!     
        if(iorien.gt.0) then
          ep(1)=ep(1)+cm1(1)*gr(1,1)+cm1(2)*gr(2,1)+cm1(4)*gr(3,1)+
     &         2.d0*(cm1(7)*gr(4,1)+cm1(11)*gr(5,1)+cm1(16)*gr(6,1))
          ep(2)=ep(2)+cm1(2)*gr(1,1)+cm1(3)*gr(2,1)+cm1(5)*gr(3,1)+
     &         2.d0*(cm1(8)*gr(4,1)+cm1(12)*gr(5,1)+cm1(17)*gr(6,1))
          ep(3)=ep(3)+cm1(4)*gr(1,1)+cm1(5)*gr(2,1)+cm1(6)*gr(3,1)+
     &         2.d0*(cm1(9)*gr(4,1)+cm1(13)*gr(5,1)+cm1(18)*gr(6,1))
          ep(4)=ep(4)+cm1(7)*gr(1,1)+cm1(8)*gr(2,1)+cm1(9)*gr(3,1)+
     &         2.d0*(cm1(10)*gr(4,1)+cm1(14)*gr(5,1)+cm1(19)*gr(6,1))
          ep(5)=ep(5)+cm1(11)*gr(1,1)+cm1(12)*gr(2,1)+cm1(13)*gr(3,1)+
     &         2.d0*(cm1(14)*gr(4,1)+cm1(15)*gr(5,1)+cm1(20)*gr(6,1))
          ep(6)=ep(6)+cm1(16)*gr(1,1)+cm1(17)*gr(2,1)+cm1(18)*gr(3,1)+
     &         2.d0*(cm1(19)*gr(4,1)+cm1(20)*gr(5,1)+cm1(21)*gr(6,1))
        else
          ep(1)=ep(1)+cm1(1)*gr(1,1)+cm1(2)*gr(2,1)+cm1(4)*gr(3,1)
          ep(2)=ep(2)+cm1(2)*gr(1,1)+cm1(3)*gr(2,1)+cm1(5)*gr(3,1)
          ep(3)=ep(3)+cm1(4)*gr(1,1)+cm1(5)*gr(2,1)+cm1(6)*gr(3,1)
          ep(4)=ep(4)+2.d0*cm1(7)*gr(4,1)
          ep(5)=ep(5)+2.d0*cm1(8)*gr(5,1)
          ep(6)=ep(6)+2.d0*cm1(9)*gr(6,1)
        endif
!     
!     update the consistency parameter
!     
        dg=dg+ddg
!     
!     end of major loop
!     
      enddo
!     
!     storing the stress
!     
      do i=1,6
        stre(i)=stri(i)
      enddo
!     
!     calculating the tangent stiffness matrix
!     
      if(icmd.ne.3) then
!     
!     determining p
!     
        gr(1,1)=1.d0 
        gr(1,2)=0.d0 
        gr(2,2)=1.d0 
        gr(1,3)=0.d0 
        gr(2,3)=0.d0 
        gr(3,3)=1.d0 
        gr(1,4)=0.d0
        gr(2,4)=0.d0 
        gr(3,4)=0.d0 
        gr(4,4)=0.5d0
        gr(1,5)=0.d0
        gr(2,5)=0.d0
        gr(3,5)=0.d0
        gr(4,5)=0.d0
        gr(5,5)=0.5d0
        gr(1,6)=0.d0
        gr(2,6)=0.d0
        gr(3,6)=0.d0
        gr(4,6)=0.d0
        gr(5,6)=0.d0
        gr(6,6)=0.5d0
        do i=1,6
          do j=1,i-1
            gr(i,j)=gr(j,i)
          enddo
        enddo
        nrhs=6
!     
        call dgetrs('No transpose',neq,nrhs,gl,lda,ipiv,gr,ldb,info)
        if(info.ne.0) then
          write(*,*) '*ERROR in sc.f: linear equation solver'
          write(*,*) '       exited with error: info = ',info
          call exit(201)
        endif
!     
        stiff(1)=gr(1,1)-gm1*Pn(1)*Pn(1)
        stiff(2)=gr(1,2)-gm1*Pn(1)*Pn(2)
        stiff(3)=gr(2,2)-gm1*Pn(2)*Pn(2)
        stiff(4)=gr(1,3)-gm1*Pn(1)*Pn(3)
        stiff(5)=gr(2,3)-gm1*Pn(2)*Pn(3)
        stiff(6)=gr(3,3)-gm1*Pn(3)*Pn(3)
        stiff(7)=gr(1,4)-gm1*Pn(1)*Pn(4)
        stiff(8)=gr(2,4)-gm1*Pn(2)*Pn(4)
        stiff(9)=gr(3,4)-gm1*Pn(3)*Pn(4)
        stiff(10)=gr(4,4)-gm1*Pn(4)*Pn(4)
        stiff(11)=gr(1,5)-gm1*Pn(1)*Pn(5)
        stiff(12)=gr(2,5)-gm1*Pn(2)*Pn(5)
        stiff(13)=gr(3,5)-gm1*Pn(3)*Pn(5)
        stiff(14)=gr(4,5)-gm1*Pn(4)*Pn(5)
        stiff(15)=gr(5,5)-gm1*Pn(5)*Pn(5)
        stiff(16)=gr(1,6)-gm1*Pn(1)*Pn(6)
        stiff(17)=gr(2,6)-gm1*Pn(2)*Pn(6)
        stiff(18)=gr(3,6)-gm1*Pn(3)*Pn(6)
        stiff(19)=gr(4,6)-gm1*Pn(4)*Pn(6)
        stiff(20)=gr(5,6)-gm1*Pn(5)*Pn(6)
        stiff(21)=gr(6,6)-gm1*Pn(6)*Pn(6)
!     
      endif
!     
!     updating the state variables
!     
      xstate(1,iint,iel)=eeq+c0*dg
      do i=1,6
        xstate(1+i,iint,iel)=ep(i)
      enddo
!     
      depvisc=0.d0
!     
      return
      end
