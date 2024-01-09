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
      subroutine umat_abaqusnl(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,kinc,pnewdt,nmethod,iperturb)
!
!     calculates stiffness and stresses for a nonlinear material
!     defined by an ABAQUS umat routine
!
!     CAVE: the thermal expansion MUST be ISOTROPIC for this
!           routine to work!!!!
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
!                        card *USER MATERIAL (actual # =
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
!     xokl(3,3)          deformation gradient at the start of the increment
!     voj                Jacobian at the start of the increment
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
!
!     This routine allows for the use of an ABAQUS umat user subroutine
!     in CalculiX. 
!
!     Note that the following fields are not supported
!     so far: sse,spd,scd,rpl,ddsddt,drplde,drpldt,predef,
!     dpred,pnewdt,celent,layer,kspt
!
!     Furthermore, the following fields have a different meaning in
!     ABAQUS and CalculiX:
!
!     temp:  in CalculiX: temperature at the end of the increment
!              in ABAQUS: temperature at the start of the increment
!     dtemp: in CalculiX: zero
!              in ABAQUS: temperature increment
!
!     Reference for the derivative of the logarithmic strain tensor
!     and the inverse right stretch tensor w.r.t. the Cauchy-Green
!     tensor: Dhondt, G., Hackenberg, H.P., Use of a rotation-invariant
!     linear strain measure for linear elastic crack propagation
!     calculations, Eng.Frac.Mech. 247(2021) 107634.
!
      implicit none
!     
      character*80 amat
!     
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),i,
     &     iorien,nmethod,iperturb(*),istep,keltot(4,36),
     &     ndi,nshr,ntens,nprops,layer,kspt,jstep(4),kinc,kal(2,6),
     &     kel(4,21),j1,j2,j3,j4,j5,j6,j7,j8,jj,n,ier,j,matz
!     
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &     vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &     time,ttime,skl(3,3),xa(3,3),ya(3,3,3,3),
     &     xstateini(nstate_,mi(1),*),w(3),fv1(3),fv2(3),d(6),c(6),
     &     v1,v2,v3,c0(6),r(3,3),r0(3,3),eln0(6),eln(6),e(3,3),tkl(3,3),
     &     u(6),c2(6),dd,um1(3,3),z(3,3),u0(3,3),ctot(3,3),expansion,
     &     c0tot(3,3),xm1(6),xm2(6),xm3(6),xstate(nstate_,mi(1),*)
!     
      real*8 ddsdde(6,6),sse,spd,scd,rpl,ddsddt(6),drplde(6),
     &     drpldt,stran(6),dstran(6),abqtime(2),predef,temp,dtemp,
     &     dpred,drot(3,3),celent,pnewdt
!     
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!     
      d=(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
!     
!     filling field jstep
!     
      jstep(1)=istep
      jstep(2)=nmethod
      jstep(3)=iperturb(2)
      if(iperturb(1).eq.1) then
        jstep(4)=1
      else
        jstep(4)=0
      endif
!     
!     calculating the logarithmic mechanical strain at the
!     start of the increment
!     
      e(1,1)=emec0(1)
      e(2,2)=emec0(2)
      e(3,3)=emec0(3)
      e(1,2)=emec0(4)
      e(1,3)=emec0(5)
      e(2,3)=emec0(6)
      e(2,1)=emec0(4)
      e(3,1)=emec0(5)
      e(3,2)=emec0(6)
!     
!     calculating the eigenvalues and eigenvectors
!     
      n=3
      matz=1
!     
      call rs(n,n,e,w,matz,z,fv1,fv2,ier)
!     
      if(ier.ne.0) then
        write(*,*) '
     &*ERROR calculating the eigenvalues/vectors in umat_abaqusnl'
        call exit(201)
      endif
!     
!     calculating the principal stretches at the start of the increment
!     
      do i=1,3
        if(w(i).lt.-0.5d0) then
          write(*,*) '*ERROR in umat_abaqusnl: negative eigenvalue'
          write(*,*) '       of the Cauchy-Green tensor;'
          call exit(201)
        else
          w(i)=dsqrt(2.d0*w(i)+1.d0)
        endif
      enddo
!     
!     calculating the invariants at the start of the increment
!     
      v1=w(1)+w(2)+w(3)
      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
      v3=w(1)*w(2)*w(3)
!     
!     calculating the right Cauchy-Green tensor at the start of the
!     increment
!     
      do i=1,3
        c0(i)=2.d0*emec0(i)+1.d0
      enddo
      do i=4,6
        c0(i)=2.d0*emec0(i)
      enddo
!     
!     calculating the square of the right Cauchy-Green tensor at the
!     start of the increment
!     
      c2(1)=c0(1)*c0(1)+c0(4)*c0(4)+c0(5)*c0(5)
      c2(2)=c0(4)*c0(4)+c0(2)*c0(2)+c0(6)*c0(6)
      c2(3)=c0(5)*c0(5)+c0(6)*c0(6)+c0(3)*c0(3)
      c2(4)=c0(1)*c0(4)+c0(4)*c0(2)+c0(5)*c0(6)
      c2(5)=c0(1)*c0(5)+c0(4)*c0(6)+c0(5)*c0(3)
      c2(6)=c0(4)*c0(5)+c0(2)*c0(6)+c0(6)*c0(3)
!     
!     calculating the right stretch tensor at the start of the increment
!     (cf. Simo and Hughes, Computational Inelasticity)
!     
!     This is the mechanical right stretch tensor, i.e. without
!     thermal expansion (since it is based on emec0)
!     
      dd=v1*v2-v3
      do i=1,6
        u(i)=(-c2(i)+(v1*v1-v2)*c0(i)+v1*v3*d(i))/dd
      enddo
!     
      u0(1,1)=u(1)
      u0(2,2)=u(2)
      u0(3,3)=u(3)
      u0(1,2)=u(4)
      u0(1,3)=u(5)
      u0(2,3)=u(6)
      u0(2,1)=u(4)
      u0(3,1)=u(5)
      u0(3,2)=u(6)
!     
!     calculating the inverse of the right stretch tensor at the start
!     of the increment
!     
      um1(1,1)=(c0(1)-v1*u(1)+v2)/v3
      um1(2,2)=(c0(2)-v1*u(2)+v2)/v3
      um1(3,3)=(c0(3)-v1*u(3)+v2)/v3
      um1(1,2)=(c0(4)-v1*u(4))/v3
      um1(1,3)=(c0(5)-v1*u(5))/v3
      um1(2,3)=(c0(6)-v1*u(6))/v3
      um1(2,1)=um1(1,2)
      um1(3,1)=um1(1,3)
      um1(3,2)=um1(2,3)
!     
!     calculation of the local rotation tensor at the start of the
!     increment
!     
      do i=1,3
        do j=1,3
          r0(i,j)=xokl(i,1)*um1(1,j)+xokl(i,2)*um1(2,j)+
     &         xokl(i,3)*um1(3,j)
        enddo
      enddo
!     
!     calculating the logarithmic strain at the start of the increment
!     
      do i=1,3
        w(i)=dlog(w(i))
      enddo
!     
!     logarithmic strain in global coordinates at the start of the
!     increment = ln(lambda) N diadic N;
!     
!     the true logarithmic strain is ln(lambda) n diadic n;
!     the rotation tensor R = n diadic N;
!     therefore, ln(lambda) N diadic N is R^T true logarithmic strain R, 
!     in other words the corotational true logarithmic strain (required
!     by the abaqus routine)
!     
      eln0(1)=z(1,1)*z(1,1)*w(1)+z(1,2)*z(1,2)*w(2)+
     &     z(1,3)*z(1,3)*w(3)          
      eln0(2)=z(2,1)*z(2,1)*w(1)+z(2,2)*z(2,2)*w(2)+
     &     z(2,3)*z(2,3)*w(3)          
      eln0(3)=z(3,1)*z(3,1)*w(1)+z(3,2)*z(3,2)*w(2)+
     &     z(3,3)*z(3,3)*w(3)          
      eln0(4)=z(1,1)*z(2,1)*w(1)+z(1,2)*z(2,2)*w(2)+
     &     z(1,3)*z(2,3)*w(3)          
      eln0(5)=z(1,1)*z(3,1)*w(1)+z(1,2)*z(3,2)*w(2)+
     &     z(1,3)*z(3,3)*w(3)          
      eln0(6)=z(2,1)*z(3,1)*w(1)+z(2,2)*z(3,2)*w(2)+
     &     z(2,3)*z(3,3)*w(3)          
!     
!     calculating the logarithmic mechanical strain at the
!     end of the increment
!     
      e(1,1)=emec(1)
      e(2,2)=emec(2)
      e(3,3)=emec(3)
      e(1,2)=emec(4)
      e(1,3)=emec(5)
      e(2,3)=emec(6)
      e(2,1)=emec(4)
      e(3,1)=emec(5)
      e(3,2)=emec(6)
!     
!     calculating the eigenvalues and eigenvectors
!     
      call rs(n,n,e,w,matz,z,fv1,fv2,ier)
!     
      if(ier.ne.0) then
        write(*,*) '
     &*ERROR calculating the eigenvalues/vectors in umat_abaqusnl'
        call exit(201)
      endif
!     
!     calculating the principal stretches at the end of the increment
!     
      do i=1,3
        if(w(i).lt.-0.5d0) then
          write(*,*) '*ERROR in umat_abaqusnl: negative eigenvalue'
          write(*,*) '         of the Cauchy-Green tensor;'
          call exit(201)
        else
          w(i)=dsqrt(2.d0*w(i)+1.d0)
        endif
      enddo
!     
!     calculating the invariants at the end of the increment
!     
      v1=w(1)+w(2)+w(3)
      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
      v3=w(1)*w(2)*w(3)
!     
!     calculating the right Cauchy-Green tensor at the end of the
!     increment
!     
      do i=1,3
        c(i)=2.d0*emec(i)+1.d0
      enddo
      do i=4,6
        c(i)=2.d0*emec(i)
      enddo
!     
!     calculating the square of the right Cauchy-Green tensor at the
!     end of the increment
!     
      c2(1)=c(1)*c(1)+c(4)*c(4)+c(5)*c(5)
      c2(2)=c(4)*c(4)+c(2)*c(2)+c(6)*c(6)
      c2(3)=c(5)*c(5)+c(6)*c(6)+c(3)*c(3)
      c2(4)=c(1)*c(4)+c(4)*c(2)+c(5)*c(6)
      c2(5)=c(1)*c(5)+c(4)*c(6)+c(5)*c(3)
      c2(6)=c(4)*c(5)+c(2)*c(6)+c(6)*c(3)
!     
!     calculating the right stretch tensor at the end of the increment
!     (cf. Simo and Hughes, Computational Inelasticity)
!     
!     This is the mechanical right stretch tensor, i.e. without
!     thermal expansion (since it is based on emec)
!     
      dd=v1*v2-v3
      do i=1,6
        u(i)=(-c2(i)+(v1*v1-v2)*c(i)+v1*v3*d(i))/dd
      enddo
!     
!     calculating the inverse of the right stretch tensor at the end
!     of the increment
!     
      um1(1,1)=(c(1)-v1*u(1)+v2)/v3
      um1(2,2)=(c(2)-v1*u(2)+v2)/v3
      um1(3,3)=(c(3)-v1*u(3)+v2)/v3
      um1(1,2)=(c(4)-v1*u(4))/v3
      um1(1,3)=(c(5)-v1*u(5))/v3
      um1(2,3)=(c(6)-v1*u(6))/v3
      um1(2,1)=um1(1,2)
      um1(3,1)=um1(1,3)
      um1(3,2)=um1(2,3)
!     
!     calculation of the local rotation tensor at the end of the
!     increment
!     
      do i=1,3
        do j=1,3
          r(i,j)=xkl(i,1)*um1(1,j)+xkl(i,2)*um1(2,j)+
     &         xkl(i,3)*um1(3,j)
        enddo
      enddo
!     
!     calculating the logarithmic strain at the end of the increment
!     Elog=Z.ln(w).Z^T
!     
      do i=1,3
        w(i)=dlog(w(i))
      enddo
!     
!     logarithmic strain in global coordinates at the end of the
!     increment = ln(lambda) N diadic N;
!     
!     the true logarithmic strain is ln(lambda) n diadic n;
!     the rotation tensor R = n diadic N;
!     therefore, ln(lambda) N diadic N is R^T true logarithmic strain R, 
!     in other words the corotational true logarithmic strain (required
!     by the abaqus routine)
!     
      xm1(1)=z(1,1)*z(1,1)
      xm1(2)=z(2,1)*z(2,1)
      xm1(3)=z(3,1)*z(3,1)
      xm1(4)=z(1,1)*z(2,1)
      xm1(5)=z(1,1)*z(3,1)
      xm1(6)=z(2,1)*z(3,1)
!     
      xm2(1)=z(1,2)*z(1,2)
      xm2(2)=z(2,2)*z(2,2)
      xm2(3)=z(3,2)*z(3,2)
      xm2(4)=z(1,2)*z(2,2)
      xm2(5)=z(1,2)*z(3,2)
      xm2(6)=z(2,2)*z(3,2)
!     
      xm3(1)=z(1,3)*z(1,3)
      xm3(2)=z(2,3)*z(2,3)
      xm3(3)=z(3,3)*z(3,3)
      xm3(4)=z(1,3)*z(2,3)
      xm3(5)=z(1,3)*z(3,3)
      xm3(6)=z(2,3)*z(3,3)
!     
      eln(1)=xm1(1)*w(1)+xm2(1)*w(2)+xm3(1)*w(3)          
      eln(2)=xm1(2)*w(1)+xm2(2)*w(2)+xm3(2)*w(3)          
      eln(3)=xm1(3)*w(1)+xm2(3)*w(2)+xm3(3)*w(3)          
      eln(4)=xm1(4)*w(1)+xm2(4)*w(2)+xm3(4)*w(3)          
      eln(5)=xm1(5)*w(1)+xm2(5)*w(2)+xm3(5)*w(3)          
      eln(6)=xm1(6)*w(1)+xm2(6)*w(2)+xm3(6)*w(3)  
c     write(*,*) 'iel', iel
c     write(*,*) 'emec',(emec(i),i=1,6)
c     write(*,*) 'eln',(eln(i),i=1,6)
c     write(*,*) 'r0',((r0(i,j),j=1,3),i=1,3)
c     write(*,*) 'r',((r(i,j),j=1,3),i=1,3)
!     
!     calculating the incremental rotation tensor
!     drot=r.r0^T
!     
      do i=1,3
        do j=1,3
          drot(i,j)=r(i,1)*r0(j,1)+r(i,2)*r0(j,2)+r(i,3)*r0(j,3)
        enddo
      enddo
!     
      ntens=6
!     
      do i=1,nstate_
        xstate(i,iint,iel)=xstateini(i,iint,iel)
      enddo
!     
      abqtime(1)=time-dtime
      abqtime(2)=ttime+time-dtime
!     
      temp=t1l
      dtemp=0.d0
!     
      ndi=3
      nshr=3
      ntens=ndi+nshr
!     
      nprops=-kode-100
!     
!     taking local material orientations into account
!     
      if(iorien.ne.0) then
        call transformatrix(orab(1,iorien),pgauss,skl)
!     
!     rotating the strain at the start of the increment
!     into the local system: Elog'=T^T.Elog.T
!     
        xa(1,1)=eln0(1)
        xa(1,2)=eln0(4)
        xa(1,3)=eln0(5)
        xa(2,1)=eln0(4)
        xa(2,2)=eln0(2)
        xa(2,3)=eln0(6)
        xa(3,1)=eln0(5)
        xa(3,2)=eln0(6)
        xa(3,3)=eln0(3)
!     
        do jj=1,6
          stran(jj)=0.d0
          j1=kal(1,jj)
          j2=kal(2,jj)
          do j3=1,3
            do j4=1,3
              stran(jj)=stran(jj)+
     &             xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
            enddo
          enddo
        enddo
!     
!     rotating the strain at the end of the increment
!     into the local system
!     
        xa(1,1)=eln(1)
        xa(1,2)=eln(4)
        xa(1,3)=eln(5)
        xa(2,1)=eln(4)
        xa(2,2)=eln(2)
        xa(2,3)=eln(6)
        xa(3,1)=eln(5)
        xa(3,2)=eln(6)
        xa(3,3)=eln(3)
!     
        do jj=1,6
          dstran(jj)=-stran(jj)
          j1=kal(1,jj)
          j2=kal(2,jj)
          do j3=1,3
            do j4=1,3
              dstran(jj)=dstran(jj)+
     &             xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
            enddo
          enddo
        enddo
      else
        do jj=1,6
          stran(jj)=eln0(jj)
          dstran(jj)=eln(jj)-eln0(jj)
        enddo
      endif
!     
!     rotating the stress into the local system
!     s'=J^(-1).U.S.U^T (no orientation card) or
!     s'=J^(-1).U.T^T.S.T.U^T (orientation card)
!     
!     Here, U is the total right stretch tensor, i.e. the
!     thermal effects have to be included.
!     
!     Isotropic thermal expansion is assumed. Indeed, only
!     in that case the rotation tensor with expansion equals
!     the rotation tensor without thermal expansion, and      
!     the right stretch tensor with expansion equals the right
!     stretch tensor without expansion times (1+alpha*Delta T).
!     
!     calculate the Cauchy-Green tensor at the beginning of 
!     the increment
!     
      do i=1,3
        c0tot(i,i)=xokl(1,i)*xokl(1,i)+xokl(2,i)*xokl(2,i)+
     &       xokl(3,i)*xokl(3,i)
      enddo
!     
!     C0tot=C0*(1+alpha*Delta T)**2
!     expansion:=1+alpha*Delta T
!     
      expansion=dsqrt((c0tot(1,1)/c0(1)+c0tot(2,2)/c0(2)+
     &     c0tot(3,3)/c0(3))/3.d0)
!     
!     calculate the total right stretch tensor at the start
!     of the increment (= mechanical right stretch tensor
!     multiplied with the expansion at the start of the increment)
!     
      do i=1,3
        do j=1,3
          u0(i,j)=u0(i,j)*expansion
        enddo
      enddo
!     
      if(iorien.ne.0) then
        do i=1,3
          do j=1,3
            tkl(i,j)=u0(i,1)*skl(j,1)+u0(i,2)*skl(j,2)+
     &           u0(i,3)*skl(j,3)
          enddo
        enddo
      else
        do i=1,3
          do j=1,3
            tkl(i,j)=u0(i,j)
          enddo
        enddo
      endif
!     
      xa(1,1)=stre(1)
      xa(1,2)=stre(4)
      xa(1,3)=stre(5)
      xa(2,1)=stre(4)
      xa(2,2)=stre(2)
      xa(2,3)=stre(6)
      xa(3,1)=stre(5)
      xa(3,2)=stre(6)
      xa(3,3)=stre(3)
!     
      do jj=1,6
        stre(jj)=0.d0
        j1=kal(1,jj)
        j2=kal(2,jj)
        do j3=1,3
          do j4=1,3
            stre(jj)=stre(jj)+
     &           xa(j3,j4)*tkl(j1,j3)*tkl(j2,j4)
          enddo
        enddo
        stre(jj)=stre(jj)/voj
      enddo
!     
!     changing physical strain into engineering strain
!     ABAQUS uses the engineering strain!
!     
      do i=4,6
        stran(i)=2.d0*stran(i)
        dstran(i)=2.d0*dstran(i)
      enddo
!     
      if(amat(1:1).eq.'@') then
!     
        call call_external_umat(stre,xstate(1,iint,iel),ddsdde,
     &       sse,spd,scd,rpl,ddsddt,drplde,drpldt,stran,dstran,
     &       abqtime,dtime,temp,dtemp ,predef,dpred,amat,ndi,nshr,
     &       ntens,nstate_,elconloc,nprops,pgauss,drot,pnewdt,
     &       celent,xokl,xkl,iel,iint,layer,kspt,jstep,kinc)
!     
      elseif(amat(1:11).eq.'JOHNSONCOOK') then
!     
        call umat_johnson_cook(stre,xstate(1,iint,iel),ddsdde,sse,
     &       spd,scd,rpl,ddsddt,
     &       drplde,drpldt,stran,dstran,abqtime,dtime,temp,dtemp,
     &       predef,dpred,amat,ndi,nshr,ntens,nstate_,elconloc,nprops,
     &       pgauss,drot,pnewdt,celent,xokl,xkl,iel,iint,layer,kspt,
     &       jstep,kinc)
!     
      else
!     
        call umat(stre,xstate(1,iint,iel),ddsdde,sse,spd,scd,rpl,ddsddt
     &       ,drplde,drpldt,stran,dstran,abqtime,dtime,temp,dtemp
     &       ,predef,dpred,amat,ndi,nshr,ntens,nstate_,elconloc,nprops
     &       ,pgauss,drot,pnewdt,celent,xokl,xkl,iel,iint,layer,kspt
     &       ,jstep,kinc)
!     
      endif
!     
!     rotating the stress into the global system
!     S=J.U^(-1).s'.U^(-T)  (no orientation card) or
!     S=J.T.U^(-1).s'.U^(-T).T^T (orientation card)
!     
!     
!     Here, U is the total right stretch tensor, i.e. the
!     thermal effects have to be included.
!     
!     Isotropic thermal expansion is assumed. Indeed, only
!     in that case the rotation tensor with expansion equals
!     the rotation tensor without thermal expansion, and      
!     the right stretch tensor with expansion equals the right
!     stretch tensor without expansion times (1+alpha*Delta T).
!     
!     calculate the Cauchy-Green tensor at the beginning of 
!     the increment
!     
      do i=1,3
        ctot(i,i)=xkl(1,i)*xkl(1,i)+xkl(2,i)*xkl(2,i)+
     &       xkl(3,i)*xkl(3,i)
      enddo
!     
!     Ctot=C*(1+alpha*Delta T)**2
!     expansion:=1+alpha*Delta T
!     
      expansion=dsqrt((ctot(1,1)/c(1)+ctot(2,2)/c(2)+
     &     ctot(3,3)/c(3))/3.d0)
!     
c!     calculate the inverse total right stretch tensor at the end
c!     of the increment (= inverse mechanical right stretch tensor
c!     divided by the expansion at the end of the increment)
c!     
c      do i=1,3
c        do j=1,3
c          um1(i,j)=um1(i,j)/expansion
c        enddo
c      enddo
!     
      if(iorien.ne.0) then
        do i=1,3
          do j=1,3
            tkl(i,j)=skl(i,1)*um1(1,j)+skl(i,2)*um1(2,j)+
     &           skl(i,3)*um1(3,j)
          enddo
        enddo
      else
        do i=1,3
          do j=1,3
            tkl(i,j)=um1(i,j)
          enddo
        enddo
      endif
!     
      xa(1,1)=stre(1)
      xa(1,2)=stre(4)
      xa(1,3)=stre(5)
      xa(2,1)=stre(4)
      xa(2,2)=stre(2)
      xa(2,3)=stre(6)
      xa(3,1)=stre(5)
      xa(3,2)=stre(6)
      xa(3,3)=stre(3)
!     
      do jj=1,6
        stre(jj)=0.d0
        j1=kal(1,jj)
        j2=kal(2,jj)
        do j3=1,3
          do j4=1,3
            stre(jj)=stre(jj)+
     &           xa(j3,j4)*tkl(j1,j3)*tkl(j2,j4)
          enddo
        enddo
        stre(jj)=stre(jj)*vj/expansion**2
      enddo
!     
!     calculate the stiffness matrix (the matrix is symmetrized)
!     
      if(icmd.ne.3) then
        stiff(1)=ddsdde(1,1)
        stiff(2)=(ddsdde(1,2)+ddsdde(2,1))/2.d0
        stiff(3)=ddsdde(2,2)
        stiff(4)=(ddsdde(1,3)+ddsdde(3,1))/2.d0
        stiff(5)=(ddsdde(2,3)+ddsdde(3,2))/2.d0
        stiff(6)=ddsdde(3,3)
        stiff(7)=(ddsdde(1,4)+ddsdde(4,1))/2.d0
        stiff(8)=(ddsdde(2,4)+ddsdde(4,2))/2.d0
        stiff(9)=(ddsdde(3,4)+ddsdde(4,3))/2.d0
        stiff(10)=ddsdde(4,4)
        stiff(11)=(ddsdde(1,5)+ddsdde(5,1))/2.d0
        stiff(12)=(ddsdde(2,5)+ddsdde(5,2))/2.d0
        stiff(13)=(ddsdde(3,5)+ddsdde(5,3))/2.d0
        stiff(14)=(ddsdde(4,5)+ddsdde(5,4))/2.d0
        stiff(15)=ddsdde(5,5)
        stiff(16)=(ddsdde(1,6)+ddsdde(6,1))/2.d0
        stiff(17)=(ddsdde(2,6)+ddsdde(6,2))/2.d0
        stiff(18)=(ddsdde(3,6)+ddsdde(6,3))/2.d0
        stiff(19)=(ddsdde(4,6)+ddsdde(6,4))/2.d0
        stiff(20)=(ddsdde(5,6)+ddsdde(6,5))/2.d0
        stiff(21)=ddsdde(6,6)
!     
!     rotating the stiffness coefficients into the global system
!     
        call anisotropic(stiff,ya)
!     
        do jj=1,21
          j1=kel(1,jj)
          j2=kel(2,jj)
          j3=kel(3,jj)
          j4=kel(4,jj)
          stiff(jj)=0.d0
          do j5=1,3
            do j6=1,3
              do j7=1,3
                do j8=1,3
                  stiff(jj)=stiff(jj)+ya(j5,j6,j7,j8)*
     &                 tkl(j1,j5)*tkl(j2,j6)*tkl(j3,j7)*tkl(j4,j8)
                enddo
              enddo
            enddo
          enddo
          stiff(jj)=stiff(jj)*vj/expansion**2
        enddo
      endif
!     
      return
      end
