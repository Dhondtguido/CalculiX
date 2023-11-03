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
      subroutine umat_abaqusnl_total(amat,iel,iint,kode,elconloc,emec,
     &     emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &     icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &     iorien,pgauss,orab,istep,kinc,pnewdt,nmethod,iperturb,
     &     plconloc,depvisc)
!
!     converting nonlinear fields used by CalculiX into nonlinear
!     fields used by Abaqus and vice versa.
!
!     before the umat call: convert the Lagrange strain into
!                           corotational logarithmic strain
!     after the umat call: convert corotational Cauchy stress
!                          into PK2 stress and convert the derivative
!                          of the corotational Cauchy stress w.r.t. the
!                          corotational logarithmic strain into the
!                          derivative of the PK2 stress w.r.t. the
!                          Lagrange strain
! 
!     the appendix "total" means that in this routine only the mechanical
!     strain at the end of the increment is converted from Lagrange to
!     logarithmic. The mechanical strain at the start of the increment is 
!     left unchanged and so is the stress. So the constitutive equation 
!     should be based solely on the mechanical strain at the end of the
!     increment, on not on its change within the increment.
! 
!     emec is the mechanical Lagrange strain (excluding thermal effects)
!     xkl is the total deformation gradient (including thermal effects)
!     vj is the total Jacobian (including thermal effects)
! 
!     Reference for the derivative of a function of C w.r.t C (C is
!     the Cauchy-Green tensor):
! 
!     Dhondt, G., Hackenberg, H.-P., Use of a rotation-invariant linear
!     strain measure for linear elastic crack propagation calculations,
!     Eng. Fract. Mech. 247(2021) 107634.
! 
      implicit none
!     
      character*80 amat,amatloc
!     
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),i,
     &     iorien,nmethod,iperturb(*),istep,nprops,jstep(4),kinc,
     &     kel(4,21),j1,j2,j3,j4,j5,j6,j7,j8,jj,n,ier,j,matz,kal(2,6),
     &     mel(4,36),jm,jn,js,jt
!     
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &     vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &     time,ttime,skl(3,3),scc(3,3),ya(3,3,3,3),z(3,3),depvisc,
     &     xstateini(nstate_,mi(1),*),w(3),fv1(3),fv2(3),d(6),c(3,3),
     &     v1,v2,v3,eln(6),e(3,3),tkl(3,3),u(6),c2(6),dd,um1(3,3),
     &     expansion,ctot(3,3),ddsdde(6,6),spd,rpl,pnewdt,stran(6),temp,
     &     xstate(nstate_,mi(1),*),xm1(6),xm2(6),plconloc(802),
     &     xm3(6),wum1(3),weln(3),um1new(3,3),xa(3,3),ym1(3,3),
     &     ym2(3,3),ym3(3,3),sccum1(3,3),cm1(3,3),um1sccum1(3,3),
     &     c4(21),d4(21),dxm1dc(21),dxm2dc(21),dxm3dc(21),
     &     dum1dc(21),delndc(21),dum1dcexp(3,3,3,3),delndcexp(3,3,3,3),
     &     stiffexp(3,3,3,3),term1(36),term24(36),term3(36),
     &     stiffasym(36),aa(21),bb(21),dwum1(3),dweln(3),ym14(21),
     &     ym24(21),ym34(21)
!     
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!     
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!     
      mel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &     1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &     3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &     1,2,2,3,1,3,2,3,2,3,2,3,2,2,1,1,3,3,1,1,3,3,2,2,
     &     1,2,1,1,1,2,2,2,1,2,3,3,1,3,1,1,1,3,2,2,1,3,3,3,
     &     1,3,1,2,2,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3/),(/4,36/))
!     
      d=(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
!     
!     fourth order identity tensor     
!     
      d4=(/1.d0,0.d0,1.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.5d0,0.d0,0.d0,
     &0.d0,0.d0,0.5d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.5d0/)
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
c      write(*,*) 'umat_abaqusnl_total ',(emec(i),i=1,6)
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
!     calculating the - principal stretches square at the end of the increment
!                       = eigenvalues of the right Cauchy-Green tensor
!                     - eigenvalues of the logarithmic strain
!                     - eigenvalues of the inverse right stretch tensor
!     
      do i=1,3
        if(w(i).lt.-0.5d0) then
          write(*,*) '*ERROR in umat_abaqusnl: negative eigenvalue'
          write(*,*) '         of the Cauchy-Green tensor;'
          call exit(201)
        else
          w(i)=2.d0*w(i)+1.d0
          weln(i)=dlog(w(i))/2.d0
          wum1(i)=1.d0/dsqrt(w(i))
        endif
      enddo
!     
!     calculating the mechanical right Cauchy-Green tensor at the end of the
!     increment
!     
      do i=1,3
        c(i,i)=2.d0*e(i,i)+1.d0
      enddo
      c(1,2)=2.d0*e(1,2)
      c(1,3)=2.d0*e(1,3)
      c(2,3)=2.d0*e(2,3)
      c(2,1)=c(1,2)
      c(3,1)=c(1,3)
      c(3,2)=c(2,3)
!     
!     calculating the structural tensors    
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
!     corotational mechanical logarithmic strain in global coordinates at the
!     end of the increment = ln(lambda) N diadic N;
!     
!     the true logarithmic strain is ln(lambda) n diadic n;
!     the rotation tensor R = n diadic N;
!     therefore, ln(lambda) N diadic N is R^T true logarithmic strain R, 
!     in other words the corotational true logarithmic strain (required
!     by the abaqus routine)
!     
      eln(1)=xm1(1)*weln(1)+xm2(1)*weln(2)+xm3(1)*weln(3)          
      eln(2)=xm1(2)*weln(1)+xm2(2)*weln(2)+xm3(2)*weln(3)          
      eln(3)=xm1(3)*weln(1)+xm2(3)*weln(2)+xm3(3)*weln(3)          
      eln(4)=xm1(4)*weln(1)+xm2(4)*weln(2)+xm3(4)*weln(3)          
      eln(5)=xm1(5)*weln(1)+xm2(5)*weln(2)+xm3(5)*weln(3)          
      eln(6)=xm1(6)*weln(1)+xm2(6)*weln(2)+xm3(6)*weln(3)  
!     
!     calculating the mechanical inverse right stretch tensor U^{-1}
!     
      um1(1,1)=xm1(1)*wum1(1)+xm2(1)*wum1(2)+xm3(1)*wum1(3)          
      um1(2,2)=xm1(2)*wum1(1)+xm2(2)*wum1(2)+xm3(2)*wum1(3)          
      um1(3,3)=xm1(3)*wum1(1)+xm2(3)*wum1(2)+xm3(3)*wum1(3)          
      um1(1,2)=xm1(4)*wum1(1)+xm2(4)*wum1(2)+xm3(4)*wum1(3)          
      um1(1,3)=xm1(5)*wum1(1)+xm2(5)*wum1(2)+xm3(5)*wum1(3)          
      um1(2,3)=xm1(6)*wum1(1)+xm2(6)*wum1(2)+xm3(6)*wum1(3)
      um1(2,1)=um1(1,2)
      um1(3,1)=um1(1,3)
      um1(3,2)=um1(2,3)
!     
      do i=1,nstate_
        xstate(i,iint,iel)=xstateini(i,iint,iel)
      enddo
!     
      temp=t1l
!     
      nprops=-kode-100
!     
!     taking local material orientations into account
!     
      if(iorien.ne.0) then
        call transformatrix(orab(1,iorien),pgauss,skl)
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
      else
        do jj=1,6
          stran(jj)=eln(jj)
        enddo
      endif
!     
      if(amat(1:11).eq.'JOHNSONCOOK') then
!     
!     changing physical strain into engineering strain
!     ABAQUS uses the engineering strain!
!     
        do i=4,6
          stran(i)=2.d0*stran(i)
        enddo
!     
        call umat_johnson_cook(stre,xstate(1,iint,iel),ddsdde,
     &       spd,rpl,stran,dtime,temp,elconloc,nprops,
     &       iel,iint,jstep,kinc,icmd)
        if(icmd.ne.3) then
!     
          stiff(1)=ddsdde(1,1)
          stiff(2)=ddsdde(2,1)
          stiff(3)=ddsdde(2,2)
          stiff(4)=ddsdde(3,1)
          stiff(5)=ddsdde(3,2)
          stiff(6)=ddsdde(3,3)
          stiff(7)=ddsdde(4,1)
          stiff(8)=ddsdde(4,2)
          stiff(9)=ddsdde(4,3)
          stiff(10)=ddsdde(4,4)
          stiff(11)=ddsdde(5,1)
          stiff(12)=ddsdde(5,2)
          stiff(13)=ddsdde(5,3)
          stiff(14)=ddsdde(5,4)
          stiff(15)=ddsdde(5,5)
          stiff(16)=ddsdde(6,1)
          stiff(17)=ddsdde(6,2)
          stiff(18)=ddsdde(6,3)
          stiff(19)=ddsdde(6,4)
          stiff(20)=ddsdde(6,5)
          stiff(21)=ddsdde(6,6)
        endif
!     
      elseif(amat(1:11).eq.'ANISO_CREEP') then
!     
!       orthotropic elasticity with isotropic creep defined by
!       a user subroutine     
!     
        amatloc(1:69)=amat(12:80)
        amatloc(70:80)='           '
        call umat_aniso_creep(amatloc,
     &       iel,iint,kode,elconloc,stran,emec0,
     &       beta,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &       icmd,ielas,mi(1),nstate_,xstateini,xstate,stre,stiff,
     &       iorien,pgauss,orab,nmethod,pnewdt,depvisc)
!          
      elseif(kode.eq.-50) then
!     
!       deformation plasticity
!     
        call umat_def_plas(elconloc,stiff,stran,icmd,stre,
     &       xstate(1,iint,iel),iel,iint)
!     
      elseif(kode.eq.-53) then
!     
!       Mohr-Coulomb plasticity
!     
        call mohrcoulomb(elconloc,plconloc,xstate,xstateini,
     &       stiff,stran,icmd,beta,stre,
     &       ielas,dtime,time,ttime,iel,iint,nstate_,mi,pnewdt)
!     
      elseif(kode.eq.-54) then
!     
!       orthotropic elasticity with isotropic plasticity
!     
        call ortho_plas(amat,iel,iint,kode,elconloc,stran,
     &       emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &       icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,iorien,
     &       pgauss,orab,nmethod,pnewdt,plconloc)
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
!     calculate the total Cauchy-Green tensor (based on the
!     total deformation gradient, i.e. including thermal effects     
!     
      do i=1,3
        do j=1,3
          ctot(i,j)=xkl(1,i)*xkl(1,i)+xkl(2,i)*xkl(2,i)+
     &         xkl(3,i)*xkl(3,i)
        enddo
      enddo
!     
!     Ctot=C*(1+alpha*Delta T)**2
!     expansion:=1+alpha*Delta T
!     
      expansion=dsqrt((ctot(1,1)/c(1,1)+ctot(2,2)/c(2,2)+
     &     ctot(3,3)/c(3,3))/3.d0)
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
!     scc is the corotational Cauchy stress as matrix    
!     
      scc(1,1)=stre(1)
      scc(1,2)=stre(4)
      scc(1,3)=stre(5)
      scc(2,1)=stre(4)
      scc(2,2)=stre(2)
      scc(2,3)=stre(6)
      scc(3,1)=stre(5)
      scc(3,2)=stre(6)
      scc(3,3)=stre(3)
!     
      do jj=1,6
        stre(jj)=0.d0
        j1=kal(1,jj)
        j2=kal(2,jj)
        do j3=1,3
          do j4=1,3
            stre(jj)=stre(jj)+
     &           scc(j3,j4)*tkl(j1,j3)*tkl(j2,j4)
          enddo
        enddo
        stre(jj)=stre(jj)*vj/expansion**2
      enddo
!     
!     calculate the stiffness matrix (the matrix is symmetrized)
!     
      if(icmd.ne.3) then
c!     
c!     rotating the stiffness coefficients into the global system
c!     
c        call anisotropic(stiff,ya)
c!     
c        do jj=1,21
c          j1=kel(1,jj)
c          j2=kel(2,jj)
c          j3=kel(3,jj)
c          j4=kel(4,jj)
c          stiff(jj)=0.d0
c          do j5=1,3
c            do j6=1,3
c              do j7=1,3
c                do j8=1,3
c                  stiff(jj)=stiff(jj)+ya(j5,j6,j7,j8)*
c     &                 tkl(j1,j5)*tkl(j2,j6)*tkl(j3,j7)*tkl(j4,j8)
c                enddo
c              enddo
c            enddo
c          enddo
c          stiff(jj)=stiff(jj)*vj/expansion**2
c        enddo
!     
!       sccum1=scc*um1
!
        do i=1,3
          do j=1,3
            sccum1(i,j)=scc(i,1)*um1(1,j)+scc(i,2)*um1(2,j)+
     &           scc(i,3)*um1(3,j)
          enddo
        enddo
!     
!       cm1=um1*um1 (symmetric): inverse mechanical right Cauchy-Green
!       tensor
!
        do i=1,3
          do j=i,3
            cm1(i,j)=um1(i,1)*um1(1,j)+um1(i,2)*um1(2,j)+
     &           um1(i,3)*um1(3,j)
          enddo
        enddo
        cm1(2,1)=cm1(1,2)
        cm1(3,1)=cm1(1,3)
        cm1(3,2)=cm1(2,3)
!     
!       um1sccum1=um1*sccum1 (symmetric)
!
        do i=1,3
          do j=i,3
            um1sccum1(i,j)=um1(i,1)*sccum1(1,j)+um1(i,2)*sccum1(2,j)+
     &           um1(i,3)*sccum1(3,j)
          enddo
        enddo
        um1sccum1(2,1)=um1sccum1(1,2)
        um1sccum1(3,1)=um1sccum1(1,3)
        um1sccum1(3,2)=um1sccum1(2,3)
!     
!       ym1, ym2 and ym3 are the structural tensors as matrices
!     
        ym1(1,1)=xm1(1)
        ym1(1,2)=xm1(4)
        ym1(1,3)=xm1(5)
        ym1(2,1)=xm1(4)
        ym1(2,2)=xm1(2)
        ym1(2,3)=xm1(6)
        ym1(3,1)=xm1(5)
        ym1(3,2)=xm1(6)
        ym1(3,3)=xm1(3)
!     
        ym2(1,1)=xm2(1)
        ym2(1,2)=xm2(4)
        ym2(1,3)=xm2(5)
        ym2(2,1)=xm2(4)
        ym2(2,2)=xm2(2)
        ym2(2,3)=xm2(6)
        ym2(3,1)=xm2(5)
        ym2(3,2)=xm2(6)
        ym2(3,3)=xm2(3)
!     
        ym3(1,1)=xm3(1)
        ym3(1,2)=xm3(4)
        ym3(1,3)=xm3(5)
        ym3(2,1)=xm3(4)
        ym3(2,2)=xm3(2)
        ym3(2,3)=xm3(6)
        ym3(3,1)=xm3(5)
        ym3(3,2)=xm3(6)
        ym3(3,3)=xm3(3)
!
!       calculating the derivative of the eigenvalue functions
!
        do i=1,3
          dweln(i)=1.d0/(2.d0*w(i))
          dwum1(i)=-wum1(i)/(2.d0*w(i))
        enddo
!
!       calculating A, B, dM_i/dC, dU^{-1}/dC and dln(e)/dC
!       (all of them symmetric, i.e. 21 constants)
!
          write(*,*) 'eigenvalues ',w(1),w(2),w(3)
        if((dabs(w(1)-w(2)).lt.1.d-10).and.
     &       (dabs(w(2)-w(3)).lt.1.d-10)) then
!
!         three equal eigenvalues
!     
        do jj=1,21
            dum1dc(jj)=d4(jj)*dwum1(1)
            delndc(jj)=d4(jj)*dweln(1)
          enddo
        elseif(dabs(w(1)-w(2)).lt.1.d-10) then
          write(*,*) 'entering the new branch'
!
!         first and second eigenvalue are equal    
!         (M3 o I + I o M3)'    
!         note: o is the dyadic product     
!     
          c4(1)=2.d0*ym3(1,1)
          c4(2)=0.d0
          c4(3)=2.d0*ym3(2,2)
          c4(4)=0.d0
          c4(5)=0.d0
          c4(6)=2.d0*ym3(3,3)
          c4(7)=ym3(2,1)
          c4(8)=ym3(1,2)
          c4(9)=0.d0
          c4(10)=(ym3(2,2)+ym3(1,1))/2.d0
          c4(11)=ym3(3,1)
          c4(12)=0.d0
          c4(13)=ym3(1,3)
          c4(14)=ym3(3,2)/2.d0
          c4(15)=(ym3(3,3)+ym3(1,1))/2.d0
          c4(16)=0.d0
          c4(17)=ym3(3,2)
          c4(18)=ym3(2,3)
          c4(19)=ym3(3,1)/2.d0
          c4(20)=ym3(2,1)/2.d0
          c4(21)=(ym3(3,3)+ym3(2,2))/2.d0
!
          do jj=1,21
            j1=kel(1,jj)
            j2=kel(2,jj)
            j3=kel(3,jj)
            j4=kel(4,jj)
            ym34(jj)=ym3(j1,j2)*ym3(j3,j4)
            dum1dc(jj)=dwum1(3)*ym34(jj)+
     &           (wum1(3)-wum1(1))/(w(3)-w(1))*(c4(jj)-2.d0*ym34(jj))+
     &           dwum1(1)*(d4(jj)+ym34(jj)-c4(jj))
            delndc(jj)=dweln(3)*ym34(jj)+
     &           (weln(3)-weln(1))/(w(3)-w(1))*(c4(jj)-2.d0*ym34(jj))+
     &           dweln(1)*(d4(jj)+ym34(jj)-c4(jj))
          enddo
         elseif(dabs(w(2)-w(3)).lt.1.d-10) then
c          write(*,*) 'entering the new branch'
!
!         second and third eigenvalue are equal    
!         (M1 o I + I o M1)'    
!         note: o is the dyadic product     
!     
          c4(1)=2.d0*ym1(1,1)
          c4(2)=0.d0
          c4(3)=2.d0*ym1(2,2)
          c4(4)=0.d0
          c4(5)=0.d0
          c4(6)=2.d0*ym1(3,3)
          c4(7)=ym1(2,1)
          c4(8)=ym1(1,2)
          c4(9)=0.d0
          c4(10)=(ym1(2,2)+ym1(1,1))/2.d0
          c4(11)=ym1(3,1)
          c4(12)=0.d0
          c4(13)=ym1(1,3)
          c4(14)=ym1(3,2)/2.d0
          c4(15)=(ym1(3,3)+ym1(1,1))/2.d0
          c4(16)=0.d0
          c4(17)=ym1(3,2)
          c4(18)=ym1(2,3)
          c4(19)=ym1(3,1)/2.d0
          c4(20)=ym1(2,1)/2.d0
          c4(21)=(ym1(3,3)+ym1(2,2))/2.d0
!
          do jj=1,21
            j1=kel(1,jj)
            j2=kel(2,jj)
            j3=kel(3,jj)
            j4=kel(4,jj)
            ym14(jj)=ym1(j1,j2)*ym1(j3,j4)
            dum1dc(jj)=dwum1(1)*ym14(jj)+
     &           (wum1(1)-wum1(2))/(w(1)-w(2))*(c4(jj)-2.d0*ym14(jj))+
     &           dwum1(2)*(d4(jj)+ym14(jj)-c4(jj))
            delndc(jj)=dweln(1)*ym14(jj)+
     &           (weln(1)-weln(2))/(w(1)-w(2))*(c4(jj)-2.d0*ym14(jj))+
     &           dweln(2)*(d4(jj)+ym14(jj)-c4(jj))
          enddo
        elseif(dabs(w(3)-w(1)).lt.1.d-10) then
c          write(*,*) 'entering the new branch'
!
!         third and first eigenvalue are equal    
!         (M2 o I + I o M2)'    
!         note: o is the dyadic product     
!     
          c4(1)=2.d0*ym2(1,1)
          c4(2)=0.d0
          c4(3)=2.d0*ym2(2,2)
          c4(4)=0.d0
          c4(5)=0.d0
          c4(6)=2.d0*ym2(3,3)
          c4(7)=ym2(2,1)
          c4(8)=ym2(1,2)
          c4(9)=0.d0
          c4(10)=(ym2(2,2)+ym2(1,1))/2.d0
          c4(11)=ym2(3,1)
          c4(12)=0.d0
          c4(13)=ym2(1,3)
          c4(14)=ym2(3,2)/2.d0
          c4(15)=(ym2(3,3)+ym2(1,1))/2.d0
          c4(16)=0.d0
          c4(17)=ym2(3,2)
          c4(18)=ym2(2,3)
          c4(19)=ym2(3,1)/2.d0
          c4(20)=ym2(2,1)/2.d0
          c4(21)=(ym2(3,3)+ym2(2,2))/2.d0
!
          do jj=1,21
            j1=kel(1,jj)
            j2=kel(2,jj)
            j3=kel(3,jj)
            j4=kel(4,jj)
            ym24(jj)=ym2(j1,j2)*ym2(j3,j4)
            dum1dc(jj)=dwum1(2)*ym24(jj)+
     &           (wum1(2)-wum1(3))/(w(2)-w(3))*(c4(jj)-2.d0*ym24(jj))+
     &           dwum1(3)*(d4(jj)+ym24(jj)-c4(jj))
            delndc(jj)=dweln(2)*ym24(jj)+
     &           (weln(2)-weln(3))/(w(2)-w(3))*(c4(jj)-2.d0*ym24(jj))+
     &           dweln(3)*(d4(jj)+ym24(jj)-c4(jj))
          enddo
       else
!
!         derivative of c*c w.r.t. c (symmetric 4th order tensor)     
!     
          c4(1)=2.d0*c(1,1)
          c4(2)=0.d0
          c4(3)=2.d0*c(2,2)
          c4(4)=0.d0
          c4(5)=0.d0
          c4(6)=2.d0*c(3,3)
          c4(7)=c(2,1)
          c4(8)=c(1,2)
          c4(9)=0.d0
          c4(10)=(c(2,2)+c(1,1))/2.d0
          c4(11)=c(3,1)
          c4(12)=0.d0
          c4(13)=c(1,3)
          c4(14)=c(3,2)/2.d0
          c4(15)=(c(3,3)+c(1,1))/2.d0
          c4(16)=0.d0
          c4(17)=c(3,2)
          c4(18)=c(2,3)
          c4(19)=c(3,1)/2.d0
          c4(20)=c(2,1)/2.d0
          c4(21)=(c(3,3)+c(2,2))/2.d0
!
          do jj=1,21
            j1=kel(1,jj)
            j2=kel(2,jj)
            j3=kel(3,jj)
            j4=kel(4,jj)
            ym14(jj)=ym1(j1,j2)*ym1(j3,j4)
            ym24(jj)=ym2(j1,j2)*ym2(j3,j4)
            ym34(jj)=ym3(j1,j2)*ym3(j3,j4)
            aa(jj)=d4(jj)-ym14(jj)-ym24(jj)-ym34(jj)
            bb(jj)=c4(jj)-2.d0*(w(1)*ym14(jj)
     &           +w(2)*ym24(jj)
     &           +w(3)*ym34(jj))
            dxm1dc(jj)=(bb(jj)-(w(2)+w(3))*aa(jj))/
     &           ((w(1)-w(2))*(w(1)-w(3)))
            dxm2dc(jj)=(bb(jj)-(w(1)+w(3))*aa(jj))/
     &           ((w(2)-w(1))*(w(2)-w(3)))
            dxm3dc(jj)=(bb(jj)-(w(1)+w(2))*aa(jj))/
     &           ((w(3)-w(1))*(w(3)-w(2)))
            dum1dc(jj)=ym14(jj)*dwum1(1)+wum1(1)*dxm1dc(jj)
     &           +ym24(jj)*dwum1(2)+wum1(2)*dxm2dc(jj)
     &           +ym34(jj)*dwum1(3)+wum1(3)*dxm3dc(jj)
            delndc(jj)=ym14(jj)*dweln(1)+weln(1)*dxm1dc(jj)
     &           +ym24(jj)*dweln(2)+weln(2)*dxm2dc(jj)
     &           +ym34(jj)*dweln(3)+weln(3)*dxm3dc(jj)
          enddo
        endif
!     
!        expanding dum1dc(21), delndc(21) and stiff(21) into     
!        dum1dcexp(3,3,3,3), delndcexp(3,3,3,3) and stiffexp(3,3,3,3) 
!     
         call anisotropic(dum1dc,dum1dcexp)
         call anisotropic(delndc,delndcexp)
         call anisotropic(stiff,stiffexp)
!     
!        collecting all terms in the final expression (possibly     
!        asymmetric, i.e. 36 terms)
!
         do jj=1,36
           j1=mel(1,jj)
           j2=mel(2,jj)
           j3=mel(3,jj)
           j4=mel(4,jj)
           term1(jj)=um1sccum1(j1,j2)*cm1(j3,j4)
           term24(jj)=dum1dcexp(j1,1,j3,j4)*sccum1(1,j2)+
     &                dum1dcexp(j1,2,j3,j4)*sccum1(2,j2)+
     &                dum1dcexp(j1,3,j3,j4)*sccum1(3,j2)+
     &                dum1dcexp(j2,1,j3,j4)*sccum1(1,j1)+
     &                dum1dcexp(j2,2,j3,j4)*sccum1(2,j1)+
     &                dum1dcexp(j2,3,j3,j4)*sccum1(3,j1)
           term3(jj)=0.d0
           do jm=1,3
             do jn=1,3
               do js=1,3
                 do jt=1,3
                   term3(jj)=term3(jj)+um1(j1,jm)*stiffexp(jm,jn,js,jt)*
     &                       delndcexp(js,jt,j3,j4)*um1(jn,j2)
                 enddo
               enddo
             enddo
           enddo
           stiffasym(jj)=vj*(term1(jj)+2.d0*(term24(jj)+term3(jj)))/
     &          expansion**2
         enddo
!
!        symmetrizing stiffasym -> stiff
!
         stiff(1)=stiffasym(1)
         stiff(2)=(stiffasym(2)+stiffasym(22))/2.d0
         stiff(3)=stiffasym(3)
         stiff(4)=(stiffasym(4)+stiffasym(23))/2.d0
         stiff(5)=(stiffasym(5)+stiffasym(24))/2.d0
         stiff(6)=stiffasym(6)
         stiff(7)=(stiffasym(7)+stiffasym(25))/2.d0
         stiff(8)=(stiffasym(8)+stiffasym(26))/2.d0
         stiff(9)=(stiffasym(9)+stiffasym(27))/2.d0
         stiff(10)=stiffasym(10)
         stiff(11)=(stiffasym(11)+stiffasym(28))/2.d0
         stiff(12)=(stiffasym(12)+stiffasym(29))/2.d0
         stiff(13)=(stiffasym(13)+stiffasym(30))/2.d0
         stiff(14)=(stiffasym(14)+stiffasym(31))/2.d0
         stiff(15)=stiffasym(15)
         stiff(16)=(stiffasym(16)+stiffasym(32))/2.d0
         stiff(17)=(stiffasym(17)+stiffasym(33))/2.d0
         stiff(18)=(stiffasym(18)+stiffasym(34))/2.d0
         stiff(19)=(stiffasym(19)+stiffasym(35))/2.d0
         stiff(20)=(stiffasym(20)+stiffasym(36))/2.d0
         stiff(21)=stiffasym(21)
c
c         write(*,*) 'umat_abaqusnl_total'
c         do i=1,21
cc     write(*,*) i,stiff(i),stiffasym(i)
c           write(*,'(i5,5(1x,e11.4))') i,stiff(i),
c     &          dabs(stiff(i)-stiffasym(i))/dabs(stiff(i)),
c     &          vj*term1(i)/expansion**2,
c     &          vj*2.d0*term24(i)/expansion**2,
c     &          vj*2.d0*term3(i)/expansion**2
c         enddo
!     
       endif
!     
      return
      end
