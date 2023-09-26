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
      implicit none
!     
      character*80 amat,amatloc
!     
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),i,
     &     iorien,nmethod,iperturb(*),istep,nprops,jstep(4),kinc,
     &     kel(4,21),j1,j2,j3,j4,j5,j6,j7,j8,jj,n,ier,j,matz,kal(2,6)
!     
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &     vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &     time,ttime,skl(3,3),xa(3,3),ya(3,3,3,3),z(3,3),depvisc,
     &     xstateini(nstate_,mi(1),*),w(3),fv1(3),fv2(3),d(6),c(6),
     &     v1,v2,v3,eln(6),e(3,3),tkl(3,3),u(6),c2(6),dd,um1(3,3),
     &     expansion,ctot(3,3),ddsdde(6,6),spd,rpl,pnewdt,stran(6),temp,
     &     xstate(nstate_,mi(1),*),xm1(6),xm2(6),plconloc(802),
     &     xm3(6),wum1(3),weln(3),um1new(3,3)
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
!     calculating the principal stretches square at the end of the increment
!                     eigenvalues of the logarithmic strain
!                     eigenvalues of the inverse right stretch tensor
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
c!     
c!     calculating the invariants at the end of the increment
c!     
c      v1=w(1)+w(2)+w(3)
c      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
c      v3=w(1)*w(2)*w(3)
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
c!     
c!     calculating the square of the right Cauchy-Green tensor at the
c!     end of the increment
c!     
c      c2(1)=c(1)*c(1)+c(4)*c(4)+c(5)*c(5)
c      c2(2)=c(4)*c(4)+c(2)*c(2)+c(6)*c(6)
c      c2(3)=c(5)*c(5)+c(6)*c(6)+c(3)*c(3)
c      c2(4)=c(1)*c(4)+c(4)*c(2)+c(5)*c(6)
c      c2(5)=c(1)*c(5)+c(4)*c(6)+c(5)*c(3)
c      c2(6)=c(4)*c(5)+c(2)*c(6)+c(6)*c(3)
c!     
c!     calculating the right stretch tensor at the end of the increment
c!     (cf. Simo and Hughes, Computational Inelasticity)
c!     
c!     This is the mechanical right stretch tensor, i.e. without
c!     thermal expansion (since it is based on emec)
c!     
c      dd=v1*v2-v3
c      do i=1,6
c        u(i)=(-c2(i)+(v1*v1-v2)*c(i)+v1*v3*d(i))/dd
c      enddo
c!     
c!     calculating the inverse of the right stretch tensor at the end
c!     of the increment
c!     
c      um1(1,1)=(c(1)-v1*u(1)+v2)/v3
c      um1(2,2)=(c(2)-v1*u(2)+v2)/v3
c      um1(3,3)=(c(3)-v1*u(3)+v2)/v3
c      um1(1,2)=(c(4)-v1*u(4))/v3
c      um1(1,3)=(c(5)-v1*u(5))/v3
c      um1(2,3)=(c(6)-v1*u(6))/v3
c      um1(2,1)=um1(1,2)
c      um1(3,1)=um1(1,3)
c      um1(3,2)=um1(2,3)
c!     
c!     calculating the logarithmic strain at the end of the increment
c!     Elog=Z.ln(w).Z^T
c!     
c      do i=1,3
c        w(i)=dlog(w(i))
c      enddo
c!     
c!     logarithmic strain in global coordinates at the end of the
c!     increment = ln(lambda) N diadic N;
c!     
c!     the true logarithmic strain is ln(lambda) n diadic n;
c!     the rotation tensor R = n diadic N;
c!     therefore, ln(lambda) N diadic N is R^T true logarithmic strain R, 
c!     in other words the corotational true logarithmic strain (required
c!     by the abaqus routine)
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
!     corotational logarithmic strain in global coordinates at the end of the
!     increment = ln(lambda) N diadic N;
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
!     calculating the inverse right stretch tensor U^{-1}
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
c      write(*,*) 'umat_abaqusnl_total',um1(1,1),um1(1,1)
c      write(*,*) 'umat_abaqusnl_total',um1(2,2),um1(2,2)
c      write(*,*) 'umat_abaqusnl_total',um1(3,3),um1(3,3)
c      write(*,*) 'umat_abaqusnl_total',um1(1,2),um1(1,2)
c      write(*,*) 'umat_abaqusnl_total',um1(1,3),um1(1,3)
c      write(*,*) 'umat_abaqusnl_total',um1(2,3),um1(2,3)
c      write(*,*)
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
c        if(amat(1:11).eq.'JOHNSONCOOK') then
c          stiff(1)=ddsdde(1,1)
c          stiff(2)=ddsdde(2,1)
c          stiff(3)=ddsdde(2,2)
c          stiff(4)=ddsdde(3,1)
c          stiff(5)=ddsdde(3,2)
c          stiff(6)=ddsdde(3,3)
c          stiff(7)=ddsdde(4,1)
c          stiff(8)=ddsdde(4,2)
c          stiff(9)=ddsdde(4,3)
c          stiff(10)=ddsdde(4,4)
c          stiff(11)=ddsdde(5,1)
c          stiff(12)=ddsdde(5,2)
c          stiff(13)=ddsdde(5,3)
c          stiff(14)=ddsdde(5,4)
c          stiff(15)=ddsdde(5,5)
c          stiff(16)=ddsdde(6,1)
c          stiff(17)=ddsdde(6,2)
c          stiff(18)=ddsdde(6,3)
c          stiff(19)=ddsdde(6,4)
c          stiff(20)=ddsdde(6,5)
c          stiff(21)=ddsdde(6,6)
c        endif
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
