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
      subroutine umat_lin_el_corot(
     &        amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,eloc,nlgeom_undo)
!
!     calculates stiffness and stresses for a user defined material
!     law
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
!     eloc(6)            Lagrange total strain tensor (component order:
!                        11,22,33,12,13,23) at the end of the increment
!                        (thermal strains are subtracted)
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
!     eloc(6)            linear total strain tensor for large
!                        rotations (component order:
!                        11,22,33,12,13,23) at the end of the increment
!                        (thermal strains are subtracted)
!     nlgeom_undo        0: Lagrange strain goes out
!                        1: linear strain for large rotations goes out
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &     kel(4,21),n,matz,ier,i,j,kal(2,6),j1,j2,j3,j4,nlgeom_undo,
     &     nconstants
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,pnewdt,xstate(nstate_,mi(1),*),e(3,3),we(3),
     &  xstateini(nstate_,mi(1),*),wc(3),z(3,3),fv1(3),fv2(3),eps,
     &  pi,young,fla(3),xm1(3,3),xm2(3,3),xm3(3,3),dfla(3),a(21),b(21),
     &  d(3,3),c(3,3),xmm1(21),xmm2(21),xmm3(21),ca(3),cb(3),u(6),
     &  dude(21),um,un,al,am1,am2,ee,eloc(6),eth(6),elin(6)
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &          1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &          3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &          1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!
      d=reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),
     &     (/3,3/))
!
      nlgeom_undo=1
!
!     corotational linear elastic material. The constants are the
!     the same as for linear elastic material (2 for isotropic,
!     9 for orthotropic and 21 for completely anisotropic material)
!
!     PK2=(Hooke matrix)*(U-I)
!
!     U is the right stretch tensor = sum sqrt(Lambda_i)*M_i
!     = square root of C
!
      pi=4.d0*datan(1.d0)
!
c      if(elconloc(1).lt.1.d-30) then
c         write(*,*) '*ERROR in umat_lin_el_corot: the Young modulus'
c         write(*,*) '       is too small'
c         call exit(201)
c      endif
c!
c      if(elconloc(2).lt.1.d-30) then
c         write(*,*) '*ERROR in umat_lin_el_corot: maximum pressure'
c         write(*,*) '       value is too small'
c         call exit(201)
c      endif
c      young=elconloc(1)
c      eps=elconloc(2)*pi/young
!
!     calculation of the eigenvalues and eigenvectors of the
!     Lagrange strain
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
c      if(iint.eq.1) write(*,*) 'emec'
c      if(iint.eq.1) write(*,100) (emec(i),i=1,6)
!
      n=3
      matz=1
      ier=0
!
      call rs(n,n,e,we,matz,z,fv1,fv2,ier)
!
      if(ier.ne.0) then
         write(*,*) '
     &  *ERROR calculating the eigenvalues/vectors in umat_tension'
         call exit(201)
      endif
!
!     calculating the eigenvalues of the Cauchy tensor
!     at the end of the increment (eigenvalues of the Cauchy tensor)
!
      do i=1,3
         wc(i)=2.d0*we(i)+1.d0
c         if(iint.eq.1)
c     &       write(*,*) 'eigenvalues',i,wc(i),z(1,i),z(2,i),z(3,i)
      enddo
!
!     check for 3 equal eigenvalues
!
      if((dabs(we(3)-we(2)).lt.1.d-10).and.
     &   (dabs(we(2)-we(1)).lt.1.d-10)) then
         fla(1)=dsqrt(wc(1))
         do i=1,3
            u(i)=fla(1)
         enddo
         do i=4,6
            u(i)=0.d0
         enddo
c         if(iint.eq.1) write(*,*) '1 diff u'
c         if(iint.eq.1) write(*,100) (u(i),i=1,6)
!
         if(icmd.ne.3) then
            dfla(1)=1.d0/(2.d0*fla(1))
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
               dude(j)=dfla(1)*(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))
            enddo
c            if(iint.eq.1) write(*,*) '1 diff dude'
c            if(iint.eq.1) write(*,100) (dude(i),i=1,21)
         endif
!
!     2 equal eigenvalues: w2=w3
!
      elseif(dabs(we(3)-we(2)).lt.1.d-10) then
!
!        calculating the structural tensor for the unequal
!        eigenvalue of the Cauchy and the Lagrange tensor
!         
         do i=1,3
            do j=1,3
               xm1(i,j)=z(i,1)*z(j,1)
            enddo
         enddo
!
!        calculating the modified Young's modulus (due to the
!        tension-only modification)
!
         do i=1,2
            fla(i)=dsqrt(wc(i))
         enddo
!
!        calculating the right stretch tensor
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            u(j)=fla(1)*xm1(j1,j2)+fla(2)*(d(j1,j2)-xm1(j1,j2))
         enddo
c         if(iint.eq.1) write(*,*) '2 diff u (w2=w3)'
c         if(iint.eq.1) write(*,100) (u(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the matrix dU/dE
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=1,2
               dfla(i)=1.d0/(2.d0*fla(i))
            enddo
!
!           matrix dU/dE
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
!
!              calculating the auxiliary field xmm1
!
               xmm1(j)=xm1(j1,j2)*xm1(j3,j4)
!
!              calculating the auxiliary fields a and b
!
               a(j)=(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))/2.d0
               b(j)=(d(j3,j1)*xm1(j4,j2)+d(j4,j1)*xm1(j3,j2)+
     &               d(j4,j2)*xm1(j1,j3)+d(j3,j2)*xm1(j1,j4))/2.d0
!
!              calculating the matrix dU/dE
!
               dude(j)=2.d0*(dfla(1)*xmm1(j)
     &                 +dfla(2)*(a(j)+xmm1(j)-b(j))
     &                 +(fla(1)-fla(2))*(b(j)-2.d0*xmm1(j))
     &                 /(2.d0*(we(1)-we(2))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff dude (w2=w3)'
c            if(iint.eq.1) write(*,100) (dude(i),i=1,21)
         endif
!
!     2 equal eigenvalues: w1=w2
!
      elseif(dabs(we(2)-we(1)).lt.1.d-10) then
!
!        calculating the structural tensor for the unequal
!        eigenvalue of the Cauchy and the Lagrange tensor
!         
         do i=1,3
            do j=1,3
               xm3(i,j)=z(i,3)*z(j,3)
            enddo
         enddo
!
!        calculating the modified Young's modulus (due to the
!        tension-only modification)
!
         do i=2,3
            fla(i)=dsqrt(wc(i))
         enddo
!
!        calculating the right stretch tensor
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            u(j)=fla(3)*xm3(j1,j2)+fla(2)*(d(j1,j2)-xm3(j1,j2))
         enddo
c         if(iint.eq.1) write(*,*) '2 diff u (w1=w2)'
c         if(iint.eq.1) write(*,100) (u(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the matrix dU/dE
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=2,3
               dfla(i)=1.d0/(2.d0*fla(i))
            enddo
!
!           matrix dU/dE
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
!
!              calculating the auxiliary field xmm3
!
               xmm3(j)=xm3(j1,j2)*xm3(j3,j4)
!
!              calculating the auxiliary fields a and b
!
               a(j)=(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))/2.d0
               b(j)=(d(j3,j1)*xm3(j4,j2)+d(j4,j1)*xm3(j3,j2)+
     &               d(j4,j2)*xm3(j1,j3)+d(j3,j2)*xm3(j1,j4))/2.d0
!
!              calculating the matrix dU/dE
!
               dude(j)=2.d0*(dfla(3)*xmm3(j)
     &                 +dfla(2)*(a(j)+xmm3(j)-b(j))
     &                 +(fla(2)-fla(3))*(b(j)-2.d0*xmm3(j))
     &                 /(2.d0*(we(2)-we(3))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff dude (w1=w2)'
c            if(iint.eq.1) write(*,100) (dude(i),i=1,21)
         endif
!
!     2 equal eigenvalues: w1=w3
!
      elseif(dabs(we(3)-we(1)).lt.1.d-10) then
!
!        calculating the structural tensor for the unequal
!        eigenvalue of the Cauchy and the Lagrange tensor
!         
         do i=1,3
            do j=1,3
               xm3(i,j)=z(i,3)*z(j,3)
            enddo
         enddo
!
!        calculating the modified Young's modulus (due to the
!        tension-only modification)
!
         do i=2,3
            fla(i)=dsqrt(wc(i))
         enddo
!
!        calculating the right stretch tensor
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            u(j)=fla(2)*xm3(j1,j2)+fla(3)*(d(j1,j2)-xm3(j1,j2))
         enddo
c         if(iint.eq.1) write(*,*) '2 diff u (w1=w2)'
c         if(iint.eq.1) write(*,100) (u(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the matrix dU/dE
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=2,3
               dfla(i)=1.d0/(2.d0*fla(i))
            enddo
!
!           matrix dU/dE
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
!
!              calculating the auxiliary field xmm3
!
               xmm3(j)=xm3(j1,j2)*xm3(j3,j4)
!
!              calculating the auxiliary fields a and b
!
               a(j)=(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))/2.d0
               b(j)=(d(j3,j1)*xm3(j4,j2)+d(j4,j1)*xm3(j3,j2)+
     &               d(j4,j2)*xm3(j1,j3)+d(j3,j2)*xm3(j1,j4))/2.d0
!
!              calculating the matrix dU/dE
!
               dude(j)=2.d0*(dfla(2)*xmm3(j)
     &                 +dfla(3)*(a(j)+xmm3(j)-b(j))
     &                 +(fla(3)-fla(2))*(b(j)-2.d0*xmm3(j))
     &                 /(2.d0*(we(3)-we(2))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff dude (w1=w2)'
c            if(iint.eq.1) write(*,100) (dude(i),i=1,21)
         endif
!
!     3 different eigenvalues
!
      else
!
!        calculating the structural tensors of the Cauchy and the
!        Lagrange tensor
!         
         do i=1,3
            do j=1,3
               xm1(i,j)=z(i,1)*z(j,1)
               xm2(i,j)=z(i,2)*z(j,2)
               xm3(i,j)=z(i,3)*z(j,3)
            enddo
         enddo
!
!        calculating the modified Young's modulus (due to the
!        tension-only modification)
!
         do i=1,3
            fla(i)=dsqrt(wc(i))
         enddo
!
!        calculating the right stretch tensor
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            u(j)=fla(1)*xm1(j1,j2)+fla(2)*xm2(j1,j2)
     &                               +fla(3)*xm3(j1,j2)
         enddo
c         if(iint.eq.1) write(*,*) '3 diff u'
c         if(iint.eq.1) write(*,100) (u(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the matrix dU/dE
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=1,3
               dfla(i)=1.d0/(2.d0*fla(i))
            enddo
!
            cb(1)=1.d0/(4.d0*(we(1)-we(2))*(we(1)-we(3)))
            ca(1)=-(wc(2)+wc(3))*cb(1)
            cb(2)=1.d0/(4.d0*(we(2)-we(3))*(we(2)-we(1)))
            ca(2)=-(wc(3)+wc(1))*cb(2)
            cb(3)=1.d0/(4.d0*(we(3)-we(1))*(we(3)-we(2)))
            ca(3)=-(wc(1)+wc(2))*cb(3)
!
!           Cauchy tensor
!
            do i=1,3
               do j=1,3
                  c(i,j)=2.d0*e(i,j)
               enddo
               c(i,i)=c(i,i)+1.d0
            enddo
!
!           matrix dU/dE
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
!
!              calculating the auxiliary fields xmm1,xmm2,xmm3
!
               xmm1(j)=xm1(j1,j2)*xm1(j3,j4)
               xmm2(j)=xm2(j1,j2)*xm2(j3,j4)
               xmm3(j)=xm3(j1,j2)*xm3(j3,j4)
!
!              calculating the auxiliary fields a and b
!              (Dhondt, G., The Finite Element Method for Three-
!              Dimensional Thermomechanical Applications, Wiley 2004,
!              formulae (4.203) and (4.204))
!
               a(j)=(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))/2.d0
     &             -xmm1(j)-xmm2(j)-xmm3(j)
               b(j)=(d(j3,j1)*c(j4,j2)+d(j4,j1)*c(j3,j2)+
     &               d(j4,j2)*c(j1,j3)+d(j3,j2)*c(j1,j4))/2.d0
     &             -2.d0*(wc(1)*xmm1(j)+wc(2)*xmm2(j)+wc(3)*xmm3(j))
!
!              calculating the matrix dU/dE
!
               dude(j)=2.d0*(dfla(1)*xmm1(j)+dfla(2)*xmm2(j)+
     &                        dfla(3)*xmm3(j)
     &                 +fla(1)*(cb(1)*b(j)+ca(1)*a(j))
     &                 +fla(2)*(cb(2)*b(j)+ca(2)*a(j))
     &                 +fla(3)*(cb(3)*b(j)+ca(3)*a(j)))
            enddo
c            if(iint.eq.1) write(*,*) '3 diff dude'
c            if(iint.eq.1) write(*,100) (dude(i),i=1,21)
         endif
      endif
!
      nconstants=-kode-100
!
      if(nconstants.eq.2) then
!
!        isotropic
!
!        calculating the elastic constants
!
         ee=elconloc(1)
         un=elconloc(2)
         um=ee/(1.d0+un)
         al=un*um/(1.d0-2.d0*un)
         um=um/2.d0
         am1=al+2.d0*um
         am2=2.d0*um
!
!        determining the thermal strain
!
c         do i=1,6
c            eth(i)=eloc(i)-emec(i)
c         enddo
!
!        calculating the strain tensor
!         
         do i=1,3
            elin(i)=u(i)-1.d0
         enddo
         do i=4,6
            elin(i)=u(i)
         enddo
!
!        determining the new total strain
!
c         do i=1,6
c            eloc(i)=elin(i)+eth(i)
c         enddo
!
!        calculating the stresses         
!     
         stre(1)=am1*elin(1)+al*(elin(2)+elin(3))-beta(1)
         stre(2)=am1*elin(2)+al*(elin(1)+elin(3))-beta(2)
         stre(3)=am1*elin(3)+al*(elin(1)+elin(2))-beta(3)
         stre(4)=am2*elin(4)-beta(4)
         stre(5)=am2*elin(5)-beta(5)
         stre(6)=am2*elin(6)-beta(6)
!         
         if(icmd.ne.3) then
c         if(iint.eq.1) then
cc     write(*,100) (stre(i),i=1,6)
c            write(*,*) 'dude'
c            write(*,100) (dude(i),i=1,6)
c            write(*,100) (dude(i),i=7,10)
c            write(*,100) (dude(i),i=11,15)
c            write(*,100) (dude(i),i=16,21)
c         endif
!     
!     calculating the stiffness matrix
!
!     stiff_klpq=lambda*(delta_kl*dude_mmpq+delta_pq*dude_mmkl)/2
!               +2*mu*dude_klpq            
!
            do i=1,21
               stiff(i)=2.d0*um*dude(i)
            enddo
!
            stiff( 1)=stiff( 1)+al*(
     &        dude( 1)+
     &        dude( 2)+
     &        dude( 4)+
     &        dude( 1)+
     &        dude( 2)+
     &        dude( 4)
     &        )/2.d0
            stiff( 2)=stiff( 2)+al*(
     &        dude( 2)+
     &        dude( 3)+
     &        dude( 5)+
     &        dude( 1)+
     &        dude( 2)+
     &        dude( 4)
     &        )/2.d0
            stiff( 3)=stiff( 3)+al*(
     &        dude( 2)+
     &        dude( 3)+
     &        dude( 5)+
     &        dude( 2)+
     &        dude( 3)+
     &        dude( 5)
     &        )/2.d0
            stiff( 4)=stiff( 4)+al*(
     &        dude( 4)+
     &        dude( 5)+
     &        dude( 6)+
     &        dude( 1)+
     &        dude( 2)+
     &        dude( 4)
     &        )/2.d0
            stiff( 5)=stiff( 5)+al*(
     &        dude( 4)+
     &        dude( 5)+
     &        dude( 6)+
     &        dude( 2)+
     &        dude( 3)+
     &        dude( 5)
     &        )/2.d0
            stiff( 6)=stiff( 6)+al*(
     &        dude( 4)+
     &        dude( 5)+
     &        dude( 6)+
     &        dude( 4)+
     &        dude( 5)+
     &        dude( 6)
     &        )/2.d0
            stiff( 7)=stiff( 7)+al*(
     &        dude( 7)+
     &        dude( 8)+
     &        dude( 9)
     &        )/2.d0
            stiff( 8)=stiff( 8)+al*(
     &        dude( 7)+
     &        dude( 8)+
     &        dude( 9)
     &        )/2.d0
            stiff( 9)=stiff( 9)+al*(
     &        dude( 7)+
     &        dude( 8)+
     &        dude( 9)
     &        )/2.d0
            stiff(11)=stiff(11)+al*(
     &        dude(11)+
     &        dude(12)+
     &        dude(13)
     &        )/2.d0
            stiff(12)=stiff(12)+al*(
     &        dude(11)+
     &        dude(12)+
     &        dude(13)
     &        )/2.d0
            stiff(13)=stiff(13)+al*(
     &        dude(11)+
     &        dude(12)+
     &        dude(13)
     &        )/2.d0
            stiff(16)=stiff(16)+al*(
     &        dude(16)+
     &        dude(17)+
     &        dude(18)
     &        )/2.d0
            stiff(17)=stiff(17)+al*(
     &        dude(16)+
     &        dude(17)+
     &        dude(18)
     &        )/2.d0
            stiff(18)=stiff(18)+al*(
     &        dude(16)+
     &        dude(17)+
     &        dude(18)
     &        )/2.d0
c            stiff(1)=stiff(1)+al*(dude(1)+dude(2)+dude(4))
c            stiff(2)=stiff(2)+al*(dude(2)+dude(3)+dude(5))
c            stiff(3)=stiff(3)+al*(dude(4)+dude(5)+dude(6))
c            stiff(4)=stiff(4)+al*(dude(2)+dude(3)+dude(5))
c            stiff(5)=stiff(5)+al*(dude(4)+dude(5)+dude(6))
c            stiff(6)=stiff(6)+al*(dude(4)+dude(5)+dude(6))
c         if(iint.eq.1) then
cc            write(*,100) (stre(i),i=1,6)
c            write(*,*) 'stiff'
c            write(*,100) (stiff(i),i=1,6)
c            write(*,100) (stiff(i),i=7,10)
c            write(*,100) (stiff(i),i=11,15)
c            write(*,100) (stiff(i),i=16,21)
c         endif
         endif
      endif
 100  format(6(1x,e15.8))
!
      return
      end


