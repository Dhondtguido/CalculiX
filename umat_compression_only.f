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
      subroutine umat_compression_only(
     &        amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt,ipkon)
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
!     elconloc(* )       user defined constants defined by the keyword
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
!                        Notice that the matrix is an integral part of the 
!                        fourth order material tensor, i.e. the Voigt notation
!                        is not used.
!     pnewdt             to be specified by the user if the material
!                        routine is unable to return the stiffness matrix
!                        and/or the stress due to divergence within the
!                        routine. pnewdt is the factor by which the time
!                        increment is to be multiplied in the next
!                        trial and should exceed zero but be less than 1.
!                        Default is -1 indicating that the user routine
!                        has converged.
!     ipkon(*)           ipkon(iel) points towards the position in field
!                        kon prior to the first node of the element's
!                        topology. If ipkon(iel) is set to -1, the 
!                        element is removed from the mesh
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &  ipkon(*),kel(4,21),n,matz,ier,i,j,kal(2,6),j1,j2,j3,j4
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,pnewdt,xstate(nstate_,mi(1),*),e(3,3),we(3),
     &  xstateini(nstate_,mi(1),*),wc(3),z(3,3),fv1(3),fv2(3),eps,
     &  pi,young,fla(3),xm1(3,3),xm2(3,3),xm3(3,3),dfla(3),a(21),b(21),
     &  d(3,3),c(3,3),xmm1(21),xmm2(21),xmm3(21),ca(3),cb(3)
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &          1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &          3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &          1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!
      d=reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),
     &       (/3,3/))
!
!     compression-only material. The constants are:
!     1. Young's modulus
!     2. maximum tension allowed (in absolute value)
!
      pi=4.d0*datan(1.d0)
!
      if(elconloc(1).lt.1.d-30) then
         write(*,*) '*ERROR in umat_compression_only: the Young modulus'
         write(*,*) '       is too small'
         call exit(201)
      endif
!
      if(elconloc(2).lt.1.d-30) then
         write(*,*) '*ERROR in umat_compression_only: maximum tension'
         write(*,*) '       value is too small'
         call exit(201)
      endif
      young=elconloc(1)
      eps=elconloc(2)*pi/young
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
     &  *ERROR calculating the eigenvalues/vectors in umat_compression'
         call exit(201)
      endif
!
!     calculating the eigenvalues of the Cauchy tensor
!     at the end of the increment (eigenvalues of the Cauchy tensor)
!
      do i=1,3
         wc(i)=2.d0*we(i)+1.d0
c         if(iint.eq.1)
c     &       write(*,*) 'eigenvalues',i,we(i),z(1,i),z(2,i),z(3,i)
      enddo
!
!     check for 3 equal eigenvalues
!
      if((dabs(we(3)-we(2)).lt.1.d-10).and.
     &   (dabs(we(2)-we(1)).lt.1.d-10)) then
         fla(1)=young*we(1)*(0.5d0+datan(-we(1)/eps)/pi)
         do i=1,3
            stre(i)=fla(1)-beta(i)
         enddo
         do i=4,6
            stre(i)=0.d0-beta(i)
         enddo
c         if(iint.eq.1) write(*,*) '1 diff stre'
c         if(iint.eq.1) write(*,100) (stre(i),i=1,6)
!
         if(icmd.ne.3) then
            dfla(1)=young*((0.5d0+datan(-we(1)/eps)/pi)/2.d0
     &             -we(1)/(2.d0*pi*eps*(1.d0+(we(1)/eps)**2)))
!
            do j=1,21
               j1=kel(1,j)
               j2=kel(2,j)
               j3=kel(3,j)
               j4=kel(4,j)
               stiff(j)=dfla(1)*(d(j3,j1)*d(j4,j2)+d(j3,j2)*d(j4,j1))
            enddo
c            if(iint.eq.1) write(*,*) '1 diff stiff'
c            if(iint.eq.1) write(*,100) (stiff(i),i=1,21)
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
!        compression-only modification)
!
         do i=1,2
            fla(i)=young*we(i)*(0.5d0+datan(-we(i)/eps)/pi)
         enddo
!
!        calculating the stresses
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            stre(j)=fla(1)*xm1(j1,j2)+fla(2)*(d(j1,j2)-xm1(j1,j2))
     &             -beta(j)
         enddo
c         if(iint.eq.1) write(*,*) '2 diff stre (w2=w3)'
c         if(iint.eq.1) write(*,100) (stre(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the stiffness matrix
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=1,2
               dfla(i)=young*((0.5d0+datan(-we(i)/eps)/pi)/2.d0
     &                -we(i)/(2.d0*pi*eps*(1.d0+(we(i)/eps)**2)))
            enddo
!
!           stiffness matrix
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
!              calculating the stiffness matrix
!
               stiff(j)=2.d0*(dfla(1)*xmm1(j)
     &                 +dfla(2)*(a(j)+xmm1(j)-b(j))
     &                 +(fla(1)-fla(2))*(b(j)-2.d0*xmm1(j))
     &                 /(2.d0*(we(1)-we(2))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff stiff (w2=w3)'
c            if(iint.eq.1) write(*,100) (stiff(i),i=1,21)
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
!        compression-only modification)
!
         do i=2,3
            fla(i)=young*we(i)*(0.5d0+datan(-we(i)/eps)/pi)
         enddo
!
!        calculating the stresses
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            stre(j)=fla(3)*xm3(j1,j2)+fla(2)*(d(j1,j2)-xm3(j1,j2))
     &             -beta(j)
         enddo
c         if(iint.eq.1) write(*,*) '2 diff stre (w1=w2)'
c         if(iint.eq.1) write(*,100) (stre(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the stiffness matrix
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=2,3
               dfla(i)=young*((0.5d0+datan(-we(i)/eps)/pi)/2.d0
     &                -we(i)/(2.d0*pi*eps*(1.d0+(we(i)/eps)**2)))
            enddo
!
!           stiffness matrix
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
!              calculating the stiffness matrix
!
               stiff(j)=2.d0*(dfla(3)*xmm3(j)
     &                 +dfla(2)*(a(j)+xmm3(j)-b(j))
     &                 +(fla(2)-fla(3))*(b(j)-2.d0*xmm3(j))
     &                 /(2.d0*(we(2)-we(3))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff stiff (w1=w2)'
c            if(iint.eq.1) write(*,100) (stiff(i),i=1,21)
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
!        compression-only modification)
!
         do i=2,3
            fla(i)=young*we(i)*(0.5d0+datan(-we(i)/eps)/pi)
         enddo
!
!        calculating the stresses
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            stre(j)=fla(2)*xm3(j1,j2)+fla(3)*(d(j1,j2)-xm3(j1,j2))
     &             -beta(j)
         enddo
c         if(iint.eq.1) write(*,*) '2 diff stre (w1=w2)'
c         if(iint.eq.1) write(*,100) (stre(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the stiffness matrix
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=2,3
               dfla(i)=young*((0.5d0+datan(-we(i)/eps)/pi)/2.d0
     &                -we(i)/(2.d0*pi*eps*(1.d0+(we(i)/eps)**2)))
            enddo
!
!           stiffness matrix
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
!              calculating the stiffness matrix
!
               stiff(j)=2.d0*(dfla(2)*xmm3(j)
     &                 +dfla(3)*(a(j)+xmm3(j)-b(j))
     &                 +(fla(3)-fla(2))*(b(j)-2.d0*xmm3(j))
     &                 /(2.d0*(we(3)-we(2))))
            enddo
c            if(iint.eq.1) write(*,*) '2 diff stiff (w1=w2)'
c            if(iint.eq.1) write(*,100) (stiff(i),i=1,21)
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
!        compression-only modification)
!
         do i=1,3
            fla(i)=young*we(i)*(0.5d0+datan(-we(i)/eps)/pi)
         enddo
!
!        calculating the stresses
!
         do j=1,6
            j1=kal(1,j)
            j2=kal(2,j)
            stre(j)=fla(1)*xm1(j1,j2)+fla(2)*xm2(j1,j2)
     &                               +fla(3)*xm3(j1,j2)-beta(j)
         enddo
c         if(iint.eq.1) write(*,*) '3 diff stre'
c         if(iint.eq.1) write(*,100) (stre(i),i=1,6)
!
         if(icmd.ne.3) then
!
!           calculating the stiffness matrix
!
!           derivative of the modified young's modulus w.r.t.
!           the Cauchy eigenvalues
!            
            do i=1,3
               dfla(i)=young*((0.5d0+datan(-we(i)/eps)/pi)/2.d0
     &                -we(i)/(2.d0*pi*eps*(1.d0+(we(i)/eps)**2)))
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
!           stiffness matrix
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
!              calculating the stiffness matrix
!
               stiff(j)=2.d0*(dfla(1)*xmm1(j)+dfla(2)*xmm2(j)+
     &                        dfla(3)*xmm3(j)
     &                 +fla(1)*(cb(1)*b(j)+ca(1)*a(j))
     &                 +fla(2)*(cb(2)*b(j)+ca(2)*a(j))
     &                 +fla(3)*(cb(3)*b(j)+ca(3)*a(j)))
            enddo
c            if(iint.eq.1) write(*,*) '3 diff stiff'
c            if(iint.eq.1) write(*,100) (stiff(i),i=1,21)
         endif
      endif
!
c 100  format(6(1x,e11.4))
!
      return
      end


