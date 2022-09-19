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
      subroutine umat_undo_nlgeom_lin_iso_el(amat,iel,iint,kode,
     &        elconloc,emec,emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,
     &        dtime,time,ttime,icmd,ielas,mi,nstate_,xstateini,xstate,
     &        stre,stiff,iorien,pgauss,orab,eloc,nlgeom_undo)
!
!     calculates stiffness and stresses for a linear elastic isotropic
!     material with special modification of the strain tensor
!
!     The strain tensor that enters the routine is the Lagrange strain
!     tensor. This tensor is modified into an infinitesimal strain tensor
!     for large rotations. Application: singular stress/strain field at a
!     crack tip
!
!     Procedure: from the deformation gradient F the rotation is 
!                removed to yield U=R^(-1).F
!                The resulting field U is linearized
!                E=((U-I)+(U-I)^T)/2 and subsequently
!                used to calculate the stresses by a
!                linear elastic isotropic relationship
!
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
!     emec(6)            linear mechanical strain tensor for large
!                        rotations (component order:
!                        11,22,33,12,13,23) at the end of the increment
!                        (thermal strains are subtracted)
!     eloc(6)            linear total strain tensor for large
!                        rotations (component order:
!                        11,22,33,12,13,23) at the end of the increment
!     nlgeom_undo        0: Lagrange strain goes out
!                        1: linear strain for large rotations goes out
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &  i,j,k,nlgeom_undo,n,matz,ier
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,ukl(3,3),ca(3,3),elag(3,3),elin(3,3),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  e,un,al,um,am1,am2,eme1(6),eloc(6),w(3),z(3,3),fv1(3),
     &  fv2(3),v1,v2,v3,c(6),c2(6),dd,d(6),u(6),ur(3,3)
!
      d=(/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
!
      nlgeom_undo=1
!
!     calculate the rotationless linearized strain
!
!     a) calculating the Cauchy-Green tensor C=F^T.F
!
      do i=1,3
         do j=1,3
            ca(i,j)=0.d0
            do k=1,3
               ca(i,j)=ca(i,j)+xkl(k,i)*xkl(k,j)
            enddo
         enddo
      enddo
c      write(*,*) 'umat_undo..xkl',xkl(1,1),xkl(1,2),xkl(1,3)
c      write(*,*) 'umat_undo..xkl',xkl(2,1),xkl(2,2),xkl(2,3)
c      write(*,*) 'umat_undo..xkl',xkl(3,1),xkl(3,2),xkl(3,3)
c      write(*,*) 'umat_undo..ca',ca(1,1)
!
!     b) calculating the principal stretches
!
      n=3
      matz=1
!
      call rs(n,n,ca,w,matz,z,fv1,fv2,ier)
!
      if(ier.ne.0) then
         write(*,*) '
     & *ERROR calculating the eigenvalues/vectors in '
         write(*,*) '       umat_undo_nlgeom_lin_iso_el'
         call exit(201)
      endif
      w(1)=dsqrt(w(1))
      w(2)=dsqrt(w(2))
      w(3)=dsqrt(w(3))
!
!     c) calculating the invariants of U
!
      v1=w(1)+w(2)+w(3)
      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
      v3=w(1)*w(2)*w(3)
c      write(*,*) 'umat_undo..v1,v2,v3',v1,v2,v3
!
!     d) storing C as vector
!
      c(1)=ca(1,1)
      c(2)=ca(2,2)
      c(3)=ca(3,3)
      c(4)=ca(1,2)
      c(5)=ca(1,3)
      c(6)=ca(2,3)
!
!     e) calculating C.C
!
      c2(1)=c(1)*c(1)+c(4)*c(4)+c(5)*c(5)
      c2(2)=c(4)*c(4)+c(2)*c(2)+c(6)*c(6)
      c2(3)=c(5)*c(5)+c(6)*c(6)+c(3)*c(3)
      c2(4)=c(1)*c(4)+c(4)*c(2)+c(5)*c(6)
      c2(5)=c(1)*c(5)+c(4)*c(6)+c(5)*c(3)
      c2(6)=c(4)*c(5)+c(2)*c(6)+c(6)*c(3)
!
!     f) calculating the right stretch tensor U
!        (cf. Simo and Hughes, Computational Inelasticity)
!
      dd=v1*v2-v3
      do i=1,6
         u(i)=(-c2(i)+(v1*v1-v2)*c(i)+v1*v3*d(i))/dd
      enddo
!
!     g) storing U as matrix (U=R^(-1).F)
!
      ur(1,1)=u(1)
      ur(2,2)=u(2)
      ur(3,3)=u(3)
      ur(1,2)=u(4)
      ur(1,3)=u(5)
      ur(2,3)=u(6)
      ur(2,1)=u(4)
      ur(3,1)=u(5)
      ur(3,2)=u(6)
c      write(*,*) 'umat_undo..ur',ur(1,1),ur(1,2),ur(1,3)
c      write(*,*) 'umat_undo..ur',ur(2,1),ur(2,2),ur(2,3)
c      write(*,*) 'umat_undo..ur',ur(3,1),ur(3,2),ur(3,3)
!
!     i) calculating the rotationless displacement gradient u_k,l
!        = (U-I)
!
      do i=1,3
         do j=1,3
            ukl(i,j)=ur(i,j)
         enddo
         ukl(i,i)=ukl(i,i)-1.d0
      enddo
!
!     j) calculating the total rotationless linearized strain
!        ((U-I)+(U-I)^T)/2
!
      do i=1,3
         do j=i,3
            elin(i,j)=(ukl(i,j)+ukl(j,i))/2.d0
         enddo
      enddo
c      write(*,*) 'umat_undo..elin',elin(1,1)
!
!     k) calculating the mechanical rotationless linearized strain =
!        total rotationless linearized strain -
!        (total Lagrange strain - mechanical Lagrange strain)
!
      emec(1)=elin(1,1)-eloc(1)+emec(1)
      emec(2)=elin(2,2)-eloc(2)+emec(2)
      emec(3)=elin(3,3)-eloc(3)+emec(3)
      emec(4)=elin(1,2)-eloc(4)+emec(4)
      emec(5)=elin(1,3)-eloc(5)+emec(5)
      emec(6)=elin(2,3)-eloc(6)+emec(6)
!
      eloc(1)=elin(1,1)
      eloc(2)=elin(2,2)
      eloc(3)=elin(3,3)
      eloc(4)=elin(1,2)
      eloc(5)=elin(1,3)
      eloc(6)=elin(2,3)
!
c      do i=1,6
c         write(*,*) 'umat...lin_el',time,iel,iint,eloc(i),stre(i)
c      enddo
c      write(*,*) 'umat_undo..eloc',eloc(1)
!
!     calculating the stresses
!
      e=elconloc(1)
      un=elconloc(2)
      al=un*e/(1.d0+un)/(1.d0-2.d0*un)
      um=e/2.d0/(1.d0+un)
      am1=al+2.d0*um
      am2=2.d0*um
!      
      stre(1)=am1*emec(1)+al*(emec(2)+emec(3))-beta(1)
      stre(2)=am1*emec(2)+al*(emec(1)+emec(3))-beta(2)
      stre(3)=am1*emec(3)+al*(emec(1)+emec(2))-beta(3)
      stre(4)=am2*emec(4)-beta(4)
      stre(5)=am2*emec(5)-beta(5)
      stre(6)=am2*emec(6)-beta(6)
c      write(*,*) 'umat_undo..stre',stre(1)
!
c      do i=1,6
c         write(*,*) 'umat...lin_el',time,iel,iint,eloc(i),stre(i)
c      enddo
!
      if(icmd.ne.3) then
!
!        calculating the stiffness matrix
!
         stiff(1)=al+2.d0*um
         stiff(2)=al
         stiff(3)=al+2.d0*um
         stiff(4)=al
         stiff(5)=al
         stiff(6)=al+2.d0*um
         stiff(7)=0.d0
         stiff(8)=0.d0
         stiff(9)=0.d0
         stiff(10)=um
         stiff(11)=0.d0
         stiff(12)=0.d0
         stiff(13)=0.d0
         stiff(14)=0.d0
         stiff(15)=um
         stiff(16)=0.d0
         stiff(17)=0.d0
         stiff(18)=0.d0
         stiff(19)=0.d0
         stiff(20)=0.d0
         stiff(21)=um
      endif
!
      return
      end
