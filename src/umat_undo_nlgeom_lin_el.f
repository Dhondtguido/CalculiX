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
      subroutine umat_undo_nlgeom_lin_el(amat,iel,iint,kode,
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
!     Procedure: The Biot strain U-I is
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
     &  i,j,k,nlgeom_undo,n,matz,ier,nconstants,mattyp
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
      write(*,*) 'umat_undo_nlgeom_lin_el ielas=',ielas
!     
      nlgeom_undo=1
!     
!     a) calculate the Cauchy-Green tensor C from the mechanical Lagrange strain
!     
      do i=1,6
        c(i)=2*emec(i)
      enddo
      do i=1,3
        c(i)=c(i)+1.d0
      enddo
!     
!     b) storing C as matrix
!     
      ca(1,1)=c(1)
      ca(2,2)=c(2)
      ca(3,3)=c(3)
      ca(1,2)=c(4)
      ca(1,3)=c(5)
      ca(2,3)=c(6)
      ca(2,1)=c(4)
      ca(3,1)=c(5)
      ca(3,2)=c(6)
!     
!     c) calculating the principal stretches
!     
      n=3
      matz=1
!     
      call rs(n,n,ca,w,matz,z,fv1,fv2,ier)
!     
      if(ier.ne.0) then
        write(*,*) '
     &*ERROR calculating the eigenvalues/vectors in '
        write(*,*) '       umat_undo_nlgeom_lin_el'
        call exit(201)
      endif
      w(1)=dsqrt(w(1))
      w(2)=dsqrt(w(2))
      w(3)=dsqrt(w(3))
!     
!     d) calculating the invariants of U
!     
      v1=w(1)+w(2)+w(3)
      v2=w(1)*w(2)+w(2)*w(3)+w(3)*w(1)
      v3=w(1)*w(2)*w(3)
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
!     (cf. Simo and Hughes, Computational Inelasticity)
!     
      dd=v1*v2-v3
      do i=1,6
        u(i)=(-c2(i)+(v1*v1-v2)*c(i)+v1*v3*d(i))/dd
      enddo
!     
!     g) calculating the Biot strain = U-I
!     
      do i=1,3
        emec(i)=u(i)-1.d0
      enddo
      do i=4,6
        emec(i)=u(i)
      enddo
!     
      nconstants=-kode-100
!
!     calculating the stress and the linear elastic material data
!
      call linel(nconstants,mattyp,beta,emec,stre,stiff,elconloc,
     &     iorien,orab,pgauss)
!
c      do i=1,6
c         write(*,*) 'umat...lin_el',time,iel,iint,elin(i),stre(i)
c      enddo
!
!     for a user material (umat) the material is considered to be
!     fully anisotropic; rewriting the isotropic or orthotropic
!     data in an anisotropic format
!
      if(icmd.ne.3) then
         if(mattyp.eq.1) then
!
!           isotropic
!        
            e=elconloc(1)
            un=elconloc(2)
            um=e/(1.d0+un)
            al=un*um/(1.d0-2.d0*un)
            um=um/2.d0
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
         elseif(mattyp.eq.2) then
!     
!           orthotropic
!     
            stiff(10)=stiff(7)
            stiff(15)=stiff(8)
            stiff(21)=stiff(9)
            do i=7,9
               stiff(i)=0.d0
            enddo
            do i=11,14
               stiff(i)=0.d0
            enddo
            do i=16,20
               stiff(i)=0.d0
            enddo
         endif
      endif
!
      return
      end
