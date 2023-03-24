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
      subroutine umat_ideal_gas(amat,iel,iint,kode,elconloc,emec,emec0,
     &     beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &     icmd,ielas,mi,
     &     nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab,
     &     physcon)
!
!     calculates stiffness and stresses for an ideal gas
!     For this material there is just one material constant equal to
!     initial density x specific gas constant
!     The user should list this constant (which does not depend on
!     the temperature)
!     underneath the *USER MATERIAL,CONSTANTS=1 card. The name of the
!     material has to start with IDEAL_GAS, e.g. IDEAL_GAS_AIR or
!     IDEAL_GAS_NITROGEN etc.
!
!     This routine should only be used with nlgeom=yes and
!     *PHYSICAL CONSTANTS must be used to define absolute zero temperature
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
!     physcon(1)         absolute zero (for the temperature)
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
     &  i,
     &  kk(84),k,l,m,n,nt
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  rho0r,c(3,3),cinv(3,3),v3,didc(3,3),physcon(*)
!
      kk=(/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &  1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,3,3,1,3,
     &  1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3,2,3,
     &  2,3,2,3/)
!
      rho0r=elconloc(1)
!
!     calculation of the Green deformation tensor from the Lagrange
!     strain
!     
      do i=1,3
         c(i,i)=emec(i)*2.d0+1.d0
      enddo
      c(1,2)=2.d0*emec(4)
      c(1,3)=2.d0*emec(5)
      c(2,3)=2.d0*emec(6)
!
!     calculation of the third invariant of C
!
      v3=c(1,1)*(c(2,2)*c(3,3)-c(2,3)*c(2,3))
     &     -c(1,2)*(c(1,2)*c(3,3)-c(1,3)*c(2,3))
     &     +c(1,3)*(c(1,2)*c(2,3)-c(1,3)*c(2,2))
!
!     inversion of c
!     
      cinv(1,1)=(c(2,2)*c(3,3)-c(2,3)*c(2,3))/v3
      cinv(2,2)=(c(1,1)*c(3,3)-c(1,3)*c(1,3))/v3
      cinv(3,3)=(c(1,1)*c(2,2)-c(1,2)*c(1,2))/v3
      cinv(1,2)=(c(1,3)*c(2,3)-c(1,2)*c(3,3))/v3
      cinv(1,3)=(c(1,2)*c(2,3)-c(2,2)*c(1,3))/v3
      cinv(2,3)=(c(1,2)*c(1,3)-c(1,1)*c(2,3))/v3
      cinv(2,1)=cinv(1,2)
      cinv(3,1)=cinv(1,3)
      cinv(3,2)=cinv(2,3)
c      write(*,*) 'cinv umat_ideal_gas1 ',cinv(1,1),cinv(2,2),cinv(3,3)
c      write(*,*) 'cinv umat_ideal_gas2 ',cinv(1,2),cinv(1,3),cinv(2,3)
!
!     changing the meaning of v3
!
      v3=v3*rho0r*(t1l-physcon(1))
!
!     stress from the  mechanical strain: 
!     sigma=-rho_0*r*T*C^{-1}
!     
      stre(1)=-v3*cinv(1,1)
      stre(2)=-v3*cinv(2,2)
      stre(3)=-v3*cinv(3,3)
      stre(4)=-v3*cinv(1,2)
      stre(5)=-v3*cinv(1,3)
      stre(6)=-v3*cinv(2,3)
      write(*,*) 'umat_ideal_gas ',stre(1),stre(2),stre(4)
!     
!     tangent=-2*rho_0*r*T*d(C^{-1})/dC
!     
      if(icmd.ne.3) then
!     
        nt=0
        do i=1,21
          k=kk(nt+1)
          l=kk(nt+2)
          m=kk(nt+3)
          n=kk(nt+4)
          nt=nt+4
! 
!         in my opinion, the "-" before v3 in the next line should be
!         deleted. However, without "-" divergence occurs;
!         
!
          stiff(i)=-v3*
     &         (cinv(k,m)*cinv(n,l)+cinv(k,n)*cinv(m,l))
c          stiff(i)=-v3*2.d0*
c     &         (cinv(k,m)*cinv(n,l))
        enddo
c      write(*,*) 'umat_ideal_gas: ',iel,iint,stiff(1),stiff(21)
      endif
!     
      return
      end
