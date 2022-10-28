!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     The program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine umat_single_crystal_creep(amat,iel,iint,kode,elconloc,
     &        emec,emec0,beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,
     &        ttime,icmd,ielas,
     &        mi,nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,
     &        orab,pnewdt)
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
      integer convergence
!
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien
!
      integer i,j,k,l,ipiv(18),info,neq,lda,ldb,
     &  nrhs,iplas,icounter
!
      real*8 ep0(6),ep(6),pnewdt,
     &  dg(18),ddg(18),xm(6,18),ck(18),cn(18),
     &  stri(6),htri(18),sg(18),r(6),xmc(6,18),
     &  t(6),gl(18,18),gr(18,18),ee(6),c1111,c1122,c1212,dd,
     &  skl(3,3),xmtran(3,3),ddsdde(6,6),xx(6,18)
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  elas(21),time,ttime
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
!     crystallographic slip planes:
!
!     1.  n=1,1,1   l=1,-1,0
!     2.  n=1,1,1   l=1,0,-1
!     3.  n=1,1,1   l=0,1,-1
!     4.  n=1,-1,1  l=0,1,1
!     5.  n=1,-1,1  l=1,0,-1
!     6.  n=1,-1,1  l=1,1,0
!     7.  n=1,-1,-1 l=0,1,-1
!     8.  n=1,-1,-1 l=1,0,1
!     9.  n=1,-1,-1 l=1,1,0
!     10. n=1,1,-1  l=0,1,1
!     11. n=1,1,-1  l=1,0,1
!     12. n=1,1,-1  l=1,-1,0
!     13. n=1,0,0   l=0,1,1
!     14. n=1,0,0   l=0,1,-1
!     15. n=0,1,0   l=1,0,1
!     16. n=0,1,0   l=1,0,-1
!     17. n=0,0,1   l=1,1,0
!     18. n=0,0,1   l=1,-1,0
!
      xm=reshape((
     &   /0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.2041241452319d+00,-0.2041241452319d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.0000000000000d+00,-0.2041241452319d+00,
     &    0.0000000000000d+00, 0.4082482904639d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00,-0.2041241452319d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.4082482904639d+00, 0.4082482904639d+00,
     &    0.2041241452319d+00, 0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &   -0.2041241452319d+00, 0.0000000000000d+00, 0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.2041241452319d+00, 0.2041241452319d+00,
     &    0.0000000000000d+00,-0.4082482904639d+00, 0.4082482904639d+00,
     &    0.2041241452319d+00,-0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &   -0.2041241452319d+00, 0.0000000000000d+00,-0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.2041241452319d+00,-0.2041241452319d+00,
     &    0.0000000000000d+00, 0.4082482904639d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.0000000000000d+00, 0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.2041241452319d+00, 0.2041241452319d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.3535533905933d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00,-0.3535533905933d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.0000000000000d+00, 0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.0000000000000d+00,-0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.3535533905933d+00, 0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.3535533905933d+00,-0.3535533905933d+00/
     &    ),(/6,18/))
!
      xx=reshape((
     &   /0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.2041241452319d+00,-0.2041241452319d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.0000000000000d+00,-0.2041241452319d+00,
     &    0.0000000000000d+00, 0.4082482904639d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00,-0.2041241452319d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.4082482904639d+00, 0.4082482904639d+00,
     &    0.2041241452319d+00, 0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &   -0.2041241452319d+00, 0.0000000000000d+00, 0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.2041241452319d+00, 0.2041241452319d+00,
     &    0.0000000000000d+00,-0.4082482904639d+00, 0.4082482904639d+00,
     &    0.2041241452319d+00,-0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &   -0.2041241452319d+00, 0.0000000000000d+00,-0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.2041241452319d+00,-0.2041241452319d+00,
     &    0.0000000000000d+00, 0.4082482904639d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.2041241452319d+00, 0.0000000000000d+00,
     &    0.4082482904639d+00, 0.0000000000000d+00,-0.4082482904639d+00,
     &    0.2041241452319d+00, 0.0000000000000d+00, 0.2041241452319d+00,
     &    0.4082482904639d+00,-0.4082482904639d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00,-0.2041241452319d+00, 0.2041241452319d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.3535533905933d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00,-0.3535533905933d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.0000000000000d+00, 0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.3535533905933d+00, 0.0000000000000d+00,-0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.3535533905933d+00, 0.3535533905933d+00,
     &    0.0000000000000d+00, 0.0000000000000d+00, 0.0000000000000d+00,
     &    0.0000000000000d+00, 0.3535533905933d+00,-0.3535533905933d+00/
     &    ),(/6,18/))
!
c      pnewdt=-1.d0
!
!     elastic constants
!
      c1111=elconloc(1)
      c1122=elconloc(2)
      c1212=elconloc(3)
!
      if(iorien.gt.0) then
         call transformatrix(orab(1,iorien),pgauss,skl)
         do k=1,18
            do i=1,3
               do j=i,3
                  xmtran(i,j)=skl(i,1)*skl(j,1)*xx(1,k)+
     &                        skl(i,2)*skl(j,2)*xx(2,k)+
     &                        skl(i,3)*skl(j,3)*xx(3,k)+
     &                        (skl(i,1)*skl(j,2)+
     &                         skl(i,2)*skl(j,1))*xx(4,k)+
     &                        (skl(i,1)*skl(j,3)+
     &                         skl(i,3)*skl(j,1))*xx(5,k)+
     &                        (skl(i,2)*skl(j,3)+
     &                         skl(i,3)*skl(j,2))*xx(6,k)
               enddo
            enddo
            xm(1,k)=xmtran(1,1)
            xm(2,k)=xmtran(2,2)
            xm(3,k)=xmtran(3,3)
            xm(4,k)=xmtran(1,2)
            xm(5,k)=xmtran(1,3)
            xm(6,k)=xmtran(2,3)
         enddo
!
         elas( 1)=
     &      skl(1,1)*skl(1,1)*skl(1,1)*skl(1,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(1,2)*skl(1,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(1,3)*skl(1,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(1,1)*skl(1,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(1,2)*skl(1,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,1)*skl(1,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,3)*skl(1,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,1)*skl(1,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,2)*skl(1,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(1,1)*skl(1,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(1,2)*skl(1,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(1,3)*skl(1,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(1,2)*skl(1,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(1,3)*skl(1,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,1)*skl(1,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,3)*skl(1,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,2)*skl(1,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,3)*skl(1,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(1,1)*skl(1,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,2)*skl(1,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,3)*skl(1,3)*c1111
         elas( 2)=
     &      skl(1,1)*skl(1,1)*skl(2,1)*skl(2,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(2,2)*skl(2,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(2,3)*skl(2,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(2,1)*skl(2,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(2,2)*skl(2,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(2,1)*skl(2,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(2,3)*skl(2,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(2,1)*skl(2,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(2,2)*skl(2,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(2,1)*skl(2,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(2,2)*skl(2,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(2,3)*skl(2,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(2,2)*skl(2,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(2,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(2,1)*skl(2,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(2,3)*skl(2,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(2,2)*skl(2,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(2,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(2,1)*skl(2,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(2,2)*skl(2,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(2,3)*skl(2,3)*c1111
         elas( 3)=
     &      skl(2,1)*skl(2,1)*skl(2,1)*skl(2,1)*c1111+
     &      skl(2,1)*skl(2,1)*skl(2,2)*skl(2,2)*c1122+
     &      skl(2,1)*skl(2,1)*skl(2,3)*skl(2,3)*c1122+
     &      skl(2,1)*skl(2,2)*skl(2,1)*skl(2,2)*c1212+
     &      skl(2,1)*skl(2,2)*skl(2,2)*skl(2,1)*c1212+
     &      skl(2,1)*skl(2,3)*skl(2,1)*skl(2,3)*c1212+
     &      skl(2,1)*skl(2,3)*skl(2,3)*skl(2,1)*c1212+
     &      skl(2,2)*skl(2,1)*skl(2,1)*skl(2,2)*c1212+
     &      skl(2,2)*skl(2,1)*skl(2,2)*skl(2,1)*c1212+
     &      skl(2,2)*skl(2,2)*skl(2,1)*skl(2,1)*c1122+
     &      skl(2,2)*skl(2,2)*skl(2,2)*skl(2,2)*c1111+
     &      skl(2,2)*skl(2,2)*skl(2,3)*skl(2,3)*c1122+
     &      skl(2,2)*skl(2,3)*skl(2,2)*skl(2,3)*c1212+
     &      skl(2,2)*skl(2,3)*skl(2,3)*skl(2,2)*c1212+
     &      skl(2,3)*skl(2,1)*skl(2,1)*skl(2,3)*c1212+
     &      skl(2,3)*skl(2,1)*skl(2,3)*skl(2,1)*c1212+
     &      skl(2,3)*skl(2,2)*skl(2,2)*skl(2,3)*c1212+
     &      skl(2,3)*skl(2,2)*skl(2,3)*skl(2,2)*c1212+
     &      skl(2,3)*skl(2,3)*skl(2,1)*skl(2,1)*c1122+
     &      skl(2,3)*skl(2,3)*skl(2,2)*skl(2,2)*c1122+
     &      skl(2,3)*skl(2,3)*skl(2,3)*skl(2,3)*c1111
         elas( 4)=
     &      skl(1,1)*skl(1,1)*skl(3,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(3,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(3,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(3,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(3,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(3,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(3,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(3,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(3,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(3,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(3,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(3,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(3,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(3,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(3,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(3,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(3,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(3,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(3,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(3,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(3,3)*skl(3,3)*c1111
         elas( 5)=
     &      skl(2,1)*skl(2,1)*skl(3,1)*skl(3,1)*c1111+
     &      skl(2,1)*skl(2,1)*skl(3,2)*skl(3,2)*c1122+
     &      skl(2,1)*skl(2,1)*skl(3,3)*skl(3,3)*c1122+
     &      skl(2,1)*skl(2,2)*skl(3,1)*skl(3,2)*c1212+
     &      skl(2,1)*skl(2,2)*skl(3,2)*skl(3,1)*c1212+
     &      skl(2,1)*skl(2,3)*skl(3,1)*skl(3,3)*c1212+
     &      skl(2,1)*skl(2,3)*skl(3,3)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,1)*skl(3,1)*skl(3,2)*c1212+
     &      skl(2,2)*skl(2,1)*skl(3,2)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,2)*skl(3,1)*skl(3,1)*c1122+
     &      skl(2,2)*skl(2,2)*skl(3,2)*skl(3,2)*c1111+
     &      skl(2,2)*skl(2,2)*skl(3,3)*skl(3,3)*c1122+
     &      skl(2,2)*skl(2,3)*skl(3,2)*skl(3,3)*c1212+
     &      skl(2,2)*skl(2,3)*skl(3,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,1)*skl(3,1)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,1)*skl(3,3)*skl(3,1)*c1212+
     &      skl(2,3)*skl(2,2)*skl(3,2)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,2)*skl(3,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,3)*skl(3,1)*skl(3,1)*c1122+
     &      skl(2,3)*skl(2,3)*skl(3,2)*skl(3,2)*c1122+
     &      skl(2,3)*skl(2,3)*skl(3,3)*skl(3,3)*c1111
         elas( 6)=
     &      skl(3,1)*skl(3,1)*skl(3,1)*skl(3,1)*c1111+
     &      skl(3,1)*skl(3,1)*skl(3,2)*skl(3,2)*c1122+
     &      skl(3,1)*skl(3,1)*skl(3,3)*skl(3,3)*c1122+
     &      skl(3,1)*skl(3,2)*skl(3,1)*skl(3,2)*c1212+
     &      skl(3,1)*skl(3,2)*skl(3,2)*skl(3,1)*c1212+
     &      skl(3,1)*skl(3,3)*skl(3,1)*skl(3,3)*c1212+
     &      skl(3,1)*skl(3,3)*skl(3,3)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,1)*skl(3,1)*skl(3,2)*c1212+
     &      skl(3,2)*skl(3,1)*skl(3,2)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,2)*skl(3,1)*skl(3,1)*c1122+
     &      skl(3,2)*skl(3,2)*skl(3,2)*skl(3,2)*c1111+
     &      skl(3,2)*skl(3,2)*skl(3,3)*skl(3,3)*c1122+
     &      skl(3,2)*skl(3,3)*skl(3,2)*skl(3,3)*c1212+
     &      skl(3,2)*skl(3,3)*skl(3,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,1)*skl(3,1)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,1)*skl(3,3)*skl(3,1)*c1212+
     &      skl(3,3)*skl(3,2)*skl(3,2)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,2)*skl(3,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,3)*skl(3,1)*skl(3,1)*c1122+
     &      skl(3,3)*skl(3,3)*skl(3,2)*skl(3,2)*c1122+
     &      skl(3,3)*skl(3,3)*skl(3,3)*skl(3,3)*c1111
         elas( 7)=
     &      skl(1,1)*skl(1,1)*skl(1,1)*skl(2,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(1,2)*skl(2,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(1,3)*skl(2,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(1,1)*skl(2,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(1,2)*skl(2,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,1)*skl(2,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,3)*skl(2,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,1)*skl(2,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,2)*skl(2,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(1,1)*skl(2,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(1,2)*skl(2,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(1,3)*skl(2,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(1,2)*skl(2,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(1,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,1)*skl(2,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,3)*skl(2,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,2)*skl(2,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(1,1)*skl(2,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,2)*skl(2,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,3)*skl(2,3)*c1111
         elas( 8)=
     &      skl(2,1)*skl(2,1)*skl(1,1)*skl(2,1)*c1111+
     &      skl(2,1)*skl(2,1)*skl(1,2)*skl(2,2)*c1122+
     &      skl(2,1)*skl(2,1)*skl(1,3)*skl(2,3)*c1122+
     &      skl(2,1)*skl(2,2)*skl(1,1)*skl(2,2)*c1212+
     &      skl(2,1)*skl(2,2)*skl(1,2)*skl(2,1)*c1212+
     &      skl(2,1)*skl(2,3)*skl(1,1)*skl(2,3)*c1212+
     &      skl(2,1)*skl(2,3)*skl(1,3)*skl(2,1)*c1212+
     &      skl(2,2)*skl(2,1)*skl(1,1)*skl(2,2)*c1212+
     &      skl(2,2)*skl(2,1)*skl(1,2)*skl(2,1)*c1212+
     &      skl(2,2)*skl(2,2)*skl(1,1)*skl(2,1)*c1122+
     &      skl(2,2)*skl(2,2)*skl(1,2)*skl(2,2)*c1111+
     &      skl(2,2)*skl(2,2)*skl(1,3)*skl(2,3)*c1122+
     &      skl(2,2)*skl(2,3)*skl(1,2)*skl(2,3)*c1212+
     &      skl(2,2)*skl(2,3)*skl(1,3)*skl(2,2)*c1212+
     &      skl(2,3)*skl(2,1)*skl(1,1)*skl(2,3)*c1212+
     &      skl(2,3)*skl(2,1)*skl(1,3)*skl(2,1)*c1212+
     &      skl(2,3)*skl(2,2)*skl(1,2)*skl(2,3)*c1212+
     &      skl(2,3)*skl(2,2)*skl(1,3)*skl(2,2)*c1212+
     &      skl(2,3)*skl(2,3)*skl(1,1)*skl(2,1)*c1122+
     &      skl(2,3)*skl(2,3)*skl(1,2)*skl(2,2)*c1122+
     &      skl(2,3)*skl(2,3)*skl(1,3)*skl(2,3)*c1111
         elas( 9)=
     &      skl(3,1)*skl(3,1)*skl(1,1)*skl(2,1)*c1111+
     &      skl(3,1)*skl(3,1)*skl(1,2)*skl(2,2)*c1122+
     &      skl(3,1)*skl(3,1)*skl(1,3)*skl(2,3)*c1122+
     &      skl(3,1)*skl(3,2)*skl(1,1)*skl(2,2)*c1212+
     &      skl(3,1)*skl(3,2)*skl(1,2)*skl(2,1)*c1212+
     &      skl(3,1)*skl(3,3)*skl(1,1)*skl(2,3)*c1212+
     &      skl(3,1)*skl(3,3)*skl(1,3)*skl(2,1)*c1212+
     &      skl(3,2)*skl(3,1)*skl(1,1)*skl(2,2)*c1212+
     &      skl(3,2)*skl(3,1)*skl(1,2)*skl(2,1)*c1212+
     &      skl(3,2)*skl(3,2)*skl(1,1)*skl(2,1)*c1122+
     &      skl(3,2)*skl(3,2)*skl(1,2)*skl(2,2)*c1111+
     &      skl(3,2)*skl(3,2)*skl(1,3)*skl(2,3)*c1122+
     &      skl(3,2)*skl(3,3)*skl(1,2)*skl(2,3)*c1212+
     &      skl(3,2)*skl(3,3)*skl(1,3)*skl(2,2)*c1212+
     &      skl(3,3)*skl(3,1)*skl(1,1)*skl(2,3)*c1212+
     &      skl(3,3)*skl(3,1)*skl(1,3)*skl(2,1)*c1212+
     &      skl(3,3)*skl(3,2)*skl(1,2)*skl(2,3)*c1212+
     &      skl(3,3)*skl(3,2)*skl(1,3)*skl(2,2)*c1212+
     &      skl(3,3)*skl(3,3)*skl(1,1)*skl(2,1)*c1122+
     &      skl(3,3)*skl(3,3)*skl(1,2)*skl(2,2)*c1122+
     &      skl(3,3)*skl(3,3)*skl(1,3)*skl(2,3)*c1111
         elas(10)=
     &      skl(1,1)*skl(2,1)*skl(1,1)*skl(2,1)*c1111+
     &      skl(1,1)*skl(2,1)*skl(1,2)*skl(2,2)*c1122+
     &      skl(1,1)*skl(2,1)*skl(1,3)*skl(2,3)*c1122+
     &      skl(1,1)*skl(2,2)*skl(1,1)*skl(2,2)*c1212+
     &      skl(1,1)*skl(2,2)*skl(1,2)*skl(2,1)*c1212+
     &      skl(1,1)*skl(2,3)*skl(1,1)*skl(2,3)*c1212+
     &      skl(1,1)*skl(2,3)*skl(1,3)*skl(2,1)*c1212+
     &      skl(1,2)*skl(2,1)*skl(1,1)*skl(2,2)*c1212+
     &      skl(1,2)*skl(2,1)*skl(1,2)*skl(2,1)*c1212+
     &      skl(1,2)*skl(2,2)*skl(1,1)*skl(2,1)*c1122+
     &      skl(1,2)*skl(2,2)*skl(1,2)*skl(2,2)*c1111+
     &      skl(1,2)*skl(2,2)*skl(1,3)*skl(2,3)*c1122+
     &      skl(1,2)*skl(2,3)*skl(1,2)*skl(2,3)*c1212+
     &      skl(1,2)*skl(2,3)*skl(1,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(2,1)*skl(1,1)*skl(2,3)*c1212+
     &      skl(1,3)*skl(2,1)*skl(1,3)*skl(2,1)*c1212+
     &      skl(1,3)*skl(2,2)*skl(1,2)*skl(2,3)*c1212+
     &      skl(1,3)*skl(2,2)*skl(1,3)*skl(2,2)*c1212+
     &      skl(1,3)*skl(2,3)*skl(1,1)*skl(2,1)*c1122+
     &      skl(1,3)*skl(2,3)*skl(1,2)*skl(2,2)*c1122+
     &      skl(1,3)*skl(2,3)*skl(1,3)*skl(2,3)*c1111
         elas(11)=
     &      skl(1,1)*skl(1,1)*skl(1,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(1,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(1,3)*skl(3,3)*c1111
         elas(12)=
     &      skl(2,1)*skl(2,1)*skl(1,1)*skl(3,1)*c1111+
     &      skl(2,1)*skl(2,1)*skl(1,2)*skl(3,2)*c1122+
     &      skl(2,1)*skl(2,1)*skl(1,3)*skl(3,3)*c1122+
     &      skl(2,1)*skl(2,2)*skl(1,1)*skl(3,2)*c1212+
     &      skl(2,1)*skl(2,2)*skl(1,2)*skl(3,1)*c1212+
     &      skl(2,1)*skl(2,3)*skl(1,1)*skl(3,3)*c1212+
     &      skl(2,1)*skl(2,3)*skl(1,3)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,1)*skl(1,1)*skl(3,2)*c1212+
     &      skl(2,2)*skl(2,1)*skl(1,2)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,2)*skl(1,1)*skl(3,1)*c1122+
     &      skl(2,2)*skl(2,2)*skl(1,2)*skl(3,2)*c1111+
     &      skl(2,2)*skl(2,2)*skl(1,3)*skl(3,3)*c1122+
     &      skl(2,2)*skl(2,3)*skl(1,2)*skl(3,3)*c1212+
     &      skl(2,2)*skl(2,3)*skl(1,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,1)*skl(1,1)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,1)*skl(1,3)*skl(3,1)*c1212+
     &      skl(2,3)*skl(2,2)*skl(1,2)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,2)*skl(1,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,3)*skl(1,1)*skl(3,1)*c1122+
     &      skl(2,3)*skl(2,3)*skl(1,2)*skl(3,2)*c1122+
     &      skl(2,3)*skl(2,3)*skl(1,3)*skl(3,3)*c1111
         elas(13)=
     &      skl(3,1)*skl(3,1)*skl(1,1)*skl(3,1)*c1111+
     &      skl(3,1)*skl(3,1)*skl(1,2)*skl(3,2)*c1122+
     &      skl(3,1)*skl(3,1)*skl(1,3)*skl(3,3)*c1122+
     &      skl(3,1)*skl(3,2)*skl(1,1)*skl(3,2)*c1212+
     &      skl(3,1)*skl(3,2)*skl(1,2)*skl(3,1)*c1212+
     &      skl(3,1)*skl(3,3)*skl(1,1)*skl(3,3)*c1212+
     &      skl(3,1)*skl(3,3)*skl(1,3)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,1)*skl(1,1)*skl(3,2)*c1212+
     &      skl(3,2)*skl(3,1)*skl(1,2)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,2)*skl(1,1)*skl(3,1)*c1122+
     &      skl(3,2)*skl(3,2)*skl(1,2)*skl(3,2)*c1111+
     &      skl(3,2)*skl(3,2)*skl(1,3)*skl(3,3)*c1122+
     &      skl(3,2)*skl(3,3)*skl(1,2)*skl(3,3)*c1212+
     &      skl(3,2)*skl(3,3)*skl(1,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,1)*skl(1,1)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,1)*skl(1,3)*skl(3,1)*c1212+
     &      skl(3,3)*skl(3,2)*skl(1,2)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,2)*skl(1,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,3)*skl(1,1)*skl(3,1)*c1122+
     &      skl(3,3)*skl(3,3)*skl(1,2)*skl(3,2)*c1122+
     &      skl(3,3)*skl(3,3)*skl(1,3)*skl(3,3)*c1111
         elas(14)=
     &      skl(1,1)*skl(2,1)*skl(1,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(2,1)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(2,1)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(2,2)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(2,2)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(2,3)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(2,3)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(2,1)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(2,1)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(2,2)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(2,2)*skl(1,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(2,2)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(2,3)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(2,3)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(2,1)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(2,1)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(2,2)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(2,2)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(2,3)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(2,3)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(2,3)*skl(1,3)*skl(3,3)*c1111
         elas(15)=
     &      skl(1,1)*skl(3,1)*skl(1,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(3,1)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(3,1)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(3,2)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(3,2)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(3,3)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(3,3)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(3,1)*skl(1,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(3,1)*skl(1,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(3,2)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(3,2)*skl(1,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(3,2)*skl(1,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(3,3)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(3,3)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(3,1)*skl(1,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(3,1)*skl(1,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(3,2)*skl(1,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(3,2)*skl(1,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(3,3)*skl(1,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(3,3)*skl(1,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(3,3)*skl(1,3)*skl(3,3)*c1111
         elas(16)=
     &      skl(1,1)*skl(1,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(1,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(1,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(1,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(1,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(1,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(1,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(1,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(1,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(1,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(1,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(1,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(1,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(1,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(1,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(1,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(1,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(1,3)*skl(2,3)*skl(3,3)*c1111
         elas(17)=
     &      skl(2,1)*skl(2,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(2,1)*skl(2,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(2,1)*skl(2,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(2,1)*skl(2,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(2,1)*skl(2,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(2,1)*skl(2,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(2,1)*skl(2,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(2,2)*skl(2,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(2,2)*skl(2,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(2,2)*skl(2,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(2,2)*skl(2,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(2,2)*skl(2,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(2,2)*skl(2,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(2,3)*skl(2,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(2,3)*skl(2,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(2,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(2,3)*skl(2,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(2,3)*skl(2,3)*skl(2,3)*skl(3,3)*c1111
         elas(18)=
     &      skl(3,1)*skl(3,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(3,1)*skl(3,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(3,1)*skl(3,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(3,1)*skl(3,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(3,1)*skl(3,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(3,1)*skl(3,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(3,1)*skl(3,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(3,2)*skl(3,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(3,2)*skl(3,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(3,2)*skl(3,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(3,2)*skl(3,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(3,2)*skl(3,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(3,2)*skl(3,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(3,3)*skl(3,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(3,3)*skl(3,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(3,3)*skl(3,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(3,3)*skl(3,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(3,3)*skl(3,3)*skl(2,3)*skl(3,3)*c1111
         elas(19)=
     &      skl(1,1)*skl(2,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(2,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(2,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(2,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(2,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(2,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(2,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(2,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(2,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(2,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(2,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(2,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(2,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(2,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(2,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(2,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(2,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(2,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(2,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(2,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(2,3)*skl(2,3)*skl(3,3)*c1111
         elas(20)=
     &      skl(1,1)*skl(3,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(1,1)*skl(3,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,1)*skl(3,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,1)*skl(3,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,1)*skl(3,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,1)*skl(3,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,1)*skl(3,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,2)*skl(3,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(1,2)*skl(3,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(1,2)*skl(3,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,2)*skl(3,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(1,2)*skl(3,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(1,2)*skl(3,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,2)*skl(3,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(3,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(1,3)*skl(3,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(1,3)*skl(3,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(1,3)*skl(3,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(1,3)*skl(3,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(1,3)*skl(3,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(1,3)*skl(3,3)*skl(2,3)*skl(3,3)*c1111
         elas(21)=
     &      skl(2,1)*skl(3,1)*skl(2,1)*skl(3,1)*c1111+
     &      skl(2,1)*skl(3,1)*skl(2,2)*skl(3,2)*c1122+
     &      skl(2,1)*skl(3,1)*skl(2,3)*skl(3,3)*c1122+
     &      skl(2,1)*skl(3,2)*skl(2,1)*skl(3,2)*c1212+
     &      skl(2,1)*skl(3,2)*skl(2,2)*skl(3,1)*c1212+
     &      skl(2,1)*skl(3,3)*skl(2,1)*skl(3,3)*c1212+
     &      skl(2,1)*skl(3,3)*skl(2,3)*skl(3,1)*c1212+
     &      skl(2,2)*skl(3,1)*skl(2,1)*skl(3,2)*c1212+
     &      skl(2,2)*skl(3,1)*skl(2,2)*skl(3,1)*c1212+
     &      skl(2,2)*skl(3,2)*skl(2,1)*skl(3,1)*c1122+
     &      skl(2,2)*skl(3,2)*skl(2,2)*skl(3,2)*c1111+
     &      skl(2,2)*skl(3,2)*skl(2,3)*skl(3,3)*c1122+
     &      skl(2,2)*skl(3,3)*skl(2,2)*skl(3,3)*c1212+
     &      skl(2,2)*skl(3,3)*skl(2,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(3,1)*skl(2,1)*skl(3,3)*c1212+
     &      skl(2,3)*skl(3,1)*skl(2,3)*skl(3,1)*c1212+
     &      skl(2,3)*skl(3,2)*skl(2,2)*skl(3,3)*c1212+
     &      skl(2,3)*skl(3,2)*skl(2,3)*skl(3,2)*c1212+
     &      skl(2,3)*skl(3,3)*skl(2,1)*skl(3,1)*c1122+
     &      skl(2,3)*skl(3,3)*skl(2,2)*skl(3,2)*c1122+
     &      skl(2,3)*skl(3,3)*skl(2,3)*skl(3,3)*c1111
      endif
!
      do i=1,6
         ep0(i)=xstateini(i,iint,iel)
      enddo
!
!     elastic strains
!
      do i=1,6
         ee(i)=emec(i)-ep0(i)
      enddo
!
!     (visco)plastic constants: octahedral slip system
!
      do i=1,12
         ck(i)=elconloc(4)
         cn(i)=elconloc(5)
      enddo
!
!     (visco)plastic constants: cubic slip system
!
      do i=13,18
         ck(i)=elconloc(6)
         cn(i)=elconloc(7)
      enddo
!
!     global trial stress tensor
!
      if(iorien.gt.0) then
         stri(1)=elas(1)*ee(1)+elas(2)*ee(2)+elas(4)*ee(3)+
     &     2.d0*(elas(7)*ee(4)+elas(11)*ee(5)+elas(16)*ee(6))
     &     -beta(1)
         stri(2)=elas(2)*ee(1)+elas(3)*ee(2)+elas(5)*ee(3)+
     &     2.d0*(elas(8)*ee(4)+elas(12)*ee(5)+elas(17)*ee(6))
     &     -beta(2)
         stri(3)=elas(4)*ee(1)+elas(5)*ee(2)+elas(6)*ee(3)+
     &     2.d0*(elas(9)*ee(4)+elas(13)*ee(5)+elas(18)*ee(6))
     &     -beta(3)
         stri(4)=elas(7)*ee(1)+elas(8)*ee(2)+elas(9)*ee(3)+
     &     2.d0*(elas(10)*ee(4)+elas(14)*ee(5)+elas(19)*ee(6))
     &     -beta(4)
         stri(5)=elas(11)*ee(1)+elas(12)*ee(2)+elas(13)*ee(3)+
     &     2.d0*(elas(14)*ee(4)+elas(15)*ee(5)+elas(20)*ee(6))
     &     -beta(5)
         stri(6)=elas(16)*ee(1)+elas(17)*ee(2)+elas(18)*ee(3)+
     &     2.d0*(elas(19)*ee(4)+elas(20)*ee(5)+elas(21)*ee(6))
     &     -beta(6)
      else
         stri(1)=c1111*ee(1)+c1122*(ee(2)+ee(3))-beta(1)
         stri(2)=c1111*ee(2)+c1122*(ee(1)+ee(3))-beta(2)
         stri(3)=c1111*ee(3)+c1122*(ee(1)+ee(2))-beta(3)
         stri(4)=2.d0*c1212*ee(4)-beta(4)
         stri(5)=2.d0*c1212*ee(5)-beta(5)
         stri(6)=2.d0*c1212*ee(6)-beta(6)
      endif
!
!     stress radius in each slip plane
!
      do i=1,18
         sg(i)=xm(1,i)*stri(1)+xm(2,i)*stri(2)+xm(3,i)*stri(3)+
     &      2.d0*(xm(4,i)*stri(4)+xm(5,i)*stri(5)+xm(6,i)*stri(6))
      enddo
!
!     evaluation of the yield surface
!
      do i=1,18
         htri(i)=dabs(sg(i))
      enddo
!
      if(ielas.eq.1) then
!
!        elastic stress
!
         do i=1,6
            stre(i)=stri(i)
         enddo
!
!        elastic stiffness
!
         if(icmd.ne.3) then
            if(iorien.gt.0) then
               do i=1,21
                  stiff(i)=elas(i)
               enddo
            else
               stiff(1)=c1111
               stiff(2)=c1122
               stiff(3)=c1111
               stiff(4)=c1122
               stiff(5)=c1122
               stiff(6)=c1111
               stiff(7)=0.d0
               stiff(8)=0.d0
               stiff(9)=0.d0
               stiff(10)=c1212
               stiff(11)=0.d0
               stiff(12)=0.d0
               stiff(13)=0.d0
               stiff(14)=0.d0
               stiff(15)=c1212
               stiff(16)=0.d0
               stiff(17)=0.d0
               stiff(18)=0.d0
               stiff(19)=0.d0
               stiff(20)=0.d0
               stiff(21)=c1212
            endif
!
!           updating the state variables
!     
            do i=1,6
               xstate(i,iint,iel)=ep0(i)
            enddo
            do i=7,24
               xstate(i,iint,iel)=xstateini(i,iint,iel)
            enddo
         endif
!
         return
      endif
!
!     plastic deformation
!
      nrhs=1
      lda=18
      ldb=18
!
!     initializing the state variables
!
      do i=1,6
         ep(i)=ep0(i)
      enddo
      do i=1,18
         dg(i)=0.d0
c         dg(i)=dtime*htri(i)**cn(i)/ck(i)**cn(i)
      enddo
!
!     major loop
!
      icounter=0
      do
         icounter=icounter+1
         if(icounter.gt.100) then
            pnewdt=0.25d0
            return
c            write(*,*) 
c     &       '*ERROR in umat_single_crystal_creep: no convergence'
c            call exit(201)
         endif
!     
!     elastic strains
!     
         do i=1,6
            ee(i)=emec(i)-ep(i)
c         write(*,*) 'umat_single_crystal_creep ee',icounter,i,ee(i)
         enddo
!     
!     global trial stress tensor
!     
         if(iorien.gt.0) then
            stri(1)=elas(1)*ee(1)+elas(2)*ee(2)+elas(4)*ee(3)+
     &           2.d0*(elas(7)*ee(4)+elas(11)*ee(5)+elas(16)*ee(6))
     &           -beta(1)
            stri(2)=elas(2)*ee(1)+elas(3)*ee(2)+elas(5)*ee(3)+
     &           2.d0*(elas(8)*ee(4)+elas(12)*ee(5)+elas(17)*ee(6))
     &           -beta(2)
            stri(3)=elas(4)*ee(1)+elas(5)*ee(2)+elas(6)*ee(3)+
     &           2.d0*(elas(9)*ee(4)+elas(13)*ee(5)+elas(18)*ee(6))
     &           -beta(3)
            stri(4)=elas(7)*ee(1)+elas(8)*ee(2)+elas(9)*ee(3)+
     &           2.d0*(elas(10)*ee(4)+elas(14)*ee(5)+elas(19)*ee(6))
     &           -beta(4)
            stri(5)=elas(11)*ee(1)+elas(12)*ee(2)+elas(13)*ee(3)+
     &           2.d0*(elas(14)*ee(4)+elas(15)*ee(5)+elas(20)*ee(6))
     &           -beta(5)
            stri(6)=elas(16)*ee(1)+elas(17)*ee(2)+elas(18)*ee(3)+
     &           2.d0*(elas(19)*ee(4)+elas(20)*ee(5)+elas(21)*ee(6))
     &           -beta(6)
         else
            stri(1)=c1111*ee(1)+c1122*(ee(2)+ee(3))-beta(1)
            stri(2)=c1111*ee(2)+c1122*(ee(1)+ee(3))-beta(2)
            stri(3)=c1111*ee(3)+c1122*(ee(1)+ee(2))-beta(3)
            stri(4)=2.d0*c1212*ee(4)-beta(4)
            stri(5)=2.d0*c1212*ee(5)-beta(5)
            stri(6)=2.d0*c1212*ee(6)-beta(6)
         endif
         do i=1,6
c         write(*,*) 'umat_single_crystal_creep stri',icounter,i,stri(i)
         enddo
!     
!     stress radius in each slip plane
!     
         do i=1,18
            sg(i)=xm(1,i)*stri(1)+xm(2,i)*stri(2)+xm(3,i)*stri(3)+
     &           2.d0*(xm(4,i)*stri(4)+xm(5,i)*stri(5)+xm(6,i)*stri(6))
c         write(*,*) 'umat_single_crystal_creep sg',icounter,i,sg(i)
         enddo
!     
!     evaluation of the yield surface
!     
         do i=1,18
c         write(*,*) 'umat_single_crystal_creep ck,dtime,cn,dg',icounter,
c     &         i,ck(i),dtime,cn(i),dg(i)
            htri(i)=dabs(sg(i))-ck(i)*(dg(i)/dtime)**(1.d0/cn(i))
c         write(*,*) 'umat_single_crystal_creep htri',icounter,i,htri(i)
         enddo
!     
!     replace sg(i) by sgn(sg(i))
!     
         do i=1,18
            if(sg(i).lt.0.d0) then
               sg(i)=-1.d0
            else
               sg(i)=1.d0
            endif
         enddo
!     
!     determining the residual matrix
!     
         do i=1,6  
            r(i)=ep0(i)-ep(i)
         enddo
         do i=1,18
            do j=1,6
               r(j)=r(j)+xm(j,i)*sg(i)*dg(i)
            enddo
         enddo
!     
!     check convergence
!  
         if(icounter.gt.1) then
            convergence=1
            do i=1,18
c               if(dabs(htri(i)).gt.1.d-10) then
               if(dabs(htri(i)).gt.1.d-3) then
c                  write(*,*) 'divergence htri dg ',
c     &                    icounter,i,htri(i),dg(i)
                  convergence=0
                  exit
               endif
            enddo
            if(convergence.eq.1) exit
         endif
c         do i=1,18
c            if(dabs(htri(i)).gt.1.d-5) then
c               write(*,*) 'htri ',i,htri(i)
c               convergence=0
c               exit
c            endif
c         enddo
c         if(convergence.eq.1) then
c            dd=0.d0
c            do i=1,6
c               dd=dd+r(i)*r(i)
c            enddo
c            dd=sqrt(dd)
c            if(dd.gt.1.d-10) then
c               convergence=0
c               write(*,*) 'dd ',dd
c            else
c               exit
c            endif
c         endif
!     
!     compute xmc=c:xm
!     
         do i=1,18
            if(iorien.gt.0) then
               xmc(1,i)=elas(1)*xm(1,i)+elas(2)*xm(2,i)+
     &              elas(4)*xm(3,i)+2.d0*(elas(7)*xm(4,i)+
     &              elas(11)*xm(5,i)+elas(16)*xm(6,i))
               xmc(2,i)=elas(2)*xm(1,i)+elas(3)*xm(2,i)+
     &              elas(5)*xm(3,i)+2.d0*(elas(8)*xm(4,i)+
     &              elas(12)*xm(5,i)+elas(17)*xm(6,i))
               xmc(3,i)=elas(4)*xm(1,i)+elas(5)*xm(2,i)+
     &              elas(6)*xm(3,i)+2.d0*(elas(9)*xm(4,i)+
     &              elas(13)*xm(5,i)+elas(18)*xm(6,i))
               xmc(4,i)=elas(7)*xm(1,i)+elas(8)*xm(2,i)+
     &              elas(9)*xm(3,i)+2.d0*(elas(10)*xm(4,i)+
     &              elas(14)*xm(5,i)+elas(19)*xm(6,i))
               xmc(5,i)=elas(11)*xm(1,i)+elas(12)*xm(2,i)+
     &              elas(13)*xm(3,i)+2.d0*(elas(14)*xm(4,i)+
     &              elas(15)*xm(5,i)+elas(20)*xm(6,i))
               xmc(6,i)=elas(16)*xm(1,i)+elas(17)*xm(2,i)+
     &              elas(18)*xm(3,i)+2.d0*(elas(19)*xm(4,i)+
     &              elas(20)*xm(5,i)+elas(21)*xm(6,i))
            else
               xmc(1,i)=c1111*xm(1,i)+c1122*(xm(2,i)+xm(3,i))
               xmc(2,i)=c1111*xm(2,i)+c1122*(xm(1,i)+xm(3,i))
               xmc(3,i)=c1111*xm(3,i)+c1122*(xm(1,i)+xm(2,i))
               xmc(4,i)=2.d0*c1212*xm(4,i)
               xmc(5,i)=2.d0*c1212*xm(5,i)
               xmc(6,i)=2.d0*c1212*xm(6,i)
            endif
         enddo
!     
         neq=18
!     
!     filling the LHS
!     
         do i=1,18
            do j=1,18
               if(i.ne.j) then
                  gl(i,j)=(xm(1,i)*xmc(1,j)+
     &                 xm(2,i)*xmc(2,j)+xm(3,i)*xmc(3,j)+2.d0*
     &                 (xm(4,i)*xmc(4,j)+xm(5,i)*xmc(5,j)+
     &                 xm(6,i)*xmc(6,j)))
     &                 *sg(i)*sg(j)
               else
                  gl(i,j)=(xm(1,i)*xmc(1,j)+
     &                 xm(2,i)*xmc(2,j)+xm(3,i)*xmc(3,j)+2.d0*
     &                 (xm(4,i)*xmc(4,j)+xm(5,i)*xmc(5,j)+
     &                 xm(6,i)*xmc(6,j)))
               endif
            enddo
            if(dg(i).gt.1.d-30) then
               gl(i,i)=gl(i,i)+
     &              (dg(i)/dtime)**(1.d0/cn(i)-1.d0)*ck(i)/
     &              (cn(i)*dtime)
            else
!     
!     for gamma ein default of 1.d-30 is taken to
!     obtain a finite gradient
!     
               gl(i,i)=gl(i,i)+
     &              (1.d-30/dtime)**(1.d0/cn(i)-1.d0)*ck(i)/
     &              (cn(i)*dtime)
            endif
         enddo
!     
!     filling the RHS
!     
         do i=1,18
!
            gr(i,1)=htri(i)
cccccc     &           -ck(i)*(dg(i)/dtime)**(1.d0/cn(i))
!
            do j=1,6
               t(j)=xmc(j,i)*sg(i)
            enddo
            do j=1,6
               gr(i,1)=gr(i,1)-t(j)*r(j)
            enddo
            gr(i,1)=gr(i,1)
     &           -t(4)*r(4)-t(5)*r(5)-t(6)*r(6)
         enddo
!     
!     solve gl*ddg=gr
!     
         call dgesv(neq,nrhs,gl,lda,ipiv,gr,ldb,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in umat_single_crystal_creep:'
            write(*,*) '       linear equation solver exited'
            write(*,*) '       with error: info = ',info
            call exit(201)
         endif
!     
         do i=1,18
            ddg(i)=gr(i,1)
         enddo
!     
!     updating the plastic strain
!     
         do j=1,6
            ep(j)=ep(j)+r(j)
            do i=1,18
               ep(j)=ep(j)+xm(j,i)*sg(i)*ddg(i)
            enddo
c         write(*,*) 'umat_single_crystal_creep ep',icounter,j,ep(j)
         enddo
!
!        updating the consistency parameter
!
         do i=1,18
            dg(i)=max(dg(i)+ddg(i),0.d0)
c         write(*,*) 'umat_single_crystal_creep ddg',icounter,i,ddg(i)
         enddo
c         write(*,*) 'umat_single_crystal_creep ',icounter,dg(2)
!     
!     end of major loop
!     
      enddo
!     
!     inversion of G
!     
      do i=1,neq
         do j=1,neq
            gr(i,j)=0.d0
         enddo
         gr(i,i)=1.d0
      enddo
      nrhs=neq
      call dgetrs('No transpose',neq,nrhs,gl,lda,ipiv,gr,ldb,info)
      if(info.ne.0) then
         write(*,*) '*ERROR in sc.f: linear equation solver'
         write(*,*) '       exited with error: info = ',info
         call exit(201)
      endif
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
         if(iorien.gt.0) then
            ddsdde(1,1)=elas(1)
            ddsdde(1,2)=elas(2)
            ddsdde(1,3)=elas(4)
            ddsdde(1,4)=elas(7)
            ddsdde(1,5)=elas(11)
            ddsdde(1,6)=elas(16)
c            ddsdde(2,1)=elas(2)
            ddsdde(2,2)=elas(3)
            ddsdde(2,3)=elas(5)
            ddsdde(2,4)=elas(8)
            ddsdde(2,5)=elas(12)
            ddsdde(2,6)=elas(17)
c            ddsdde(3,1)=elas(4)
c            ddsdde(3,2)=elas(5)
            ddsdde(3,3)=elas(6)
            ddsdde(3,4)=elas(9)
            ddsdde(3,5)=elas(13)
            ddsdde(3,6)=elas(18)
c            ddsdde(4,1)=elas(7)
c            ddsdde(4,2)=elas(8)
c            ddsdde(4,3)=elas(9)
            ddsdde(4,4)=elas(10)
            ddsdde(4,5)=elas(14)
            ddsdde(4,6)=elas(19)
c            ddsdde(5,1)=elas(11)
c            ddsdde(5,2)=elas(12)
c            ddsdde(5,3)=elas(13)
c            ddsdde(5,4)=elas(14)
            ddsdde(5,5)=elas(15)
            ddsdde(5,6)=elas(20)
c            ddsdde(6,1)=elas(16)
c            ddsdde(6,2)=elas(17)
c            ddsdde(6,3)=elas(18)
c            ddsdde(6,4)=elas(19)
c            ddsdde(6,5)=elas(20)
            ddsdde(6,6)=elas(21)
         else
            do i=1,6
c               do j=1,6
               do j=i,6
                  ddsdde(i,j)=0.d0
               enddo
            enddo
            do i=1,3
               ddsdde(i,i)=c1111
            enddo
            do i=1,3
               do j=i+1,3
                  ddsdde(i,j)=c1122
               enddo
c               do j=1,i-1
c                  ddsdde(i,j)=c1122
c               enddo
               ddsdde(i+3,i+3)=c1212
            enddo
         endif
         do i=1,18
            do j=1,18
               do k=1,6
c                  do l=1,6
                  do l=k,6
                     ddsdde(k,l)=ddsdde(k,l)-
     &               gr(i,j)*xmc(k,i)*sg(i)*xmc(l,j)*sg(j)
                  enddo
               enddo
            enddo
         enddo
!
!        symmatrizing the stiffness matrix: to be changed: the
!        matrix is symmetric!
!
         stiff(1)=ddsdde(1,1)
         stiff(2)=ddsdde(1,2)
         stiff(3)=ddsdde(2,2)
         stiff(4)=ddsdde(1,3)
         stiff(5)=ddsdde(2,3)
         stiff(6)=ddsdde(3,3)
         stiff(7)=ddsdde(1,4)
         stiff(8)=ddsdde(2,4)
         stiff(9)=ddsdde(3,4)
         stiff(10)=ddsdde(4,4)
         stiff(11)=ddsdde(1,5)
         stiff(12)=ddsdde(2,5)
         stiff(13)=ddsdde(3,5)
         stiff(14)=ddsdde(4,5)
         stiff(15)=ddsdde(5,5)
         stiff(16)=ddsdde(1,6)
         stiff(17)=ddsdde(2,6)
         stiff(18)=ddsdde(3,6)
         stiff(19)=ddsdde(4,6)
         stiff(20)=ddsdde(5,6)
         stiff(21)=ddsdde(6,6)
!
      endif
!
!     updating the state variables
!
      do i=1,6
         xstate(i,iint,iel)=ep(i)
      enddo
      do i=7,24
         xstate(i,iint,iel)=xstateini(i,iint,iel)+dg(i-6)
      enddo
!
      return
      end
