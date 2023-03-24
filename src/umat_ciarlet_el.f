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
!     Author: Sven Kassbohm
!     Permission to use this routine was given by the author on
!     January 17, 2015.
!
      subroutine umat_ciarlet_el(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
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
!     ttime              total time at the start of the current increment
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
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!     
!     Ela-Pla-modifications:
!
      real*8 C1(3,3), C1i(3,3), detC1, dS2dE(3,3,3,3), S2PK1(3,3)
      real*8 lamda, mu,un,e
!
!     change on 190315: input are Young's modulus and Poisson's
!     coefficient instead of the Lame constants
!
      e=elconloc(1)
      un=elconloc(2)
      mu=e/(1.d0+un)
      lamda=mu*un/(1.d0-2.d0*un)
      mu=mu/2.d0

c      lamda=elconloc(1)
c      mu=elconloc(2)

C     ciarlet_doc
C
C     Compute pullback C1(3,3) of metric tensor g.
c     C1  = matmul(transpose(xkl),xkl)


!
!     change on 190315: the Cauchy tensor should be
!     computed based on the mechanical Lagrange strain in 
!     order to discard the thermal strains; the deformation
!     gradient still contains the thermal strain effect
!      
!     right Cauchy-Green tensor (eloc contains the Lagrange strain,
!     including thermal strain)
!
      C1(1,1)=2.d0*emec(1)+1.d0
      C1(2,2)=2.d0*emec(2)+1.d0
      C1(3,3)=2.d0*emec(3)+1.d0
      C1(1,2)=2.d0*emec(4)
      C1(2,1)=C1(1,2)
      C1(1,3)=2.d0*emec(5)
      C1(3,1)=C1(1,3)
      C1(2,3)=2.d0*emec(6)
      C1(3,2)=C1(2,3)
C
C     2PK-Stress at C:
      call S2_Ciarlet(C1,lamda,mu, S2PK1,C1i,detC1)
C     Tangent/linearization:
      call dS2_Ciarlet_dE(lamda,mu,detC1,C1i,dS2dE)

C     Dhondt-Voigt-Order, see umat_lin_iso_el.f
      call voigt_33mat_to_6vec(S2PK1,stre)
      call voigt_tetrad_to_matrix(dS2dE,stiff)

C     ciarlet_cod

      return
      end



      subroutine voigt_33mat_to_6vec(mat,vec)
      implicit none
C     in:
      double precision mat(3,3)
C     out:
      double precision vec(6)
C     11,22,33,12,13,23
      vec(1)=mat(1,1)
      vec(2)=mat(2,2)
      vec(3)=mat(3,3)
      vec(4)=mat(1,2)
      vec(5)=mat(1,3)
      vec(6)=mat(2,3)

      return
      end


      subroutine voigt_tetrad_to_matrix(C,stiff)
      implicit none
C     in:
      double precision C(3,3,3,3)
C     out:
      double precision stiff(21)
      
C     indices according to 
C     1) file anisotropic.f and 
C     2) html-documentation
C     -> Making C symmetric in the **last two** indices:
C     
C     stiff stores:
C     1111, 1122, 2222, 1133
C     2233, 3333, 1112, 2212
C
C     1111
      stiff(1)  = C(1,1,1,1)
C     1122 = 2211
      stiff(2)  = C(1,1,2,2)
C     2222
      stiff(3)  = C(2,2,2,2)
C     1133 = 3311
      stiff(4)  = C(1,1,3,3)
C     2233 = 3322
      stiff(5)  = C(2,2,3,3)
C     3333
      stiff(6)  = C(3,3,3,3)
C     1112 = 1121 = 1211 = 2111
      stiff(7)  = 0.5d0*(C(1,1,1,2)+C(1,1,2,1))
C     2212 = 1222 = 2122 = 2221
      stiff(8)  = 0.5d0*(C(2,2,1,2)+C(2,2,2,1))
C
C     stiff stores:
C     3312, 1212, 1113, 2213
C     3313, 1213, 1313, 1123
C     
C     3312 = 1233 = 2133 = 3321
      stiff(9)  = 0.5d0*(C(3,3,1,2)+C(3,3,2,1))
C     1212 = 1221 = 2112 = 2121
      stiff(10) = 0.5d0*(C(1,2,1,2)+C(1,2,2,1))
C     1113 = 1131 = 1311 = 3111
      stiff(11) = 0.5d0*(C(1,1,1,3)+C(1,1,3,1))
C     2213 = 1322 = 2231 = 3122
      stiff(12) = 0.5d0*(C(2,2,1,3)+C(2,2,3,1))
C     
C     3313 = 1333 = 3133 = 3331
      stiff(13) = 0.5d0*(C(3,3,1,3)+C(3,3,3,1))     
C     1213 = 1231 = 1312 = 1321 = 2113 = 2131 = 3112 = 3121:
      stiff(14) = 0.5d0*(C(1,2,1,3)+C(1,2,3,1))
C     1313 = 1331 = 3113 = 3131
      stiff(15) = 0.5d0*(C(1,3,1,3)+C(1,3,3,1))
C     1123 = 1132 = 2311 = 3211
      stiff(16) = 0.5d0*(C(1,1,2,3)+C(1,1,3,2))
C
C     stiff stores:
C     2223, 3323, 1223, 1323
C     2323
C
C     2223 = 2232 = 2322 = 3222
      stiff(17) = 0.5d0*(C(2,2,2,3)+C(2,2,3,2)) 
C     3323 = 2333 = 3233 = 3332
      stiff(18) = 0.5d0*(C(3,3,2,3)+C(3,3,3,2))
C     1223 = 1232 = 2123 = 2132 = 2312 = 2321 = 3212 = 3221
      stiff(19) = 0.5d0*(C(1,2,2,3)+C(1,2,3,2))
C     1323 = 1332 = 2313 = 2331 = 3123 = 3132 = 3213 = 3231
      stiff(20) = 0.5d0*(C(1,3,2,3)+C(1,3,3,2))
C     2323 = 2332
      stiff(21) = 0.5d0*(C(2,3,2,3)+C(2,3,3,2))

      return
      end
      
      subroutine S2_Ciarlet(Cb,lamda,mu, S2,Cbi,detC)
C     computes 2PK stresses for Ciarlet-model
C     
      implicit none
C     input:
      double precision Cb(3,3),lamda,mu
C     output:
      double precision S2(3,3), Cbi(3,3), detC
C     local:
      double precision id(3,3)
      double precision f1
      integer I,J
      logical OK_FLAG

C     Set id:
      do I=1,3
         do J=1,3
            id(I,J)=0.d0
         enddo
         id(I,I)=1.d0
      enddo

C     Inverse of Cb and det(C):
      call M33INV_DET(Cb, Cbi, detC, OK_FLAG)
C
C
      f1 = 0.5d0*lamda*(detC-1.d0)
      do I=1,3
         do J=1,3
            S2(I,J)=f1*Cbi(I,J) + mu * ( id(I,J)-Cbi(I,J) )
         enddo
      enddo
      return
      end

      subroutine dS2_Ciarlet_dE(lamda,mu,detC,Cbi,A)
C     computes linearization of 2PK stresses with E
C     for Ciarlet-model
      implicit none
C     input:
      double precision lamda,mu,detC,Cbi(3,3)
C     output:
      double precision A(3,3,3,3)
C     local:
      double precision f1,f2
      integer I,J,K,L
      
      f1 = lamda*detC
      f2 = 2.d0*mu - lamda*(detC-1.d0)
      do I=1,3
         do J=1,3
            do K=1,3
               do L=1,3
                  A(I,J,K,L) =  
     &                 f1 * Cbi(K,L)*Cbi(I,J) +
     &                 f2 * Cbi(I,K)*Cbi(L,J)     
               enddo
            enddo
         enddo
      enddo
      return
      end

      subroutine M33INV_DET(A, AINV, DET, OK_FLAG)
C     
C
C     !******************************************************************
C     !  M33INV_DET  -  Compute the inverse of a 3x3 matrix.
C     !
C     !  A       = input 3x3 matrix to be inverted
C     !  AINV    = output 3x3 inverse of matrix A
C     !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, 
C     !            and .FALSE. if the input matrix is singular.
C     !******************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)  
     &      - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) 
     &      + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN
      END


