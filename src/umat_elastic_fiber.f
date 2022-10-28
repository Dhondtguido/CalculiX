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
      subroutine umat_elastic_fiber
     &       (amat,iel,iint,kode,elconloc,emec,emec0,
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
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),nfiber,
     &  i,
     &  j,k,l,m,n,ioffset,nt,kk(84),iorien
!
      real*8 elconloc(*),stiff(21),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,c(3,3),a(3),pgauss(3),
     &  orab(7,*),skl(3,3),aa(3),emec(6),time,ttime
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &   constant(21),dd,dm(3,3,4),djdc(3,3,4),d2jdc2(3,3,3,3,4),
     &   v1,v1b,v3,v3bi,v4(4),v4br(4),djbdc(3,3,4),d2jbdc2(3,3,3,3,4),
     &   didc(3,3,3),d2idc2(3,3,3,3,3),dibdc(3,3,3),d2ibdc2(3,3,3,3,3),
     &   dudc(3,3),d2udc2(3,3,3,3),v33,cinv(3,3),xk1,xk2,d(3,3),term
!
      kk=(/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &  1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,3,3,1,3,
     &  1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3,2,3,
     &  2,3,2,3/)
!
!     calculating the transformation matrix
!
      if(iorien.gt.0) then
         call transformatrix(orab(1,iorien),pgauss,skl)
      endif
!
!     # of fibers
!
      nfiber=(-kode-102)/4
      do i=1,-kode-100
         constant(i)=elconloc(i)
      enddo
      if(dabs(constant(2)).lt.1.d-10) then
         constant(2)=1.d0/(20.d0*constant(1))
      endif
!
!     calculation of the Green deformation tensor for the
!     mechanical strain
!
      do i=1,3
         c(i,i)=emec(i)*2.d0+1.d0
      enddo
      c(1,2)=2.d0*emec(4)
      c(1,3)=2.d0*emec(5)
      c(2,3)=2.d0*emec(6)
!
!     creation of the delta Dirac matrix d
!
      do i=1,3
         d(i,i)=1.d0
      enddo
      d(1,2)=0.d0
      d(1,3)=0.d0
      d(2,3)=0.d0
!
!     calculation of the structural tensors
!
      do k=1,nfiber
         ioffset=4*k-1
         a(1)=constant(ioffset)
         a(2)=constant(ioffset+1)
         dd=a(1)*a(1)+a(2)*a(2)
         if(dd.gt.1.d0) then
            write(*,*) '*ERROR in umat_el_fiber: components of'
            write(*,*) '       direction vector ',k,' are too big'
            call exit(201)
         endif
         a(3)=dsqrt(1.d0-dd)
!
!        check for local coordinate systems
!
         if(iorien.gt.0) then
            do j=1,3
               aa(j)=a(j)
            enddo
            do j=1,3
               a(j)=skl(j,1)*aa(1)+skl(j,2)*aa(2)+skl(j,3)*aa(3)
            enddo
         endif
!
         do j=1,3
            do i=1,j
               dm(i,j,k)=a(i)*a(j)
            enddo
         enddo
      enddo
!
!     calculation of the invariants 
!
      v1=c(1,1)+c(2,2)+c(3,3)
      v3=c(1,1)*(c(2,2)*c(3,3)-c(2,3)*c(2,3))
     &     -c(1,2)*(c(1,2)*c(3,3)-c(1,3)*c(2,3))
     &     +c(1,3)*(c(1,2)*c(2,3)-c(1,3)*c(2,2))
      do j=1,nfiber
         v4(j)=dm(1,1,j)*c(1,1)+dm(2,2,j)*c(2,2)+dm(3,3,j)*c(3,3)+
     &       2.d0*(dm(1,2,j)*c(1,2)+dm(1,3,j)*c(1,3)+dm(2,3,j)*c(2,3))
      enddo
!
      v33=v3**(-1.d0/3.d0)
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
!
!     first derivative of the invariants with respect to c(k,l)
!
      do l=1,3
         do k=1,l
            didc(k,l,1)=d(k,l)
            didc(k,l,3)=v3*cinv(k,l)
            do j=1,nfiber
               djdc(k,l,j)=dm(k,l,j)
            enddo
         enddo
      enddo
!
!     second derivative of the invariants w.r.t. c(k,l) 
!     and c(m,n)
!
      if(icmd.ne.3) then
         nt=0
         do i=1,21
            k=kk(nt+1)
            l=kk(nt+2)
            m=kk(nt+3)
            n=kk(nt+4)
            nt=nt+4
            d2idc2(k,l,m,n,1)=0.d0
            d2idc2(k,l,m,n,3)=v3*(cinv(m,n)*cinv(k,l)-
     &           (cinv(k,m)*cinv(n,l)+cinv(k,n)*cinv(m,l))/2.d0)
            do j=1,nfiber
               d2jdc2(k,l,m,n,j)=0.d0
            enddo
         enddo
      endif
!
!     derivatives for the reduced invariants
!
      v1b=v1*v33
      v3bi=1.d0/dsqrt(v3)
      do j=1,nfiber
         v4br(j)=v4(j)*v33-1.d0
      enddo
!
!     first derivative of the reduced c-invariants w.r.t. c(k,l)
!
      do l=1,3
         do k=1,l
            dibdc(k,l,1)=-v33**4*v1*didc(k,l,3)/3.d0
     &           +v33*didc(k,l,1)
            do j=1,nfiber
               djbdc(k,l,j)=-v33**4*v4(j)*didc(k,l,3)/3.d0
     &           +v33*djdc(k,l,j)
            enddo
         enddo
      enddo
!
!     second derivative of the reduced c-invariants w.r.t. c(k,l)
!     and c(m,n)
!
      if(icmd.ne.3) then
         nt=0
         do i=1,21
            k=kk(nt+1)
            l=kk(nt+2)
            m=kk(nt+3)
            n=kk(nt+4)
            nt=nt+4
            d2ibdc2(k,l,m,n,1)=4.d0/9.d0*v33**7*v1*didc(k,l,3)
     &           *didc(m,n,3)-v33**4/3.d0*(didc(m,n,1)*didc(k,l,3)
     &           +didc(k,l,1)*didc(m,n,3))-v33**4/3.d0*v1*
     &           d2idc2(k,l,m,n,3)+v33*d2idc2(k,l,m,n,1)
            do j=1,nfiber
               d2jbdc2(k,l,m,n,j)=4.d0/9.d0*v33**7*v4(j)*didc(k,l,3)
     &              *didc(m,n,3)-v33**4/3.d0*(djdc(m,n,j)*didc(k,l,3)
     &              +djdc(k,l,j)*didc(m,n,3))-v33**4/3.d0*v4(j)*
     &              d2idc2(k,l,m,n,3)+v33*d2jdc2(k,l,m,n,j)
            enddo
         enddo
      endif
!
!     calculation of the stress
!     the anisotropy is only taken into account for v4br(j)>=0
!
      do l=1,3
         do k=1,l
            dudc(k,l)=constant(1)*dibdc(k,l,1)+
     &        (1.d0-v3bi)*didc(k,l,3)/constant(2)
            do j=1,nfiber
               if(v4br(j).lt.0.d0) cycle
               if(xk2*v4br(j)**2.gt.227.d0) then
                  write(*,*) '*ERROR in umat_elastic_fiber'
                  write(*,*) '       fiber extension is too large'
                  write(*,*) '       for exponential function'
                  call exit(201)
               endif
               ioffset=4*j
               xk1=constant(ioffset+1)
               xk2=constant(ioffset+2)
               dudc(k,l)=dudc(k,l)+xk1*v4br(j)*
     &           dexp(xk2*v4br(j)**2)*djbdc(k,l,j)
            enddo
         enddo
      enddo
!
!     calculation of the stiffness matrix
!     the anisotropy is only taken into account for v4br(j)>=0
!
      if(icmd.ne.3) then
         nt=0
         do i=1,21
            k=kk(nt+1)
            l=kk(nt+2)
            m=kk(nt+3)
            n=kk(nt+4)
            nt=nt+4
            term=constant(1)*d2ibdc2(k,l,m,n,1)+
     &           v3bi**3*didc(k,l,3)*didc(m,n,3)/(2.d0*constant(2))
     &           +(1.d0-v3bi)*d2idc2(k,l,m,n,3)/constant(2)
            do j=1,nfiber
               if(v4br(j).lt.0.d0) cycle
               ioffset=4*j
               xk1=constant(ioffset+1)
               xk2=constant(ioffset+2)
               term=term+xk1*dexp(xk2*v4br(j)**2)*
     &            (djbdc(k,l,j)*djbdc(m,n,j)*(1.d0+2.d0*xk2*v4br(j)**2)+
     &             v4br(j)*d2jbdc2(k,l,m,n,j))
            enddo
            d2udc2(k,l,m,n)=term
         enddo
      endif
!
!     storing the stiffness matrix and/or the stress
!
      if(icmd.ne.3) then
!
!        storing the stiffness matrix
!
         nt=0
         do i=1,21
            k=kk(nt+1)
            l=kk(nt+2)
            m=kk(nt+3)
            n=kk(nt+4)
            nt=nt+4
            stiff(i)=4.d0*d2udc2(k,l,m,n)
         enddo
      endif
!
!     store the stress at mechanical strain
!
      stre(1)=2.d0*dudc(1,1)
      stre(2)=2.d0*dudc(2,2)
      stre(3)=2.d0*dudc(3,3)
      stre(4)=2.d0*dudc(1,2)
      stre(5)=2.d0*dudc(1,3)
      stre(6)=2.d0*dudc(2,3)
!
      return
      end





