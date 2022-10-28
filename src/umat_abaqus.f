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
      subroutine umat_abaqus(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,kinc,pnewdt,nmethod,iperturb)
!
!     calculates stiffness and stresses for a nonlinear material
!     defined by an ABAQUS umat routine
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
!
!     This routine allows for the use of an ABAQUS umat user subroutine
!     in CalculiX. 
!
!     Note that the following fields are not supported
!     so far: sse,spd,scd,rpl,ddsddt,drplde,drpldt,predef,
!     dpred,drot,pnewdt,celent,layer,kspt
!
!     Furthermore, the following fields have a different meaning in
!     ABAQUS and CalculiX:
!
!     stran: in CalculiX: Lagrangian strain tensor
!              in ABAQUS: logarithmic strain tensor
!     dstran: in CalculiX: Lagrangian strain increment tensor
!              in ABAQUS: logarithmic strain increment tensor
!     temp:  in CalculiX: temperature at the end of the increment
!              in ABAQUS: temperature at the start of the increment
!     dtemp: in CalculiX: zero
!              in ABAQUS: temperature increment
!
!     Because of this, this routine should only be used for small
!     deformations and small rotations (in that case all strain 
!     measures basically reduce to the infinitesimal strain). 
!
      implicit none
!
      character*80 amat
!
      integer ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),i,
     &  iorien,nmethod,iperturb(*),istep,
     &  ndi,nshr,ntens,nprops,layer,kspt,jstep(4),kinc,kal(2,6),
     &  kel(4,21),j1,j2,j3,j4,j5,j6,j7,j8,jj
!
      real*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,skl(3,3),xa(3,3),ya(3,3,3,3),xstate(nstate_,mi(1),*),
     &  xstateini(nstate_,mi(1),*)
!
      real*8 ddsdde(6,6),sse,spd,scd,rpl,ddsddt(6),drplde(6),
     &  drpldt,stran(6),dstran(6),abqtime(2),predef,temp,dtemp,
     &  dpred,drot(3,3),celent,pnewdt
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      kel=reshape((/1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &          1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,
     &          3,3,1,3,1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,
     &          1,2,2,3,1,3,2,3,2,3,2,3/),(/4,21/))
!
      drot=reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),
     &            (/3,3/))
!
!     filling field jstep
!
      jstep(1)=istep
      jstep(2)=nmethod
      jstep(3)=iperturb(2)
      if(iperturb(1).eq.1) then
         jstep(4)=1
      else
         jstep(4)=0
      endif
!
!     calculating the mechanical strain
!
      do i=1,6
         stran(i)=emec0(i)
         dstran(i)=emec(i)-emec0(i)
      enddo
!      
      ntens=6
!
      do i=1,nstate_
         xstate(i,iint,iel)=xstateini(i,iint,iel)
      enddo
!
      abqtime(1)=time-dtime
      abqtime(2)=ttime+time-dtime
!
      temp=t1l
      dtemp=0.d0
!
      ndi=3
      nshr=3
      ntens=ndi+nshr
!
      nprops=-kode-100
c      nprops=21
!
!     taking local material orientations into account
!
      if(iorien.ne.0) then
         call transformatrix(orab(1,iorien),pgauss,skl)
!
!        rotating the stress into the local system
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
     &                 xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
               enddo
            enddo
         enddo
!
!        rotating the strain into the local system
!
         xa(1,1)=stran(1)
         xa(1,2)=stran(4)
         xa(1,3)=stran(5)
         xa(2,1)=stran(4)
         xa(2,2)=stran(2)
         xa(2,3)=stran(6)
         xa(3,1)=stran(5)
         xa(3,2)=stran(6)
         xa(3,3)=stran(3)
!
         do jj=1,6
            stran(jj)=0.d0
            j1=kal(1,jj)
            j2=kal(2,jj)
            do j3=1,3
               do j4=1,3
                  stran(jj)=stran(jj)+
     &                 xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
               enddo
            enddo
         enddo
!
!        rotating the strain increment into the local system
!
         xa(1,1)=dstran(1)
         xa(1,2)=dstran(4)
         xa(1,3)=dstran(5)
         xa(2,1)=dstran(4)
         xa(2,2)=dstran(2)
         xa(2,3)=dstran(6)
         xa(3,1)=dstran(5)
         xa(3,2)=dstran(6)
         xa(3,3)=dstran(3)
!
         do jj=1,6
            dstran(jj)=0.d0
            j1=kal(1,jj)
            j2=kal(2,jj)
            do j3=1,3
               do j4=1,3
                  dstran(jj)=dstran(jj)+
     &                 xa(j3,j4)*skl(j3,j1)*skl(j4,j2)
               enddo
            enddo
         enddo
      endif
!
!     changing physical strain into engineering strain
!
      do i=4,6
         stran(i)=2.d0*stran(i)
         dstran(i)=2.d0*dstran(i)
      enddo
!
      if(amat(1:1).eq.'@') then
!
         call call_external_umat(stre,xstate(1,iint,iel),ddsdde,
     &        sse,spd,scd,rpl,ddsddt,drplde,drpldt,stran,dstran,
     &        abqtime,dtime,temp,dtemp ,predef,dpred,amat,ndi,nshr,
     &        ntens,nstate_,elconloc,nprops,pgauss,drot,pnewdt,
     &        celent,xokl,xkl,iel,iint,layer,kspt,jstep,kinc)
!
      else
!
         call umat(stre,xstate(1,iint,iel),ddsdde,sse,spd,scd,rpl,
     &        ddsddt,drplde,drpldt,stran,dstran,abqtime,dtime,temp,
     &        dtemp,predef,dpred,amat,ndi,nshr,ntens,nstate_,elconloc,
     &        nprops,pgauss,drot,pnewdt,celent,xokl,xkl,iel,iint,layer,
     &        kspt,jstep,kinc)
!
      endif
!
!     taking local material orientations into account
!
      if(iorien.ne.0) then
!
!        rotating the stress into the global system
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
     &                 xa(j3,j4)*skl(j1,j3)*skl(j2,j4)
               enddo
            enddo
         enddo
      endif
!
!     calculate the stiffness matrix (the matrix is symmetrized)
!
      if(icmd.ne.3) then
         stiff(1)=ddsdde(1,1)
         stiff(2)=(ddsdde(1,2)+ddsdde(2,1))/2.d0
         stiff(3)=ddsdde(2,2)
         stiff(4)=(ddsdde(1,3)+ddsdde(3,1))/2.d0
         stiff(5)=(ddsdde(2,3)+ddsdde(3,2))/2.d0
         stiff(6)=ddsdde(3,3)
         stiff(7)=(ddsdde(1,4)+ddsdde(4,1))/2.d0
         stiff(8)=(ddsdde(2,4)+ddsdde(4,2))/2.d0
         stiff(9)=(ddsdde(3,4)+ddsdde(4,3))/2.d0
         stiff(10)=ddsdde(4,4)
         stiff(11)=(ddsdde(1,5)+ddsdde(5,1))/2.d0
         stiff(12)=(ddsdde(2,5)+ddsdde(5,2))/2.d0
         stiff(13)=(ddsdde(3,5)+ddsdde(5,3))/2.d0
         stiff(14)=(ddsdde(4,5)+ddsdde(5,4))/2.d0
         stiff(15)=ddsdde(5,5)
         stiff(16)=(ddsdde(1,6)+ddsdde(6,1))/2.d0
         stiff(17)=(ddsdde(2,6)+ddsdde(6,2))/2.d0
         stiff(18)=(ddsdde(3,6)+ddsdde(6,3))/2.d0
         stiff(19)=(ddsdde(4,6)+ddsdde(6,4))/2.d0
         stiff(20)=(ddsdde(5,6)+ddsdde(6,5))/2.d0
         stiff(21)=ddsdde(6,6)
!
         if(iorien.ne.0) then
!
!        rotating the stiffness coefficients into the global system
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
     &                       skl(j1,j5)*skl(j2,j6)*skl(j3,j7)*skl(j4,j8)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
!
      return
      end
