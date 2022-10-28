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
      subroutine e_c3d_u1(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &  ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,
     &  t0,t1,ithermal,vold,iperturb,nelemload,
     &  sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &  nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &  matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &  ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &  ne0,ipkon,thicke,
     &  integerglob,doubleglob,tieset,istartset,iendset,ialset,ntie,
     &  nasym,ielprop,prop)
!
!     computation of the element matrix and rhs for the user element
!     of type u1
!
!     This is a beam type element. Reference:
!     Yunhua Luo, An Efficient 3D Timoshenko Beam Element with
!     Consistent Shape Functions, Adv. Theor. Appl. Mech., Vol. 1,
!     2008, no. 3, 95-106
!
!     special case for which the beam axis goes through the
!     center of gravity of the cross section and the 1-direction
!     corresponds with a principal axis
!
!
!     INPUT:
!
!     co(1..3,i)         coordinates of node i
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     lakonl             label of current element nelem (character*8)
!     p1(1..3)           coordinates of a point on the rotation axis
!                        (if applicable)
!     p2(1..3)           unit vector on rotation axis (if applicable)
!     omx                rotational speed square
!     bodyfx(1..3)       acceleration vector
!     nbody              number of body loads
!     nelem              element number
!     nmethod            procedure:
!                        1: static analysis
!                        2: frequency analysis  
!                        3: buckling analysis 
!                        4: (linear or nonlinear) dynamic analysis 
!                        5: steady state dynamics analysis 
!                        6: Coriolis frequency calculation 
!                        7: flutter frequency calculation 
!                        8:  magnetostatics 
!                        9:  magnetodynamics 
!                        10: electromagnetic eigenvalue problems 
!                        11: superelement creation or Green function 
!                            calculation 
!                        12: sensitivity analysis  
!     elcon              elastic constants (cf. List of variables and
!                        their meaning in the User's Manual)
!     nelcon             integers describing the elastic constant fields
!                        (cf. User's Manual)
!     rhcon              density constants (cf. User's Manual)
!     nrhcon             integers describing the density constant fields
!                        (cf. User's Manual)
!     alcon              thermal expansion constants (cf. User's Manual)
!     nalcon             integers describing the thermal expansion 
!                        constants (cf. User's Manual)
!     alzero             thermal expansion reference values (cf. User's Manual)
!     ielmat(i)          material for element i
!     ielorien(i)        orientation for element i
!     norien             number of orientations
!     orab(7,*)          description of all local coordinate systems.
!                        (cf. List of variables and their meaning in the
!                        User's manual)
!     ntmat_             maximum number of material temperature points
!     t0(i)              temperature in node i at start of calculation
!     t1(i)              temperature in node i at the end of the current
!                        increment
!     ithermal           cf. ithermal(1) in the List of variables and
!                        their meaning in the User's Manual
!     vold(0..mi(2),i)   value of the variables in node i at the end
!                        of the previous iteration
!     iperturb(*)        describes the kind of nonlinearity of the
!                        calculation, cf. User's Manual
!     nelemload          field describing facial distributed load (cf.
!                        User's Manual)
!     sideload           field describing facial distributed load (cf.
!                        User's Manual)
!     xload              facial distributed load values (cf.
!                        User's Manual)
!     nload              number of facial distributed loads
!     idist              0: no distributed forces in the model
!                        1: distributed forces are present in the model
!                        (body forces or thermal loads or residual 
!                        stresses or distributed facial loads)
!     sti(i,j,k)         stress component i in integration point j
!                        of element k at the end of the previous step
!     stx(i,j,k)         stress component i in integration point j
!                        of element k at the end of a static step
!                        describing a reference buckling load
!     iexpl              0: structure: implicit dynamics,
!                           fluid: incompressible
!     iexpl              1: structure: implicit dynamics,
!                           fluid: compressible
!     iexpl              2: structure: explicit dynamics,
!                           fluid: incompressible
!     iexpl              3: structure: explicit dynamics,
!                           fluid: compressible
!     plicon,nplicon     fields describing isotropic hardening of
!                        a plastic material or spring constants of
!                        a nonlinear spring (cf. User's Manual)
!     plkcon,nplkcon     fields describing kinematic hardening of
!                        a plastic material or gap conductance
!                        constants (cf. User's Manual)
!     xstiff(i,j,k)      stiffness (i=1...21) and conductivity
!                        (i=22..27) constants in integration point j
!                        of element k
!     npmat_             maximum number of plastic constants
!     dtime              length of present time increment
!     matname(i)         name of material i
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     mi(3)              max number of layers in the structure
!     ncmat_             max number of elastic constants
!     mass(1)            if 1: mass matrix needed, else not
!     mass(2)            if 1: heat capacity matrix needed, else not
!     stiffness          if 1: stiffness matrix needed, else not
!     buckling           if 1: linear buckling calculation, else not
!     rhsi               if 1: right hand side needed, else not
!     intscheme          1: integration point scheme for the calculation
!                           of the initial acceleration in a dynamic
!                           step
!                        0: any other case
!     ttime              total time at the start of the present increment
!     time               step time at the end of the present increment
!     istep              current step number
!     iinc               current increment number within the actual step
!     coriolis           if 1: Coriolis forces are taken into account,
!                              else not
!     xloadold           load at the end of the previous increment
!                        (cf. User's Manual)
!     reltime            relative step time (between 0 and 1)
!     ipompc(i)          pointer into field nodempc and coefmpc for
!                        MPC i
!     nodempc            field containing the nodes and dof's of all 
!                        MPC's (cf. User's Manual)
!     coefmpc            field containing the coefficients of all
!                        MPC's  (cf. User's Manual)
!     nmpc               number of MPC's in the model
!     ikmpc(i)           field containing an integer code number 
!                        describing the dependent dof of some MPC:
!                        ikmpc(i)=8*(node-1)+dir, where the dependent
!                        dof is applied in direction dir of node "node"
!                        ikmpc is sorted in ascending order
!     ilmpc(i)           MPC number corresponding to ikmpc(i) (cf.
!                        User's manual)
!     veold(j,i)         time rate of variable j in node i at the end
!                        of the previous iteration
!     ne0                element number before the first contact element;
!                        contact elements are appended after ne0
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     thicke(j,i)        thickness of layer j in node i
!     integerglob        integer field needed for determining the
!                        interface loading of a submodel
!     doubleglob         real field needed for determining the
!                        interface loading of a submodel
!     tieset(1..3,*)     tie character information for tie i (cf.
!                        User's Manual)
!     istartset          integer describing set information (cf.
!                        User's Manual)
!     iendset            integer describing set information (cf.
!                        User's Manual)
!     ialset             integer describing set information (cf.
!                        User's Manual)
!     ntie               number of ties in the model
!     nasym              1: asymmetric contributions, else purely
!                        symmetric
!     ielprop(i)         points to the location in field prop preceding
!                        the properties of element i
!     prop(*)            contains the properties and some beam 
!                        elements (cf. User's Manual)
!
!
!     OUTPUT:
!
!     s                  element stiffness matrix
!     sm                 element mass matrix
!     ff                 element load vector
!     xload              facial distributed load values (cf.
!                        User's Manual)
!     nmethod            may be changed by the user, e.g. a
!                        change to 0 for negative Jacobians leads 
!                        to a program exit in the calling program
!                        for the other value: cf. above
!
      implicit none
!
      integer mass,stiffness,buckling,rhsi,coriolis
!
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
      character*81 tieset(3,*)
!
      integer konl(26),nelemload(2,*),nbody,nelem,mi(*),kon(*),
     &  ielprop(*),index,mattyp,ithermal(*),iperturb(*),nload,idist,
     &  i,j,k,i1,nmethod,kk,nelcon(2,*),nrhcon(*),nalcon(2,*),
     &  ielmat(mi(3),*),ielorien(mi(3),*),ipkon(*),indexe,
     &  ntmat_,nope,norien,ihyper,iexpl,kode,imat,iorien,istiff,
     &  ncmat_,intscheme,istep,iinc,ipompc(*),nodempc(3,*),
     &  nmpc,ikmpc(*),ilmpc(*),ne0,ndof,istartset(*),iendset(*),
     &  ialset(*),ntie,integerglob(*),nasym,nplicon(0:ntmat_,*),
     &  nplkcon(0:ntmat_,*),npmat_
!
      real*8 co(3,*),xl(3,26),veold(0:mi(2),*),rho,s(60,60),bodyfx(3),
     &  ff(60),elconloc(ncmat_),coords(3),p1(3),
     &  p2(3),eth(6),rhcon(0:1,ntmat_,*),reltime,prop(*),tm(3,3),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &  xloadold(2,*),vold(0:mi(2),*),xload(2,*),omx,e,un,um,tt,
     &  sm(60,60),sti(6,mi(1),*),stx(6,mi(1),*),t0l,t1l,coefmpc(*),
     &  elas(21),thicke(mi(3),*),doubleglob(*),dl,e2(3),e3(3),
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,tmg(12,12),
     &  a,xi11,xi12,xi22,xk,e1(3),offset1,offset2,y1,y2,y3,z1,z2,z3,
     &  sg(12,12),elcon(0:ncmat_,ntmat_,*)
!
      indexe=ipkon(nelem)
!
!     material and orientation
!
      imat=ielmat(1,nelem)
      amat=matname(imat)
      if(norien.gt.0) then
         iorien=max(0,ielorien(1,nelem))
      else
         iorien=0
      endif
!     
      if(nelcon(1,imat).lt.0) then
         ihyper=1
      else
         ihyper=0
      endif
!
!     properties of the cross section
!
      index=ielprop(nelem)
      a=prop(index+1)
      xi11=prop(index+2)
      xi12=prop(index+3)
      xi22=prop(index+4)
      xk=prop(index+5)
      e2(1)=prop(index+6)
      e2(2)=prop(index+7)
      e2(3)=prop(index+8)
      offset1=prop(index+9)
      offset2=prop(index+10)
!
      nope=2
      ndof=6
!
      if(dabs(xi12).gt.0.d0) then
         write(*,*) '*ERROR in e_c3d_u1: no nonzero cross moment'
         write(*,*) '       of inertia for this type of element'
         call exit(201)
      endif
!
      if(dabs(offset1).gt.0.d0) then
         write(*,*) '*ERROR in e_c3d_u1: no offset in 1-direction'
         write(*,*) '       for this type of element'
         call exit(201)
      endif
!
      if(dabs(offset2).gt.0.d0) then
         write(*,*) '*ERROR in e_c3d_u1: no offset in 2-direction'
         write(*,*) '       for this type of element'
         call exit(201)
      endif
!
!     computation of the coordinates of the local nodes
!
      do i=1,nope
         konl(i)=kon(indexe+i)
         do j=1,3
            xl(j,i)=co(j,konl(i))
         enddo
      enddo
!
!     local axes e1-e2-e3  (e2 is given by the user (in the input deck
!     this is called e1), e1 is parallel to the beam axis)
!
      do j=1,3
         e1(j)=xl(j,2)-xl(j,1)
      enddo
!
      dl=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
!
      do j=1,3
         e1(j)=e1(j)/dl
      enddo
!
!     e3 = e1 x e2
!
      e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
      e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
      e3(3)=e1(1)*e2(2)-e1(2)*e2(1)
c      write(*,*) 'u1 trans1 ',e1(1),e1(2),e1(3)
c      write(*,*) 'u1 trans2 ',e2(1),e2(2),e2(3)
c      write(*,*) 'u1 trans3 ',e3(1),e3(2),e3(3)
!
!     transformation matrix from the global into the local system
!
      do j=1,3
         tm(1,j)=e1(j)
         tm(2,j)=e2(j)
         tm(3,j)=e3(j)
      enddo
!
!     initialisation for distributed forces
!
      if(rhsi.eq.1) then
        if(idist.ne.0) then
          do i=1,ndof*nope
            ff(i)=0.d0
          enddo
        endif
      endif
!
!     displacements for 2nd order static and modal theory
!
      if(((iperturb(1).eq.1).or.(iperturb(2).eq.1)).and.
     &          (stiffness.eq.1).and.(buckling.eq.0)) then
         write(*,*) '*ERROR in e_c3d_u1: no second order'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
!
!     initialisation of sm
!
      if((mass.eq.1).or.(buckling.eq.1).or.(coriolis.eq.1)) then
         write(*,*) '*ERROR in e_c3d_u1: no dynamic or buckling'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
!
!     initialisation of s
!
      do i=1,ndof*nope
        do j=1,ndof*nope
          s(i,j)=0.d0
        enddo
      enddo
!
!     computation of the matrix
!
!     calculating the temperature in the integration
!     point
!     
      t0l=0.d0
      t1l=0.d0
      if(ithermal(1).eq.1) then
         do i1=1,nope
            t0l=t0l+t0(konl(i1))/2.d0
            t1l=t1l+t1(konl(i1))/2.d0
         enddo
      elseif(ithermal(1).ge.2) then
         write(*,*) '*ERROR in e_c3d_u1: no thermal'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
      tt=t1l-t0l
!     
!     calculating the coordinates of the integration point
!     for material orientation purposes (for cylindrical
!     coordinate systems)
!     
      if(iorien.gt.0) then
         write(*,*) '*ERROR in e_c3d_u1: no orientation'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
!     
!     for deformation plasticity: calculating the Jacobian
!     and the inverse of the deformation gradient
!     needed to convert the stiffness matrix in the spatial
!     frame of reference to the material frame
!     
      kode=nelcon(1,imat)
!
      if(kode.eq.2) then
         mattyp=1
      else
         mattyp=0
      endif
!     
!     material data and local stiffness matrix
!     
      istiff=0
      call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     imat,amat,iorien,coords,orab,ntmat_,elas,rho,
     &     nelem,ithermal,alzero,mattyp,t0l,t1l,
     &     ihyper,istiff,elconloc,eth,kode,plicon,
     &     nplicon,plkcon,nplkcon,npmat_,
     &     plconloc,mi(1),dtime,kk,
     &     xstiff,ncmat_)
!     
      if(mattyp.eq.1) then
         e=elconloc(1)
         un=elconloc(2)
         um=e/(1.d0+un)
         um=um/2.d0
      elseif(mattyp.eq.2) then
         write(*,*) '*ERROR in e_c3d_u1: no orthotropic material'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      else
         write(*,*) '*ERROR in e_c3d_u1: no anisotropic material'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
!     
!     initialization for the body forces
!     
      if(rhsi.eq.1) then
         if(nbody.ne.0) then
            write(*,*) '*ERROR in e_c3d_u1: no body forces'
            write(*,*) '       for this type of element'
            call exit(201)
         endif
      endif
!     
      if(buckling.eq.1) then
         write(*,*) '*ERROR in e_c3d_u1: no buckling '
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
!     
!     determination of the stiffness, and/or mass and/or
!     buckling matrix
!     
      if((stiffness.eq.1).or.(mass.eq.1).or.(buckling.eq.1).or.
     &     (coriolis.eq.1)) then
!     
         if(((iperturb(1).ne.1).and.(iperturb(2).ne.1)).or.
     &        (buckling.eq.1)) then
!     
!     stiffness matrix
!   
            y1=xk*um*a*e*xi11*(12.d0*e*xi11+xk*um*a*dl*dl)
            y2=(12.d0*e*xi11-xk*um*a*dl*dl)**2
            y3=4.d0*e*xi11*((xk*um*a)**2*dl**4+
     &         3.d0*xk*um*a*dl*dl*e*xi11+
     &         36.d0*(e*xi11)**2)
            z1=xk*um*a*e*xi22*(12.d0*e*xi22+xk*um*a*dl*dl)
            z2=(12.d0*e*xi22-xk*um*a*dl*dl)**2
            z3=4.d0*e*xi22*((xk*um*a)**2*dl**4+
     &         3.d0*xk*um*a*dl*dl*e*xi22+
     &         36.d0*(e*xi22)**2)
!  
!     stiffness matrix S' in local coordinates
!
!     this if from Appendix A of the article by Yunhua Luo;
!     the following corrections were made to the equations in the
!     article: 1) in line 4: k(3,3)=k(3,9) should be replaced by
!                 k(3,3)=-k(3,9)
!              2) in line 6: k(4,4)=k(7,7) should be replaced by
!                 k(4,4)=k(10,10)
!
            s(1,1)=e*a/dl
            s(1,7)=-s(1,1)
            s(2,2)=12.d0*y1/(dl*y2)
            s(2,6)=6.d0*y1/y2
            s(2,8)=-s(2,2)
            s(2,12)=s(2,6)
            s(3,3)=12.d0*z1/(dl*z2)
            s(3,5)=-6.d0*z1/z2
            s(3,9)=-s(3,3)
            s(3,11)=s(3,5)
            s(4,4)=um*(xi11+xi22)/dl
            s(4,10)=-s(4,4)
            s(5,5)=z3/(dl*z2)
            s(5,9)=6.d0*z1/z2
            s(5,11)=-2.d0*e*xi22*(72.d0*(e*xi22)**2-
     &              (xk*um*a)**2*dl**4-
     &              30.d0*xk*um*a*dl*dl*e*xi22)/(dl*z2)
            s(6,6)=y3/(dl*y2)
            s(6,8)=-6.d0*y1/y2
            s(6,12)=-2.d0*e*xi11*(-(xk*um*a)**2*dl**4-
     &              30.d0*xk*um*a*dl*dl*e*xi11+
     &              72.d0*(e*xi11)**2)/(dl*y2)
            s(1,7)=-s(1,1)
            s(7,7)=s(1,1)
            s(8,8)=12.d0*y1/(dl*y2)
            s(8,12)=-6.d0*y1/y2
            s(9,9)=12.d0*z1/(dl*z2)
            s(9,11)=6.d0*z1/z2
            s(10,10)=s(4,4)
            s(11,11)=z3/(dl*z2)
            s(12,12)=y3/(dl*y2)
!
!           completing the symmetric part
! 
            do i=1,12
               do j=1,i
                  s(i,j)=s(j,i)
               enddo
            enddo
!
!           12 x 12 transformation matrix
!
            do i=1,12
               do j=1,12
                  tmg(i,j)=0.d0
               enddo
            enddo
            do i=1,3
               do j=1,3
                  tmg(i,j)=tm(i,j)
                  tmg(i+3,j+3)=tm(i,j)
                  tmg(i+6,j+6)=tm(i,j)
                  tmg(i+9,j+9)=tm(i,j)
               enddo
            enddo
!
!           stiffness matrix in global coordinates: S = T^T*S'*T
!
c            do i=1,12
c              do j=1,12
c                write(*,*) 'u1 lstiffness',i,j,s(i,j)
c              enddo
c            enddo
!                
            do i=1,12
               do j=1,12
                  sg(i,j)=0.d0
                  do k=1,12
                     sg(i,j)=sg(i,j)+s(i,k)*tmg(k,j)
                  enddo
               enddo
            enddo
!
!           only upper triangular matrix
!
            do i=1,12
               do j=i,12
                  s(i,j)=0.d0
                  do k=1,12
                     s(i,j)=s(i,j)+tmg(k,i)*sg(k,j)
                   enddo
c                   write(*,*) 'u1 stiffness',i,j,s(i,j)
               enddo
            enddo
!
         endif
!     
      endif
!     
      return
      end

