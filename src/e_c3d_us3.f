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
      subroutine e_c3d_us3(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
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
     &  nplkcon(0:ntmat_,*),npmat_,jjj,id,ig
!
      real*8 co(3,*),xl(3,26),veold(0:mi(2),*),rho,s(60,60),bodyfx(3),
     &  ff(60),elconloc(ncmat_),coords(3),p1(3),
     &  p2(3),eth(6),rhcon(0:1,ntmat_,*),reltime,prop(*),tm(3,3),
     &  alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &  xloadold(2,*),vold(0:mi(2),*),xload(2,*),omx,e,un,um,tt,
     &  sm(60,60),sti(6,mi(1),*),stx(6,mi(1),*),t0l,t1l,coefmpc(*),
     &  elas(21),thicke(mi(3),*),doubleglob(*),dl,e2(3),e3(3),h,
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,xg(3,3),
     &  x(3,3),Ds(2,2),Qs(2,2),Qin(3,3),Dm(3,3),Db(3,3),Dc(3,3),
     &  Kp(18,18),Km(18,18),tmg(18,18),Kshell(18,18),dett,dettt,
     &  Mshell(18,18),Dmi(3,3),Dbi(3,3),Dsi(2,2),di,gpthick(3,2),
     &  pres,fpres(18),Ae,elcon(0:ncmat_,ntmat_,*)
      !      
      gpthick(1,1) = +1.d0
      gpthick(2,1) =  0.d0
      gpthick(3,1) = -1.d0
      gpthick(1,2) = 2.d0/6.d0
      gpthick(2,2) = 8.d0/6.d0
      gpthick(3,2) = 2.d0/6.d0 
      !
      di = 1.d0/3.d0
      !
      indexe = ipkon(nelem)
      index=ielprop(nelem)
      !
      nope = 3
      !
      !   properties of the cross section
      !
      h     = prop(index+1)	  ! user section porp.
      dett  = h/2.d0
      dettt = h**3/8.d0
      !
      !   coordinates
      !
      do i = 1,nope
         konl(i) = kon(indexe+i)
         do j = 1,3
            xl(i,j) = co(j,konl(i))
            xg(i,j) = co(j,konl(i)) ! xg(node,xi)
         enddo
      enddo  
      !
      call us3_csys(xg,tm,tmg) ! e0-frame
      !call us3_csys_cr(xg,tm,tmg) !e0-frame
      !
      ! nodal coordinates in e0-frame
      !
      x(:,:) = 0.d0
      x(1,:) = matmul(tm,xg(1,:))
      x(2,:) = matmul(tm,xg(2,:))
      x(3,:) = matmul(tm,xg(3,:))
      !
      !     displacements for 2nd order static and modal theory
      !
      if(((iperturb(1).eq.1).or.(iperturb(2).eq.1)).and.
     &          (stiffness.eq.1).and.(buckling.eq.0)) then
         write(*,*) '*ERROR in e_c3d_US3: no second order'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif
      !
      if((buckling.eq.1).or.(coriolis.eq.1)) then
         write(*,*) '*ERROR in e_c3d_US3: no buckling'
         write(*,*) '       calculation for this type of element'
         call exit(201)
      endif      
      !
      !     material and orientation
      !
      imat = ielmat(1,nelem)
      amat = matname(imat)
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
      kode = nelcon(1,imat)
      !
      if(kode.eq.2) then
         mattyp=1
      else
         mattyp=0
      endif
      !
      Dm(:,:) = 0.d0
      Db(:,:) = 0.d0
      Ds(:,:) = 0.d0
      do jjj = 1,3  ! workaround for thermals
        istiff = 1
!
!     in materialdata_me the following data are calculated:
!     - the tangent stiffness, taken from resultsmech.f
!     - the density. For solid calculations this is the density in the undeformed
!       state. Any temperature effect due to loading is taken into account by
!       the Jacobian, and not by temperature interpolation. Therefore, the
!       temperature can be set to zero here.        
!     
        t0l=0.d0
        t1l=0.d0
        call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     imat,amat,iorien,coords,orab,ntmat_,elas,rho,
     &     nelem,ithermal,alzero,mattyp,t0l,t1l,
     &     ihyper,istiff,elconloc,eth,kode,plicon,
     &     nplicon,plkcon,nplkcon,npmat_,
     &     plconloc,mi(1),dtime,jjj,
     &     xstiff,ncmat_)
        !
        e   = elas(1)
        un  = elas(2)
        !
        call us3_linel_Qi(e,un,Qin,Qs) 
        Dm = Dm+Qin*gpthick(jjj,2)*dett
        Db = Db+Qin*((gpthick(jjj,1))**2)*gpthick(jjj,2)*dettt 
        Ds = Ds+Qs*gpthick(jjj,2)*dett 
      enddo   
      !
      !   stiffness matirx (6 dofs 3 nodes -> 18x18)  
      !
      call us3_Kp(x,Db,Ds,Kp) ! plate part (CS-DSG)   - in e0-frame
      call us3_Km(x,Km,Qin,h) ! membrane part (ANDES) - in e0-frame
      Kshell = Km + Kp
      !
      call us3_M(x,h,rho,Mshell)
      !
      Kshell = matmul(matmul(transpose(tmg),Kshell),tmg)
      Mshell = matmul(matmul(transpose(tmg),Mshell),tmg)
      !
      do i=1,18
        do j=1,18
            sm(i,j) = Mshell(i,j) 
            s(i,j)  = Kshell(i,j)
        enddo
      enddo 
!
! pressure loads defined y/n
!
      if(rhsi.eq.1) then
        if(idist.ne.0) then
          ff(:)=0.d0
        endif
      endif
!      
      if(nload.gt.0) then ! idlist=1 distributed loads defined
        !
        call nident2(nelemload,nelem,nload,id)
        ! 
        ! *SURFACE arg (sideload 1 NEG, 2 POS):
        !
c        ff(:) = 0.d0
        if(sideload(id)(1:2).eq.'P2') then
          pres = -xload(1,id)
        elseif(sideload(id)(1:2).eq.'P1') then
          pres = +xload(1,id) 
        endif
        !
        call us3_ae(x,Ae) ! element area
        ff(3)  = di*Ae*pres
        ff(9)  = di*Ae*pres
        ff(15) = di*Ae*pres
        !
      endif
      !
      return
      end

