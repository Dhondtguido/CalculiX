!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine resultsmech_us3(co,kon,ipkon,lakon,ne,v,
     &     stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &     ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &     veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &     xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &     ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &     reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &     ikin,nal,ne0,thicke,emeini,i,ielprop,prop,t0g,t1g)
!     
!     calculates nal,qa,fn,xstiff,ener,eme,eei,stx for user element 1
!     
!     This is a beam type element. Reference:
!     Yunhua Luo, An Efficient 3D Timoshenko Beam Element with
!     Consistent Shape Functions, Adv. Theor. Appl. Mech., Vol. 1,
!     2008, no. 3, 95-106
!     
!     
!     special case for which the beam axis goes through the
!     center of gravity of the cross section and the 1-direction
!     corresponds with a principal axis
!     
!     note that the strain components are Lagrange strain components,
!     the stress components are Piola-Kirchhoff components of the
!     second kind. For linear geometric calculations the strains
!     reduce to the infinitesimal strains
!     
!     
!     INPUT:
!     
!     co(1..3,i)         coordinates of node i
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     lakon(i)           label of element i (character*8)
!     ne                 highest element number in the mesh
!     v(0..mi(2),i)      value of the variables in node i at the end
!                        of the present iteration
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
!     ithermal(1..2)     cf. List of variables and
!                        their meaning in the User's Manual
!     prestr(i,j,k)      residual stress component i in integration point j
!                        of element k 
!     iprestr            if 0: no residual stresses
!                        else: residual stresses
!     iperturb(*)        describes the kind of nonlinearity of the
!                        calculation, cf. User's Manual
!     iout               if -2: v is assumed to be known and is used to
!                               calculate strains, stresses..., no result output
!                               corresponds to iout=-1 with in addition the
!                               calculation of the internal energy density
!                        if -1: v is assumed to be known and is used to
!                               calculate strains, stresses..., no result 
!                               output; is used to take changes in SPC's and 
!                               MPC's at the start of a new increment or 
!                               iteration into account
!                        if 0: v is calculated from the system solution
!                              and strains, stresses.. are calculated, 
!                              no result output
!                        if 1: v is calculated from the system solution and 
!                              strains, stresses.. are calculated, requested 
!                              results output
!                        if 2: v is assumed to be known and is used to 
!                              calculate strains, stresses..., requested 
!                              results output
!     vold(0..mi(2),i)   value of the variables in node i at the end
!                        of the previous iteration
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
!     veold(j,i)         time rate of variable j in node i at the end
!                        of the previous iteration
!     dtime              length of present time increment
!     time               step time at the end of the present increment
!     ttime              total time at the start of the present increment
!     plicon,nplicon     fields describing isotropic hardening of
!                        a plastic material or spring constants of
!                        a nonlinear spring (cf. User's Manual)
!     plkcon,nplkcon     fields describing kinematic hardening of
!                        a plastic material or gap conductance
!                        constants (cf. User's Manual)
!     xstateini(i,j,k)   state variable component i in integration point j
!                        of element k at the start of the present increment
!     xstate(i,j,k)      state variable component i in integration point j
!                        of element k at the end of the present increment
!     npmat_             maximum number of plastic constants
!     matname(i)         name of material i
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     mi(3)              max number of layers in the structure
!     ielas              0: no elastic iteration: irreversible effects
!                           are allowed
!                        1: elastic iteration, i.e. no irreversible
!                           deformation allowed
!     icmd               not equal to 3: calculate stress and stiffness
!                        3: calculate only stress
!     ncmat_             max number of elastic constants
!     nstate_            max number of state variables in any integration
!                        point
!     stiini(i,j,k)      stress component i in integration point j
!                        of element k at the start of the present
!                        increment (= end of last increment)
!     vini(0..mi(2),i)   value of the variables in node i at the start
!                        of the present increment
!     enerini(1,j,k)     internal energy density at integration point j
!                        of element k at the start of the present increment
!     istep              current step number
!     iinc               current increment number within the actual step
!     reltime            relative step time (between 0 and 1)
!     calcul_fn          if 0: no nodal forces have to be calculated
!                        else: nodal forces are required on output
!     calcul_qa          if 0: no mean forces have to be calculated
!                        else: mean forces are required on output
!     calcul_cauchy      if 0: no Cauchy stresses are required
!                        else: Cauchy stresses are required on output: have
!                        to be calculated from the PK2 stresses
!     nener              if 0: internal energy calculation is not required
!                        else: internal energy is required on output
!     ikin               if 0: kinetic energy calculation is not requred
!                        else: kinetic energy is required on output
!     ne0                largest element number without contact elements (are
!                        stored after all other elements)
!     thicke(j,i)        layer thickness for layer j in element i
!     emeini(i,j,k)      mechanical strain component i in integration point j
!                        of element k at the start of the present increment
!     i                  actual element at stake
!     ielprop(i)         points to the location in field prop preceding
!                        the properties of element i
!     prop(*)            contains the properties and some beam 
!                        elements (cf. User's Manual)
!     
!     
!     OUTPUT:
!     
!     stx(i,j,k)         stress component i in integration point j
!                        of element k at the end of the present
!                        iteration
!     eme(i,j,k)         mechanical strain component i in integration point j
!                        of element k at the end of the present iteration
!     fn(j,i)            nodal force component j in node i
!     qa(1..4)           qa(1): average force
!                        qa(2): average flux
!                        ... cf. User's Manual
!     xstiff(i,j,k)      stiffness (i=1...21) and conductivity
!                        (i=22..27) constants in integration point j
!                        of element k
!     ener(1,j,k)        internal energy density at integration point j
!                        of element k at the end of the present increment
!     eei(i,j,k)         total strain component i in integration point j
!                        of element k at the end of the present iteration
!                        (only for iout>0)
!     nal                number of nodal force contributions
!     
      implicit none
!
      character*8 lakon(*)
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),mi(*),
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &     ielorien(mi(3),*),ntmat_,ipkon(*),ne0,nelem,
     &     istep,iinc,ne,mattyp,ithermal(*),iprestr,i,j,k,m1,m2,jj,
     &     i1,kk,nener,indexe,nope,norien,iperturb(*),iout,
     &     nal,icmd,ihyper,nmethod,kode,imat,iorien,ielas,
     &     istiff,ncmat_,nstate_,ikin,ielprop(*),
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &     calcul_cauchy,calcul_qa,index,node,jjj,j1,j2
     !     
      real*8 co(3,*),v(0:mi(2),*),stiini(6,mi(1),*),
     &     stx(6,mi(1),*),xl(3,26),vl(0:mi(2),26),stre(6),prop(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &     alcon(0:6,ntmat_,*),vini(0:mi(2),*),q1,
     &     alzero(*),orab(7,*),elas(21),rho,fn(0:mi(2),*),
     &     q(0:mi(2),26),t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &     vold(0:mi(2),*),eloc(9),elconloc(ncmat_),eth(6),coords(3),
     &     ener(2,mi(1),*),emec(6),eei(6,mi(1),*),enerini(2,mi(1),*),
     &     veold(0:mi(2),*),e,un,um,tt,dl,qa(*),t0l,t1l,dtime,time,
     &     ttime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &     xstateini(nstate_,mi(1),*),tm(3,3),a,reltime,kap,
     &     thicke(mi(3),*),emeini(6,mi(1),*),aly,alz,bey,bez,xi(2),
     &     vlp(6,2),aa,Bs(2,18),Bb(3,18),Bm(3,18),tau(2),Smf(3),
     &     h,xg(3,3),x(3,3),ushell(18),fintg(18),ueg(18),pgauss(3),
     &     Ds(2,2),Qs(2,2),Qin(3,3),Dm(3,3),Db(3,3),temp_grad0,
     &     Kp(18,18),Km(18,18),tmg(18,18),Kshell(18,18),
     &     alp(3),temp_grad,thick_gp(2),ftherm(18),beta(6),
     &     d3,Ae,tmt(3,3),Ys(2),Emf(3),gpthick(3,2),
     &     Ethm(3),Ethf(3),Em(3),Ef(3),alph6(6),t0g(2,*),t1g(2,*),
     &     Dsi(2,2),Dmi(3,3),Dbi(3,3),di,dett,dettt
!
      nope = 3 
!
      Bs(:,:) = 0.d0
      Bb(:,:) = 0.d0
      Bm(:,:) = 0.d0
      beta(:) = 0.d0
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
      indexe=ipkon(i)
!
! properties of the shell section
!
      index = ielprop(i)
      h = prop(index+1)
      dett  = h/2.d0
      dettt = h**3/8.d0
!
! material and orientation
!
      imat   = ielmat(1,i)
      amat   = matname(imat)
!
      if(norien.gt.0) then
        iorien = max(0,ielorien(1,i))
      else
        iorien = 0
      endif
!      
      if(nelcon(1,imat).lt.0) then
        ihyper = 1
      else
        ihyper = 0
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
!     computation of the coordinates and displacements of the local nodes
!
! vl: u,v,w,rx,ry,rz
      do k = 1,nope
        konl(k) = kon(indexe+k)
        do j = 1,3
          xg(k,j)   = co(j,konl(k)) ! global coordinates element nodes
          vl(j,k)   = v(j,konl(k)) ! uvw
          vl(j+3,k) = v(j+3,konl(k)) ! rxyz	    
        enddo
      enddo
!     
!     gp temp gradient based nodal gradients 
      if(ithermal(1).ne.0) then
        temp_grad  = (t1g(1,konl(1))+t1g(1,konl(2))+t1g(1,konl(3)))*di
        temp_grad0 = (t0g(1,konl(1))+t0g(1,konl(2))+t0g(1,konl(3)))*di
      endif
!
! e0-frame base tm (3x3) -> tmg (18x18)
      call us3_csys(xg,tm,tmg)  ! | call us3_csys_cr(xg,tm,tmg)
      call us3_xu(x,ushell,ueg,xg,tm,vl(1:6,1:3)) ! coords,displ in e0-frame & global displ
!
!  stress and strains
!
      call us3_Bs(x,Bs)         ! constant
      call us3_Bb(x(:,1),x(:,2),Bb) ! constant
!
      Dm(:,:) = 0.d0
      Db(:,:) = 0.d0
      Ds(:,:) = 0.d0
!
      do jjj = 1,3              ! ----  start loop gauß points ----
!
! natural coordinates z-dirc. -1,0,1
!
        if(jjj.eq.1) then
          aa = -0.5d0*h
        elseif(jjj.eq.2) then
          aa = 0.d0
        else
          aa = 0.5d0*h
        endif
!
        if(ithermal(1).ne.0) then
          t0l = (t0(konl(1))+t0(konl(2))+t0(konl(3)))/3.d0
     &         +temp_grad0*aa
          t1l = (t1(konl(1))+t1(konl(2))+t1(konl(3)))/3.d0
     &         +temp_grad*aa
        endif
!
        istiff = 0
        call us3_materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &       nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &       elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &       istiff,elconloc,eth,kode,plicon,nplicon,
     &       plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jjj,
     &       xstiff,ncmat_)     
!
        if(mattyp.eq.1) then
          e   = elconloc(1)
          un  = elconloc(2)
          rho = rhcon(1,1,imat)
          alp(1) = eth(1)       !alcon(1,1,imat)    
          alp(2) = eth(1)       !alcon(1,1,imat)    
          alp(3) = 0.d0 
        elseif(mattyp.eq.2) then
          write(*,*) '*ERROR in e_c3d_US3: no orthotropic material'
          write(*,*) '       calculation for this type of element'
          call exit(201)
        else
          write(*,*) '*ERROR in e_c3d_US3: no anisotropic material'
          write(*,*) '       calculation for this type of element'
          call exit(201)
        endif
!
        call us3_linel_Qi(e,un,Qin,Qs) 
        Dm = Dm + Qin*gpthick(jjj,2)*dett
        Db = Db + Qin*((gpthick(jjj,1))**2)*gpthick(jjj,2)*dettt
        Ds = Ds + Qs*gpthick(jjj,2)*dett 
!
!  save elastic constants in xstff -> e_c3d_us3
!
! material data; for linear elastic materials
! this includes the calculation of the stiffness
! matrix 
        kode = 2
        call linel(kode,mattyp,beta,eme,stre,elas,elconloc,
     &       iorien,orab,pgauss)
        do m1=1,21
          xstiff(m1,jjj,i) = elas(m1) ! elas for each gp saved in xstiff    
        enddo
!
        call bm_ANDES(x,Bm,Qin,h,di,di)
!
! Total strains
!
        Ys  = matmul(Bs,ushell)
        Em  = matmul(Bm,ushell)
        Ef  = aa*(matmul(Bb,ushell))
        Emf = Em + Ef
!
! -> eei total strains out
!
        if(iout.gt.0) then 
          eei(1,jjj,i)  = Emf(1)
          eei(2,jjj,i) = Emf(2)
          eei(3,jjj,i) = 0.d0
          eei(4,jjj,i) = Emf(3)
          eei(5,jjj,i) = Ys(2)
          eei(6,jjj,i) = Ys(1) 
        endif   
!
        if(ithermal(1).ne.0) then
          t0l = (t0(konl(1))+t0(konl(2))+t0(konl(3)))/3.d0+temp_grad0*aa
          t1l = (t1(konl(1))+t1(konl(2))+t1(konl(3)))/3.d0+temp_grad*aa
!
          Ethm = alp*(t1l-t0l)  ! mem
          Ethf = alp*(temp_grad-temp_grad0)*aa ! bend.
          Emf = Emf-(Ethf+Ethm) 
        endif
!
! mechanical strains
!
        eme(1,jjj,i) = Emf(1)
        eme(2,jjj,i) = Emf(2)
        eme(3,jjj,i) = 0.d0
        eme(4,jjj,i) = Emf(3)
        eme(5,jjj,i) = Ys(2)
        eme(6,jjj,i) = Ys(1)        
!
        Smf = matmul(Qin,Emf) 
        tau = matmul(Qs,Ys)
!
        stx(1,jjj,i) = Smf(1)   ! sxx
        stx(2,jjj,i) = Smf(2)   ! syy
        stx(3,jjj,i) = 0.d0     ! szz
        stx(4,jjj,i) = Smf(3)   ! txy
        stx(5,jjj,i) = tau(1)   ! tyz
        stx(6,jjj,i) = tau(2)   ! txz
!
        stre(1) = Smf(1)        ! sxx
        stre(2) = Smf(2)        ! syy
        stre(3) = 0.d0          ! szz
        stre(4) = Smf(3)        ! txy
        stre(5) = tau(1)        ! tyz
        stre(6) = tau(2)        ! txz
!
      enddo                     ! ----  end loop gauß points ----
!
! Thermal loads
!   
      if(ithermal(1).ne.0) then
!
        call us3_ae(x,Ae)       ! element area
        ftherm(:) = 0.d0
!
        do k = 1,3
!
          t0l = ((t0(1)+t0(2)+t0(3))/3.d0) !+temp_grad0*0.5*h*gpthick(k,1)
          t1l = ((t1(1)+t1(2)+t1(3))/3.d0) !+temp_grad*0.5*h*gpthick(k,1)
!
          istiff=1
          call us3_materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &         nalcon,imat,amat,iorien,coords,orab,ntmat_,elas,rho,
     &         i,ithermal,alzero,mattyp,t0l,t1l,
     &         ihyper,istiff,elconloc,eth,kode,plicon,
     &         nplicon,plkcon,nplkcon,npmat_,
     &         plconloc,mi(1),dtime,k,
     &         xstiff,ncmat_)
          e     = elconloc(1)
          un    = elconloc(2)
          rho   = rhcon(1,1,imat)        
          alp(1) = eth(1)       !alcon(1,1,imat)    
          alp(2) = eth(2)       !alcon(1,1,imat)    
          alp(3) = 0.d0 
!
          call us3_linel_Qi(e,un,Qin,Qs)
!
          ftherm=ftherm-matmul(matmul(transpose(Bm),Qin),alp)
     &         *(t1l-t0l)*Ae*gpthick(k,2)*dett
          ftherm=ftherm-matmul(matmul(transpose(Bb),Qin) ,alp)*
     &         (temp_grad-temp_grad0)*Ae*(gpthick(k,1)**2)*dettt*
     &         gpthick(k,2)
        enddo
!     t0l = ((t0(1)-t0(2)-t0(3))/3.d0)!+temp_grad0*0.5*h*g_p(k,1)
!     t1l = ((t1(1)+t1(2)+t1(3))/3.d0)!+temp_grad*0.5*h*g_p(k,1)
!     ftherm=ftherm-matmul(matmul(transpose(Bm),Dm),alp)
!     &     *(t1l-t0l)*Ae
!     ftherm=ftherm-matmul(matmul(transpose(Bb),Db)
!     &     ,alp)*(temp_grad-temp_grad0)*Ae*h
        do k = 1,nope
          konl(k) = kon(indexe+k)
          do j = 1,3
            fn(j,konl(k))   = fn(j,konl(k))  +ftherm(j+(k-1)*6) ! f_th
            fn(j+3,konl(k)) = fn(j+3,konl(k))+ftherm(j+3+(k-1)*6) ! r_th	    
          enddo         
        enddo
      endif       
!
! Internal forces based on displacments
!
      if(calcul_fn.eq.1)then     
!
!   stiffness matirx (6 dofs 3 nodes -> 18x18)  
!
        call us3_Kp(x,Db,Ds,Kp) ! plate part (CS-DSG) - in e0-frame
        call us3_Km(x,Km,Dm/h,h) ! membrane part (ANDES) - in e0-frame
!
        Kshell = matmul(matmul(transpose(tmg),(Km+Kp)),tmg)      
        fintg  = matmul(Kshell,ueg)
        do k = 1,nope
          konl(k) = kon(indexe+k)
          do j = 1,3
            fn(j,konl(k))   = fn(j,konl(k))  +fintg(j+(k-1)*6) ! fi
            fn(j+3,konl(k)) = fn(j+3,konl(k))+fintg(j+3+(k-1)*6) ! ri	    
          enddo         
        enddo
      endif
!
!     q contains the contributions to the nodal force in the nodes
!     belonging to the element at stake from other elements (elements
!     already treated). These contributions have to be
!     subtracted to get the contributions attributable to the element
!     at stake only
!
      if(calcul_qa.eq.1) then
        do m1=1,nope
          do m2=1,3
            qa(1)=qa(1)+dabs(fn(m2,konl(m1))-q(m2,m1))
          enddo
        enddo
        nal=nal+3*nope
      endif
!
      return
      end
