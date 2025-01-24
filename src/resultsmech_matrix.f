!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2022 Guido Dhondt
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
      subroutine resultsmech_matrix(co,kon,ipkon,lakon,ne,v,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &  reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,nal,ne0,thicke,emeini,i,ielprop,prop,t0g,t1g)
!
!     calculates fn for a substructure (also called superelement,
!                provided using *MATRIX ASSEMBLE
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
!                        are allowed
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
!     enerini(j,k)       internal energy density at integration point j
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
!                              to be calculated from the PK2 stresses
!     nener              if 0: internal energy calculation is not required
!                        else: internal energy is required on output
!     ikin               if 0: kinetic energy calculation is not required
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
!     t0g(1..2,i)        temperature gradient in node i at start of calculation
!     t1g(1..2,i)        temperature gradient in node i at the end of the 
!                        current increment
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
!     ener(j,k)          internal energy density at integration point j
!                        of element k at the end of the present increment
!     eei(i,j,k)         total strain component i in integration point j
!                        of element k at the end of the present iteration
!                        (only for iout>0)
!     nal                number of nodal force contributions
!
      implicit none
!     
      character*8 lakon(*)
      character*80 amat,matname(*),filestiff
!     
      integer kon(*),konl(26),mi(*),kal(2,6),j1,j2,j3,j4,id1,i2,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),node1,
     &     ielorien(mi(3),*),ntmat_,ipkon(*),ne0,node2,idof1,idof2,
     &     istep,iinc,ne,mattyp,ithermal(*),iprestr,i,j,k,m1,m2,jj,
     &     i1,kk,nener,indexe,nope,norien,iperturb(*),iout,irow,
     &     nal,icmd,ihyper,nmethod,kode,imat,iorien,ielas,icol,
     &     istiff,ncmat_,nstate_,ikin,ielprop(*),
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &     calcul_cauchy,calcul_qa,index,node
!     
      real*8 co(3,*),v(0:mi(2),*),stiini(6,mi(1),*),xa(3,3),
     &     stx(6,mi(1),*),xl(3,20),vl(0:mi(2),20),stre(6),prop(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &     alcon(0:6,ntmat_,*),vini(0:mi(2),*),value,
     &     alzero(*),orab(7,*),rho,fn(0:mi(2),*),
     &     q(0:mi(2),20),t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &     vold(0:mi(2),*),eloc(9),elconloc(ncmat_),eth(6),coords(3),
     &     ener(mi(1),*),emec(6),eei(6,mi(1),*),enerini(mi(1),*),
     &     veold(0:mi(2),*),e,un,um,tt,dl,qa(*),t0l,t1l,dtime,time,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &     xstateini(nstate_,mi(1),*),tm(3,3),a,reltime,ttime,
     &     thicke(mi(3),*),emeini(6,mi(1),*),aly,alz,bey,bez,xi(2),
     &     vlp(6,2),xi11,xi12,xi22,xk,offset1,offset2,e1(3),e2(3),e3(3),
     &     t0g(2,*),t1g(2,*)
!     
      real*8,dimension(:),allocatable::fnl
!     
      indexe=ipkon(i)
!     
!     material name = filename for the stiffness matrix
!
      imat=ielmat(1,i)
      filestiff=matname(imat)
!
      nope=ichar(lakon(i)(8:8))
      allocate(fnl(3*nope))
!
      open(20,file=filestiff,status='old',err=2)
!
      do
        read(20,*,end=1) node1,idof1,node2,idof2,value
        call nident(kon(indexe+1),node1,nope,id1)
        irow=(id1-1)*3+idof1
        fnl(irow)=fnl(irow)+value*vold(idof2,node2)
        fn(idof1,node1)=fn(idof1,node1)+value*vold(idof2,node2)
      enddo
!    
 1    close(20)
!    
      if(calcul_qa.eq.1) then
        do m1=1,3*nope
          qa(1)=qa(1)+dabs(fnl(m1))
        enddo
        nal=nal+3*nope
      endif
      deallocate(fnl)
!
      return
 2    write(*,*) '*ERROR reading stiffness matrix file ',filestiff
      write(*,*) '       for element with label ',lakon(i)(1:5)
      write(*,*) '       possible reasons:'
      write(*,*) '       1) incorrect format of the data'
      write(*,*) '       2) user element was not properly coded'
      call exit(201)
      end
