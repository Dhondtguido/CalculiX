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
      subroutine film(h,sink,temp,kstep,kinc,time,noel,npt,
     &  coords,jltyp,field,nfield,loadtype,node,area,vold,mi,
     &  ipkon,kon,lakon,iponoel,inoel,ielprop,prop,ielmat,
     &  shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon,
     &  ipobody,xbody,ibody,heatnod,heatfac)
!
!     user subroutine film
!
!
!     INPUT:
!
!     sink               most recent sink temperature
!     temp               current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     field              currently not used
!     nfield             currently not used (value = 1)
!     loadtype           load type label
!     node               network node (only for forced convection)
!     area               area covered by the integration point
!     vold(0..4,1..nk)   actual solution field in all nodes; 
!                        for structural nodes:
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!                        for network nodes
!                        0: total temperature (at end nodes)
!                           = static temperature for liquids
!                        1: mass flow (at middle nodes)
!                        2: total pressure (at end nodes)
!                           = static pressure for liquids
!                        3: static temperature (at end nodes; only for gas)
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     mi(3)              max # of layers in any element
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     lakon(i)           contains the label of element i
!     iponoel(i)         the network elements to which node i belongs
!                        are stored in inoel(1,iponoel(i)),
!                        inoel(1,inoel(2,iponoel(i)))...... until
!                        inoel(2,inoel(2,inoel(2......)=0
!     inoel(1..2,*)      field containing the network elements
!     ielprop(i)         points to the location in field prop preceding
!                        the properties of element i
!     prop(*)            contains the properties of all network elements. The
!                        properties of element i start at prop(ielprop(i)+1)
!                        and continues until all properties are covered. The
!                        appropriate amount of properties depends on the
!                        element label. The kind of properties, their
!                        number and their order corresponds
!                        to the description in the user's manual,
!                        cf. the sections "Fluid Section Types"
!     ielmat(j,i)        contains the material number for element i
!                        and layer j
!     shcon(0,j,i)       temperature at temperature point j of material i
!     shcon(1,j,i)       specific heat at constant pressure at the
!                        temperature point j of material i
!     shcon(2,j,i)       dynamic viscosity at the temperature point j of
!                        material i
!     shcon(3,1,i)       specific gas constant of material i
!     nshcon(i)          number of temperature data points for the specific
!                        heat of material i
!     rhcon(0,j,i)       temperature at density temperature point j of 
!                        material i
!     rhcon(1,j,i)       density at the density temperature point j of
!                        material i
!     nrhcon(i)          number of temperature data points for the density
!                        of material i
!     ntmat_             maximum number of temperature data points for 
!                        any material property for any material
!     ncocon(1,i)        number of conductivity constants for material i
!     ncocon(2,i)        number of temperature data points for the 
!                        conductivity coefficients of material i
!     cocon(0,j,i)       temperature at conductivity temperature point
!                        j of material i
!     cocon(k,j,i)       conductivity coefficient k at conductivity
!                        temperature point j of material i
!     ipobody(1,i)       points to an entry in fields ibody and xbody 
!                        containing the body load applied to element i, 
!                        if any, else 0
!     ipobody(2,i)       index referring to the line in field ipobody
!                        containing a pointer to the next body load
!                        applied to element i, else 0
!     ibody(1,i)         code identifying the kind of body load i:
!                        1=centrifugal, 2=gravity, 3=generalized gravity
!     ibody(2,i)         amplitude number for load i
!     ibody(3,i)         load case number for load i
!     xbody(1,i)         size of body load i
!     xbody(2..4,i)      for centrifugal loading: point on the axis,
!                        for gravity loading with known gravity vector:
!                          normalized gravity vector
!     xbody(5..7,i)      for centrifugal loading: normalized vector on the
!                          rotation axis
!
!     OUTPUT:
!
!     h(1)               magnitude of the film coefficient
!     h(2)               not used; please do NOT assign any value
!     sink               (updated) sink temperature (need not be
!                        defined for forced convection)
!     heatnod            extra heat flow going to the network node
!                        (zero if not specified)
!     heatfac            extra heat flow going to the structural face
!                        (zero if not specified)
!           
      implicit none
!
      character*8 lakon(*)
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,nfield,node,mi(*),ipkon(*),
     &  kon(*),iponoel(*),inoel(2,*),ielprop(*),ielmat(mi(3),*),ntmat_,
     &  nshcon(*),nrhcon(*),ncocon(2,*),nodem,indexprop,indexe,
     &  iel1,iel2,ielup,iit,imat,icase,itherm,ipobody(2,*),
     &  ibody(3,*)
!
      real*8 h(2),sink,time(2),coords(3),temp,field(nfield),area,
     &  vold(0:mi(2),*),prop(*),shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  cocon(0:6,ntmat_,*),rho,r,pt,Pr,xl,Tt,Ts,xflow,Tsold,Re,um,
     &  xks,xkappa,xlambda,f,cp,A,D,form_fact,xbody(7,*),heatnod,
     &  heatfac
!
      if((node.eq.0).or.(iponoel(node).eq.0)) then
!
!     simple example: constant film coefficient
!
         h(1)=200.d0
      else
!
!     complicated example: forced convection with a film coefficient
!     determined by the Gnielinski equation (pipe application for
!     turbulent flow),
!     which runs like:
!
!     NuD=f/8*(ReD-1000)*pr/(1+12.7*sqrt(f/8)*(Pr**(2/3)-1))
!
!     NuD=h*D/xlambda: Nusselt number
!     ReD=massflow*D/(um*A): Reynolds number
!     Pr=um*cp/xlambda: Prandl number
!     h: film coefficient
!     D: hydraulic diameter
!     A: area of the pipe cross section
!     f: friction coefficient (cf. User's manual for formula)
!     xlambda: heat conduction coefficient
!     um: dynamic viscosity
!     cp: heat capacity
!
!     The convection node comes from the calling program: "node"
!
!     The elements connected to this node are determined by:
!     iel1=inoel(1,iponoel(node)),
!     iel2=inoel(1,inoel(2,iponoel(node))),
!     .... until
!     inoel(2,inoel(2,inoel(2,.......inoel(2,iponoel(node)))))=0
!
!     Let us assumed that "node" belongs to exactly two network
!     elements, i.e. inoel(2,inoel(2,iponoel(node)))=0
!     
!     Let us look at the upstream element. To determine the
!     upstream element one looks at the sign of the mass flow
!     in the middle node of the elements.
!
!     The middle node nodem of element iel1 is given by
!     kon(ipkon(iel1)+2), the end nodes by 
!     kon(ipkon(iel1)+1) and kon(ipkon(iel1)+3). The mass
!     flow in the element is given by vold(1,nodem). If this value 
!     is positive and kon(ipkon(iel1)+3)=node or if is negative
!     and kon(ipkon(iel1)+1)=node iel1 is the upstream element
!
!     Let us assume this element is a gas pipe element, i.e.
!     lakon(iel1)='DGAPFA  '. The properties of such an element
!     are the cross sectional area, the hydraulic diameter.... (in 
!     that order), i.e. prop(ielprop(iel1)+1)=A,
!     prop(ielprop(iel1)+2)=D,...
!
!     The material number of element iel1 is given by ielmat(iel1)
!     with this information and the gas temperature (sink=Tfluid) the
!     gas material properties can be determined: xlambda by calling
!     materialdata_cond, um by calling materialdata_dvi and cp by
!     calling materialdata_cp
!
!     This yields all data needed to determine the film coefficient
!     with the Gnielinski equation
!
!     For laminar flow (Re < 3000) the following laminar flow
!     equation is used:
!
!     h=4.36*xlambda/D
!
!        elements belonging to the network node
!
         iel1=inoel(1,iponoel(node))
         if(inoel(2,iponoel(node)).ne.0) then
            iel2=inoel(1,inoel(2,iponoel(node)))
!
!           check whether the node belongs to maximum two elements
!
            if(inoel(2,inoel(2,iponoel(node))).ne.0) then
               write(*,*) 'ERROR in film: the network node'
               write(*,*) '      belongs to more than 2 elements'
               call exit(201)
            endif
         else
!
!           node (must be an end node) belongs only to one
!           element (pure thermal network without entry or
!           exit element)
!
            iel2=iel1
         endif
!
!        only adiabatic gas pipes considered (heat exchange only
!        takes place at the end nodes, not within the pipe)
!
         icase=0
!
!        check whether iel1 is the upstream element
!
         indexe=ipkon(iel1)
         nodem=kon(indexe+2)
         if(((vold(1,nodem).ge.0.d0).and.(kon(indexe+3).eq.node)).or.
     &      ((vold(1,nodem).le.0.d0).and.(kon(indexe+1).eq.node))) then
            ielup=iel1
         else
            ielup=iel2
         endif
!
!        check whether the upstream element is an adiabatic gas
!        pipe element or a White-Colebrook liquid pipe
!
         if((lakon(ielup)(1:6).ne.'DGAPFA').and.
     &      (lakon(ielup)(1:7).ne.'DLIPIWC')) then
            write(*,*) 'ERROR in film: upstream element',ielup
            write(*,*) '      is no adiabatic'
            write(*,*) '      gas pipe nor White-Colebrook'
            write(*,*) '      liquid pipe'
            call exit(201)
         endif
!
         indexprop=ielprop(ielup)
         A=prop(indexprop+1)
         D=prop(indexprop+2)
         xl=prop(indexprop+3)
         xks=prop(indexprop+4)
         form_fact=prop(indexprop+5)
!
!        material of upstream element
!
         imat=ielmat(1,ielup)
!
         xflow=dabs(vold(1,nodem))
         Tt=vold(0,node)
         pt=vold(2,node)
!
         if(lakon(ielup)(2:2).eq.'G') then
!
!           compressible flow
!
            iit=0
            Ts=Tt
            do
               iit=iit+1
               Tsold=Ts
               call materialdata_tg(imat,ntmat_,Ts,shcon,nshcon,cp,r,
     &              um,rhcon,nrhcon,rho)
               call ts_calc(xflow,Tt,pt,xkappa,r,A,Ts,icase)
               if((dabs(Ts-Tsold).le.1.d-5*Ts).or.(iit.eq.10)) exit
            enddo
!
            call materialdata_tg(imat,ntmat_,Ts,shcon,nshcon,cp,r,
     &           um,rhcon,nrhcon,rho)
         else
!
!           incompressible flow
!
            Ts=Tt
            call materialdata_cp(imat,ntmat_,Ts,shcon,nshcon,cp)
!
            itherm=2
            call materialdata_dvi(shcon,nshcon,imat,um,Ts,ntmat_,
     &           itherm)
!
         endif
!
         call materialdata_cond(imat,ntmat_,Ts,cocon,ncocon,xlambda)
!
         Re=xflow*D/(um*A)
!
         if(Re.lt.3000.d0) then
!
!           laminar flow
!
            h(1)=4.36d0*xlambda/D
         else
!
!           turbulent flow
!
            if(Re.gt.5.e6) then
               write(*,*) '*ERROR in film: Reynolds number ',Re
               write(*,*) '       is outside valid range'
               call exit(201)
            endif
!
            Pr=um*cp/xlambda
!
            if((Pr.lt.0.5d0).or.(Pr.gt.2000.d0)) then
               write(*,*) '*ERROR in film: Prandl number ',Pr
               write(*,*) '       is outside valid range'
               call exit(201)
            endif
!
            call friction_coefficient(xl,D,xks,Re,form_fact,f)
!     
            h(1)=f/8.d0*(Re-1000.d0)*Pr/
     &           (1.d0+12.7d0*dsqrt(f/8.d0)*(Pr**(2.d0/3.d0)-1.d0))
     &           *xlambda/D
         endif
      endif
!
      return
      end

