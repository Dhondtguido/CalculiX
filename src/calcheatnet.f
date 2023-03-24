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
      subroutine calcheatnet(nelem,lakon,ipkon,kon,v,ielprop,prop,
     &  ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,ipobody,ibody,xbody,
     &  mi,nacteq,bc,qat,nalt)
!
!     user subroutine film
!
!
!     INPUT:
!
!     nelem              actual element number
!     lakon(i)           contains the label of element i
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     v(0..4,1..nk)      actual solution field in all nodes; 
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
!     ntmat_             maximum number of temperature data points for 
!                        any material property for any material
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
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     mi(3)              max # of layers in any element
!     nacteq(j,i)      contains the number of an equation expressed at
!                        network node i (i=1..,ntg, ntg is the number of
!                        network nodes)
!                        j=0: energy conservation
!                        j=1: mass conservation
!                        j=2: element equation
!                        j=3: geometric equation
!
!     OUTPUT:
!
!     bc(*)              right hand side of the network equation system
!     qat                sum of energy contributions in the network
!     nalt               number of energy contributions: used to calculate
!                        a typical energy flow (used in the convergence
!                        criteria)
!
      implicit none
!
      character*8 lakon(*)
!
      integer nelem,ipkon(*),kon(*),nshcon(*),nrhcon(*),ipobody(2,*),
     &  ibody(3,*),mi(*),ntmat_,ielprop(*),ielmat(mi(3),*),node1,
     &  node2,nodem,imat,nalt,ieq,nacteq(0:3,*),index
!
      real*8 shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),xbody(7,*),heat,
     &  v(0:mi(2),*),prop(*),xflow,Uout,Tg1,Tg2,Rin,Rout,R1,R2,r,cp,rho,
     &  dvi,Uin,gastemp,cp_cor,bc(*),qat,om
!
!
!
      heat=0.d0
!     
      nodem=kon(ipkon(nelem)+2)
      xflow=v(1,nodem)
      node1=kon(ipkon(nelem)+1)
      node2=kon(ipkon(nelem)+3)
      if((node1.eq.0).or.(node2.eq.0)) return
!
      if(lakon(nelem)(2:3).eq.'VO') then
!
         index=ielprop(nelem)
!         
         if(xflow.gt.0d0) then
            R1=prop(index+2)
            R2=prop(index+1)
            Rout=R2
            Rin=R1
         else
            R1=prop(index+2)
            R2=prop(index+1)
            Rout=R1
            Rin=R2
         endif
!     
!     computing temperature corrected Cp=Cp(T) coefficient
!     
         Tg1=v(0,node1)
         Tg2=v(0,node2)
         if((lakon(nelem)(2:3).ne.'LP').and.
     &      (lakon(nelem)(2:3).ne.'LI')) then
            gastemp=(Tg1+Tg2)/2.d0
         else
            if(xflow.gt.0) then
               gastemp=Tg1
            else
               gastemp=Tg2
            endif
         endif
!
         imat=ielmat(1,nelem)
         call materialdata_tg(imat,ntmat_,gastemp,
     &        shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
!
         call cp_corrected(cp,Tg1,Tg2,cp_cor)
!     
c         Uout=prop(index+5)*Rout
c         Uin=prop(index+5)*Rin
!     
!     free and forced vortices with temperature 
!     change in the relative system of coordinates 
!     
         if((lakon(nelem)(2:5).eq.'VOFR') .and.
     &        (nint(prop(index+8)).eq.(-1))) then
!     
            Uout=prop(index+7)*Rout
            Uin=prop(index+7)*Rin
!     
            heat=0.5d0*Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
         elseif (((lakon(nelem)(2:5).eq.'VOFO')
     &           .and.(nint(prop(index+6)).eq.(-1)))) then
!     
            Uout=prop(index+5)*Rout
            Uin=prop(index+5)*Rin
!     
            heat=0.5d0*Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
!     forced vortices with temperature change in the absolute system
!     
         elseif((lakon(nelem)(2:5).eq.'VOFO')
     &           .and.((nint(prop(index+6)).eq.1))) then
!     
            Uout=prop(index+5)*Rout
            Uin=prop(index+5)*Rin
            heat=Cp/Cp_cor*(Uout**2-Uin**2)*xflow
!     
         endif
      elseif(lakon(nelem)(2:5).eq.'GAPR') then
!
!        heat production in a rotating pipe (in the relative
!        system)         
!         
         index=ielprop(nelem)
!         
         R1=prop(index+8)
         R2=prop(index+9)
         om=prop(index+10)
         if(xflow.gt.0.d0) then
            Rin=R1
            Rout=R2
         else
            Rin=R2
            Rout=R1
         endif
         Uin=Rin*om
         Uout=Rout*om
!     
         heat=(Uout**2-Uin**2)*xflow/2.d0
!     
      elseif(lakon(nelem)(2:2).eq.'U') then
!
!        insert here the heat generated in user defined network elements
!
!        START insert
!
         heat=0.d0
!
!        END insert
! 
      else
         return
      endif
!     
!     including the resulting additional heat flux in the energy equation
!     
      if(xflow.gt.0d0)then
         ieq=nacteq(0,node2)
         if(ieq.ne.0) then
            bc(ieq)=bc(ieq)+heat
            qat=qat+dabs(heat)
            nalt=nalt+1
         endif
      else
         ieq=nacteq(0,node1)
         if(ieq.ne.0) then
            bc(ieq)=bc(ieq)+heat
            qat=qat+dabs(heat)
            nalt=nalt+1
         endif
      endif
!     
      return
      end

