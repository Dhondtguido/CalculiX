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
      subroutine user_network_element_p0(node1,node2,nodem,nelem,lakon,
     &     kon,ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
!     
!     user subroutine user_network_element
!
!     skeleton file
!
!     INPUT:
!
!     node1              first node in element topology
!     node2              third node in element topology
!     nodem              second node in element topology (middle node)
!     nelem              element number
!     lakon(i)           label of element i
!     kon                connectivity list of all elements; the topology
!                        of element i starts at kon(ipkon(i))
!     ipkon(i)           pointer of element i into list kon
!     nactdog(j,i)       global degree of freedom in the network
!                        equation system of local dof j (0,1 or 2 for
!                        networks) of node i. If nactdog(j,i)=0 the
!                        variable is known
!     ielprop(i)         pointer for element i into field prop. The 
!                        properties for element i start at 
!                        prop(ielprop(i)+1,...)
!     prop               field of all properties
!     kflag              indicates what information should be returned
!                        by the routine:
!                        0: identity
!                        1: xflow
!                        2: numf, nodef, idirf, f, df
!                        3: none
!     v(0..mi(2),i)      values at node i in the current network
!                        iteration (0=total temperature,
!                        1=mass flow, 2=total pressure for network nodes,
!                        0=temperature, 1..3=displacements for structural
!                        nodes)
!     cp                 specific heat at constant pressure corresponding
!                        to a mean static temperature across the element
!     R                  specific gas constant
!     physcon(1..)       physical constants (e.g. physcon(1) is absolute
!                        zero in the unit systemof the user; cf. the
!                        user's manual for the other entries)
!     dvi                dynamical viscosity corresponding
!                        to a mean static temperature across the element
!     set(i)             set name corresponding to set i
!     co(1..3,i)         coordinates of node i in the global system
!     vold(0..mi(2),i)   values at node i at the start of the current network
!                        iterations (0=total temperature,
!                        1=mass flow, 2=total pressure for network nodes,
!                        0=temperature, 1..3=displacements for structural
!                        nodes)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...; the other
!                        values of mi are not relevant here
!     ttime              total time at the end of the current
!                        thermo-mechanical increment. To reach the end
!                        of this increment several thermo-mechanical
!                        iterations are performed. For each of these
!                        iterations a loop of network iterations is 
!                        performed
!     time               step time a the end of the current thermo-
!                        mechanical increment
!     iaxial             number of times the current structure fits into
!                        360 degrees
!     iplausi            flag telling whether any plausibility checks
!                        have been violated up to entry in this routine
!                        0: plausibility checks not satisfied
!                        1: plausibility checks (if any) are satisfied
!     
!
!     OUTPUT:
!
!     identity           if .true. the user_network_element routine is
!                        not needed (all variables known)
!     xflow              mass flow
!     f                  value of the element equation
!     nodef              nodes corresponding to the variables in the 
!                        element equation
!     idirf              degrees of freedom corresponding to the variables
!                        in the element equation
!     df                 derivatives of the element equation w.r.t. its
!                        variables
!     numf               number of variables in the element equation
!     iplausi            flag telling whether any plausibility checks
!                        were violated at return time from this routine
!                        0: plausibility checks not satisfied
!                        1: plausibility checks (if any) are satisfied
!                        only feasible change within this routine is from
!                        1 to 0.
!     
!     NOTE: to convert total temperature into static temperatures
!           subroutine 
!           call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
!           may be used (cf. user_netowrk_element_p1.f for an example
!           of its use).
!
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),kflag,ipkon(*),kon(*),
     &     iaxial,mi(*),inv,index,icase,iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),R,cp,physcon(*),dvi,
     &     co(3,*),vold(0:mi(2),*),ttime,time,xmach,xflow_oil,Tt1,Tt2,
     &     reynolds,pt1,pt2,kappa,c,a,d,pi,T1,T2
!
!
!
      if(kflag.eq.0) then
!
!        called by envtemp.f:
!
!        check whether element equation is needed (this is the
!        case if identity=.false.
!
!         identity=?
!     
      elseif(kflag.eq.1)then
         if(v(1,nodem).ne.0.d0) then
            xflow=v(1,nodem)
            return
         endif
!
!        called by initialnet.f: 
!
!        calculation of the mass flow if everything else is known
!     
!         xflow=?
!     
      elseif(kflag.eq.2)then
!
!        called by resultnet.f and mafillnet.f
!     
!         numf=?
!         nodef(1...numf)=?
!         idirf(1...numf)=?
!         f=?
!         df(1...numf)=?
!     
      elseif(kflag.eq.3) then
!
!        called by flowoutput.f
!
!        storage in the .net-file
!        this is an example
!     
         write(1,*) ''
         write(1,55) ' from node ',node1,
     &   ' to node ', node2,' :   air massflow rate = '
     &,inv*xflow,' ',
     &   ', oil massflow rate = ',xflow_oil,' '
!          
         if(inv.eq.1) then
            write(1,56)'       Inlet node ',node1,' :   Tt1 = ',Tt1,
     &           '  , Ts1 = ',T1,'  , Pt1 = ',pt1, ' '
!             
            write(1,*)'             Element ',nelem,lakon(nelem)
            write(1,57)'             dyn.visc = ',dvi,'  , Re = '
     &           ,reynolds,', M = ',xmach
!
            write(1,56)'      Outlet node ',node2,' :   Tt2 = ',Tt2,
     &           '  , Ts2 = ',T2,'  , Pt2 = ',pt2,' '
!     
         else if(inv.eq.-1) then
            write(1,56)'       Inlet node ',node2,':    Tt1 = ',Tt1,
     &           '  , Ts1 = ',T1,' , Pt1 = ',pt1, ' '
     &          
            write(1,*)'             element R    ',set(numf)(1:30)
            write(1,57)'             dyn.visc. = ',dvi,' , Re ='
     &           ,reynolds,', M = ',xmach
! 
            write(1,56)'      Outlet node ',node1,':    Tt2 = ',Tt2,
     &           '  , Ts2 = ',T2,'  , Pt2 = ',pt2, ' '
         endif
!
      endif
!     
 55   format(1x,a,i6,a,i6,a,e11.4,a,a,e11.4,a)
 56   format(1x,a,i6,a,e11.4,a,e11.4,a,e11.4,a)
 57   format(1x,a,e11.4,a,e11.4,a,e11.4)
!  
!     following lines are needed because the element equations
!     usually specify the mass flow for the complete cross section
!
!     in CalculiX, axisymmetric elements are expanded into 3D
!     using a sector of 360°/iaxial. Therefore, the mass flow and
!     the derivative of f w.r.t. the mass flow have to be adjusted
!     appropriately.     
!
      xflow=xflow/iaxial
!
!     only if the mass flow is an active degree of freedom:
!
!     df(mass_flow_dof)=df(mass flow dof)*iaxial
!     
      return
      end
