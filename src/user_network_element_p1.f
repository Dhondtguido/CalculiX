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
      subroutine user_network_element_p1(node1,node2,nodem,nelem,lakon,
     &     kon,ipkon,nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
!     
!     UP1 element: mass flow = c * dsqrt(pt2-pt1)
!
!     f:= c * dsqrt(pt2-pt1) - mass flow
!
!     for this element it is known that the flow direction 
!     has to correspond to the pressure drop
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
      pi=4.d0*datan(1.d0)   
      if(kflag.eq.0) then
!
!        called by envtemp.f:
!
!        check whether element equation is needed (this is the
!        case if identity=.false.
!
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
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
         index=ielprop(nelem)
         c=prop(index+2)
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
         if(pt1.ge.pt2) then
            inv=1
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
         endif
!     
         xflow=c*dsqrt(pt1-pt2)
!     
      elseif(kflag.eq.2)then
!
!        called by resultnet.f and mafillnet.f
!     
         numf=3
         index=ielprop(nelem)
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
         if(pt1.ge.pt2) then
!
!           flow direction corresponds to element orientation
!
            inv=1
            xflow=v(1,nodem)*iaxial
            nodef(1)=node1
            nodef(2)=nodem
            nodef(3)=node2
         else
!
!           flow direction does not correspond to element orientation
!
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            nodef(1)=node2
            nodef(2)=nodem
            nodef(3)=node1
         endif
!     
         idirf(1)=2
         idirf(2)=1
         idirf(3)=2
!
!        nodef and idirf show that the variables are (if inv=1):
!        1. total pressure in node 1
!        2. mass flow
!        3. total pressure in node 2
!
         c=prop(index+2)
!
         f=c*dsqrt(pt1-pt2)-xflow
         df(1)=c/(2.d0*dsqrt(pt1-pt2))
         df(2)=-1.d0*inv
         df(3)=-df(1)
!     
!     output
!     
      elseif(kflag.eq.3) then
!
!        called by flowoutput.f
!
!        storage in the .net-file
!     
         index=ielprop(nelem)
         a=prop(index+1)
         kappa=(cp/(cp-R))
         icase=1
!
         pt1=v(2,node1)
         pt2=v(2,node2)
         if(pt1.ge.pt2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)-physcon(1)
            Tt2=v(0,node2)-physcon(1)
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial 
            Tt1=v(0,node2)-physcon(1)
            Tt2=v(0,node1)-physcon(1)
         endif
!
!        calculation of the static temperatures
!
         if(lakon(nelem)(3:3).eq.'P') then
!
!           "pipe"-like element: total and static temperatures
!           differ
!
            call ts_calc(xflow,Tt1,pt1,kappa,r,A,T1,icase)
            call ts_calc(xflow,Tt2,pt2,kappa,r,A,T2,icase)
         else
!
!           "chamber-connecting"-like element: total and static
!           temperatures are equal
!
            T1=Tt1
            T2=Tt2
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1d-30) then
            write(*,*) '*ERROR in orifice: '
            write(*,*) '       no dynamic viscosity defined'
            write(*,*) '       dvi= ',dvi
            call exit(201)
         endif 
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         a=prop(index+1)
         d=dsqrt(a*4/Pi)           
         reynolds=dabs(xflow)*d/(dvi*a)
         xmach=dsqrt(((pt1/pt2)**((kappa-1.d0)/kappa)-1.d0)*2.d0/
     &          (kappa-1.d0))
!
         xflow_oil=0.d0
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
!     degree of freedom 2 is the mass flow dof (cf. field idirf)
!
      xflow=xflow/iaxial
      df(2)=df(2)*iaxial
!     
      return
      end
