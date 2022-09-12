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
      subroutine calcgeomelemnet(vold,co,prop,lakon,nelem,ttime,
     &  time,ielprop,mi,A,A2,d,l,s)
!     
!     user subroutine
!
!     calculate the cross section for flexible network
!     elements and user defined network elements
!
!     INPUT:
!
!     vold(0..4,1..nk)   solution field in all nodes; 
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
!     co(3,1..nk)        coordinates of all nodes
!                        1: coordinate in global x-direction
!                        2: coordinate in global y-direction
!                        3: coordinate in global z-direction
!     prop(*)            contains the properties of all network elements. The
!                        properties of element i start at prop(ielprop(i)+1)
!                        and continues until all properties are covered. The
!                        appropriate amount of properties depends on the
!                        element label. The kind of properties, their
!                        number and their order corresponds
!                        to the description in the user's manual,
!                        cf. the sections "Fluid Section Types"
!     lakon(i)           contains the label of element i
!     nelem              actual element number
!     ttime              total time at the start of actual thermomechanical
!                        increment
!     time               step time at the end of the actual thermomechanical
!                        increment
!     ielprop(i)         points to the location in field prop preceding
!                        the properties of element i
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedom per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     A2                 presently not used
!     d                  presently not used
!     l                  presently not used
!     s                  presently not used
!
!     OUTPUT:
!
!     A                  cross section
!
      implicit none
!     
      character*8 lakon(*)
!     
      integer nelem,ielprop(*),mi(*),nodea,nodeb,nodec,index
!
      real*8 vold(0:mi(2),*),co(3,*),prop(*),ttime,time,pi,radius,
     &  A,A2,d,l,s
!
!
!
      index=ielprop(nelem)
      pi=4.d0*datan(1.d0)
!
      A=0.d0
!
      if((lakon(nelem)(2:6).eq.'GAPFA').or.
     &   (lakon(nelem)(2:6).eq.'GAPFI')) then
!     
         if(lakon(nelem)(7:8).eq.'FR') then
!     
!     flexible radius
!     
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
!     
            A=pi*radius**2
!     
         elseif(lakon(nelem)(7:8).eq.'RL') then
!     
!     flexible radius and length
!     
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            nodec=nint(prop(index+3))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            A=pi*radius**2
         endif
      elseif(lakon(nelem)(2:3).eq.'UP') then
!
!        determine the relevant cross section (used in the
!        calculation of the static temperature from the
!        total temperature) for user network elements
!
!        START insert
!
         A=prop(index+1)
!
!        END insert
!
      endif
!     
      return
      end
