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
      subroutine cload(xload,kstep,kinc,time,node,idof,coords,vold,
     &  mi,ntrans,trab,inotr,veold)
!
!     user subroutine cload
!
!
!     INPUT:
!
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     node               node number
!     idof               degree of freedom
!     coords(1..3)       global coordinates of the node
!     vold(0..mi(2)
!              ,1..nk)   solution field in all nodes (for modal
!                        dynamics: in all nodes for which output
!                        was requested or a force was applied)
!                        (not available for CFD-calculations)
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: not used
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     veold(0..mi(2)
!               ,1..nk)  For non-CFD-calculations:
!                        derivative of the solution field w.r.t.
!                        time in all nodes(for modal
!                        dynamics: in all nodes for which output
!                        was requested or a force was applied)
!                        0: temperature rate
!                        1: velocity in global x-direction
!                        2: velocity in global y-direction
!                        3: velocity in global z-direction
!
!                        For CFD-calculations:
!                        0: temperature
!                        1: velocity in global x-direction
!                        2: velocity in global y-direction
!                        3: velocity in global z-direction
!                        4: static pressure
!
!     ntrans             number of transform definitions
!     trab(1..6,i)       coordinates of two points defining transform i
!     trab(7,i)          -1: cylindrical transformation
!                         1: rectangular transformation 
!     inotr(1,j)         transformation number applied to node j
!     inotr(2,j)         a SPC in a node j in which a transformation
!                        applied corresponds to a MPC. inotr(2,j) 
!                        contains the number of a new node generated
!                        for the inhomogeneous part of the MPC
!
!     OUTPUT:
!
!     xload              concentrated load in direction idof of node
!                        "node" (global coordinates)
!           
      implicit none
!
      integer kstep,kinc,node,idof,mi(*),ntrans,inotr(2,*),itr 
!
      real*8 xload,time(2),coords(3),vold(0:mi(2),*),trab(7,*),
     &  veold(0:mi(2),*),a(3,3),ve1,ve2,ve3,f1,f2,f3
!
!     displacements vold and velocities veold are given in
!     the global system
!
!     example how to transform the velocity into the local system
!     defined in node "node"
!
      if(ntrans.eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
      if(itr.ne.0) then
         call transformatrix(trab(1,itr),coords,a)
         ve1=veold(1,node)*a(1,1)+veold(2,node)*a(2,1)
     &     +veold(3,node)*a(3,1)
         ve2=veold(1,node)*a(1,2)+veold(2,node)*a(2,2)
     &     +veold(3,node)*a(3,2)
         ve3=veold(1,node)*a(1,3)+veold(2,node)*a(2,3)
     &     +veold(3,node)*a(3,3)
!
!     suppose you know the size of the force in local coordinates:
!     f1, f2 and f3. Calculating the size of the force in 
!     direction idof in global coordinates is done in the following
!     way:
!
         xload=f1*a(idof,1)+f2*a(idof,2)+f3*a(idof,3)
      else
!
!        no local system defined in node "node"; suppose the force in
!        global coordinates has components f1, f2 and f3
!
         if(idof.eq.1) then
            xload=f1
         elseif(idof.eq.2) then
            xload=f2
         elseif(idof.eq.3) then
            xload=f3
         endif
      endif
!
      return
      end

