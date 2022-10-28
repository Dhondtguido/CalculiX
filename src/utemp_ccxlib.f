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
      subroutine utemp_ccxlib(temp,msecpt,kstep,kinc,time,node,coords,
     &  vold,mi,iponoel,inoel,ipobody,xbody,ibody)
!
!     user subroutine utemp
!
!
!     INPUT:
!
!     msecpt             number of temperature values (for volume elements:1)
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     node               node number
!     coords(1..3)       global coordinates of the node
!     vold(0..4,1..nk)   solution field in all nodes 
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
!     iponoel(i)         the network elements to which node i belongs
!                        are stored in inoel(1,iponoel(i)),
!                        inoel(1,inoel(2,iponoel(i)))...... until
!                        inoel(2,inoel(2,inoel(2......)=0
!     inoel(1..2,*)      field containing the network elements
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
!     temp(1..msecpt)    temperature in the node
!           
      implicit none
!
      integer msecpt,kstep,kinc,node,mi(*),iponoel(*),inoel(2,*),
     &  ipobody(2,*),ibody(3,*)
! 
      real*8 temp(msecpt),time(2),coords(3),vold(0:mi(2),*),xbody(7,*)
!
!
!
      temp(1)=293.d0
!
      return
      end

