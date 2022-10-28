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
      subroutine umpc_dist(x,u,f,a,jdof,n,force,iit,idiscon)
!
!     updates the coefficients in a dist mpc (name DIST)
!
!     a dist mpc specifies that the distance between two nodes
!     a and b must not exceed value d
!
!     input nodes: a,a,a,b,b,b,c
!
!     node c is a fictitious node. The value d must be assigned
!     to the first coordinate of node c by means of a *NODE card;
!     the other coordinates of the node can be arbitrary.
! 
!     A value of zero must be assigned to the first DOF of node c by using
!     a *BOUNDARY card. The second DOF of node c is not constrained and is 
!     used when the distance between nodes a and b is less than d: in
!     that case there is no constraint at all.
!
!     INPUT:
!
!     x(3,n)             Carthesian coordinates of the nodes in the
!                        user mpc.
!     u(3,n)             Actual displacements of the nodes in the
!                        user mpc.     
!     jdof               Actual degrees of freedom of the mpc terms
!     n                  number of terms in the user mpc
!     force              Actual value of the mpc force
!     iit                iteration number
!
!     OUTPUT:
!
!     f                  Actual value of the mpc. If the mpc is
!                        exactly satisfied, this value is zero
!     a(n)               coefficients of the linearized mpc
!     jdof               Corrected degrees of freedom of the mpc terms
!     idiscon            0: no discontinuity
!                        1: discontinuity
!                        If a discontinuity arises the previous
!                        results are not extrapolated at the start of
!                        a new increment
!
      implicit none
!
      integer jdof(*),n,iit,ifix,idiscon
!
      real*8 x(3,*),u(3,*),f,a(*),dist(3),force
!
c      write(*,*) (jdof(i),i=1,7)
      if(jdof(7).eq.1) then
         ifix=1
      else
         ifix=0
         jdof(7)=2
      endif
!
      dist(1)=x(1,1)+u(1,1)-x(1,4)-u(1,4)
      dist(2)=x(2,1)+u(2,1)-x(2,4)-u(2,4)
      dist(3)=x(3,1)+u(3,1)-x(3,4)-u(3,4)
!
      f=dist(1)**2+dist(2)**2+dist(3)**2-x(1,7)**2
!
c      write(*,*) 'mpcforc=, f= ',force,f
!
      a(7)=-1.d0
!
!     only one change per increment is allowed
!        (change= from free to linked or vice versa)
!     ifix=0: free
!     ifix=1: linked
!
      if(ifix.eq.0) then
!
!        previous state: free
!
         if(f.lt.0) then
!
!           new state: free
!
            f=0.d0
         elseif(iit.le.1) then
!
!           new state: linked
!
            write(*,*) 'switch to linked'
            write(*,*)
            jdof(7)=1
            idiscon=1
         else
!
!           new state: free
!
            f=0.d0
         endif
      else
!
!        previous state: linked
!
         if(force.le.0.d0) then
!
!           new state: linked
!
         elseif(iit.le.1) then
!
!           new state: free
!
            write(*,*) 'switch to free'
            write(*,*)
            jdof(7)=2
            f=0.d0
            idiscon=1
         else
!
!           new state: linked
!
         endif
      endif
!
      if(dabs(dist(jdof(1))).gt.1.d-10) then
         a(1)=2.d0*dist(jdof(1))
         if(jdof(1).eq.1) then
            jdof(2)=2
            jdof(3)=3
         elseif(jdof(1).eq.2) then
            jdof(2)=3
            jdof(3)=1
         else
            jdof(2)=1
            jdof(3)=2
         endif
         a(2)=2.d0*dist(jdof(2))
         a(3)=2.d0*dist(jdof(3))
      else
         if(jdof(1).eq.3) then
            jdof(1)=1
         else
            jdof(1)=jdof(1)+1
         endif
         if(dabs(dist(jdof(1))).gt.1.d-10) then
            a(1)=2.d0*dist(jdof(1))
            if(jdof(1).eq.1) then
               jdof(2)=2
               jdof(3)=3
            elseif(jdof(1).eq.2) then
               jdof(2)=3
               jdof(3)=1
            else
               jdof(2)=1
               jdof(3)=2
            endif
            a(2)=2.d0*dist(jdof(2))
            a(3)=2.d0*dist(jdof(3))
         else
            if(jdof(1).eq.3) then
               jdof(1)=1
            else
               jdof(1)=jdof(1)+1
            endif
            if(dabs(dist(jdof(1))).gt.1.d-10) then
               a(1)=2.d0*dist(jdof(1))
               if(jdof(1).eq.1) then
                  jdof(2)=2
                  jdof(3)=3
               elseif(jdof(1).eq.2) then
                  jdof(2)=3
                  jdof(3)=1
               else
                  jdof(2)=1
                  jdof(3)=2
               endif
               a(2)=2.d0*dist(jdof(2))
               a(3)=2.d0*dist(jdof(3))
            endif
         endif
      endif
!
      a(4)=-2.d0*dist(1)
      a(5)=-2.d0*dist(2)
      a(6)=-2.d0*dist(3)
      jdof(4)=1
      jdof(5)=2
      jdof(6)=3
!
c      write(*,*) 'jdof,a'
c      do i=1,7
c         write(*,*) jdof(i),a(i)
c      enddo
c      write(*,*) 'f ',f
!
      return
      end













