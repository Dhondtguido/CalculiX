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
      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi,
     &     sti,xstateini,xstate,nstate_,dtime)
!
!     user subroutine dflux
!
!
!     INPUT:
!
!     sol                current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        1  = body flux
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     temp               currently not used
!     press              currently not used
!     loadtype           load type label
!     area               for surface flux: area covered by the
!                            integration point
!                        for body flux: volume covered by the
!                            integration point
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     co(3,1..nk)        coordinates of all nodes
!                        1: coordinate in global x-direction
!                        2: coordinate in global y-direction
!                        3: coordinate in global z-direction
!     lakonl             element label
!     konl(1..20)        nodes belonging to the element
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     nmpc               number of MPC's
!     ikmpc(1..nmpc)     ordered global degrees of freedom of the MPC's
!                        the global degree of freedom is
!                        8*(node-1)+direction of the dependent term of
!                        the MPC (direction = 0: temperature;
!                        1-3: displacements; 4: static pressure;
!                        5-7: rotations)
!     ilmpc(1..nmpc)     ilmpc(i) is the MPC number corresponding
!                        to the reference number in ikmpc(i)   
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     sti(i,j,k)         actual Cauchy stress component i at integration
!                        point j in element k. The components are
!                        in the order xx,yy,zz,xy,xz,yz
!     xstateini(i,j,k)   value of the state variable i at integration
!                        point j in element k at the beginning of the
!                        present increment
!     xstate(i,j,k)   value of the state variable i at integration
!                        point j in element k at the end of the
!                        present increment
!     nstate_            number of state variables
!     dtime              time length of the increment
!
!
!     OUTPUT:
!
!     flux(1)            magnitude of the flux
!     flux(2)            not used; please do NOT assign any value
!     iscale             determines whether the flux has to be
!                        scaled for increments smaller than the 
!                        step time in static calculations
!                        0: no scaling
!                        1: scaling (default)
!           
      implicit none
!
      character*8 lakonl
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,konl(20),ipompc(*),nstate_,i,
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),node,idof,id,iscale,mi(*)
!
      real*8 flux(2),time(2),coords(3),sol,temp,press,vold(0:mi(2),*),
     &  area,co(3,*),coefmpc(*),sti(6,mi(1),*),xstate(nstate_,mi(1),*),
     &  xstateini(nstate_,mi(1),*),dtime
!
!
!
!     Start of your own code.
!
!     example: heat generation due to plasticity for an isotropic
!              material
!
!     the plastic strain tensor is stored in state variables 2...7
!
!
      flux(1)=(sti(1,npt,noel)*
     &            (xstate(2,npt,noel)-xstateini(2,npt,noel))+
     &         sti(2,npt,noel)*
     &            (xstate(3,npt,noel)-xstateini(3,npt,noel))+
     &         sti(3,npt,noel)*
     &            (xstate(4,npt,noel)-xstateini(4,npt,noel))+
     &   2.d0*(sti(4,npt,noel)*
     &            (xstate(5,npt,noel)-xstateini(5,npt,noel))+
     &         sti(5,npt,noel)*
     &            (xstate(6,npt,noel)-xstateini(6,npt,noel))+
     &         sti(6,npt,noel)*
     &            (xstate(7,npt,noel)-xstateini(7,npt,noel))))/
     &         dtime
!
      return
      end

