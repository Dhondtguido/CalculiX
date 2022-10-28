!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2017 Guido Dhondt
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
      subroutine dload(f,kstep,kinc,time,noel,npt,layer,kspt,
     &     coords,jltyp,loadtype,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,veold,
     &     rho,amat,mi)
!
!     user subroutine dload
!
!
!     INPUT:
!
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     layer              currently not used
!     kspt               currently not used
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        21 = face 1 
!                        22 = face 2 
!                        23 = face 3 
!                        24 = face 4 
!                        25 = face 5 
!                        26 = face 6
!     loadtype           load type label
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     veold(0..3,1..nk)  derivative of the solution field w.r.t.
!                        time in all nodes
!                        0: temperature rate
!                        1: velocity in global x-direction
!                        2: velocity in global y-direction
!                        3: velocity in global z-direction
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
!     rho                local density 
!     amat               material name
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!
!     f                  magnitude of the distributed load
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
      character*80 amat
!
      integer kstep,kinc,noel,npt,jltyp,layer,kspt,konl(20),iscale,
     &  mi(*)
c
      integer ifaceq(8,6),ifacet(6,4),ifacew(8,5),ig,nelem,nopes,
     &  iflag,i,j,nope,ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),
     &  node,idof,id
!
      real*8 f,time(2),coords(3),vold(0:mi(2),*),co(3,*),rho
      real*8  coefmpc(*),veold(0:mi(2),*)
!
      intent(in) kstep,kinc,time,noel,npt,layer,kspt,
     &     coords,jltyp,loadtype,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,
     &     rho,amat,mi
!
      intent(inout) f,iscale
!
!     the code starting here up to the end of the file serves as
!     an example for combined mechanical-lubrication problems. 
!     Please replace it by your own code for your concrete application.
!
c
      real*8 c,p1,p2,p3,x1,x2,x3,x,xs1,xs2,xs3,p
c
c     Geschwindigkeit der Explosionsfront
      c=1.0e6
c
c     Druck vor der Exfront
      p1=0.0
c
c     Spitzendruck
      p2=100.0
c
c     Druck am Ende des Peaks, Enddruck
      P3=3.0
c
c     Mitbewegte Koordinaten der Explosionsfront
      x3=-100.0
      x2=-50.0
      x1=0.0
c
c     Koordinate des Integrationspunktes
      x=coords(1)
c
c     Zeitabhängige Koordinaten
      xs1=x1+c*time(2)
      xs2=x2+c*time(2)
      xs3=x3+c*time(2)
c
c     Druckverlauf
      if(x.le.xs3) then
      p=p3
      elseif(x.le.xs2) then
      p=p3+(p2-p3)/(xs2-xs3)*(x-xs3)
      elseif(x.le.xs1) then
      p=p2+(p1-p2)/(xs1-xs2)*(x-xs2)
      else
      p=p1
      endif
c
c     return value
      f=p

      iscale=0
!
      return
      end

