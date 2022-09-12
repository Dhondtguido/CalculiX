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
      subroutine cflux(flux,msecpt,kstep,kinc,time,node,coords,vold,
     &  mi,ipkon,kon,lakon,iponoel,inoel,ielprop,prop,ielmat,
     &  shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
!
!     user subroutine cflux
!
!
!     INPUT:
!
!     msecpt             number of flux values (for volume elements:1)
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
!
!     OUTPUT:
!
!     flux(1..msecpt)   concentrated flux in the node
!           
      implicit none
!
      character*8 lakon(*)
!
      integer msecpt,kstep,kinc,node,mi(*),ipkon(*),kon(*),
     &  iponoel(*),inoel(2,*),ielprop(*),ielmat(mi(3),*),nshcon(*),
     &  nrhcon(*),ncocon(2,*),ntmat_
! 
      real*8 flux(msecpt),time(2),coords(3),vold(0:mi(2),*),
     &  prop(*),shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  cocon(0:6,ntmat_,*)
!
!
!
      flux(1)=10.d0
!
      return
      end

