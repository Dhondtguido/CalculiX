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
!     construction of the B matrix      
!     
      subroutine networkmpc_lhs(i,ipompc,nodempc,coefmpc,
     &  labmpc,v,nactdog,ac,j,mi,nteq,ipkon,kon,lakon,iponoel,
     &  inoel,ielprop,prop,ielmat,
     &  shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
!
!     user defined network mpc: calculation of the left hand
!     side
!
!     INPUT:
!
!     i                  MPC number
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     labmpc(*)          label of the MPC. For user-defined
!                        network MPC's it starts with NETWORK;
!                        the remaining 13 characters can be used
!                        to distinguish between different kinds of 
!                        network user MPC's
!     v(0..mi(2),1..nk)  actual solution field in all nodes
!                        0: total temperature
!                        1: mass flow
!                        2: total pressure
!     nactdog(j,i)       determines the network equation corresponding
!                        to degree of freedom j in node i;
!                        if zero the degree of freedom is not active
!     j                  network equation corresponding to the
!                        present MPC (i.e. MPC i)
!     mi(*)              field with global information; mi(2) is the
!                        highest variable number
!     nteq               number of network equations
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
!     ac(*)              left hand side of the system of network
!                        equations; this routines should return the
!                        derivative of the network MPC at stake w.r.t.
!                        all active degrees of freedom occurring in the
!                        MPC and store them in row j of matrix ac.
!
      implicit none
!     
      character*8 lakon(*)
      character*20 labmpc(*)
!     
      integer mi(*),i,ipompc(*),nodempc(3,*),j,index,node,idir,
     &  nactdog(0:3,*),nteq,ipkon(*),kon(*),
     &  iponoel(*),inoel(2,*),ielprop(*),ielmat(mi(3),*),nshcon(*),
     &  nrhcon(*),ncocon(2,*),ntmat_
!     
      real*8 coefmpc(*),v(0:mi(2),*),ac(nteq,*),
     &  prop(*),shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  cocon(0:6,ntmat_,*)
!
!
!     
      if(labmpc(i)(8:16).eq.'QUADRATIC') then
!
!        example equation of the form
!        f:=a*v(idir1,node1)+b*v(idir2,node2)**2=0
!
!        a,idir1,node1,b,idir2,node2 are given in the input deck
!        using the *NETWORK MPC keyword
!        to be calculated: a*df/d(v(idir1,node1))
!                          b*df/d(v(idir2,node2))
!
         index=ipompc(i)
         node=nodempc(1,index)
         idir=nodempc(2,index)
         ac(j,nactdog(idir,node))=coefmpc(index)
!
         index=nodempc(3,index)
         node=nodempc(1,index)
         idir=nodempc(2,index)
!
!        if nactdog(idir,node) is zero the degree of freedom is
!        not active
!
         if(nactdog(idir,node).ne.0) then
            ac(j,nactdog(idir,node))=2.d0*coefmpc(index)*v(idir,node)
         endif
      else
         write(*,*) '*ERROR in networkmpc_lhs:'
         write(*,*) '       unknown MPC: ',labmpc(i)
         call exit(201)
      endif
!         
      return
      end
