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
      subroutine extrapolateshell_us3(yi,yn,ipkon,inum,kon,lakon,
     &  nfield,nk,ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  ielmat,thicke,ielprop,prop,i,iflag)
!
!     extrapolates field values at the integration points to the 
!     nodes for user element i of type u1
!
!     the present routine is called for each user element of type us45;
!     the field yn(j,i) contains at entry the sum of all extrapolations
!     done for component j of the variable to node i (since a node can
!     belong to n elements (not necessarily all of type user element
!     u1) at the end n extrapolated values are available
!     for this node). They are accumulated and at the end divided by n. The
!     value of n is stored in inum(i): it starts at zero and is incremented
!     by one for each extrapolation done to node i.
!       
!     the present routine cannot be used for network elements, only for
!     structural elements
!
!     INPUT:
!
!     yi(ndim,mi(1),*)   value of the variables at the integration
!                        points
!     yn(nfield,*)       sum of the extrapolated variables at the nodes
!                        from all elements treated previously
!     ipkon(i)           points to the location in field kon preceding
!                        the topology of element i
!     inum(i)            < 0: node i is a network node
!                        > 0: node i is a structural node; its value is
!                             number of extrapolations performed to this
!                             node so far
!                        =0: node i is not used
!     kon(*)             contains the topology of all elements. The
!                        topology of element i starts at kon(ipkon(i)+1)
!                        and continues until all nodes are covered. The
!                        number of nodes depends on the element label
!     lakon(i)           contains the label of element i
!     nfield             number of variables to be extrapolated
!     nk                 maximum node number in the mesh
!     ne                 maximum element number in the mesh
!     ndim               number of variables in the integration point
!                        field to be extrapolated
!     orab(7,*)          description of all local coordinate systems.
!                        (cf. List of variables and their meaning in the
!                        User's manual)
!     ielorien(i)        orientation in element i
!     co(1..3,i)         global coordinates of node i
!     iorienloc          0: extrapolated variables requested in global 
!                           coordinates
!                        1: extrapolated variables requested in local
!                           coordinates
!     cflag (char*1)     I: interpolate 3D results onto 1D/2D
!                        E: store extrapolated 1D/2D results
!                        M: store 1D section forces
!                        blank: any other case
!     ielmat(i)          material of element i
!     thicke(j,i)        thickness of layer j in node i
!     ielprop(i)         properties for element i are stored in
!                        prop(ielprop(i)+1),prop(ielprop(i)+2),....
!                        (number of properties depends on the type of
!                        element)
!     prop               property field
!     i                  number of the element for which the extrapolation
!                        is to be performed
!     iflag              -1: values for lower surface of shell requested
!                         0: values for midsurface of shell requested
!                        +1: values for upper surface of shell requested
!      
!
!     OUTPUT:
!
!     yn(nfield,*)       value of the variables at the nodes
!     inum(i)            < 0: node i is a network node
!                        > 0: node i is a structural node; inum(i)
!                             should be incremented by 1 if in the
!                             call of this routine an extrapolated value
!                             was stored for this node      
!                        =0: node i is not used
!
      implicit none
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,nfield,nk,i,ndim,
     &  iorienloc,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*),iflag,
     &  jj
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  thicke(mi(3),*)
!
!     START OF THIS SUBROUTINE
!
      integer indexe,j,k,node


      !
      if(iorienloc.ne.0) then
         write(*,*) '*ERROR in extrapolate_us3'
         write(*,*) '       no local orientation for variables'
         write(*,*) '       belonging to this type of element'
         write(*,*) '       allowed'
         call exit(201)
      endif
      !
      if(nfield.eq.6) then
         indexe=ipkon(i)
         do j=1,3 ! number of nodes
            node = kon(indexe+j)
            do k=1,nfield
            !
            ! BOT
            !
              if(iflag.eq.-1) then  ! gauss point 1
                yn(k,node) = yn(k,node) + yi(k,1,i)
              endif
            !
            ! MID
            !
              if(iflag.eq.0) then  ! gauss point 2
                yn(k,node) = yn(k,node) + yi(k,2,i)
              endif
            !
            ! TOP
            !
              if(iflag.eq.1) then  ! gauss point 3
                yn(k,node) = yn(k,node) + yi(k,3,i)
              endif              
            enddo
            inum(node)=inum(node)+1  ! add divider
         enddo
      else
         write(*,*) '*ERROR in extrapolate_us3'
         write(*,*) '       extropolation for element of type US3'
         write(*,*) '       is only coded for fields with 6'
         write(*,*) '       entries'
         call exit(201)
      endif
!
!     END OF THIS SUBROUTINE
!
      return
      end
