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
      subroutine extrapolate_us45(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  vold,force,ielmat,thicke,ielprop,prop,i)
!
!     extrapolates field values at the integration points to the 
!     nodes for user element i of type us45
!
!     the present routine is called for each user element of type us45;
!     routine to calculate inum(i) for nodes i belonging to shell
!     elements of type us45. This guarantees that the coordinates and
!     fields such as U, S or SNEG are stored for these nodes in the
!     frd-file. 
!   
!     In this routine no extrapolation is done. For 'true' shell elements
!     the extrapolation is done in extrapolateshell_us45.f
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
!     vold(j,i)          value of variable j in node i at the end
!                        of the previous iteration   
!     force              logical variable; if true the values to
!                        be extrapolated are force values; important for
!                        interpolation from 3D expanded structures on the
!                        original 1D/2D structure: forces across the
!                        expansion have to be summed, not interpolated
!     ielmat(i)          material of element i
!     thicke(j,i)        thickness of layer j in node i
!     ielprop(i)         properties for element i are stored in
!                        prop(ielprop(i)+1),prop(ielprop(i)+2),....
!                        (number of properties depends on the type of
!                        element)
!     prop               property field
!     i                  number of the element for which the extrapolation
!                        is to be performed
!
!     OUTPUT:
!
!     inum(i)            < 0: node i is a network node
!                        > 0: node i is a structural node; inum(i)
!                             should be incremented by 1 if in the
!                             call of this routine an extrapolated value
!                             was stored for this node      
!                        =0: node i is not used
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,nfield,nk,i,ndim,
     &  iorienloc,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*)
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  vold(0:mi(2),*),thicke(mi(3),*)
!
      integer indexe,j,k,node
!     
      indexe=ipkon(i)
      do j=1,3
        node=kon(indexe+j)
        inum(node)=inum(node)+1
      enddo
!
      return
      end
