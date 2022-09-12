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
      subroutine generateeminterfaces(istartset,iendset,
     &  ialset,iactive,ipkon,lakon,kon,ikmpc,nmpc,nafaces)
!
!     determines the interfaces in between the phi-domain and any
!     other domain, i.e.
!
!     faces belonging to the phi domain and adjacent to the A- or
!           A-V-domain 
!     faces belonging to the A-domain adjacent to the phi-domain
!     faces belonging to the A-V-domain adjacent to the phi-domain
!
      implicit none
!
      character*8 lakon(*)
!
      integer istartset(*),iendset(*),ialset(*),iactive(3),
     &  ipkon(*),kon(*),iset,ntotal,nope,nopes,indexe,nodef(9),
     &  ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),m,
     &  node,idof,jdof(3),idummy,nlength,kflag,ikmpc(*),nmpc,id,
     &  iface,nelem,jface,i,j,k,nafaces
!
      data kflag /1/
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
!     The sets iactive(1),iactive(2) and iactive(3) contain all
!     external surfaces of the phi-domain, the A-V-domain and the 
!     A-domain, respectively. In what follows these surfaces are
!     shrunk to those faces in common between the phi-domain and the
!     others, i.e. iactive(1) contains the faces of the phi-domain
!     adjacent to the A-V-domain and A-domain, iactive(2) contains the
!     faces of the A-V-domain adjacent to the phi-domain and iactive(3)
!     contains the faces of the A-domain adjacent to the phi-domain
!
!     1) determining the faces of the phi-domain in all nodes of which a
!        MPC for the x-component of A exists.
!     2) determining the faces of the A-V-domain in all nodes of which
!        a MPC for phi exists (dof 5)
!     3) determining the faces of the A-domain in all nodes of which 
!        a MPC for phi exists (dof 5)
!
      jdof(1)=1
      jdof(2)=5
      jdof(3)=5
!
      do m=1,3
         if(iactive(m).gt.0) then
            iset=iactive(m)
            ntotal=0
            loop1: do j=istartset(iset),iendset(iset)
               iface=ialset(j)
               nelem=int(iface/10.d0)
               jface=iface-10*nelem
!     
!              nodes belonging to the element (nope) and to the
!              face (nopes) 
!     
               if(lakon(nelem)(4:5).eq.'8R') then
                  nopes=4
                  nope=8
               elseif(lakon(nelem)(4:4).eq.'8') then
                  nopes=4
                  nope=8
               elseif(lakon(nelem)(4:6).eq.'20R') then
                  nopes=8
                  nope=20
               elseif(lakon(nelem)(4:5).eq.'20') then
                  nopes=8
                  nope=20
               elseif(lakon(nelem)(4:5).eq.'10') then
                  nopes=6
                  nope=10
               elseif(lakon(nelem)(4:4).eq.'4') then
                  nopes=3
                  nope=4
!     
!     treatment of wedge faces
!     
               elseif(lakon(nelem)(4:4).eq.'6') then
                  nope=6
                  if(jface.le.2) then
                     nopes=3
                  else
                     nopes=4
                  endif
               elseif(lakon(nelem)(4:5).eq.'15') then
                  nope=15
                  if(jface.le.2) then
                     nopes=6
                  else
                     nopes=8
                  endif
               endif
!     
!     nodes belonging to the face
!     
               indexe=ipkon(nelem)
               if((nope.eq.20).or.(nope.eq.8)) then
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifaceq(k,jface))
                  enddo
               elseif((nope.eq.10).or.(nope.eq.4)) then
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifacet(k,jface))
                  enddo
               elseif(nope.eq.15) then
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifacew2(k,jface))
                  enddo
               else
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifacew1(k,jface))
                  enddo
               endif
!     
               do i=1,nopes
                  node=nodef(i)
                  idof=8*(node-1)+jdof(m)
                  call nident(ikmpc,idof,nmpc,id)
                  if(id.gt.0) then
                     if(ikmpc(id).eq.idof) cycle
                  endif
                  cycle loop1
               enddo
!
               ialset(istartset(iset)+ntotal)=iface
               ntotal=ntotal+1
               cycle
            enddo loop1
            iendset(iset)=istartset(iset)+ntotal-1
         endif
      enddo
!
!     sorting the sets
!
      do m=1,3
         if(iactive(m).eq.0) cycle
         nlength=iendset(iactive(m))-istartset(iactive(m))+1
         call isortii(ialset(istartset(iactive(m))),idummy,
     &                 nlength,kflag)
      enddo
!
!     nafaces is the total number of faces belonging to the A-V-domain
!     or A-domain and adjacent to the phi-domain. The nodes on these
!     faces are subject to A.n=0
!
      if((iactive(2).eq.0).and.(iactive(3).eq.0)) then
         nafaces=0
      elseif(iactive(2).eq.0) then
         nafaces=iendset(iactive(3))-istartset(iactive(3))+1
      elseif(iactive(3).eq.0) then
         nafaces=iendset(iactive(2))-istartset(iactive(2))+1
      else
         nafaces=iendset(iactive(2))-istartset(iactive(2))+1
     &        +iendset(iactive(3))-istartset(iactive(3))+1
      endif
!     
      return
      end
