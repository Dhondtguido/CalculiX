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
      subroutine interpolsubmodel(integerglob,doubleglob,value,
     &     coo,iselect,nselect,nodeface,tieset,istartset,iendset,
     &     ialset,ntie,entity,set,nset)
!
!     interpolates for a node with coordinates in "coo" the
!     "nselect" values with relative positions in "iselect" within the
!     global mesh stored in integerglob and doubleglob. The fields
!     integerglob and doubleglob are created and filled in
!     getglobalresults.c. 

!     The domain of the global model within which 
!     the interpolation takes place can be limited to an element
!     subset. To this end the submodel to which node "node" belongs is
!     determined. The submodels are stored in tieset(1..3,i), i=1..ntie.
!     tieset(1,i)(81:81)='S' if the tie is a submodel. In that case the
!     set corresponding to the submodel boundary is stored in 
!     tieset(2,i) and the set corresponding to the global element
!     model in tieset(3,i). Notice that the submodel boundary can be
!     a nodal set or an element face set (the actual node and the actual
!     element face are stored in nodeface, respecively).
!
      implicit none
!
      character*1 entity
      character*81 tieset(3,*),globset,submset,set(*)
!
      integer integerglob(*),nselect,iselect(nselect),nodeface,
     &  istartset(*),iendset(*),ialset(*),ntie,i,islavset,iset,
     &  nlength,id,jfaces,nelems,nktet,netet,ne,nkon,nfaces,nfield,
     &  imastset,nterms,konl(20),nelem,loopa,nset,ipos
!
      real*8 doubleglob(*),value(*),coo(3),ratio(20),dist,coords(3)
!
!     if no global file was read, set results to zero
!
      if(integerglob(1).eq.0) then
         do i=1,nselect
            value(i)=0.d0
         enddo
         return
      endif
!
!     determining the submodel to which the entity "nodeface" belongs
!
      islavset=0
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'S') cycle
!
!        check whether submodel is of the right kind (nodal or
!        element face)
!
         ipos=index(tieset(2,i),' ')-1
         if(tieset(2,i)(ipos:ipos).ne.entity) cycle
         read(tieset(2,i),'(a81)') submset
         call cident81(set,submset,nset,iset)
         nlength=iendset(iset)-istartset(iset)+1
         call nident(ialset(istartset(iset)),nodeface,nlength,id)
         if(id.le.0) cycle
         if(ialset(istartset(iset)+id-1).ne.nodeface) cycle
!
!        check whether slave set is of the right 
!
         islavset=iset
         exit
      enddo
!
!     check whether a submodel was found
!
      if(islavset.eq.0) then
         if(entity.eq.'N') then
            write(*,*) '*ERROR in interpolsubmodel: node',nodeface
            write(*,*) '       does not belong to any submodel'
            call exit(201)
         else
            nelems=int(nodeface/10)
            jfaces=nodeface-nelems*10
            write(*,*) '*ERROR in interpolsubmodel: face',jfaces
            write(*,*) '       of element',nelems
            write(*,*) '       does not belong to any submodel'
            call exit(201)
         endif
      endif
!
!     determining the global element set (if zero: all global elements
!     are taken)
!
      read(tieset(3,i),'(a81)') globset
      if(globset(1:1).eq.' ') then
        imastset=0
      else
        call cident81(set,globset,nset,imastset)
      endif
!
!     reading the number of nodes, tetrahedral interpolation elements,
!     global elements, connectivity size and number of faces
!
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=13
!
!     perform the interpolation
!
      coords(1)=coo(1)
      coords(2)=coo(2)
      coords(3)=coo(3)
      loopa=8
      call basis(doubleglob(1),doubleglob(netet+1),
     &     doubleglob(2*netet+1),
     &     doubleglob(3*netet+1),doubleglob(4*netet+1),
     &     doubleglob(5*netet+1),integerglob(6),integerglob(netet+6),
     &     integerglob(2*netet+6),doubleglob(6*netet+1),
     &     integerglob(3*netet+6),nktet,netet,
     &     doubleglob(4*nfaces+6*netet+1),nfield,
     &     doubleglob(13*nktet+4*nfaces+6*netet+1),
     &     integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &     integerglob(2*ne+7*netet+6),integerglob(nkon+2*ne+7*netet+6),
     &     coords(1),coords(2),coords(3),value,ratio,iselect,nselect,
     &     istartset,iendset,ialset,imastset,
     &     integerglob(nkon+2*ne+8*netet+6),nterms,konl,nelem,loopa,
     &     dist)
!
      return
      end
