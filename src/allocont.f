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
      subroutine allocont(ncont,ntie,tieset,nset,set,istartset,
     &  iendset,ialset,lakon,ncone,tietol,ismallsliding,kind1,kind2,
     &  mortar,istep)
!
!     counting the number of triangles needed for the 
!     triangulation of the contact master surfaces
!
!     ismallsliding = 0: large sliding
!                   = 1: small sliding
!
      implicit none
!
      logical nodeslavsurf
!
      character*1 kind1,kind2
      character*8 lakon(*)
      character*81 tieset(3,*),mastset,set(*),slavset
!
      integer ncont,ntie,i,j,k,nset,istartset(*),iendset(*),ialset(*),
     &  imast,nelem,jface,ncone,islav,ismallsliding,ipos,mortar,istep,
     &  kflag,idummy,jact,id
!
      real*8 tietol(4,*)
!
!     number of master triangles
!
      ncont=0
!
!     number of slave entities (nodes for nodal surfaces and
!     faces for facial surfaces)
!
      ncone=0
!
      do i=1,ntie
!
!        check for contact conditions
!
         if((tieset(1,i)(81:81).eq.kind1).or.
     &      (tieset(1,i)(81:81).eq.kind2)) then
            if(tietol(1,i).lt.0.d0) then
               ismallsliding=1
            else
               ismallsliding=0
            endif
            mastset=tieset(3,i)
!
!           determining the master surface
!
c            do j=1,nset
c               if(set(j).eq.mastset) exit
c            enddo
            call cident81(set,mastset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(mastset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
               ipos=index(mastset,' ')
               write(*,*) '*ERROR in allocont: master surface ',
     &               mastset(1:ipos-2)
               write(*,*) '       does not exist or does not contain'
               write(*,*) '       element faces'
               call exit(201)
            endif
            imast=j
!
!           deleting identical entries in the master facial surface
!           definition (leads otherwise to problems in the triangulation
!           and the creation of imastop)
!
            if(istep.eq.1) then
               kflag=1
               call isortii(ialset(istartset(imast)),idummy,
     &            iendset(imast)-istartset(imast)+1,kflag)
               jact=istartset(imast)
               do j=istartset(imast)+1,iendset(imast)
                  if(ialset(j).eq.ialset(j-1)) cycle
                  jact=jact+1
                  ialset(jact)=ialset(j)
               enddo
               iendset(imast)=jact
            endif
!
            do j=istartset(imast),iendset(imast)
!     
               nelem=int(ialset(j)/10.d0)
               jface=ialset(j)-10*nelem
!     
               if(lakon(nelem)(4:5).eq.'20') then
                  ncont=ncont+6
               elseif(lakon(nelem)(4:4).eq.'2') then
                  ncont=ncont+8
               elseif(lakon(nelem)(4:4).eq.'8') then
                  ncont=ncont+2
               elseif(lakon(nelem)(4:5).eq.'10') then
                  ncont=ncont+4
               elseif(lakon(nelem)(4:4).eq.'4') then
                  ncont=ncont+1
               elseif(lakon(nelem)(4:5).eq.'15') then
                  if(jface.le.2) then
                     ncont=ncont+4
                  else
                     ncont=ncont+6
                  endif
               elseif(lakon(nelem)(4:4).eq.'6') then
                  if(jface.le.2) then
                     ncont=ncont+1
                  else
                     ncont=ncont+2
                  endif
               endif
            enddo
!
!           counting the slave nodes
!
            slavset=tieset(2,i)
            ipos=index(slavset,' ')-1
            if(slavset(ipos:ipos).eq.'T') then
!
!              face-to-face penalty contact (facial slave surface)
!
c               mortar=1
               nodeslavsurf=.false.
            elseif(slavset(ipos:ipos).eq.'M') then
!
!              quad-quad Mortar contact (facial slave surface)
!
c               mortar=2
               nodeslavsurf=.false.
            elseif(slavset(ipos:ipos).eq.'P') then
!
!              quad-lin Petrov Galerkin Mortar contact (facial slave surface)
!
c               mortar=4
               nodeslavsurf=.false.
            elseif(slavset(ipos:ipos).eq.'G') then
!
!              quad-quad Petrov Galerkin Mortar contact (facial slave surface)
!
c               mortar=5
               nodeslavsurf=.false.
            elseif(slavset(ipos:ipos).eq.'O') then
!
!              quad-lin Mortar contact (facial slave surface)
!
c               mortar=3
               nodeslavsurf=.false.
            else
!
!              node-to-face contact
!              default is a nodal slave surface
!
               nodeslavsurf=.true.
            endif
!
!           determining the slave surface
!
c            do j=1,nset
c               if(set(j).eq.slavset) exit
c            enddo
            call cident81(set,slavset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(slavset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
               if(mortar.eq.1) then
                  write(*,*) 
     &              '*ERROR in allocont: element slave surface ',
     &              slavset(1:ipos-1)
                  write(*,*) '       does not exist'
                  call exit(201)
               endif
               do j=1,nset
                  if((set(j)(1:ipos-1).eq.slavset(1:ipos-1)).and.
     &                 (set(j)(ipos:ipos).eq.'T')) then
                     nodeslavsurf=.false.
                     exit
                  endif
               enddo
               if(j.gt.nset) then
                  write(*,*) '*ERROR in allocont: slave surface ',
     &                 slavset(1:ipos-1)
                  write(*,*) '       does not exist'
                  call exit(201)
               endif
            endif
!
            islav=j
!
!           deleting identical entries in the slave facial surface
!           definition (leads otherwise to problems in the calculation
!           of the are corresponding to the slave nodes)
!
            if((istep.eq.1).and.((mortar.eq.1).or.(.not.nodeslavsurf)))
     &         then
               kflag=1
               call isortii(ialset(istartset(islav)),idummy,
     &            iendset(islav)-istartset(islav)+1,kflag)
               jact=istartset(islav)
               do j=istartset(islav)+1,iendset(islav)
                  if(ialset(j).eq.ialset(j-1)) cycle
                  jact=jact+1
                  ialset(jact)=ialset(j)
               enddo
               iendset(islav)=jact
            endif
!
!           counting the entities (nodes or faces) in the slave
!           surface
!
            do j=istartset(islav),iendset(islav)
               if(ialset(j).gt.0) then
                  ncone=ncone+1
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     ncone=ncone+1
                  enddo
               endif
            enddo
!
         endif
      enddo
!
      return
      end

