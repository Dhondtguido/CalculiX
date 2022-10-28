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
      subroutine createmdelem(imdnode,nmdnode,
     &              ikmpc,ilmpc,ipompc,nodempc,nmpc,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun,ithermal,imdelem,nmdelem,
     &              iponoel,inoel,prlab,prset,nprint,lakon,set,nset,
     &              ialset,ipkon,kon,istartset,iendset,nforc,
     &              ikforc,ilforc)
!
!     stores the elements 
!     1) for which results are requested in at least one node
!     2) for which there are *EL PRINT requests
!
!     stores the nodes, dofs, spcs and mpcs in the elements
!     for which there are *EL PRINT requests
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*81 prset(*),noset,set(*)
!
      integer iforc,node,imdnode(*),nmdnode,ikmpc(*),
     &  ilmpc(*),ipompc(*),nodempc(3,*),nmpc,imddof(*),nmddof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*),ithermal(*),imdelem(*),nmdelem,
     &  iponoel(*),inoel(2,*),index,id,nprint,i,j,k,l,indexe,
     &  nope,nset,nrset,ialset(*),ipkon(*),kon(*),istartset(*),
     &  iendset(*),idof,m,ikforc(*),ilforc(*),nforc
!
!     storing all elements to which nodes in imdnode belong
!     in imdelem
!
      do m=1,nmdnode
         node=imdnode(m)
!
         index=iponoel(node)
         do
            if(index.eq.0) exit
            i=inoel(1,index)
            call addimd(imdelem,nmdelem,i)
!            
            index=inoel(2,index)
         enddo
      enddo
!
!     storing the elements for which *EL PRINT was selected
!
      do m=1,nprint
         if((prlab(m)(1:4).eq.'S   ').or.
     &        (prlab(m)(1:4).eq.'E   ').or.
     &        (prlab(m)(1:4).eq.'PEEQ').or.
     &        (prlab(m)(1:4).eq.'ENER').or.
     &        (prlab(m)(1:4).eq.'SDV ').or.
     &        (prlab(m)(1:4).eq.'ELSE').or.
     &        (prlab(m)(1:4).eq.'ELKE').or.
     &        (prlab(m)(1:4).eq.'EVOL').or.
     &        (prlab(m)(1:4).eq.'HFL ')) then
            noset=prset(m)
c            nrset=0
c            do k=1,nset
c               if(set(k).eq.noset) then
c                  nrset=k
c                  exit
c               endif
c            enddo
            call cident81(set,noset,nset,id)
            nrset=0
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                nrset=id
              endif
            endif
!
!           adding the elements belonging to nrset
!
            do j=istartset(nrset),iendset(nrset)
               if(ialset(j).gt.0) then
                  i=ialset(j)
                  call addimd(imdelem,nmdelem,i)
!
!                 in order to calculate results at the integration
!                 point of an element the results must have been
!                 determined at the nodes of this element
!
                  indexe=ipkon(i)
c     Bernhardi start
                  if(lakon(i)(1:5).eq.'C3D8I') then
                     nope=11
                  elseif(lakon(i)(4:4).eq.'2') then
c     Bernhardi end
                     nope=20
                  elseif(lakon(i)(4:4).eq.'8') then
                     nope=8
                  elseif(lakon(i)(4:5).eq.'10') then
                     nope=10
                  elseif(lakon(i)(4:4).eq.'4') then
                     nope=4
                  elseif(lakon(i)(4:5).eq.'15') then
                     nope=15
                  elseif(lakon(i)(4:4).eq.'6') then
                     nope=6
                  elseif(lakon(i)(1:1).eq.'E') then
                     nope=ichar(lakon(i)(8:8))-47
                  else
                     cycle
                  endif
!     
                  do l=1,nope
                     node=kon(indexe+l)
                     call addimd(imdnode,nmdnode,node)
                     if(ithermal(1).ne.2) then
                        do k=1,3
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                        enddo
                     else
                        k=0
                        call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                       nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                       nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                       ikboun,nboun,ilboun)
                     endif
                  enddo
               else
                  i=ialset(j-2)
                  do
                     i=i-ialset(j)
                     if(i.ge.ialset(j-1)) exit
                     call addimd(imdelem,nmdelem,i)
!
!                 in order to calculate results at the integration
!                 point of an element the results must have been
!                 determined at the nodes of this element
!
                     indexe=ipkon(i)
c     Bernhardi start
                     if(lakon(i)(1:5).eq.'C3D8I') then
                        nope=11
                     elseif(lakon(i)(4:4).eq.'2') then
c     Bernhardi end
                        nope=20
                     elseif(lakon(i)(4:4).eq.'8') then
                        nope=8
                     elseif(lakon(i)(4:5).eq.'10') then
                        nope=10
                     elseif(lakon(i)(4:4).eq.'4') then
                        nope=4
                     elseif(lakon(i)(4:5).eq.'15') then
                        nope=15
                     elseif(lakon(i)(4:4).eq.'6') then
                        nope=6
                     elseif(lakon(i)(1:1).eq.'E') then
                        nope=ichar(lakon(i)(8:8))-47
                     else
                        cycle
                     endif
!     
                     do l=1,nope
                        node=kon(indexe+l)
                        call addimd(imdnode,nmdnode,node)
                        if(ithermal(1).ne.2) then
                           do k=1,3
                              call addimdnodedof(node,k,ikmpc,ilmpc,
     &                             ipompc,nodempc,nmpc,imdnode,nmdnode,
     &                             imddof,nmddof,nactdof,mi,imdmpc,
     &                             nmdmpc,imdboun,nmdboun,
     &                             ikboun,nboun,ilboun)
                           enddo
                        else
                           k=0
                           call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                          nodempc,nmpc,imdnode,nmdnode,imddof,
     &                          nmddof,nactdof,mi,imdmpc,nmdmpc,imdboun,
     &                          nmdboun,ikboun,
     &                          nboun,ilboun)
                        endif
                     enddo
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

