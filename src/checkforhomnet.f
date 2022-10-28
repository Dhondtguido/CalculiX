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
      subroutine checkforhomnet(ieg,nflow,lakon,ipkon,kon,itg,ntg,
     &     iponoel,inoel)
!     
!     check whether (simply or not) connected networks are
!     inhomogeneous
!     
      implicit none
!     
      logical gas,liquid,channel,negative,untreated
!     
      character*8 lakon(*)
!     
      integer i,ieg(*),nflow,nelem,indexe,ipkon(*),kon(*),itg(*),
     &     ntg,iponoel(*),inoel(2,*),index,ier,j,newnode,node,id
!     
      ier=0
!     
      do
        gas=.false.
        liquid=.false.
        channel=.false.
        untreated=.false.
!     
        loop1: do i=1,nflow
          nelem=ieg(i)
          if(ipkon(nelem).lt.-1) cycle
!     
!     look for a loose end: start the analysis of a new branch
!     
          if((lakon(nelem)(2:3).eq.'IO').or.
     &       (lakon(nelem)(2:5).eq.'INLT').or.
     &       (lakon(nelem)(2:6).eq.'OUTLT'))then
!
!           chamber-like entrance/exit of gas network
!
            gas=.true.
          elseif(lakon(nelem)(2:6).eq.'LIPIO') then
!
!           entrance/exit of liquid network
!
            liquid=.true.
          elseif(lakon(nelem)(2:7).eq.'LICHIO') then
!
!           entrance/exit of channel network
!
            channel=.true.
          elseif(((kon(ipkon(nelem)+1).eq.0).or.
     &            (kon(ipkon(nelem)+3).eq.0)).and.
     &           (lakon(nelem)(2:3).ne.'LP').and.
     &           (lakon(nelem)(2:3).ne.'LI')) then
!
!           pipe-connection-like entrance/exit of gas network
!
            gas=.true.
          else
            untreated=.true.
            cycle
          endif
!     
!     identify the nonzero node in the input-output element
!     
          indexe=ipkon(nelem)
          if(kon(indexe+1).ne.0) then
            node=kon(indexe+1)
          else
            node=kon(indexe+3)
          endif
!     
!     remove the element
!     
          ipkon(nelem)=-2-ipkon(nelem)
!     
!     mark the node by a negative sign
!     
          call nident(itg,node,ntg,id)
          itg(id)=-itg(id)
!     
          do
            negative=.false.
            loop2: do j=1,ntg
!     
!     look for an active node in the branch at
!     stake (marked by a negative sign in field itg)
!     
              if(itg(j).gt.0) cycle
              negative=.true.
              itg(j)=-itg(j)
              node=itg(j)
!     
!     look for a non-treated element connected to the node
!     
              index=iponoel(node)
              do
                nelem=inoel(1,index)
                if(ipkon(nelem).lt.-1) then
                  index=inoel(2,index)
                  if(index.eq.0) then
                    exit loop2
                  endif
                  cycle
                endif
!     
!     check whether the element label fits the type of
!     branch
!     
                if(gas) then
                  if((lakon(nelem)(2:3).eq.'LP').or.
     &                 (lakon(nelem)(2:3).eq.'LI')) then
                    write(*,*) '*ERROR in checkforhomnet:'
                    write(*,*) '       a branch of the'
                    write(*,*) '       network seems to be a'
                    write(*,*) '       gas branch, however,'
                    write(*,*) '       element',nelem,' has label',
     &                   lakon(nelem)
                    write(*,*)
     &                   '       and is a liquid or channel element'
                    ier=1
                  endif
                elseif(liquid) then
                  if((lakon(nelem)(2:3).ne.'LP').and.
     &                 (lakon(nelem)(2:5).ne.'LIPI').and.
     &                 (lakon(nelem)(2:5).ne.'LIPU')) then
                    write(*,*) '*ERROR in checkforhomnet:'
                    write(*,*) '       a branch of the'
                    write(*,*) '       network seems to be a'
                    write(*,*) '       liquid branch, however,'
                    write(*,*) '       element',nelem,' has label',
     &                   lakon(nelem)
                    write(*,*)
     &                   '       and is a gas or channel element'
                    ier=1
                  endif
                else
                  if(lakon(nelem)(2:5).ne.'LICH') then
                    write(*,*) '*ERROR in checkforhomnet:'
                    write(*,*) '       a branch of the'
                    write(*,*) '       network seems to be a'
                    write(*,*) '       channel branch, however,'
                    write(*,*) '       element',nelem,' has label',
     &                   lakon(nelem)
                    write(*,*)
     &                   '       and is a gas or liquid element'
                    ier=1
                  endif
                endif
!     
!     look for the other end node of the element
!     
                indexe=ipkon(nelem)
                ipkon(nelem)=-2-ipkon(nelem)
                if(kon(indexe+1).ne.node) then
                  newnode=kon(indexe+1)
                else
                  newnode=kon(indexe+3)
                endif
!     
!     check whether node exists; if so, mark the node
!     as active by giving it a negative sign
!     
                if(newnode.ne.0) then
                  call nident(itg,newnode,ntg,id)
                  itg(id)=-itg(id)
                endif
                index=inoel(2,index)
                if(index.eq.0) exit
              enddo
            enddo loop2
!     
!     if no more negative nodes: branch is finished
!     
            if(.not.negative) exit loop1
          enddo
        enddo loop1
!     
!     if no more untreated IO-elements: finished
!     
        if((.not.gas).and.(.not.liquid).and.(.not.channel)) then
          if(untreated) then
            write(*,*) '*ERROR in checkforhomnet: there seem to be'
            write(*,*) '       network elements not connected'
            write(*,*) '       to an input/output element'
            ier=1
          endif
          exit
        endif
      enddo
!     
!     reactivating the network elements
!     
      do i=1,nflow
        nelem=ieg(i)
        if(ipkon(nelem).lt.-1) ipkon(nelem)=-2-ipkon(nelem)
      enddo
!     
      if(ier.eq.1) call exit(201)
!     
      return
      end
