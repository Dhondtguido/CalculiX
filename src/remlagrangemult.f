!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
!     all slave nodes
!     - belonging to more than 1 slave surface
!     - belonging to SPCs or MPCs (as dependent are independent node)
!     are set to no-LM nodes
!
!     Author: Saskia Sitzmann
!     
!     islavact(i): -3: no slave node
!                  -2: no LM-node
!                  -1: no gap-node
!                   0: inactive node
!                   1: stick-node
!                   2: slip/active node
!     
      subroutine remlagrangemult(ntie,tieset,islavnode,imastnode,
     &     nslavnode,nmastnode,islavact,nodempc,nmpc,ipompc)
!     
!     check whether SPC's and MPC's in salve nodes are compatible
!     with mortar contact    
!     
!     author: Sitzmann,Saskia
!     
      implicit none
!
      character*81 tieset(3,*)
!     
      integer ntie,i,j,l,id,node,islavnode(*),imastnode(*),
     &     nslavnode(ntie+1),nmastnode(ntie+1),islavact(*),
     &     nodempc(3,*),index,nmpc,ipompc(*),ist,node2
!     
!     remove Lagrange Multiplier contribution for nodes which are
!     in more than one contact tie
!     
      if(ntie.gt.1) then
        do i=1,ntie
          if(tieset(1,i)(81:81).ne.'C') cycle
          do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
            if(islavact(l).gt.-1) then
              do j=1,ntie
                if(j.ne.i) then
                  if(tieset(1,j)(81:81).ne.'C') cycle
                  call nident(islavnode(nslavnode(j)+1),node,
     &                 nslavnode(j+1)-nslavnode(j),id)
                  if(id>0) then
                    if(islavnode(nslavnode(j)+id).eq.node) then
                      islavact(l)=-2
                      write(*,*)'checkspcmpc: node',node,
     &                     'tie1s',i,'tie2s',j
                      write(*,*)'in more than one contact',
     &                     'tie and set NoLM!'
                    endif
                  endif                   
                  call nident(imastnode(nmastnode(j)+1),node,
     &                 nmastnode(j+1)-nmastnode(j),id)
                  if(id>0) then
                    if(imastnode(nmastnode(j)+id).eq.node) then
                      islavact(l)=-2
                      write(*,*)'checkspcmpc: node',node,
     &                     'tie1s',i,'tie2m',j
                      write(*,*)'in more than one',
     &                     ' contact tie and set NoLM!'
                    endif
                  endif                   
                endif
              enddo
            endif
          enddo
        enddo
!     
      endif
!     
!     remove Lagrange Multiplier contribution from all slave nodes
!     involved in MPCs;
!     needed for quadratic elements
!     attention: 2D calculation are not possible right now
!     
      do i=1,nmpc
        ist=ipompc(i)
        node=nodempc(1,ist)
        do j=1,ntie
          call nident(islavnode(nslavnode(j)+1),node,
     &         nslavnode(j+1)-nslavnode(j),id)
          if(id.gt.0) then
            if(islavnode(nslavnode(j)+id).eq.node) then
              islavact(nslavnode(j)+id)=-2
            endif
          endif
        enddo 
        index=nodempc(3,ist)
!     
        if(index.ne.0) then
          do
            node2=nodempc(1,index)
            do j=1,ntie
              call nident(islavnode(nslavnode(j)+1),node2,
     &             nslavnode(j+1)-nslavnode(j),id)
              if(id.gt.0) then
                if(islavnode(nslavnode(j)+id).eq.node2) then
                  islavact(nslavnode(j)+id)=-2
                endif
              endif
            enddo
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        endif
      enddo
!     
      return
      end
      
