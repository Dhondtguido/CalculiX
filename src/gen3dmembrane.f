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
      subroutine gen3dmembrane(ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,ikmpc,ilmpc,labmpc,nk,ithermal,i)
!
!     connects nodes of 1-D and 2-D elements, for which SPC's were
!     defined, to the nodes of their expanded counterparts
!
      implicit none
!
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  ikmpc(*),ilmpc(*),i,j,idir,nk,newnode,idof,id,mpcfreenew,
     &  ithermal(*),jstart,jend
!
      real*8 coefmpc(*)
!
!     generating a hinge at a node of a membrane element
!
!     u(n_1)+u(n_3)=2*u(n)
!
      newnode=nk-1
!     
      if(ithermal(2).le.1) then
         jstart=1
         jend=3
      elseif(ithermal(2).eq.2) then
         jstart=0
         jend=0
      else
         jstart=0
         jend=3
      endif
!     
      do idir=jstart,jend
         idof=8*(newnode-1)+idir
         call nident(ikmpc,idof,nmpc,id)
         if((id.le.0).or.(ikmpc(id).ne.idof)) then
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) 
     &              '*ERROR in gen3dmembrane: increase nmpc_'
               call exit(201)
            endif
            labmpc(nmpc)='                    '
            ipompc(nmpc)=mpcfree
            do j=nmpc,id+2,-1
               ikmpc(j)=ikmpc(j-1)
               ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
            nodempc(1,mpcfree)=newnode
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) 
     &              '*ERROR in gen3dmembrane: increase memmpc_'
               call exit(201)
            endif
            nodempc(1,mpcfree)=nk+1
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) 
     &              '*ERROR in gen3dmembrane: increase memmpc_'
               call exit(201)
            endif
            nodempc(1,mpcfree)=i
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=-2.d0
            mpcfreenew=nodempc(3,mpcfree)
            if(mpcfreenew.eq.0) then
               write(*,*) 
     &              '*ERROR in gen3dmembrane: increase memmpc_'
               call exit(201)
            endif
            nodempc(3,mpcfree)=0
            mpcfree=mpcfreenew
         endif
      enddo
!     
      return
      end
