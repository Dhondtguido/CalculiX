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
      subroutine straightmpc(ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,xboun,inode,node,co,typeboun)
!
!     generates MPC's for nodes staying on a straight line defined
!     by two nodes a and b. The parameter inode indicates how many
!     times the present routine was called within the same *MPC 
!     definition. For inode=1 "node" is node a, for inode=2 "node"
!     is node b. Starting with inode=3 MPC's are defined.
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,
     &  ikmpc(*),
     &  ilmpc(*),node,id,mpcfreeold,j,idof,l,nodeboun(*),nodea,nodeb,
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,inode,jmax,k
!
      real*8 coefmpc(3,*),co(3,*),dd,dmax,xboun(*)
!
      save nodea,nodeb,jmax
!
      if(inode.eq.1) then
         nodea=node
         return
      elseif(inode.eq.2) then
         nodeb=node
         dmax=0.d0
         do k=1,3
            dd=abs((co(k,nodea)-co(k,nodeb)))
            if(dd.gt.dmax) then
               dmax=dd
               jmax=k
            endif
         enddo
         return
      endif
!
      nk=nk+1
      if(nk.gt.nk_) then
         write(*,*) '*ERROR in straightmpc: increase nk_'
         call exit(201)
      endif
      do j=1,3
         if(j.eq.jmax) cycle
         idof=8*(node-1)+j
         call nident(ikmpc,idof,nmpc,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
             write(*,*) '*WARNING in straightmpc: DOF for node ',node
             write(*,*) '         in direction ',j,' has been used'
             write(*,*) '         on the dependent side of another MPC'
             write(*,*) '         STRAIGHT constraint cannot be applied'
             cycle
            endif
         endif
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR in straightmpc: increase nmpc_'
            call exit(201)
         endif
!
         ipompc(nmpc)=mpcfree
         labmpc(nmpc)='STRAIGHT            '
!
         do l=nmpc,id+2,-1
            ikmpc(l)=ikmpc(l-1)
            ilmpc(l)=ilmpc(l-1)
         enddo
         ikmpc(id+1)=idof
         ilmpc(id+1)=nmpc
!
         nodempc(1,mpcfree)=node
         nodempc(2,mpcfree)=j
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=node
         nodempc(2,mpcfree)=jmax
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nodea
         nodempc(2,mpcfree)=j
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nodea
         nodempc(2,mpcfree)=jmax
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nodeb
         nodempc(2,mpcfree)=j
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nodeb
         nodempc(2,mpcfree)=jmax
         mpcfree=nodempc(3,mpcfree)
         nodempc(1,mpcfree)=nk
         nodempc(2,mpcfree)=j
         mpcfreeold=mpcfree
         mpcfree=nodempc(3,mpcfree)
         nodempc(3,mpcfreeold)=0
         idof=8*(nk-1)+j
         call nident(ikboun,idof,nboun,id)
         nboun=nboun+1
         if(nboun.gt.nboun_) then
            write(*,*) '*ERROR in straightmpc: increase nboun_'
            call exit(201)
         endif
         nodeboun(nboun)=nk
         ndirboun(nboun)=j
         typeboun(nboun)='U'
         do l=nboun,id+2,-1
            ikboun(l)=ikboun(l-1)
            ilboun(l)=ilboun(l-1)
         enddo
         ikboun(id+1)=idof
         ilboun(id+1)=nboun
       enddo
!
       return
       end
