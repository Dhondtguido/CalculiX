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
      subroutine planempc(ipompc,nodempc,coefmpc,
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
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,inode,jmax,k,nodec,
     &  m
!
      real*8 coefmpc(3,*),co(3,*),dd,dmax,pac(3),pbc(3),xboun(*)
!
      save nodea,nodeb,nodec,jmax
!
      if(inode.eq.1) then
         nodea=node
         return
      elseif(inode.eq.2) then
         nodeb=node
         return
      elseif(inode.eq.3) then
         nodec=node
         do j=1,3
            pac(j)=co(j,nodea)-co(j,nodec)
            pbc(j)=co(j,nodeb)-co(j,nodec)
         enddo
         dmax=abs(pac(2)*pbc(3)-pac(3)*pbc(2))
         jmax=1
         dd=abs(pac(1)*pbc(3)-pac(3)*pbc(1))
         if(dd.gt.dmax) then
            dmax=dd
            jmax=2
         endif
         dd=abs(pac(1)*pbc(2)-pac(2)*pbc(1))
         if(dd.gt.dmax) then
            dmax=dd
            jmax=3
         endif
         return
      endif
!
      nk=nk+1
      if(nk.gt.nk_) then
         write(*,*) '*ERROR in planempc: increase nk_'
         call exit(201)
      endif
!
      j=jmax
      k=j+1
      if(k.gt.3) k=1
      l=k+1
      if(l.gt.3) l=1
!
      idof=8*(node-1)+j
      call nident(ikmpc,idof,nmpc,id)
      if(id.gt.0) then
         if(ikmpc(id).eq.idof) then
            write(*,*) '*WARNING in planempc: DOF for node ',node
            write(*,*) '         in direction ',j,' has been used'
            write(*,*) '         on the dependent side of another MPC'
            write(*,*) '         PLANE constraint cannot be applied'
            return
         endif
      endif
      nmpc=nmpc+1
      if(nmpc.gt.nmpc_) then
         write(*,*) '*ERROR in planempc: increase nmpc_'
         call exit(201)
      endif
!
      ipompc(nmpc)=mpcfree
      labmpc(nmpc)='PLANE               '
!
      do m=nmpc,id+2,-1
         ikmpc(m)=ikmpc(m-1)
         ilmpc(m)=ilmpc(m-1)
      enddo
      ikmpc(id+1)=idof
      ilmpc(id+1)=nmpc
!
      nodempc(1,mpcfree)=node
      nodempc(2,mpcfree)=j
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=node
      nodempc(2,mpcfree)=k
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=node
      nodempc(2,mpcfree)=l
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodea
      nodempc(2,mpcfree)=j
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodea
      nodempc(2,mpcfree)=k
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodea
      nodempc(2,mpcfree)=l
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodeb
      nodempc(2,mpcfree)=j
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodeb
      nodempc(2,mpcfree)=k
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodeb
      nodempc(2,mpcfree)=l
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodec
      nodempc(2,mpcfree)=j
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodec
      nodempc(2,mpcfree)=k
      mpcfree=nodempc(3,mpcfree)
      nodempc(1,mpcfree)=nodec
      nodempc(2,mpcfree)=l
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
         write(*,*) '*ERROR in planempc: increase nboun_'
         call exit(201)
      endif
      nodeboun(nboun)=nk
      ndirboun(nboun)=j
      typeboun(nboun)='U'
      do m=nboun,id+2,-1
         ikboun(m)=ikboun(m-1)
         ilboun(m)=ilboun(m-1)
      enddo
      ikboun(id+1)=idof
      ilboun(id+1)=nboun
!
       return
       end
