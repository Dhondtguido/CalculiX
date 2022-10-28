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
      subroutine addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &  nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,nactdof,mi,
     &  imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun,ilboun)
!
!     node was kept by the user in a modal dynamics calculation;
!     the present routine checks DOF k of node; if this DOF belongs
!     to a MPC all independent nodes and DOF's of the MPC have to be kept
!
      implicit none
!
      integer node,k,idof,ikmpc(*),ilmpc(*),ipompc(*),nodempc(3,*),
     &  nmpc,imdnode(*),nmdnode,imddof(*),nmddof,id,ist,index,jdof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*)
!
      idof=nactdof(k,node)
c      write(*,*) 'addimdnodedof ',node,k,idof
      if(idof.le.0) then
         idof=(node-1)*8+k
!
!        checking for mpc's
!
         call nident(ikmpc,idof,nmpc,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
               call addimd(imdmpc,nmdmpc,ilmpc(id))
               id=ilmpc(id)
               ist=ipompc(id)
               index=nodempc(3,ist)
               if(index.ne.0) then
                  do
                     call addimd(imdnode,nmdnode,nodempc(1,index))
                     jdof=nactdof(nodempc(2,index),nodempc(1,index))
                     if(jdof.gt.0) call addimd(imddof,nmddof,jdof)
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
            endif
         endif
!
!        checking for spc's
!
         call nident(ikboun,idof,nboun,id)
         if(id.gt.0) then
            if(ikboun(id).eq.idof) then
               call addimd(imdboun,nmdboun,ilboun(id))
            endif
         endif
      else
         call addimd(imddof,nmddof,idof)
      endif
!
      return
      end
