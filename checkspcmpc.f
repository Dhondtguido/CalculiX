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
!     check whether SPC's and MPC's in salve nodes are compatible
!     with mortar contact 
!     for 3D calculations all slave nodes involved in SPCs/MPCs (dependent and independent) are set to noLM nodes
!     
!     Author: Saskia Sitzmann
!     
!     [in] islavnode	field storing the nodes of the slave surface
!     [in] imastnode	field storing the nodes of the master surfaces
!     [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
!     [in] nmastnode	(i)pointer into field imastnode for contact tie i
!     [in] slavnor		slave normals
!     [in,out] islavact	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
!     [in] nslavspc		(2*i) pointer to islavspc...
!     [in] islavspc         ... which stores SPCs for slave node i
!     [in] nsspc            number of SPC for slave nodes
!     [in] nslavmpc		(2*i) pointer to islavmpc...
!     [in,out] islavmpc		... which stores MPCs for slave node i
!     [in] nsmpc		number of MPC for slave nodes
!     [in] nmspc            number of SPC for master nodes
!     [in] nmastmpc		(2*i) pointer to imastmpc...
!     [in,out] imastmpc		... which stores MPCs for master node i
!     [in] nmmpc		number of MPC for master nodes
!     
      subroutine checkspcmpc(ntie,tieset,islavnode,imastnode,nslavnode,
     &     nmastnode,slavnor,islavact,nboun,ndirboun,xboun,
     &     nodempc,coefmpc,ikboun,ilboun,nmpc2,ipompc2,nodempc2,
     &     nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
     &     nmspc,nmastmpc,imastmpc,nmmpc)
!     
!     check whether SPC's and MPC's in salve nodes are compatible
!     with mortar contact    
!     
!     author: Sitzmann,Saskia
!     
!     islavmpc(2,j)=1  directional blocking
!     islavmpc(2,j)=2  cyclic symmetry
!     
!     imastmpc(2,j)=1  directional blocking
!     imastmpc(2,j)=2  cyclic symmetry
!     imastmpc(2,j)=3  spc with displacement
!     
      implicit none
!     
      logical debug,incompatible,nogap,twod
!     
      character*81 tieset(3,*)
!     
      integer ntie,i,j,k,l,dir,dirind,dirdep,id,node,
     &     islavnode(*),imastnode(*),nslavnode(ntie+1),
     &     nmastnode(ntie+1),islavact(*),nboun,ndirboun(*),
     &     nodempc(3,*),index,nmpc2,ipompc2(*),nodempc2(3,*),
     &     ikboun(*),ilboun(*),nslavspc(2,*),islavspc(2,*),nsspc,
     &     nslavmpc(2,*),islavmpc(2,*),nsmpc,nmspc,nmastmpc(2,*),
     &     imastmpc(2,*),nmmpc,ist,zs(3),dof,node2,nsl,nc
!     
      real*8  xboun(*),coefmpc(*),nn,n(3),fixed_disp,coefdep,
     &     slavnor(3,*),v(3),sp
!     
      debug=.false.
!     
!     remove Lagrange Multiplier contributino for nodes which are
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
      do i=1,nmpc2
        ist=ipompc2(i)
        node=nodempc2(1,ist)
        do j=1,ntie
          call nident(islavnode(nslavnode(j)+1),node,
     &         nslavnode(j+1)-nslavnode(j),id)
          if(id.gt.0) then
            if(islavnode(nslavnode(j)+id).eq.node) then
              islavact(nslavnode(j)+id)=-2
            endif
          endif
        enddo 
        index=nodempc2(3,ist)
!     
        if(index.ne.0) then
          do
            node2=nodempc2(1,index)
            do j=1,ntie
              call nident(islavnode(nslavnode(j)+1),node2,
     &             nslavnode(j+1)-nslavnode(j),id)
              if(id.gt.0) then
                if(islavnode(nslavnode(j)+id).eq.node2) then
                  islavact(nslavnode(j)+id)=-2
                endif
              endif
            enddo
            index=nodempc2(3,index)
            if(index.eq.0) exit
          enddo
        endif
      enddo
!     
      return
      end
      
