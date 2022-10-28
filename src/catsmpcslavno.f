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
!     Author: Saskia Sitzmann
!     
!     [in] nboun2           number of transformed SPCs
!     [in] nmpc2		number of transformed mpcs
!     [in] ipompc2           (i) pointer to nodempc and coeffmpc for transformed MPC i
!     [in] nodempc2          nodes and directions of transformed MPCs
!     [in] ikboun2           sorted dofs idof=8*(node-1)+dir for transformed SPCs
!     [in] ilboun2           transformed SPC numbers for sorted dofs
!     [in] ikmpc2 		sorted dofs idof=8*(node-1)+dir for transformed MPCs
!     [in] ilmpc2		transformed SPC numbers for sorted dofs
!     [out] nslavspc	(2*i) pointer to islavspc...
!     [out] islavspc         ... which stores SPCs for slave node i
!     [out] nsspc            number of SPC for slave nodes
!     [out] nslavmpc	(2*i) pointer to islavmpc...
!     [out] islavmpc	... which stores MPCs for slave node i
!     [out] nsmpc		number of MPC for slave nodes
!     [out] nslavspc2	(2*i) pointer to islavspc2...
!     [out] islavspc2       ... which stores transformed SPCs for slave node i
!     [out] nsspc2          number of transformed SPC for slave nodes
!     [out] nslavmpc2	(2*i) pointer to islavmpc2...
!     [out] islavmpc2	... which stores transformed MPCs for slave node i
!     [out] nsmpc2		number of transformed MPC for slave nodes 
!     [out] nmastspc	(2*i) pointer to imastspc...
!     [out] imastspc        ... which stores SPCs for master node i
!     [out] nmspc           number of SPC for master nodes
!     [out] nmastmpc	(2*i) pointer to imastmpc...
!     [out] imastmpc	... which stores MPCs for master node i
!     [out] nmmpc		number of MPC for master nodes
!     
      subroutine catsmpcslavno(ntie,islavnode,imastnode,
     &     nslavnode,nmastnode,nboun,ndirboun,nodeboun,nmpc,ipompc,
     &     nodempc,ikboun,ilboun,ikmpc,ilmpc,nboun2,nmpc2,ipompc2,
     &     nodempc2,ikboun2,ilboun2,ikmpc2,ilmpc2,
     &     nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
     &     nslavspc2,islavspc2,nsspc2,nslavmpc2,islavmpc2,nsmpc2,
     &     nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
     &     nmastmpc2,imastmpc2,nmmpc2,jobnamef)
!     
!     subroutine to catalogue SPC'S and MPC'S for master and slave node    
!     needed for mortar contact
!     
!     correlate the transformed and untransformed spc's and mpc's to the
!     slave and master nodes
!     
!     author: Sitzmann,Saskia
!     
      implicit none
!     
      logical debug,testm,isslavenode,ismastnode,
     &     test1to1
!
      character*132 jobnamef(*)
      character*256 fn
!     
      integer ntie,i,j,k,l,iwrite,ilen,
     &     id,node,islavnode(*),imastnode(*),nslavnode(ntie+1),
     &     nmastnode(ntie+1),nmmpc2,index1,nboun,ndirboun(*),
     &     nmpc,ipompc(*),nodempc(3,*),idof,nboun2,nmpc2,ipompc2(*),
     &     nodempc2(3,*),ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &     ikboun2(*),ilboun2(*),ikmpc2(*),ilmpc2(*),nslavspc(2,*),
     &     islavspc(2,*),nsspc,nslavmpc(2,*),islavmpc(2,*),nsmpc,
     &     nslavspc2(2,*),islavspc2(2,*),nsspc2,nslavmpc2(2,*),
     &     islavmpc2(2,*),nsmpc2,nmastspc(2,*),imastspc(2,*),nmspc,
     &     nmastmpc(2,*),imastmpc(2,*),nmmpc,isspc,imspc,ismpc,immpc,
     &     ist,nmastmpc2(2,*),imastmpc2(2,*),secondnode,node2,itie,
     &     nodeboun(*)
!     
      debug=.false.
!     
!     slave surfaces
!     
      isspc=0
      ismpc=0
      do i=1,ntie
        do l=nslavnode(i)+1,nslavnode(i+1)
          node=islavnode(l)
!     
!     check for SPCs
!     
          nslavspc(1,l)=isspc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
              if(idof.eq.ikboun(id)) then
                isspc=isspc+1
                islavspc(1,isspc)=ilboun(id)
              endif
            endif
          enddo
          nslavspc(2,l)=isspc
!     
!     check for MPCs
!     
          nslavmpc(1,l)=ismpc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
              if(idof.eq.ikmpc(id)) then
                ismpc=ismpc+1
                islavmpc(1,ismpc)=ipompc(ilmpc(id))
              endif
            endif
          enddo
          nslavmpc(2,l)=ismpc
        enddo
      enddo
      nsspc=isspc
      nsmpc=ismpc
!     
      isspc=0
      ismpc=0
      i=0
      do i=1,ntie
        do l=nslavnode(i)+1,nslavnode(i+1)
          node=islavnode(l)
!     
!     check for SPCs (2)
!     
          nslavspc2(1,l)=isspc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikboun2,idof,nboun2,id)
            if(id.gt.0) then
              if(idof.eq.ikboun2(id)) then
                isspc=isspc+1
                islavspc2(1,isspc)=ilboun2(id)
              endif
            endif
          enddo
          nslavspc2(2,l)=isspc
!     
!     check for MPCs (2)
!     
          nslavmpc2(1,l)=ismpc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikmpc2,idof,nmpc2,id)
            if(id.gt.0) then
              if(idof.eq.ikmpc2(id)) then
                ismpc=ismpc+1
                islavmpc2(1,ismpc)=ipompc2(ilmpc2(id))
              endif
            endif
          enddo
          nslavmpc2(2,l)=ismpc
        enddo
      enddo
      nsspc2=isspc
      nsmpc2=ismpc
!     
!     master surfaces
!     
      imspc=0
      immpc=0
      do i=1,ntie
        do l=nmastnode(i)+1,nmastnode(i+1)
          node=imastnode(l)
!     
!     check for SPCs
!     
          nmastspc(1,l)=imspc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
              if(idof.eq.ikboun(id)) then
                imspc=imspc+1
                imastspc(1,imspc)=ilboun(id)
              endif
            endif
          enddo
          nmastspc(2,l)=imspc
!     
!     check for MPCs
!     
          nmastmpc(1,l)=immpc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
              if(idof.eq.ikmpc(id)) then
                immpc=immpc+1
                imastmpc(1,immpc)=ipompc(ilmpc(id))
              endif
            endif
          enddo
          nmastmpc(2,l)=immpc
        enddo
      enddo
      nmspc=imspc
      nmmpc=immpc
!     
      immpc=0
      do i=1,ntie
        do l=nmastnode(i)+1,nmastnode(i+1)
          node=imastnode(l)
!     
!     check for MPCs (2)
!     
          nmastmpc2(1,l)=immpc
          do k=1,3
            idof=8*(node-1)+k
            call nident(ikmpc2,idof,nmpc2,id)
            if(id.gt.0) then
              if(idof.eq.ikmpc2(id)) then
                immpc=immpc+1
                imastmpc2(1,immpc)=ipompc2(ilmpc2(id))
              endif
            endif
          enddo
          nmastmpc2(2,l)=immpc
        enddo
      enddo
      nmmpc2=immpc
!
!     check for not supported mpc's with slave/master nodes involved
!
!     The following rules apply:
!
!     If the dependent node in a MPC is a slave node then:
!     - none of the independent nodes is allowed to be a master node
!     - the MPC is allowed to connect at most 2 different nodes
!
!     If the dependent node is not a slave node, then none of the
!     dependent nodes is allowed to be a slave node
!
!     This allows e.g.:
!     - a one-to-one connection of two slave or two master nodes in a
!       cyclic symmetry MPC
!     - a (non)homogeneous SPC in a local system on the slave side or
!       master side
!
!     This does not allow:
!     - a slave node connected to one or more master nodes
!     - a master node connected to one or more slave nodes
!
!     
!     opening a file to store the nodes which are not connected
!
      iwrite=0
      ilen=index(jobnamef(1),' ')-1
      fn=jobnamef(1)(1:ilen)//'_WarnSlaveNodeUnallowedMpc.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnSlaveNodeUnallowedMpc'
!
      do i=1,nmpc
        ist=ipompc(i)
        node=nodempc(1,ist)
!
!       is node slave node?
!
        isslavenode=.false.
        do j=1,ntie
          call nident(islavnode(nslavnode(j)+1),node,
     &         nslavnode(j+1)-nslavnode(j),id)
          if(id.gt.0) then
            if(islavnode(nslavnode(j)+id).eq.node) then
              isslavenode=.true.
              itie=j
            endif
          endif
        enddo
!
        if(isslavenode) then
!
!     test for (in)homogeneous SPC in local coordinates or 1-to-1 cyclic symmetry
!
          secondnode=node
          test1to1=.true.
          testm=.false.
          index1=nodempc(3,ist)
          do
            if(index1.eq.0) exit
            node2=nodempc(1,index1)
            if((node.eq.secondnode).and.(node2.ne.secondnode)) then
              secondnode=node2
            endif
            if(node2.ne.node) then
              do j=1,ntie
                call nident(imastnode(nmastnode(j)+1),node2,
     &               nmastnode(j+1)-nmastnode(j),id)
                if(id.gt.0) then
                  if(imastnode(nmastnode(j)+id).eq.node2) then
                    ismastnode=.true.
                  endif
                endif
              enddo
              if((ismastnode).and.(.not.testm)) then
                testm=.true.
                write(*,*) '*ERROR in catsmpcslavno: slave node',
     &               node,',is connected als dependent node in'
                write(*,*) '       a MPC to master node ',node2
                write(40,*) node
                iwrite=1
              endif
            endif
            if((node2.ne.node).and.(node2.ne.secondnode)) then
              test1to1=.false.
            endif
            index1=nodempc(3,index1)
          enddo
          if(.not.test1to1) then
            write(*,*) '*ERROR in catsmpcslavno: slave node',
     &           node,', is connected as dependent node in a'
            write(*,*) '       one-to-m (m>1) mpc !'
            write(40,*) node
            iwrite=1
          endif
        else
!
!         test if one of the independent nodes is a slave node
!
          index1=nodempc(3,ist)
          loop: do
            if(index1.eq.0) exit
            node2=nodempc(1,index1)
!     
!           is node2 slavenode?
!
            if(node2.ne.node) then
              do j=1,ntie
                call nident(islavnode(nslavnode(j)+1),node2,
     &               nslavnode(j+1)-nslavnode(j),id)
                if(id.gt.0) then
                  if(islavnode(nslavnode(j)+id).eq.node2) then
                    write(*,*) '*ERROR in catsmpcslavno: ',
     &                   ', invalid mpc found! ',
     &                   'slave node',node2,' is used as ',
     &                   'independent variable in a MPC with ',
     &                   'the non-slave node',node,
     &                   ' as dependent variable'
                    write(40,*) node2
                    iwrite=1
                  endif
                endif
              enddo
            endif
            index1=nodempc(3,index1)
          enddo loop
        endif
      enddo
!     
      if(debug) then
        do i=1,ntie
          do l=nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(l)
            write(*,*)'***node',node,'***'
            write(*,*) 'nodes-spc',node,nslavspc(1,l),nslavspc(2,l)
            do j=nslavspc(1,l)+1,nslavspc(2,l)
              ist=islavspc(1,j)
              write(*,*)' spc',ist,nodeboun(ist),ndirboun(ist)
            enddo
            write(*,*) 'nodes-mpc',node,nslavmpc(1,l),nslavmpc(2,l)
            do j=nslavmpc(1,l)+1,nslavmpc(2,l)
              ist=islavmpc(1,j)
              index1=nodempc(3,ist)
              write(*,*)' mpc',ist,nodempc(2,ist),nodempc(1,index1)
            enddo
            write(*,*)'nodes-mpc2',node,nslavmpc2(1,l),nslavmpc2(2,l)
            do j=nslavmpc2(1,l)+1,nslavmpc2(2,l)
              ist=islavmpc2(1,j)
              index1=nodempc2(3,ist)
              write(*,*)' mpc',ist,nodempc2(2,ist),nodempc2(1,index1)
            enddo
          enddo
        enddo
        do i=1,ntie
          do l=nmastnode(i)+1,nmastnode(i+1)
            node=imastnode(l)
            write(*,*)'***node',node,'***'
            write(*,*) 'nodem-spc',node,nmastspc(1,l),nmastspc(2,l)
            do j=nmastspc(1,l)+1,nmastspc(2,l)
              ist=imastspc(1,j)
              write(*,*)'spc',ist,nodeboun(ist),ndirboun(ist)
            enddo
            write(*,*) 'nodem-mpc',node,nmastmpc(1,l),nmastmpc(2,l)
            do j=nmastmpc(1,l)+1,nmastmpc(2,l)
              ist=imastmpc(1,j)
              index1=nodempc(3,ist)
              write(*,*)'mpc',ist,nodempc(2,ist),nodempc(1,index1)
            enddo
            write(*,*) 'nodem-mpc2',node,nmastmpc2(1,l),
     &           nmastmpc2(2,l)
            do j=nmastmpc2(1,l)+1,nmastmpc2(2,l)
              ist=imastmpc2(1,j)
              index1=nodempc2(3,ist)
              write(*,*)'mpc',ist,nodempc2(2,ist),
     &             nodempc2(1,index1)
            enddo
          enddo
        enddo
      endif
!
      if(iwrite.eq.1) then
        write(*,*) '*ERROR in catsmpcslavno:'
        write(*,*) '       slavenodes belonging to unallowed MPCs'
        write(*,*) '       are stored in file'
        write(*,*) '       ',fn(1:ilen+30)
        write(*,*) '       This file can be loaded into'
        write(*,*) '       an active cgx-session by typing'
        write(*,*) 
     &       '       read ',fn(1:ilen+30),' inp'
        write(*,*) '       Remove the faces to which these'
        write(*,*) '       nodes belong from the slave face'
        write(*,*) '       definition'
        write(*,*)
        close(40)
        call exit(201)
      else
        close(40,status='delete')
      endif
!     
      return
      end
      
