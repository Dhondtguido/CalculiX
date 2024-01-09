!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2023 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     i
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
      subroutine getdesiinfo2d(set,istartset,iendset,ialset,nset,
     &     mi,nactdof,ndesi,nodedesi,ntie,tieset,nodedesiinv,lakon,
     &     ipkon,kon,iponoelfa,nod2nd3rd,iponor2d,knor2d,
     &     iponoel2d,inoel2d,nobject,objectset,nod1st,ne,
     &     jobnamef,rig)    
!     
!     purpose of this routine:     
!     
!     1) replace the nodes in nodal design response sets which
!        belong to 2D-elements by the lowest node number in the     
!        3D-expansion of these elements
!     
!     2) if a node belonging to the design variable set 
!        - belongs to a 2D-element     
!        - is not a dependent MPC-node AND
!        - has not all its dofs restricted by SPC's
!        then it is declared as a design node. Its number is replaced   
!        in the field of design nodes "nodedesi" by the lowest node
!        number in the 3D-expansion of the element
!     
!     3) for such a design node i field node2nd3rd is created pointing    
!        to the second (node2nd3rd(1,*)) and third (node2nd3rd(2,*))
!        expansion node
!     
!     4) field nodedesiinv(1..nk) is created with nodedesinv(i)=1      
!        is node i is a design variable, else =0.
!     
!     5) field nod1st(1..nk) is created pointing to the first node 
!        in the expansion for every 2D design node
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 setname
      character*81 set(*)
      character*81 tieset(3,*)
      character*81 objectset(5,*)
      character*132 jobnamef
      character*256 fn
!     
      integer mi(*),istartset(*),iendset(*),ialset(*),ndesi,ilen,
     &     node,nodedesi(*),nset,ntie,i,j,k,l,nactdof(0:mi(2),*),index1,
     &     nodedesiinv(*),ipkon(*),kon(*),iaux,kflag,nod2nd3rd(2,*),
     &     iponoelfa(*),inoel2d(3,*),iset,iponoel2d(*),nodeold,ipos1,
     &     ipos2,ielem,iponor2d(2,*),num,knor2d(*),inode,nodenew,nope2d,
     &     ishift,nobject,iobject,numtest,nod1st(*),ne,id,iwrite,
     &     index2d,rig(*)
!
      integer,dimension(:),allocatable::itreated
!     
      setname(1:1)=' '
      ndesi=0
!     
!     Search for the name of the set with the design variables
!     
      do iset=1,ntie
        if(tieset(1,iset)(81:81).eq.'D') then
          setname=tieset(2,iset)
        endif
      enddo 
      call cident81(set,setname,nset,id)
      iset=nset+1
      if(id.gt.0) then
        if(setname.eq.set(id)) then
          iset=id
        endif
      endif
!     
!     Check for the existence of the name
!     
      if(setname(1:1).eq.' ') then
        write(*,*) '*ERROR in getdesiinfo2d: name of node set '
        write(*,*) '  has not yet been defined. '
        call exit(201)
      endif
!     
!     Change the node numbers in the sets for the objective and constraint
!     function
!
      allocate(itreated(nset))
      do i=1,nset
        itreated(i)=0
      enddo
!      
      do iobject=1,nobject
!     
!     only node-based design responses are treated
!     
        if((objectset(1,iobject)(1:12).ne."STRAINENERGY").and.
     &       (objectset(1,iobject)(1:4).ne."MASS").and.
     &       (objectset(1,iobject)(1:4).ne."EIGENFREQUENCY")) then
          call cident81(set,objectset(3,iobject),nset,i)
          if(i.gt.0) then
            if(objectset(3,iobject).eq.set(i)) then
!     
!     design variables are treated later
!     (set of design and objective variables may coincide)
!     
              if((objectset(3,iobject).ne.setname).and.
     &             (itreated(i).ne.1)) then
                itreated(i)=1
                do inode=istartset(i),iendset(i)
                  nodeold=ialset(inode)
c                  write(*,*) 'getdesiinfo2d',nodeold,iponoel2d(nodeold)
                  index2d=iponoel2d(nodeold)
                  if(index2d.eq.0) cycle
                  ielem=inoel2d(1,index2d)
!     
!                 replace 2D-element nodes belonging to the
!                 design response set by the lowest node number
!                 of the node's expansion
!     
                  if((lakon(ielem)(7:7).eq.'A').or.
     &                 (lakon(ielem)(7:7).eq.'E').or.
     &                 (lakon(ielem)(7:7).eq.'L').or.
     &                 (lakon(ielem)(7:7).eq.'S')) then
!     
                    k=inoel2d(2,index2d)
                    ipos2=iponor2d(2,ipkon(ielem)+k)
                    nodenew=knor2d(ipos2+1)
                    ialset(inode)=nodenew
                  endif
!     
                enddo
              endif
            endif
          endif
        endif
      enddo
      deallocate(itreated)
!     
!     opening a file to store the nodes which are rejected as
!     design variables
!
      iwrite=0
      ilen=index(jobnamef,' ')-1
      fn=jobnamef(1:ilen)//'_WarnNodeDesignReject.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeDesignReject'
!     
!     Rename the design variables and save the 2 additional expanded
!     nodes in nod2nd3rd:
!     if node i is a design variable --> 1st neighbor=nod2nd3rd(1,i)
!                                    --> 2nd neighbor=nod2nd3rd(2,i)
!     
      loop1: do inode=istartset(iset),iendset(iset)
        nodeold=ialset(inode)
!
!       check that de node is no expandable rigid body was defined
!       in the node
!
        if(rig(nodeold).ne.0) then
          write(*,*) '*WARNING in getdesiinfo2d:'
          write(*,*) '       in node ',nodeold,' an expandable'
          write(*,*) '       rigid body is defined and therefore'
          write(*,*) '       it is removed from the set'
          write(*,*) '       of design variables'
          write(40,*) nodeold
          iwrite=1
          cycle loop1
        endif
!     
        index2d=iponoel2d(nodeold)
        if(index2d.eq.0) cycle
        ielem=inoel2d(1,index2d)
!     
!       Determine element formulation
!     
        if((lakon(ielem)(7:7).eq.'A').or.
     &       (lakon(ielem)(7:7).eq.'E').or.
     &       (lakon(ielem)(7:7).eq.'L').or.
     &       (lakon(ielem)(7:7).eq.'S')) then
!
          k=inoel2d(2,index2d)
          ipos2=iponor2d(2,ipkon(ielem)+k)
          nodenew=knor2d(ipos2+1)
          ialset(inode)=nodenew
!     
!         check for the existence of a MPC in the node
!     
          do l=1,3
!     
!           comment for next line: nactdof can be zero since nodeold
!           is a 2D-node which is not necessarily used in the equation     
!           system (only if e.g. a SPC or MPC was defined in the 2D    
!           node)
!     
            if(nactdof(l,nodeold).ge.0) exit
!     
!           check if its an MPC(odd) or SPC(even)
!     
            num=nactdof(l,nodeold)
            numtest=num/2*2
            if(num.ne.numtest) then
              write(*,*) '*WARNING in getdesiinfo2d:'
              write(*,*) '       node ',nodeold,' is a'
              write(*,*) '       dependent dof in a MPC and'
              write(*,*) '       is removed from the set'
              write(*,*) '       of design variables'
              write(40,*) nodeold
              iwrite=1
              cycle loop1
            endif
          enddo
!     
!     check whether not all dofs are removed by SPCs
!     
          do l=1,3
!     
!           comment for next line: nactdof can be zero since nodeold
!           is a 2D-node which is not necessarily used in the equation     
!           system (only if e.g. a SPC or MPC was defined in the 2D    
!           node)
!     
            if(nactdof(l,nodeold).ge.0) exit
            if(l.eq.3) then
              write(*,*) '*WARNING in getdesiinfo2d:'
              write(*,*) '       node ',node,' has no'
              write(*,*) '       active dofs and'
              write(*,*) '       is removed from the set'
              write(*,*) '       of design variables'
              write(40,*) node
              iwrite=1
              cycle loop1
            endif
          enddo
          ialset(inode)=nodenew
          nod2nd3rd(1,nodenew)=knor2d(ipos2+2)
          nod2nd3rd(2,nodenew)=knor2d(ipos2+3)
!     
!         the lowest expanded node number is the design variable     
!     
          ndesi=ndesi+1
          nodedesi(ndesi)=nodenew
        endif
      enddo loop1
!     
!     creating field nodedesiinv indicating for each expanded node whether
!     it is a design variable or not
!     
      do i=1,ndesi
        node=nodedesi(i)
        nodedesiinv(node)=1
        node=nod2nd3rd(1,nodedesi(i))
        nodedesiinv(node)=1
        node=nod2nd3rd(2,nodedesi(i))
        nodedesiinv(node)=1
      enddo
!     
      kflag=1
      call isortii(nodedesi,iaux,ndesi,kflag)
!     
      do ielem=1,ne
        if(ipkon(ielem).lt.0) cycle   
!     
!     Determine element formulation
!     
        if((lakon(ielem)(7:7).eq.'A').or.
     &       (lakon(ielem)(7:7).eq.'E').or.
     &       (lakon(ielem)(7:7).eq.'L').or.
     &       (lakon(ielem)(7:7).eq.'S')) then
          if(lakon(ielem)(4:5).eq.'20') then
            nope2d=8
          elseif(lakon(ielem)(4:5).eq.'8R') then
            nope2d=4
          elseif(lakon(ielem)(4:5).eq.'8 ') then
            nope2d=4
          elseif(lakon(ielem)(4:5).eq.'8I') then
            nope2d=4
           elseif(lakon(ielem)(4:5).eq.'15') then
            nope2d=6
          elseif(lakon(ielem)(4:4).eq.'6') then
            nope2d=3
          else
            cycle
          endif
!     
!     of all 3D-nodes in which a 2D design variable node
!     is expanded only the first 3D-node in the expansion is
!     considered to be a design variable. All other 3D-nodes in
!     the expansion are not design nodes. For these nodes (let us
!     call them i) ipnk2dto3d(i) points to the first 3D-node in
!     the expansion (i.e the node which is taken as design node)
!     
          do k=1,nope2d
            ipos1=ipkon(ielem)+k    
            ipos2=iponor2d(2,ipos1)
            nod1st(knor2d(ipos2+1))=knor2d(ipos2+1)
            nod1st(knor2d(ipos2+2))=knor2d(ipos2+1)
            nod1st(knor2d(ipos2+3))=knor2d(ipos2+1)
          enddo
        endif
      enddo
!     
      if(iwrite.eq.1) then
        write(*,*) '*INFO in getdesiinfo2d:'
        write(*,*) '      rejected design nodes (if any) are stored in'
        write(*,*) '      file ',fn(1:ilen+25)
        write(*,*) '      This file can be loaded into'
        write(*,*) '      an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+25),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!     
      return
      end
