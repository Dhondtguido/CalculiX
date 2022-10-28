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
      subroutine getdesiinfo3d(set,istartset,iendset,ialset,nset,
     &     mi,nactdof,ndesi,nodedesi,ntie,tieset,itmp,nmpc,nodempc,
     &     ipompc,nodedesiinv,iponoel,inoel,lakon,ipkon,kon,iregion,
     &     ipoface,nodface,nk,jobnamef)    
!     
!     storing the design variables in nodedesi
!     marking which nodes are design variables in nodedesiinv
!     
!     a node is a design variable if:
!     1) it belongs to the design variable set AND
!     2) not all dofs in the node are defined by SPC's AND
!     3) no MPC is applied to any of its dofs AND
!     4) it belongs to at least one face whose number of
!     design variables exceeds half its nodes
!     
      implicit none
!     
      character*8 lakon(*)
!     
      character*81 setname
      character*81 set(*)
      character*81 tieset(3,*)
      character*132 jobnamef
      character*256 fn
!     
      integer mi(*),istartset(*),iendset(*),ialset(*),ndesi,
     &     node,nodedesi(*),nset,ntie,i,j,k,l,m,nmpc,nodempc(3,*),
     &     nactdof(0:mi(2),*),itmp(*),ntmp,index1,id,ipompc(*),
     &     nodedesiinv(*),iponoel(*),inoel(2,*),nelem,nope,nopedesi,
     &     ipkon(*),nnodes,kon(*),iregion,konl(26),iaux,kflag,
     &     ipoface(*),nodface(5,*),jfacem,nopesurf(9),ifaceq(8,6),
     &     ifacet(6,4),ifacew1(4,5),ifacew2(8,5),nopem,nk,iwrite,
     &     ilen
!     
      setname(1:1)=' '
      ndesi=0
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
!     Search for the set name of the set with the design variables
!     
      do i=1,ntie
        if(tieset(1,i)(81:81).eq.'D') then
          setname=tieset(2,i)
        endif
      enddo 
!     
!     Check for the existence of the name
!     
      if(setname(1:1).eq.' ') then
        write(*,*) '*ERROR in getdesiinfo: name of node set '
        write(*,*) '  has not yet been defined. '
        call exit(201)
      endif
!     
!     catalogue all nodes (dependent and independent) which
!     belong to MPC's and sort them in increasing order
!     
      ntmp=0
      do i=1,nmpc
        index1=ipompc(i)
        do
          if(index1.eq.0) exit
          node=nodempc(1,index1)
          call nident(itmp,node,ntmp,id)
          if(id.gt.0) then
            if(itmp(id).eq.node) then
              index1=nodempc(3,index1)
              cycle
            endif
          endif
          ntmp=ntmp+1
          do j=ntmp,id+2,-1
            itmp(j)=itmp(j-1)
          enddo
          itmp(id+1)=node
          index1=nodempc(3,index1)
        enddo
      enddo
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
!     Search the name of the node set in "set(i)" and
!     assign the nodes of the set to the appropriate variables
!     
c      do i=1,nset
c        if(setname.eq.set(i)) then  
      call cident81(set,setname,nset,i)
      if(i.gt.0) then
        if(setname.eq.set(i)) then
          loop1: do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
              node=ialset(j)
!     
!     check for SPC-constraints: if a node is constrained in 
!     all dofs it is removed from the design node set
!     
              do l=1,3
                if(nactdof(l,node).gt.0) exit
                if(l.eq.3) then
                  write(*,*) '*WARNING in getdesiinfo:'
                  write(*,*) '         node ',node,' has no'
                  write(*,*) '         active dofs and'
                  write(*,*) '         is removed from the set'
                  write(*,*) '         of design variables'
                  write(40,*) node
                  iwrite=1
                  cycle loop1
                endif
              enddo
!     
!     check for MPC-constraints
!     
              call nident(itmp,node,ntmp,id)
              if(id.gt.0) then
                if(itmp(id).eq.node) then
                  write(*,*) '*WARNING in getdesiinfo:'
                  write(*,*) '       node ',node,' is subject'
                  write(*,*) '       to MPC-constraints and'
                  write(*,*) '       is removed from the set'
                  write(*,*) '       of design variables'
                  write(40,*) node
                  iwrite=1
                  cycle loop1
                endif
              endif
!     
              ndesi=ndesi+1
              nodedesi(ndesi)=node
            else
              k=ialset(j-2)
              loop2: do
                k=k-ialset(j)
                if(k.ge.ialset(j-1)) exit
!     
!     check for SPC-constraints: if a node is constrained in 
!     all dofs it is removed from the design node set
!     
                do l=1,3
                  if(nactdof(l,k).gt.0) exit
                  if(l.eq.3) then
                    write(*,*) '*WARNING in getdesiinfo:'
                    write(*,*) '         node ',k,' has no'
                    write(*,*) '         active dofs and'
                    write(*,*) '         is removed from the set'
                    write(*,*) '         of design variables'
                    write(40,*) k
                    iwrite=1
                    cycle loop2
                  endif
                enddo
!     
!     check for MPC-constraints
!     
                call nident(itmp,k,ntmp,id)
                if(id.gt.0) then
                  if(itmp(id).eq.k) then
                    write(*,*) '*WARNING in getdesiinfo:'
                    write(*,*) '   node ',k,' is subject'
                    write(*,*) '   to MPC-constraints and'
                    write(*,*) '   is removed from the set'
                    write(*,*) '   of design variables'
                    write(40,*) k
                    iwrite=1
                    cycle loop2
                  endif
                endif
!     
                ndesi=ndesi+1
                nodedesi(ndesi)=k
              enddo loop2
            endif
          enddo loop1
        endif
      endif
c        endif
c      enddo 
!     
!     creating field nodedesiinv indicating for each node whether
!     it is a design variable or not
!     
      do i=1,ndesi
        index1=nodedesi(i)
        nodedesiinv(index1)=-1
      enddo
!     
      kflag=1
      call isortii(nodedesi,iaux,ndesi,kflag)
!     
!     A design node is also removed from nodedesi if it does not
!     belong to a face whose number of design variables exceeds half
!     of its nodes
!     
      do i=1,nk  
        node=i
        if(ipoface(node).eq.0) cycle
        index1=ipoface(node)
        do
          nelem=nodface(3,index1)
          jfacem=nodface(4,index1)
!     
          if(lakon(nelem)(4:4).eq.'8') then
            nope=8
            nopedesi=3
            nopem=4
          elseif(lakon(nelem)(4:5).eq.'20') then
            nope=20
            nopedesi=5
            nopem=8
          elseif(lakon(nelem)(4:5).eq.'10') then
            nope=10
            nopedesi=4
            nopem=6
          elseif(lakon(nelem)(4:4).eq.'4') then
            nope=4
            nopedesi=3
            nopem=3
          elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
            if(jfacem.le.2) then
              nopem=3
              nopedesi=3
            else
              nopem=4
              nopedesi=3
            endif
          elseif(lakon(nelem)(4:5).eq.'15') then
            nope=15
            if(jfacem.le.2) then
              nopem=6
              nopedesi=4
            else
              nopem=8
              nopedesi=5
            endif
          endif
          if(iregion.eq.0) nopedesi=0
!     
!     actual position of the nodes belonging to the
!     master surface
!     
          do k=1,nope
            konl(k)=kon(ipkon(nelem)+k)
          enddo
!     
          if((nope.eq.20).or.(nope.eq.8)) then
            do m=1,nopem
              nopesurf(m)=konl(ifaceq(m,jfacem))
            enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
            do m=1,nopem
              nopesurf(m)=konl(ifacet(m,jfacem))
            enddo
          elseif(nope.eq.15) then
            do m=1,nopem
              nopesurf(m)=konl(ifacew2(m,jfacem))
            enddo
          else
            do m=1,nopem
              nopesurf(m)=konl(ifacew1(m,jfacem))
            enddo
          endif
!     
!     sum up how many designvariables are on that surface
!     
          nnodes=0
          do m=1,nopem
            if(nodedesiinv(nopesurf(m)).ne.0) then
              nnodes=nnodes+1
            endif
          enddo
!     
          if(nnodes.ge.nopedesi) then
            do m=1,nopem
              if(nodedesiinv(nopesurf(m)).eq.-1) then
                nodedesiinv(nopesurf(m))=1
              endif
            enddo
          endif
          index1=nodface(5,index1)
          if(index1.eq.0) exit      
        enddo
      enddo
!     
!     if node i in nodedesi(i) is -1 --> delete node i from 
!     set of designvariables
!     
      do i=1,nk
        if(nodedesiinv(i).eq.-1) then
!     
          write(*,*) '*WARNING in getdesiinfo:'
          write(*,*) '          node ',i,' is removed'
          write(*,*) '          from the set of design'
          write(*,*) '          variables as not sufficient '
          write(*,*) '          other variables are on the  '
          write(*,*) '          surrounding element faces  '
          write(40,*) i
          iwrite=1
!     
          nodedesiinv(i)=0
          call nident(nodedesi,i,ndesi,id)
          do k=id+1,ndesi
            nodedesi(k-1)=nodedesi(k)
          enddo
          ndesi=ndesi-1    
        endif
      enddo 
!
      if(iwrite.eq.1) then
        write(*,*) '*INFO in getdesiinfo:'
        write(*,*) '      rejected design nodes are stored in'
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
