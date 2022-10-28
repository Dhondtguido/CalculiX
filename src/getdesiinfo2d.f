!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
     &     ipkon,kon,iponoelfa,iponod2dto3d,iponor2d,knor2d,
     &     iponoel2d,inoel2d,nobject,objectset,iponk2dto3d,ne,
     &     jobnamef)    
!     
!     storing the design variables in nodedesi
!     marking which nodes are design variables in nodedesiinv
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
     &     nodedesiinv(*),ipkon(*),kon(*),iaux,kflag,iponod2dto3d(2,*),
     &     iponoelfa(*),inoel2d(3,*),iset,iponoel2d(*),nodeold,ipos1,
     &     ipos2,ielem,iponor2d(2,*),num,knor2d(*),inode,nodenew,nope2d,
     &     ishift,nobject,iobject,numtest,iponk2dto3d(*),ne,id,iwrite
!     
      setname(1:1)=' '
      ndesi=0
!     
!     Search for the set name of the set with the design variables
!     
      do iset=1,ntie
        if(tieset(1,iset)(81:81).eq.'D') then
          setname=tieset(2,iset)
        endif
      enddo 
c     do iset=1,nset
c     if(setname.eq.set(iset)) then 
c     exit
c     endif
c     enddo
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
      do iobject=1,nobject
!     
!     only node-based objective functions are treated
!     
        if((objectset(1,iobject)(1:12).ne."STRAINENERGY").and.
     &       (objectset(1,iobject)(1:4).ne."MASS")) then
c     do i=1,nset
c     if(objectset(3,iobject).eq.set(i)) then
          call cident81(set,objectset(3,iobject),nset,i)
          if(i.gt.0) then
            if(objectset(3,iobject).eq.set(i)) then
!     
!     design variables are treated later
!     (set of design and objective variables may coincide)
!     
              if(objectset(3,iobject).ne.setname) then
                do inode=istartset(i),iendset(i)
                  nodeold=ialset(inode)
                  if(iponoel2d(nodeold).eq.0) cycle
                  ielem=inoel2d(1,iponoel2d(nodeold))
!     
!     Determine element formulation
!     
                  if((lakon(ielem)(7:7).eq.'A').or.
     &                 (lakon(ielem)(7:7).eq.'E').or.
     &                 (lakon(ielem)(7:7).eq.'L').or.
     &                 (lakon(ielem)(7:7).eq.'S')) then
                    if(lakon(ielem)(4:5).eq.'20') then
                      nope2d=8
                      ishift=20
                    elseif(lakon(ielem)(4:5).eq.'8R') then
                      nope2d=4
                      ishift=8
                    elseif(lakon(ielem)(4:5).eq.'8 ') then
                      nope2d=4
                      ishift=8
                    elseif(lakon(ielem)(4:5).eq.'8I') then
                      nope2d=4
                      ishift=11
                    elseif(lakon(ielem)(4:5).eq.'15') then
                      nope2d=6
                      ishift=15
                    elseif(lakon(ielem)(4:4).eq.'6') then
                      nope2d=3
                      ishift=6
                    else
                      cycle
                    endif
!     
                    do k=1,nope2d
                      ipos1=ipkon(ielem)+ishift+k     
                      if(nodeold.eq.kon(ipos1)) then
                        ipos2=iponor2d(2,ipos1-ishift)
                        nodenew=knor2d(ipos2+1)
                        ialset(inode)=nodenew
                        exit
                      endif
                    enddo
                  endif
                enddo
              endif
            endif
          endif
c     endif
c     enddo
        endif
        if(objectset(1,iobject)(4:13).eq."MEMBERSIZE") then
c     do i=1,nset
c     if(objectset(4,iobject).eq.set(i)) then
          call cident81(set,objectset(4,iobject),nset,i)
          if(i.gt.0) then
            if(objectset(4,iobject).eq.set(i)) then
!     
!     check whether design variable set is changed
!     
              if(objectset(4,iobject).ne.setname) then
                do inode=istartset(i),iendset(i)
                  nodeold=ialset(inode)
                  if(iponoel2d(nodeold).eq.0) cycle
                  ielem=inoel2d(1,iponoel2d(nodeold))
!     
!     Determine element formulation
!     
                  if((lakon(ielem)(7:7).eq.'A').or.
     &                 (lakon(ielem)(7:7).eq.'E').or.
     &                 (lakon(ielem)(7:7).eq.'L').or.
     &                 (lakon(ielem)(7:7).eq.'S')) then
                    if(lakon(ielem)(4:5).eq.'20') then
                      nope2d=8
                      ishift=20
                    elseif(lakon(ielem)(4:5).eq.'8R') then
                      nope2d=4
                      ishift=8
                    elseif(lakon(ielem)(4:5).eq.'8 ') then
                      nope2d=4
                      ishift=8
                    elseif(lakon(ielem)(4:5).eq.'8I') then
                      nope2d=4
                      ishift=11
                    elseif(lakon(ielem)(4:5).eq.'15') then
                      nope2d=6
                      ishift=15
                    elseif(lakon(ielem)(4:4).eq.'6') then
                      nope2d=3
                      ishift=6
                    else
                      cycle
                    endif
!     
                    do k=1,nope2d
                      ipos1=ipkon(ielem)+ishift+k     
                      if(nodeold.eq.kon(ipos1)) then
                        ipos2=iponor2d(2,ipos1-ishift)
                        nodenew=knor2d(ipos2+1)
                        ialset(inode)=nodenew
                        exit
                      endif
                    enddo
                  endif
                enddo
              endif
            endif
          endif
        endif 
      enddo     
!     
!     Rename the designvariables and the save the 2 additional expanded
!     nodes in iponod2dto3d:
!     for designvariables i --> 1st neighbor=iponod2dto3d(1,i)
!     2nd neighbor=iponod2dto3d(2,i)                           
!     
      loop1: do inode=istartset(iset),iendset(iset)
      nodeold=ialset(inode)
      if(iponoel2d(nodeold).eq.0) cycle
      ielem=inoel2d(1,iponoel2d(nodeold))
!     
!     Determine element formulation
!     
      if((lakon(ielem)(7:7).eq.'A').or.
     &     (lakon(ielem)(7:7).eq.'E').or.
     &     (lakon(ielem)(7:7).eq.'L').or.
     &     (lakon(ielem)(7:7).eq.'S')) then
        if(lakon(ielem)(4:5).eq.'20') then
          nope2d=8
          ishift=20
        elseif(lakon(ielem)(4:5).eq.'8R') then
          nope2d=4
          ishift=8
        elseif(lakon(ielem)(4:5).eq.'8 ') then
          nope2d=4
          ishift=8
        elseif(lakon(ielem)(4:5).eq.'8I') then
          nope2d=4
          ishift=11
        elseif(lakon(ielem)(4:5).eq.'15') then
          nope2d=6
          ishift=15
        elseif(lakon(ielem)(4:4).eq.'6') then
          nope2d=3
          ishift=6
        else
          cycle
        endif
!     
        do k=1,nope2d
          ipos1=ipkon(ielem)+ishift+k    
          if(nodeold.eq.kon(ipos1)) then
            ipos2=iponor2d(2,ipos1-ishift)
            nodenew=knor2d(ipos2+1)
!     
!     check for the existence of a MPC in the node
!     
            do l=1,3
              if(nactdof(l,nodeold).ge.0) exit
!     check if its an MPC(odd) or SPC(even)
              num=nactdof(l,nodeold)
              numtest=num/2*2
              if(num.ne.numtest) then
                write(*,*) '*WARNING in getdesiinfo2d:'
                write(*,*) '       node ',node,' is a'
                write(*,*) '       dependent dof in a MPC and'
                write(*,*) '       is removed from the set'
                write(*,*) '       of design variables'
                write(40,*) node
                iwrite=1
                cycle loop1
              endif
            enddo
!     
!     check whether not all dofs are removed by SPCs
!     
            do l=1,3
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
            iponod2dto3d(1,nodenew)=knor2d(ipos2+2)
            iponod2dto3d(2,nodenew)=knor2d(ipos2+3)
            ndesi=ndesi+1
            nodedesi(ndesi)=nodenew
            exit
          endif
        enddo
      endif
      enddo loop1
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
!     creating field nodedesiinv indicating for each node whether
!     it is a design variable or not
!     
      do i=1,ndesi
        index1=nodedesi(i)
        nodedesiinv(index1)=1
        index1=iponod2dto3d(1,nodedesi(i))
        nodedesiinv(index1)=1
        index1=iponod2dto3d(2,nodedesi(i))
        nodedesiinv(index1)=1
      enddo
!     
      kflag=1
      call isortii(nodedesi,iaux,ndesi,kflag)
!     
!     save the corresponding midnode for every node i in nk
!     --> midnode=iponod2dto3d(i)                          
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
            ishift=20
          elseif(lakon(ielem)(4:5).eq.'8R') then
            nope2d=4
            ishift=8
          elseif(lakon(ielem)(4:5).eq.'8 ') then
            nope2d=4
            ishift=8
          elseif(lakon(ielem)(4:5).eq.'8I') then
            nope2d=4
            ishift=11
          elseif(lakon(ielem)(4:5).eq.'15') then
            nope2d=6
            ishift=15
          elseif(lakon(ielem)(4:4).eq.'6') then
            nope2d=3
            ishift=6
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
            iponk2dto3d(knor2d(ipos2+1))=knor2d(ipos2+2)
!**** 
            iponk2dto3d(knor2d(ipos2+2))=knor2d(ipos2+2)
!**** 
            iponk2dto3d(knor2d(ipos2+3))=knor2d(ipos2+2)
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
