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
      subroutine rigidbodys(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,nset_,nalset,nalset_,ipompc,nodempc,coefmpc,
     &     labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,lakon,ipkon,kon,nk,nk_,
     &     nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,iperturb,ne_,
     &     ctrl,typeboun,istep,istat,n,iline,ipol,inl,ipoinp,inp,co,
     &     ipoinpc,ier)
!     
!     reading the input deck: *RIGID BODY
!     
      implicit none
!     
      character*1 typeboun(*),inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 set(*),elset,noset
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),ipompc(*),
     &     nodempc(3,*),
     &     nset,nset_,nalset,nalset_,nmpc,nmpc_,mpcfree,nk,nk_,ikmpc(*),
     &     ilmpc(*),ipkon(*),kon(*),inoset,ielset,i,node,ielement,id,
     &     indexe,nope,istep,istat,n,irefnode,irotnode,ne_,
     &     j,idof,k,nodeboun(*),ndirboun(*),ikboun(*),ilboun(*),
     &     nboun,nboun_,key,iperturb(*),ipos,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),ipoinpc(0:*),jmin,jmax,ier
!     
      real*8 coefmpc(3,*),ctrl(*),co(3,*)
!     
      jmin=1
      jmax=3
!     
      if(istep.gt.0) then
        write(*,*) 
     &       '*ERROR reading *RIGID BODY: *RIGID BODY should be placed'
        write(*,*) '  before all step definitions'
        ier=1
        return
      endif
!     
!     the *RIGID BODY option implies a nonlinear geometric 
!     calculation
!     
      if(iperturb(1).eq.1) then
        write(*,*) '*ERROR reading *RIGID BODY: the *RIGID BODY option'
        write(*,*) '       cannot be used in a perturbation step'
        ier=1
        return
      endif
!     
      elset='
     &'
      noset='
     &'
      irefnode=0
      irotnode=0
!     
      do i=2,n
        if(textpart(i)(1:6).eq.'ELSET=') then
          if(noset(1:1).eq.' ') then
            elset(1:80)=textpart(i)(7:86)
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
          else
            write(*,*) '*ERROR reading *RIGID BODY: either NSET or'
            write(*,*) '       ELSET can be specified, not both'
            ier=1
            return
          endif
        elseif(textpart(i)(1:8).eq.'PINNSET=') then
          if(elset(1:1).eq.' ') then
            noset(1:80)=textpart(i)(9:88)
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
          else
            write(*,*) '*ERROR reading *RIGID BODY: either NSET or'
            write(*,*) '       ELSET can be specified, not both'
            ier=1
            return
          endif
        elseif(textpart(i)(1:5).eq.'NSET=') then
          if(elset(1:1).eq.' ') then
            noset(1:80)=textpart(i)(6:85)
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
          else
            write(*,*) '*ERROR reading *RIGID BODY: either NSET or'
            write(*,*) '       ELSET can be specified, not both'
            ier=1
            return
          endif
        elseif(textpart(i)(1:8).eq.'REFNODE=') then
          read(textpart(i)(9:18),'(i10)',iostat=istat) irefnode
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RIGID BODY%",ier)
            return
          endif
          if(irefnode.gt.nk) then
            write(*,*) '*ERROR reading *RIGID BODY: ref node',
     &           irefnode
            write(*,*) '       has not been defined'
            ier=1
            return
          endif
        elseif(textpart(i)(1:8).eq.'ROTNODE=') then
          read(textpart(i)(9:18),'(i10)',iostat=istat) irotnode
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RIGID BODY%",ier)
            return
          endif
          if(irotnode.gt.nk) then
            write(*,*) '*ERROR reading *RIGID BODY: rot node',
     &           irotnode
            write(*,*) '       has not been defined'
            ier=1
            return
          endif
        else
          write(*,*) 
     &         '*WARNING reading *RIGID BODY: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*RIGID BODY%")
        endif
      enddo
!     
!     check whether a set was defined
!     
      if((elset(1:1).eq.' ').and.
     &     (noset(1:1).eq.' ')) then
        write(*,*) '*WARNING reading *RIGID BODY: no set defined'
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        return
      endif
!     
      inoset=0
      ielset=0
!     
!     checking whether the set exists
!     
      if(noset(1:1).ne.' ') then
c     do i=1,nset
c     if(set(i).eq.noset) then
c     inoset=i
c     exit
c     endif
c     enddo
        call cident81(set,noset,nset,id)
        if(id.gt.0) then
          if(noset.eq.set(id)) then
            inoset=id
          endif
        endif
        if(inoset.eq.0) then
          write(*,*) '*WARNING reading *RIGID BODY: node set ',noset
          write(*,*) '         does not exist'
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          return
        endif
      endif
!     
      if(elset(1:1).ne.' ') then
c     do i=1,nset
c     if(set(i).eq.elset) then
c     ielset=i
c     exit
c     endif
c     enddo
        call cident81(set,elset,nset,id)
        if(id.gt.0) then
          if(elset.eq.set(id)) then
            ielset=id
          endif
        endif
        if(ielset.eq.0) then
          write(*,*) '*WARNING reading *RIGID BODY: element set ',
     &         elset
          write(*,*) '         does not exist'
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          return
        endif
      endif
!     
!     check for the existence of irefnode and irotnode; if none were
!     defined, new nodes are generated
!     
      if(irefnode.eq.0) then
        nk=nk+1
        if(nk.gt.nk_) then
          write(*,*) '*ERROR reading *RIGID BODY: increase nk_'
          ier=1
          return
        endif
        irefnode=nk
!     
!     default position of the reference node is the origin
!     
        co(1,nk)=0.d0
        co(2,nk)=0.d0
        co(3,nk)=0.d0
      endif
!     
      if(irotnode.eq.0) then
        nk=nk+1
        if(nk.gt.nk_) then
          write(*,*) '*ERROR reading *RIGID BODY: increase nk_'
          ier=1
          return
        endif
        irotnode=nk
      endif
!     
!     check whether other equations apply to the dependent nodes
!     
      if(inoset.ne.0) then
        do i=istartset(inoset),iendset(inoset)
          node=ialset(i)
          if(node.gt.nk_) then
            write(*,*) '*ERROR reading *RIGID BODY: node ',node
            write(*,*) '       belonging to set ',noset
            write(*,*) '       has not been defined'
            ier=1
            return
          endif
          if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
          do j=1,3
            idof=8*(node-1)+j
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
              if(ikmpc(id).eq.idof) then
                write(*,*) '*WARNING reading *RIGID BODY: dof ',j
                write(*,*) '         of node ',node,' belonging'
                write(*,*) '         to a rigid body is detected'
                write(*,*) '         on the dependent side of '
                write(*,*) '         another equation; no rigid'
                write(*,*) '         body constrained applied'
              endif
            endif
          enddo
        enddo
      endif
!     
      if(ielset.ne.0) then
        do i=istartset(ielset),iendset(ielset)
          ielement=ialset(i)
          if(ielement.gt.ne_) then
            write(*,*) '*ERROR reading *RIGID BODY: element ',
     &           ielement
            write(*,*) '       belonging to set ',elset
            write(*,*) '       has not been defined'
            ier=1
            return
          endif
          if(ipkon(ielement).lt.0) cycle
          indexe=ipkon(ielement)
          if(lakon(ielement)(4:4).eq.'2') then
            nope=20
          elseif(lakon(ielement)(4:4).eq.'8') then
            nope=8
          elseif(lakon(ielement)(4:5).eq.'10') then
            nope=10
          elseif(lakon(ielement)(4:4).eq.'4') then
            nope=4
          elseif(lakon(ielement)(4:5).eq.'15') then
            nope=15
          else
            nope=6
          endif
          do k=indexe+1,indexe+nope
            node=kon(k)
            if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
            do j=1,3
              idof=8*(node-1)+j
              call nident(ikmpc,idof,nmpc,id)
              if(id.gt.0) then
                if(ikmpc(id).eq.idof) then
                  write(*,*)'*WARNING reading *RIGID BODY: dof ',
     &                 j,'of node ',node,' belonging to a'
                  write(*,*)'         rigid body is detected on th
     &e dependent side of another'
                  write(*,*)'         equation; no rigid body cons
     &trained applied'
                endif
              endif
            enddo
          enddo
        enddo
      endif
!     
!     generating the equations in basis form
!     
!     node set
!     
      if(inoset.ne.0) then
        do i=istartset(inoset),iendset(inoset)
          node=ialset(i)
          if(node.gt.0) then
            if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
            call rigidmpc(ipompc,nodempc,coefmpc,irefnode,irotnode,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &           nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,node,
     &           typeboun,co,jmin,jmax)
          else
            node=ialset(i-2)
            do
              node=node-ialset(i)
              if(node.ge.ialset(i-1)) exit
              if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
              call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &             irotnode,labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &             nk,nk_,nodeboun,ndirboun,ikboun,ilboun,nboun,
     &             nboun_,node,typeboun,co,jmin,jmax)
            enddo
          endif
        enddo
      endif
!     
!     element set
!     
      if(ielset.ne.0) then
        do i=istartset(ielset),iendset(ielset)
          ielement=ialset(i)
          if(ielement.gt.0) then
            if(ipkon(ielement).lt.0) cycle
            indexe=ipkon(ielement)
            if(lakon(ielement)(4:4).eq.'2') then
              nope=20
            elseif(lakon(ielement)(4:4).eq.'8') then
              nope=8
            elseif(lakon(ielement)(4:5).eq.'10') then
              nope=10
            elseif(lakon(ielement)(4:4).eq.'4') then
              nope=4
            elseif(lakon(ielement)(4:5).eq.'15') then
              nope=15
            else
              nope=6
            endif
            do k=indexe+1,indexe+nope
              node=kon(k)
              if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
              call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &             irotnode,labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &             nk,nk_,nodeboun,ndirboun,ikboun,ilboun,nboun,
     &             nboun_,node,typeboun,co,jmin,jmax)
            enddo
          else
            ielement=ialset(i-2)
            do
              ielement=ielement-ialset(i)
              if(ielement.ge.ialset(i-1)) exit
              if(ipkon(ielement).lt.0) cycle
              indexe=ipkon(ielement)
              if(lakon(ielement)(4:4).eq.'2') then
                nope=20
              elseif(lakon(ielement)(4:4).eq.'8') then
                nope=8
              elseif(lakon(ielement)(4:5).eq.'10') then
                nope=10
              elseif(lakon(ielement)(4:4).eq.'4') then
                nope=4
              elseif(lakon(ielement)(4:5).eq.'15') then
                nope=15
              else
                nope=6
              endif
              do k=indexe+1,indexe+nope
                node=kon(k)
                if((node.eq.irefnode).or.(node.eq.irotnode)) cycle
                call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &               irotnode,labmpc,nmpc,nmpc_,mpcfree,ikmpc,
     &               ilmpc,nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &               nboun,nboun_,node,typeboun,co,jmin,jmax)
              enddo
            enddo
          endif
        enddo
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

