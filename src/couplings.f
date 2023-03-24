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
      subroutine couplings(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,nboun,nk,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikboun,ikmpc,ilmpc,co,labmpc,istat,n,iline,ipol,
     &     inl,ipoinp,inp,ipoinpc,norien,orname,orab,irstrt,ipkon,
     &     kon,lakon,istep,ics,dcs,nk_,nboun_,nodeboun,ndirboun,
     &     typeboun,ilboun,ier,nfc,nfc_,coeffc,ikdc,ndc,ndc_,
     &     edc)
!     
!     reading the input deck: *COUPLING in combination with
!     *KINEMATIC or *DISTRIBUTING
!     
      implicit none
!
      logical cyclicsymmetry
!     
      character*1 inpc(*),surfkind,typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*80 orname(*),orientation
      character*81 set(*),surfset
      character*132 textpart(16),name
!     
      integer istartset(*),iendset(*),ialset(*),norien,irstrt(*),nface,
     &     iorientation,iface,jface,nset,nboun,istat,n,i,j,k,ibounstart,
     &     ibounend,key,nk,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &     ikboun(*),ikmpc(*),ilmpc(*),ipos,m,node,iline,ipol,inl,nope,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),jsurf,irefnode,indexe,
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),ipkon(*),
     &     kon(*),istep,nelem,ics(2,*),nodef(8),nboun_,nodeboun(*),
     &     npt,mint,mpcfreeold,id,idof,iflag,nk_,l,ndirboun(*),
     &     irotnode_kin,idupnode,ier,nopes,nfcstart,nfc,ikdc(*),ndc,
     &     nfc_,ndc_,ilboun(*),matz,ierrs
!     
      real*8 coefmpc(*),co(3,*),orab(7,*),dcs(*),areanodal(8),xl2(3,8),
     &     shp2(7,8),xsj2(3),xsj,xi,et,weight,xs2(3,2),area,sum1,sum2,
     &     pcg(3),cg(3),sum3,coeffc(0:6,*),dd,r(3),e1(3),e2(3),s(3,3),
     &     e3(3),rp1(3),rp2(3),rp3(3),edc(12,*),w(3),z(3,3),fv1(3),
     &     fv2(3),pcl(3)
!     
      include "gauss.f"
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
!     flag for shape functions
!     
      data iflag /2/
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) '*ERROR reading *COUPLING: *COUPLING'
        write(*,*)'       should be placed before all step definitions'
        ier=1
        return
      endif
!     
      label='                    '
      orientation='
     &'
      do i=1,81
        surfset(i:i)=' '
      enddo
!     
      name(1:1)=' '
      do i=2,n
        if(textpart(i)(1:8).eq.'REFNODE=') then
          read(textpart(i)(9:18),'(i10)',iostat=istat) irefnode
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*COUPLING%",ier)
            return
          endif
          if(irefnode.gt.nk) then
            write(*,*) '*ERROR reading *COUPLING: ref node',irefnode
            write(*,*) '       has not been defined'
            ier=1
            return
          endif
        else if(textpart(i)(1:8).eq.'SURFACE=') then
          surfset(1:80)=textpart(i)(9:88)
!     
          ipos=index(surfset,' ')
          surfkind='S'
          surfset(ipos:ipos)=surfkind
          call cident81(set,surfset,nset,id)
          j=nset+1
          if(id.gt.0) then
            if(surfset.eq.set(id)) then
              j=id
            endif
          endif
          if(j.gt.nset) then
            surfkind='T'
            surfset(ipos:ipos)=surfkind
            call cident81(set,surfset,nset,id)
            j=nset+1
            if(id.gt.0) then
              if(surfset.eq.set(id)) then
                j=id
              endif
            endif
            if(j.gt.nset) then
              write(*,*) '*ERROR reading *COUPLING:'
              write(*,*) '       surface ',surfset
              write(*,*) '       has not yet been defined.' 
              ier=1
              return
            endif
          endif
          jsurf=j
        elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
          orientation=textpart(i)(13:92)
        elseif(textpart(i)(1:15).eq.'CONSTRAINTNAME=') then
          name(1:117)=textpart(i)(16:132)
        else
          write(*,*) 
     &         '*WARNING reading *COUPLING: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*COUPLING%")
        endif
      enddo
!     
      if(name(1:1).eq.' ') then
        write(*,*)
     &       '*ERROR reading *COUPLING: no CONTRAINT NAME given'
        write(*,*) '  '
        call inputerror(inpc,ipoinpc,iline,
     &       "*COUPLING%",ier)
        return
      endif
!     
      if(orientation.eq.'                    ') then
        iorientation=0
      else
        do i=1,norien
          if(orname(i).eq.orientation) exit
        enddo
        if(i.gt.norien) then
          write(*,*)
     &         '*ERROR reading *COUPLING: nonexistent orientation'
          write(*,*) '  '
          call inputerror(inpc,ipoinpc,iline,
     &         "*COUPLING%",ier)
          return
        endif
        iorientation=i
      endif
!     
!     next keyword should be *KINEMATIC or *DISTRIBUTING
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if(textpart(1)(2:10).eq.'KINEMATIC') then
!     
!     catalogueing the nodes
!     
        npt=0
!     
        do j=istartset(jsurf),iendset(jsurf)
          if(ialset(j).gt.0) then
            if(surfkind.eq.'T') then
!     
!     facial surface
!     
              iface=ialset(j)
              nelem=int(iface/10)
              jface=iface-nelem*10
              indexe=ipkon(nelem)
!     
              if(lakon(nelem)(4:5).eq.'20') then
                nopes=8
              elseif(lakon(nelem)(4:4).eq.'2') then
                nopes=9
              elseif(lakon(nelem)(4:4).eq.'8') then
                nopes=4
              elseif(lakon(nelem)(4:5).eq.'10') then
                nopes=6
              elseif(lakon(nelem)(4:4).eq.'4') then
                nopes=3
              endif
!     
              if(lakon(nelem)(4:4).eq.'6') then
                if(jface.le.2) then
                  nopes=3
                else
                  nopes=4
                endif
              endif
              if(lakon(nelem)(4:5).eq.'15') then
                if(jface.le.2) then
                  nopes=6
                else
                  nopes=8
                endif
              endif   
            else
!     
!     nodal surface
!     
              nopes=1
            endif
!     
            do m=1,nopes
              if(surfkind.eq.'T') then
                if((lakon(nelem)(4:4).eq.'2').or.
     &               (lakon(nelem)(4:4).eq.'8')) then
                  node=kon(indexe+ifaceq(m,jface))
                elseif((lakon(nelem)(4:4).eq.'4').or.
     &                 (lakon(nelem)(4:5).eq.'10')) then
                  node=kon(indexe+ifacet(m,jface))
                elseif(lakon(nelem)(4:4).eq.'6') then
                  node=kon(indexe+ifacew1(m,jface))
                elseif(lakon(nelem)(4:5).eq.'15') then
                  node=kon(indexe+ifacew2(m,jface))
                endif
              else
                node =ialset(j)
              endif
!     
              call nident2(ics,node,npt,id)
              if(id.gt.0) then
                if(ics(1,id).eq.node) then
                  cycle
                endif
              endif
!     
!     updating ics
!     
              npt=npt+1
              do l=npt,id+2,-1
                ics(1,l)=ics(1,l-1)
              enddo
              ics(1,id+1)=node
            enddo
          else
!     
!     if a negative value occurs the surface has to be
!     nodal
!     
            k=ialset(j-2)
            do
              k=k-ialset(j)
              if(k.ge.ialset(j-1)) exit
              node=k
              call nident2(ics,node,npt,id)
              if(id.gt.0) then
                if(ics(1,id).eq.node) then
                  cycle
                endif
              endif
!     
!     updating ics
!     
              npt=npt+1
              do l=npt,id+2,-1
                ics(1,l)=ics(1,l-1)
              enddo
              ics(1,id+1)=node
            enddo
          endif
        enddo
!     
!     generating a rotational node and connecting the
!     rotational dofs of the reference node with the
!     translational dofs of the rotational node
!     
!     generating a rotational reference node
!     
        nk=nk+1
        if(nk.gt.nk_) then
          write(*,*) 
     &         '*ERROR reading *KINEMATIC: increase nk_'
          ier=1
          return
        endif
        irotnode_kin=nk
        do l=1,3
          co(l,nk)=co(l,irefnode)
        enddo
!     
!     generating connecting MPCs between the rotational
!     dofs of irefnode and the translational dofs of
!     irotnode_kin
!     
        do k=1,3
!     
          nmpc=nmpc+1
          if(nmpc.gt.nmpc_) then
            write(*,*) 
     &           '*ERROR reading *KINEMATIC: increase nmpc_'
            ier=1
            return
          endif
!     
!     the internal dofs for rotation are 4, 5 and 6
!     
          ipompc(nmpc)=mpcfree
          labmpc(nmpc)='ROTTRACOUPLING      '
          idof=8*(irefnode-1)+k+3
          call nident(ikmpc,idof,nmpc-1,id)
          do l=nmpc,id+2,-1
            ikmpc(l)=ikmpc(l-1)
            ilmpc(l)=ilmpc(l-1)
          enddo
          ikmpc(id+1)=idof
          ilmpc(id+1)=nmpc
!     
          nodempc(1,mpcfree)=irefnode
          nodempc(2,mpcfree)=k+4
          coefmpc(mpcfree)=1.d0
          mpcfree=nodempc(3,mpcfree)
!     
          nodempc(1,mpcfree)=irotnode_kin
          nodempc(2,mpcfree)=k
          coefmpc(mpcfree)=-1.d0
          mpcfreeold=mpcfree
          mpcfree=nodempc(3,mpcfree)
          nodempc(3,mpcfreeold)=0
        enddo
!     
        if(iorientation.gt.0) then
!     
!     duplicating the nodes
!     generating rigid body MPC's for all dofs in the new nodes
!     
          do m=1,npt
            nk=nk+1
            if(nk.gt.nk_) then
              write(*,*) 
     &             '*ERROR reading *KINEMATIC: increase nk_'
              ier=1
              return
            endif
            ics(2,m)=nk
            do k=1,3
              co(k,nk)=co(k,ics(1,m))
            enddo
            node=nk
            ibounstart=1
            ibounend=3
            call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &           irotnode_kin,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &           nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,node,
     &           typeboun,co,ibounstart,ibounend)
          enddo
        endif
!     
!     reading the degrees of freedom
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) return
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*KINEMATIC%",ier)
            return
          endif
          if(ibounstart.lt.1) then
            write(*,*) '*ERROR reading *KINEMATIC'
            write(*,*) '       starting degree of freedom cannot'
            write(*,*) '       be less than 1'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*KINEMATIC%",ier)
            return
          endif
!     
          if(textpart(2)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*BOUNDARY%",ier)
              return
            endif
          endif
          if(ibounend.gt.3) then
            write(*,*) '*ERROR reading *KINEMATIC'
            write(*,*) '       final degree of freedom cannot'
            write(*,*) '       exceed 3'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*KINEMATIC%",ier)
            return
          elseif(ibounend.lt.ibounstart) then
            write(*,*) '*ERROR reading *KINEMATIC'
            write(*,*) '       initial degree of freedom cannot'
            write(*,*) '       exceed final degree of freedom'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*KINEMATIC%",ier)
            return
          endif
!     
!     generating the MPCs
!     
          if(iorientation.eq.0) then
!     
!     generating rigid body MPC's for the appropriate dofs
!     
            do j=1,npt
              node=ics(1,j)
              call rigidmpc(ipompc,nodempc,coefmpc,irefnode,
     &             irotnode_kin,
     &             labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,
     &             nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &             node,typeboun,co,ibounstart,ibounend)
            enddo
          else
!     
!     connecting the original nodes with the duplicated nodes
!     for the appropriate dofs
!     
            do j=1,npt
              node=ics(1,j)
              idupnode=ics(2,j)
              call mpcadd(node,ibounstart,ibounend,nboun,ipompc,
     &             nodempc,coefmpc,nmpc,nmpc_,mpcfree,orab,ikboun,
     &             ikmpc,ilmpc,co,labmpc,label,idupnode,
     &             iorientation)
            enddo
          endif
        enddo
      elseif(textpart(1)(2:13).eq.'DISTRIBUTING') then
        if(surfkind.eq.'S') then
          write(*,*) '*ERROR reading *DISTRIBUTING'
          write(*,*) '       a nodal surface is not allowed'
          write(*,*) '       please use a facial surface on'
          write(*,*) '       the *COUPLING card'
          ier=1
          return
        endif
!     
!     check whether cyclic symmetric
!
        cyclicsymmetry=.false.
        do i=2,n
          if(textpart(i)(1:14).eq.'CYCLICSYMMETRY') then
            cyclicsymmetry=.true.
          else
            write(*,*) 
     &       '*WARNING reading *DISTRIBUTING: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &           "*DISTRIBUTING%")
          endif
        enddo
!     
        npt=0
        area=0.d0
!     
!     catalogueing the nodes belonging to the surface (ics(1,*))
!     catalogueing the area (dcs(*))
!     
        do k=istartset(jsurf),iendset(jsurf)
!     
!     facial surface
!     
          iface=ialset(k)
          nelem=int(iface/10)
          jface=iface-nelem*10
          indexe=ipkon(nelem)
!     
!     nodes: #nodes in the face
!     the nodes are stored in nodef(*)
!     
          if(lakon(nelem)(4:4).eq.'2') then
            nopes=8
            nface=6
          elseif(lakon(nelem)(3:4).eq.'D8') then
            nopes=4
            nface=6
          elseif(lakon(nelem)(4:5).eq.'10') then
            nopes=6
            nface=4
            nope=10
          elseif(lakon(nelem)(4:4).eq.'4') then
            nopes=3
            nface=4
            nope=4
          elseif(lakon(nelem)(4:5).eq.'15') then
            if(jface.le.2) then
              nopes=6
            else
              nopes=8
            endif
            nface=5
            nope=15
          elseif(lakon(nelem)(3:4).eq.'D6') then
            if(jface.le.2) then
              nopes=3
            else
              nopes=4
            endif
            nface=5
            nope=6
          else
            cycle
          endif
!     
!     determining the nodes of the face
!     
          if(nface.eq.4) then
            do i=1,nopes
              nodef(i)=kon(indexe+ifacet(i,jface))
            enddo
          elseif(nface.eq.5) then
            if(nope.eq.6) then
              do i=1,nopes
                nodef(i)=kon(indexe+ifacew1(i,jface))
              enddo
            elseif(nope.eq.15) then
              do i=1,nopes
                nodef(i)=kon(indexe+ifacew2(i,jface))
              enddo
            endif
          elseif(nface.eq.6) then
            do i=1,nopes
              nodef(i)=kon(indexe+ifaceq(i,jface))
            enddo
          endif
!     
!     loop over the nodes belonging to the face   
!     ics(1,*): surface node
!     dcs(*): area corresponding to this node   
!     
          do i=1,nopes
            node=nodef(i)
            call nident2(ics,node,npt,id)
            if(id.gt.0) then
              if(ics(1,id).eq.node) then
                cycle
              endif
            endif
!     
!     updating ics
!     
            npt=npt+1
            do j=npt,id+2,-1
              ics(1,j)=ics(1,j-1)
              dcs(j)=dcs(j-1)
            enddo
            ics(1,id+1)=node
            dcs(id+1)=0.d0
          enddo
!     
!     calculating the area of the face and its contributions
!     to the facial nodes
!     
!     number of integration points
!     
          if(lakon(nelem)(3:5).eq.'D8R') then
            mint=1
          elseif(lakon(nelem)(3:4).eq.'D8') then
            mint=4
          elseif(lakon(nelem)(4:6).eq.'20R') then
            mint=4
          elseif(lakon(nelem)(4:4).eq.'2') then
            mint=9
          elseif(lakon(nelem)(4:5).eq.'10') then
            mint=3
          elseif(lakon(nelem)(4:4).eq.'4') then
            mint=1
          elseif(lakon(nelem)(3:4).eq.'D6') then
            mint=1
          elseif(lakon(nelem)(4:5).eq.'15') then
            if(jface.le.2) then
              mint=3
            else
              mint=4
            endif
          endif
!     
          do i=1,nopes
            areanodal(i)=0.d0
            do j=1,3
              xl2(j,i)=co(j,nodef(i))
            enddo
          enddo
!     
          do m=1,mint
            if((lakon(nelem)(3:5).eq.'D8R').or.
     &           ((lakon(nelem)(3:4).eq.'D6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,m)
              et=gauss2d1(2,m)
              weight=weight2d1(m)
            elseif((lakon(nelem)(3:4).eq.'D8').or.
     &             (lakon(nelem)(4:6).eq.'20R').or.
     &             ((lakon(nelem)(4:5).eq.'15').and.
     &             (nopes.eq.8))) then
              xi=gauss2d2(1,m)
              et=gauss2d2(2,m)
              weight=weight2d2(m)
            elseif(lakon(nelem)(4:4).eq.'2') then
              xi=gauss2d3(1,m)
              et=gauss2d3(2,m)
              weight=weight2d3(m)
            elseif((lakon(nelem)(4:5).eq.'10').or.
     &             ((lakon(nelem)(4:5).eq.'15').and.
     &             (nopes.eq.6))) then
              xi=gauss2d5(1,m)
              et=gauss2d5(2,m)
              weight=weight2d5(m)
            elseif((lakon(nelem)(4:4).eq.'4').or.
     &             ((lakon(nelem)(3:4).eq.'D6').and.
     &             (nopes.eq.3))) then
              xi=gauss2d4(1,m)
              et=gauss2d4(2,m)
              weight=weight2d4(m)
            endif
!     
            if(nopes.eq.8) then
              call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
              call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
              call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.3) then
              call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
!     
!     calculating the total area and nodal area
!
            xsj=weight*dsqrt(xsj2(1)**2+xsj2(2)**2+xsj2(3)**2)
            area=area+xsj
            do i=1,nopes
              areanodal(i)=areanodal(i)+xsj*shp2(4,i)
            enddo
!     
          enddo
!     
!     inserting the nodal area into field dcs
!     
          do i=1,nopes
            node=nodef(i)
            call nident2(ics,node,npt,id)
            dcs(id)=dcs(id)+areanodal(i)
          enddo
!     
        enddo
!     
!     storing the weight in field dcs
!     
        do i=1,npt
          dcs(i)=dcs(i)/area
        enddo
!     
!     determining the weighted center of gravity (only if not
!     cyclic symmetric)
!
        if(.not.cyclicsymmetry) then
          do j=1,3
            cg(j)=0.d0
            do i=1,npt
              cg(j)=cg(j)+co(j,ics(1,i))*dcs(i)
            enddo
          enddo
        else
          do j=1,3
            cg(j)=co(j,irefnode)
          enddo
        endif
!
!     determine unit vectors along the principal axes of inertia
!
!     setting up the matrix with the moments of inertia
!
        do j=1,3
          do i=1,j
            s(i,j)=0.d0
          enddo
        enddo
        do i=1,npt
          do j=1,3
            r(j)=co(j,ics(1,i))-cg(j)
          enddo
          do k=1,3
            do j=1,k
              s(j,k)=s(j,k)+r(j)*r(k)*dcs(i)
            enddo
          enddo
        enddo
        s(2,1)=s(1,2)
        s(3,1)=s(1,3)
        s(3,2)=s(2,3)
c        write(*,*) 'couplings s ',s(1,1),s(1,2),s(1,3)
c        write(*,*) 'couplings s ',s(2,1),s(2,2),s(2,3)
c        write(*,*) 'couplings s ',s(3,1),s(3,2),s(3,3)
!
!       calculating the eigenvalues and eigenvectors
!
        n=3
        matz=1
        ierrs=0
        call rs(n,n,s,w,matz,z,fv1,fv2,ierrs)
        if(ierrs.ne.0) then
          write(*,*) '*ERROR in couplings while calculating the'
          write(*,*) '       principal moments of inertia'
          ier=1
          return
        endif
!
        do i=1,3
          e1(i)=z(i,3)
          e2(i)=z(i,2)
        enddo
        e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
        e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
        e3(3)=e1(1)*e2(2)-e1(2)*e2(1)
c        write(*,*) 'couplings e1 ',e1(1),e1(2),e1(3)
c        write(*,*) 'couplings e2 ',e2(1),e2(2),e2(3)
c        write(*,*) 'couplings e3 ',e3(1),e3(2),e3(3)
c        write(*,*) e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)
!
!     determine sum_i ||r_i'||**2*w_i
!     
        sum1=0.d0
        sum2=0.d0
        sum3=0.d0
        do i=1,npt
          do j=1,3
            r(j)=co(j,ics(1,i))-cg(j)
          enddo
!
!         rp1 is the projection of r on a plane orthogonal to e1 
!
          dd=e1(1)*r(1)+e1(2)*r(2)+e1(3)*r(3)
          do j=1,3
            rp1(j)=r(j)-dd*e1(j)
          enddo
          sum1=sum1+(rp1(1)**2+rp1(2)**2+rp1(3)**2)*dcs(i)
!
!         rp2 is the projection of r on a plane orthogonal to e2 
!
          dd=e2(1)*r(1)+e2(2)*r(2)+e2(3)*r(3)
          do j=1,3
            rp2(j)=r(j)-dd*e2(j)
          enddo
          sum2=sum2+(rp2(1)**2+rp2(2)**2+rp2(3)**2)*dcs(i)
!
!         rp3 is the projection of r on a plane orthogonal to e3 
!
          dd=e3(1)*r(1)+e3(2)*r(2)+e3(3)*r(3)
          do j=1,3
            rp3(j)=r(j)-dd*e3(j)
          enddo
          sum3=sum3+(rp3(1)**2+rp3(2)**2+rp3(3)**2)*dcs(i)
        enddo
!     
!       determine the distance vector between irefnode and cg
!     
        do i=1,3
          pcg(i)=co(i,irefnode)-cg(i)
        enddo
!
!       pcl contains the components of the distance vector in
!       the local coupling surface coordinate system
!
        pcl(1)=pcg(1)*e1(1)+pcg(2)*e1(2)+pcg(3)*e1(3)
        pcl(2)=pcg(1)*e2(1)+pcg(2)*e2(2)+pcg(3)*e2(3)
        pcl(3)=pcg(1)*e3(1)+pcg(2)*e3(2)+pcg(3)*e3(3)
!
!     start of the force constraints for this coupling+distributed
!     definition
!
        nfcstart=nfc+1
!     
!     generating 3 equations per node defining the transfer
!     of the distributing force and moment to this node
!
        do i=1,npt
!
!         relative position of the node
!
          do j=1,3
            r(j)=co(j,ics(1,i))-cg(j)
          enddo
!
!         rp1 is the projection of r on a plane orthogonal to e1 
!
          dd=e1(1)*r(1)+e1(2)*r(2)+e1(3)*r(3)
          do j=1,3
            rp1(j)=r(j)-dd*e1(j)
          enddo
!
!         rp2 is the projection of r on a plane orthogonal to e2 
!
          dd=e2(1)*r(1)+e2(2)*r(2)+e2(3)*r(3)
          do j=1,3
            rp2(j)=r(j)-dd*e2(j)
          enddo
!
!         rp3 is the projection of r on a plane orthogonal to e3 
!
          dd=e3(1)*r(1)+e3(2)*r(2)+e3(3)*r(3)
          do j=1,3
            rp3(j)=r(j)-dd*e3(j)
          enddo
!
!         equation in global x
!
          nfc=nfc+1
          if(nfc.gt.nfc_) then
            write(*,*) '*ERROR reading *DISTRIBUTING: increase nfc_'
            call exit(201)
          endif
!
          coeffc(0,nfc)=10*ics(1,i)+1.5
          coeffc(1,nfc)=e1(1)*dcs(i)
          coeffc(2,nfc)=e2(1)*dcs(i)
          coeffc(3,nfc)=e3(1)*dcs(i)
          coeffc(4,nfc)=(e1(2)*rp1(3)-e1(3)*rp1(2))*dcs(i)/sum1
          coeffc(5,nfc)=(e2(2)*rp2(3)-e2(3)*rp2(2))*dcs(i)/sum2
          coeffc(6,nfc)=(e3(2)*rp3(3)-e3(3)*rp3(2))*dcs(i)/sum3
!
!         correction for force not in center of gravity
!
          coeffc(1,nfc)=coeffc(1,nfc)+coeffc(5,nfc)*pcl(3)
     &                               -coeffc(6,nfc)*pcl(2)
          coeffc(2,nfc)=coeffc(2,nfc)+coeffc(6,nfc)*pcl(1)
     &                               -coeffc(4,nfc)*pcl(3)
          coeffc(3,nfc)=coeffc(3,nfc)+coeffc(4,nfc)*pcl(2)
     &                               -coeffc(5,nfc)*pcl(1)
!
!         equation in global y
!
          nfc=nfc+1
          if(nfc.gt.nfc_) then
            write(*,*) '*ERROR reading *DISTRIBUTING: increase nfc_'
            call exit(201)
          endif
!
          coeffc(0,nfc)=10*ics(1,i)+2.5
          coeffc(1,nfc)=e1(2)*dcs(i)
          coeffc(2,nfc)=e2(2)*dcs(i)
          coeffc(3,nfc)=e3(2)*dcs(i)
          coeffc(4,nfc)=(e1(3)*rp1(1)-e1(1)*rp1(3))*dcs(i)/sum1
          coeffc(5,nfc)=(e2(3)*rp2(1)-e2(1)*rp2(3))*dcs(i)/sum2
          coeffc(6,nfc)=(e3(3)*rp3(1)-e3(1)*rp3(3))*dcs(i)/sum3
!
!         correction for force not in center of gravity
!
          coeffc(1,nfc)=coeffc(1,nfc)+coeffc(5,nfc)*pcl(3)
     &                               -coeffc(6,nfc)*pcl(2)
          coeffc(2,nfc)=coeffc(2,nfc)+coeffc(6,nfc)*pcl(1)
     &                               -coeffc(4,nfc)*pcl(3)
          coeffc(3,nfc)=coeffc(3,nfc)+coeffc(4,nfc)*pcl(2)
     &                               -coeffc(5,nfc)*pcl(1)
!
!         equation in global z
!
          nfc=nfc+1
          if(nfc.gt.nfc_) then
            write(*,*) '*ERROR reading *DISTRIBUTING: increase nfc_'
            call exit(201)
          endif
!
          coeffc(0,nfc)=10*ics(1,i)+3.5
          coeffc(1,nfc)=e1(3)*dcs(i)
          coeffc(2,nfc)=e2(3)*dcs(i)
          coeffc(3,nfc)=e3(3)*dcs(i)
          coeffc(4,nfc)=(e1(1)*rp1(2)-e1(2)*rp1(1))*dcs(i)/sum1
          coeffc(5,nfc)=(e2(1)*rp2(2)-e2(2)*rp2(1))*dcs(i)/sum2
          coeffc(6,nfc)=(e3(1)*rp3(2)-e3(2)*rp3(1))*dcs(i)/sum3
!
!         correction for force not in center of gravity
!
          coeffc(1,nfc)=coeffc(1,nfc)+coeffc(5,nfc)*pcl(3)
     &                               -coeffc(6,nfc)*pcl(2)
          coeffc(2,nfc)=coeffc(2,nfc)+coeffc(6,nfc)*pcl(1)
     &                               -coeffc(4,nfc)*pcl(3)
          coeffc(3,nfc)=coeffc(3,nfc)+coeffc(4,nfc)*pcl(2)
     &                               -coeffc(5,nfc)*pcl(1)
        enddo
!
!       treating the translational degrees of freedom (default)
!
        do i=1,3
          idof=8*(irefnode-1)+i
          call nident(ikdc,idof,ndc,id)
          if(id.gt.0) then
            if(ikdc(id).eq.idof) cycle
          endif
          ndc=ndc+1
          if(ndc.gt.ndc_) then
            write(*,*) '*ERROR reading *DISTRIBUTING: increase ndc_'
            call exit(201)
          endif
          do j=ndc,id+2,-1
            ikdc(j)=ikdc(j-1)
            do k=1,12
              edc(k,j)=edc(k,j-1)
            enddo
          enddo
          ikdc(id+1)=idof
          edc(1,id+1)=nfcstart+0.5d0
          edc(2,id+1)=nfc+0.5d0
          edc(3,id+1)=iorientation+0.5d0
          edc(4,id+1)=e1(1)
          edc(5,id+1)=e1(2)
          edc(6,id+1)=e1(3)
          edc(7,id+1)=e2(1)
          edc(8,id+1)=e2(2)
          edc(9,id+1)=e2(3)
          edc(10,id+1)=e3(1)
          edc(11,id+1)=e3(2)
          edc(12,id+1)=e3(3)
        enddo
!
!     reading the degrees of freedom
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) return
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat) ibounstart
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DISTRIBUTING%",ier)
            return
          endif
          if(ibounstart.lt.1) then
            write(*,*) '*ERROR reading *DISTRIBUTING'
            write(*,*) '       starting degree of freedom cannot'
            write(*,*) '       be less than 1'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DISTRIBUTING%",ier)
            return
          endif
!     
          if(textpart(2)(1:1).eq.' ') then
            ibounend=ibounstart
          else
            read(textpart(2)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DISTRIBUTING%",ier)
              return
            endif
          endif
!     
          ibounstart=max(4,ibounstart)
          if(ibounend.gt.6) then
            write(*,*) '*ERROR reading *DISTRIBUTING'
            write(*,*) '       final degree of freedom cannot'
            write(*,*) '       exceed 6'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DISTRIBUTING%",ier)
            return
          elseif(ibounend.lt.ibounstart) then
            cycle
          endif
!     
!     treating the selected rotational degrees of freedom
!     
          do i=ibounstart,ibounend
            idof=8*(irefnode-1)+i
            call nident(ikdc,idof,ndc,id)
            if(id.gt.0) then
              if(ikdc(id).eq.idof) cycle
            endif
            ndc=ndc+1
            if(ndc.gt.ndc_) then
              write(*,*) '*ERROR reading *DISTRIBUTING: increase ndc_'
              call exit(201)
            endif
            do j=ndc,id+2,-1
              ikdc(j)=ikdc(j-1)
              do k=1,12
                edc(k,j)=edc(k,j-1)
              enddo
            enddo
            ikdc(id+1)=idof
            edc(1,id+1)=nfcstart+0.5d0
            edc(2,id+1)=nfc+0.5d0
            edc(3,id+1)=iorientation+0.5d0
            edc(4,id+1)=e1(1)
            edc(5,id+1)=e1(2)
            edc(6,id+1)=e1(3)
            edc(7,id+1)=e2(1)
            edc(8,id+1)=e2(2)
            edc(9,id+1)=e2(3)
            edc(10,id+1)=e3(1)
            edc(11,id+1)=e3(2)
            edc(12,id+1)=e3(3)
          enddo
        enddo
      else
        write(*,*)
     &       '*ERROR reading *COUPLING: the line following'
        write(*,*) '       *COUPLING must contain the'
        write(*,*) '       *KINEMATIC keyword or the'
        write(*,*) '       *DISTRIBUTING keyword'
        write(*,*) '  '
        call inputerror(inpc,ipoinpc,iline,
     &       "*COUPLING%",ier)
        ier=2
        return
      endif
!     
      return
      end

