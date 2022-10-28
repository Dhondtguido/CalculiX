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
      subroutine envtemp(itg,ieg,ntg,ntr,sideload,nelemload,
     &     ipkon,kon,lakon,ielmat,ne,nload,kontri,ntri,nloadtr,
     &     nflow,ndirboun,nactdog,nodeboun,nacteq,nboun,
     &     ielprop,prop,nteq,v,network,physcon,shcon,ntmat_,
     &     co,vold,set,nshcon,rhcon,nrhcon,mi,
     &     nmpc,nodempc,ipompc,labmpc,ikboun,nasym,ttime,time,
     &     iaxial)
!     
!     determines the number of gas temperatures and radiation
!     temperatures
!     
      implicit none
!     
      logical identity,walltemp,temperaturebc,pressurebc,massflowbcall,
     &     pressurebcall,channel
!     
      character*8 lakon(*)
      character*20 labmpc(*)
      character*20 sideload(*)
      character*81 set(*)
!     
      integer itg(*),ntg,ntr,nelemload(2,*),ipkon(*),network,mi(*),
     &     kon(*),ielmat(mi(3),*),ne,i,j,k,l,index,id,node,nload,
     &     ifaceq(8,6),ider,nasym,indexe,iaxial,networkmpcs,
     &     ifacet(6,4),ifacew(8,5),kontri3(3,1),kontri4(3,2),
     &     kontri6(3,4),kontri8(3,6),kontri(4,*),ntri,
     &     konf(8),nloadtr(*),nelem,nope,nopes,ig,nflow,ieg(*),
     &     ndirboun(*),nactdog(0:3,*),nboun,nodeboun(*),ntmat_,
     &     idir,ntq,nteq,nacteq(0:3,*),node1,node2,nodem,
     &     ielprop(*),idirf(8),iflag,imat,numf,nrhcon(*),nshcon(*),
     &     nmpc,nodempc(3,*),ipompc(*),ikboun(*),iplausi
!     
      real*8 prop(*),f,xflow,nodef(8),df(8),v(0:mi(2),*),g(3),
     &     cp,r,physcon(*),shcon(0:3,ntmat_,*),rho,ttime,time,
     &     co(3,*),dvi,vold(0:mi(2),*),rhcon(*)
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data kontri3 /1,2,3/
      data kontri4 /1,2,4,2,3,4/
      data kontri6 /1,4,6,4,5,6,4,2,5,6,5,3/
      data kontri8 /1,5,8,8,5,7,8,7,4,5,2,6,5,6,7,7,6,3/
!     
      ntg=0
      ntr=0
      ntri=0
!
      channel=.false.
      walltemp=.false.
      temperaturebc=.false.
      pressurebc=.false.
      massflowbcall=.true.
      pressurebcall=.true.
!     
      networkmpcs=0
!     
!     1) catalogueing nodes used as boundary condition for force
!        convection and network pressure 
!     2) triangulating the cavity radiation faces
!     
      do i=1,nload
        if(sideload(i)(3:4).eq.'FC') then
!
!         boundary condition node for forced convection
!
          walltemp=.true.
          call nident(itg,nelemload(2,i),ntg,id)
          if(id.gt.0) then
            if(itg(id).eq.nelemload(2,i)) then
              nactdog(0,nelemload(2,i))=1
              cycle
            endif
          endif
          ntg=ntg+1
          do j=ntg,id+2,-1
            itg(j)=itg(j-1)
          enddo
          itg(id+1)=nelemload(2,i)
          nactdog(0,nelemload(2,i))=1
!     
        elseif(sideload(i)(3:4).eq.'NP') then
!
!         boundary condition node for network pressure
!
          call nident(itg,nelemload(2,i),ntg,id)
          if(id.gt.0) then
            if(itg(id).eq.nelemload(2,i)) then
              nactdog(2,nelemload(2,i))=1
              cycle
            endif
          endif
          ntg=ntg+1
          do j=ntg,id+2,-1
            itg(j)=itg(j-1)
          enddo
          itg(id+1)=nelemload(2,i)
          nactdog(2,nelemload(2,i))=1
!     
        elseif(sideload(i)(3:4).eq.'CR') then
!
!         cavity radiation faces
!
          ntr=ntr+1
          nelem=nelemload(1,i)
          read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes in the face
!     
          if(lakon(nelem)(4:4).eq.'2') then
            nope=20
            nopes=8
          elseif(lakon(nelem)(4:4).eq.'8') then
            nope=8
            nopes=4
          elseif(lakon(nelem)(4:5).eq.'10') then
            nope=10
            nopes=6
          elseif(lakon(nelem)(4:4).eq.'4') then
            nope=4
            nopes=3
          elseif(lakon(nelem)(4:4).eq.'6') then
            nope=6
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          elseif(lakon(nelem)(4:5).eq.'15') then
            nope=15
            if(ig.le.2) then
              nopes=6
            else
              nopes=8
            endif
          endif
!     
!     nodes in the face
!     
          if((nope.eq.20).or.(nope.eq.8)) then
            do k=1,nopes
              konf(k)=kon(ipkon(nelem)+ifaceq(k,ig))
            enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
            do k=1,nopes
              konf(k)=kon(ipkon(nelem)+ifacet(k,ig))
            enddo
          else
            do k=1,nopes
              konf(k)=kon(ipkon(nelem)+ifacew(k,ig))
            enddo
          endif
!     
!     triangulation of the face
!     
          nloadtr(ntr)=i
          if((lakon(nelem)(4:4).eq.'2').or.
     &         ((lakon(nelem)(4:5).eq.'15').and.(ig.gt.2))) then
            do k=1,6
              ntri=ntri+1
              do l=1,3
                kontri(l,ntri)=konf(kontri8(l,k))
              enddo
              kontri(4,ntri)=ntr
            enddo
          elseif((lakon(nelem)(4:4).eq.'8').or.
     &           ((lakon(nelem)(4:4).eq.'6').and.(ig.gt.2))) then
            do k=1,2
              ntri=ntri+1
              do l=1,3
                kontri(l,ntri)=konf(kontri4(l,k))
              enddo
              kontri(4,ntri)=ntr
            enddo
          elseif((lakon(nelem)(4:5).eq.'10').or.
     &           ((lakon(nelem)(4:5).eq.'15').and.(ig.le.2))) then
            do k=1,4
              ntri=ntri+1
              do l=1,3
                kontri(l,ntri)=konf(kontri6(l,k))
              enddo
              kontri(4,ntri)=ntr
            enddo
          elseif((lakon(nelem)(4:4).eq.'4').or.
     &           ((lakon(nelem)(4:4).eq.'6').and.(ig.le.2))) then
            do k=1,1
              ntri=ntri+1
              do l=1,3
                kontri(l,ntri)=konf(kontri3(l,k))
              enddo
              kontri(4,ntri)=ntr
            enddo
          endif   
        endif
      enddo
!     
!     storing the gas elements in ieg
!     
      nflow=0
!     
      do i=1,ne
        if(lakon(i)(1:1).eq.'D') then
          if((lakon(i)(2:2).ne.' ').and.(network.ne.1)) then
            nflow=nflow+1
            ieg(nflow)=i
          else
            nasym=1
!     
!     removing gas nodes from the itg vector
!     belonging to 'D '-elements
!     in which a 'FC'-film condition was applied
!     
            indexe=ipkon(i)
            do j=1,3,2
              node=kon(indexe+j)
              call nident(itg,node,ntg,id)
              if(id.gt.0) then
                if(itg(id).eq.node) then
                  ntg=ntg-1
                  do k=id,ntg
                    itg(k)=itg(k+1)
                  enddo
                  nactdog(0,node)=0
                  nactdog(2,node)=0
                endif
              endif
            enddo
!     
          endif
        endif
!     
!     removing gas nodes belonging to advective elements
!     (last node in the element topology)
!     these are usually gas nodes not belonging to any
!     "D"-element
!     
        if(lakon(i)(1:7).eq.'ESPRNGF') then
          read(lakon(i)(8:8),'(i1)') nopes
          nope=nopes+1
          node=kon(ipkon(i)+nope)
!     
          call nident(itg,node,ntg,id)
          if(id.gt.0) then
            if(itg(id).eq.node) then
              ntg=ntg-1
              do k=id,ntg
                itg(k)=itg(k+1)
              enddo
              nactdog(0,node)=0
              nactdog(2,node)=0
            endif
          endif
        endif
      enddo
!     
!     mass flux nodes are also taken as unknowns in the
!     gas temperature system; determining the active 
!     degrees of freedom
!     
!     first node of the flow element
!     
      do i=1,nflow
!
        if(lakon(ieg(i))(2:5).eq.'LICH') channel=.true.
!        
        index=ipkon(ieg(i))
        node=kon(index+1)
        if (node.eq.0) cycle
        call nident(itg,node,ntg,id)
        if(id.gt.0) then
          if(itg(id).eq.node) then
            nactdog(0,node)=1
            nactdog(2,node)=1
            cycle 
          endif
        endif
        ntg=ntg+1
        do j=ntg,id+2,-1
          itg(j)=itg(j-1)
        enddo
        itg(id+1)=node
        nactdog(0,node)=1
        nactdog(2,node)=1
      enddo
!     
!     middle node of the flow element :flux
!     
      do i=1,nflow
        index=ipkon(ieg(i))
        node=kon(index+2)
        call nident(itg,node,ntg,id)
        if(id.gt.0) then
          if(itg(id).eq.node) cycle 
        endif
        ntg=ntg+1
        do j=ntg,id+2,-1
          itg(j)=itg(j-1)
        enddo
        itg(id+1)=node
        nactdog(1,node)=1
!     
!     variable geometric property
!     
        if(lakon(ieg(i))(6:7).eq.'GV') then
          index=ielprop(ieg(i))
          if(prop(index+2).le.0.d0) nactdog(3,node)=1
        elseif(lakon(ieg(i))(2:7).eq.'ACCTUB') then
!     
          index=ielprop(ieg(i))
          if(prop(index+1).eq.2) then
!     Interval factor unknown,setting a DOF for geometry
            nactdog(3,node)=1
          elseif(prop(index+1).eq.3) then
!     Hole diameter factor unknown,setting a DOF for geometry
            nactdog(3,node)=1
          endif
        endif
      enddo
!     
!     third node of the flow element
!     
      do i=1,nflow
        index=ipkon(ieg(i))
        node=kon(index+3)
        if (node.eq.0) cycle
        call nident(itg,node,ntg,id)
        if(id.gt.0) then
          if(itg(id).eq.node) then
            nactdog(0,node)=1
            nactdog(2,node)=1
            cycle
          endif
        endif
        ntg=ntg+1
        do j=ntg,id+2,-1
          itg(j)=itg(j-1)
        enddo
        itg(id+1)=node
        nactdog(0,node)=1
        nactdog(2,node)=1
      enddo
!
      if(channel) then
        nteq=0
        return
      endif
!     
!     tagging the network MPC's
!     
      do i=1,nmpc
!     
!     check whether network MPC
!     
        index=ipompc(i)
        do
          node=nodempc(1,index)
          call nident(itg,node,ntg,id)
          if(id.gt.0) then
            if(itg(id).eq.node) then
              labmpc(i)(1:7)='NETWORK'
              networkmpcs=1
              exit
            endif
          endif
          index=nodempc(3,index)
          if(index.eq.0) exit
        enddo
      enddo
!     
!     taking the MPC network nodes into account
!     
      if(networkmpcs.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:7).ne.'NETWORK') cycle
!     
          index=ipompc(i)
          do
            node=nodempc(1,index)
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                nactdog(nodempc(2,index),node)=1
                index=nodempc(3,index)
                if(index.eq.0) then
                  exit
                else
                  cycle
                endif
              endif
            endif
!     
!     adding a node to itg
!     
            ntg=ntg+1
            do j=ntg,id+2,-1
              itg(j)=itg(j-1)
            enddo
            itg(id+1)=node
            nactdog(nodempc(2,index),node)=1
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        enddo
      endif
!     
!     from this point on ntg,itg, nflow and ieg do not change
!     any more; check for homogeneous networks; a (simply or not)
!     connected network should be either
!     - a gas network  
!     - a liquid network  
!     - a channel network  
!     but no mix 
!     
!     subtracting the SPC conditions
!     
      do i=1,nboun
        node=nodeboun(i)
        call nident(itg,node,ntg,id)
        if (id.gt.0) then
          if (itg(id).eq.node) then
            idir=ndirboun(i)
            nactdog(idir,node)=0
            if(idir.eq.0) then
              temperaturebc=.true.
            elseif(idir.eq.2) then
              pressurebc=.true.
            endif
          endif
        endif
      enddo
!     
!     temporarily removing the dependent nodes of the MPC's
!     only for mass flow and pressure
!     
!     these are the only dofs for which the corresponding equations
!     (mass equilibrium in end node and element equation in middle
!     node) are not applied in the nodes where the dofs are lacking 
!     (mass flow in middle node and pressure in end node)
!     
      if(networkmpcs.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:7).ne.'NETWORK') cycle
          index=ipompc(i)
          idir=nodempc(2,index)
          if((idir.eq.1).or.(idir.eq.2)) then
            node=nodempc(1,index)
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                nactdog(idir,node)=0
              endif
            endif
          endif
        enddo
      endif
!     
!     determining the active equations
!     
!     element contributions
!     
      iflag=0
      do i=1,nflow
        nelem=ieg(i)
        index=ipkon(nelem)
        node1=kon(index+1)
        nodem=kon(index+2)
        node2=kon(index+3)
!     
!     "end of network" element ---X---O
!     1   m   2
!     
        if(node1.eq.0)then
          if ((nactdog(1,nodem).ne.0).or.(nactdog(3,nodem).ne.0))then
            nacteq(1,node2)=1   ! mass equation
          endif
          if (nactdog(0,node2).ne.0)then
            nacteq(0,node2)=1   ! energy equation
          endif
!     
!     "end of network" element node O---X---
!     1   m   2
        elseif (node2.eq.0) then
          if ((nactdog(1,nodem).ne.0).or.(nactdog(3,nodem).ne.0))then
            nacteq(1,node1)=1   ! mass equation
          endif
          if (nactdog(0,node1).ne.0)then
            nacteq(0,node1)=1   ! energy equation
          endif
!     
!     "flow element" O---X---O
!     1   m   2     
!     
        else  
          if((nactdog(1,nodem).ne.0).or.(nactdog(3,nodem).ne.0)) then
            nacteq(1,node2)=1   ! mass equation
            nacteq(1,node1)=1   ! mass equation
          endif
!     
          ider=0
          call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &         nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &         nodef,idirf,df,cp,r,rho,physcon,g,co,dvi,numf,
     &         vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &         ttime,time,iaxial,iplausi)
!     
          if (.not.identity) then
            nacteq(2,nodem)=1   ! momentum equation
          endif
!     
          if (nactdog(0,node1).ne.0)then
            nacteq(0,node1)=1   ! energy equation
          endif     
!     
          if (nactdog(0,node2).ne.0)then 
            nacteq(0,node2)=1   ! energy equation
          endif
        endif 
      enddo
!     
!     wall convective contributions
!     
      do i=1, nload
        if((sideload(i)(3:4)).eq.'FC')then
          node=nelemload(2,i)
          if (nactdog(0,node).ne.0) then
            nacteq(0,node)=1
          endif
        endif
      enddo
!     
!     restoring the dependent nodes of the MPC's
!     only for mass flow and pressure
!     
      if(networkmpcs.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:7).ne.'NETWORK') cycle
          index=ipompc(i)
          idir=nodempc(2,index)
          if((idir.eq.1).or.(idir.eq.2)) then
            node=nodempc(1,index)
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                nactdog(idir,node)=1
              endif
            endif
          endif
        enddo
      endif
!     
!     removing the energy equation from those end nodes for which
!     the temperature constitutes the first term in a network MPC
!     
      if(networkmpcs.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:7).ne.'NETWORK') cycle
          index=ipompc(i)
          idir=nodempc(2,index)
          if(idir.eq.0) then
            node=nodempc(1,index)
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                nacteq(0,node)=0
              endif
            endif
          endif
        enddo
      endif
!     
!     check whether all mass flow is known
!     
      do i=1,ntg
        node=itg(i)
        if((nactdog(1,node).ne.0).or.(nactdog(3,node).ne.0)) then
          massflowbcall=.false.
          exit
        endif
      enddo
!     
!     check whether all pressures are known
!     
      do i=1,ntg
        node=itg(i)
        if(nactdog(2,node).ne.0) then
          pressurebcall=.false.
          exit
        endif
      enddo
!     
!     check for special cases
!     
      if(massflowbcall.and.((.not.pressurebc).or.(pressurebcall))) then
!     
!     purely thermal (only set to 1 if D-type elements are present)
!     
        if(network.gt.1) network=2
        do i=1,nflow
          nelem=ieg(i)
          index=ipkon(nelem)
          node1=kon(index+1)
          nodem=kon(index+2)
          node2=kon(index+3)
          nacteq(2,nodem)=0
          if(node1.ne.0) nactdog(2,node1)=0
          if(node2.ne.0) nactdog(2,node2)=0
        enddo
      elseif((.not.temperaturebc).and.(.not.walltemp)) then
!     
!     pure liquid dynamics
!     
        write(*,*) '*INFO in envtemp: no thermal boundary conditions'
        write(*,*) '      detected; the network is considered to be'
        write(*,*) '      athermal and no gas temperatures will be'
        write(*,*) '      calculated'
        network=4
        do i=1,ntg
          node=itg(i)
          nactdog(0,node)=0
          nacteq(0,node)=0
        enddo
      elseif((.not.temperaturebc).and.walltemp) then
        write(*,*) '*ERROR in envtemp: at least one temperature'
        write(*,*) '       boundary condition must be given'
        call exit(201)
      elseif(.not.pressurebc) then
        write(*,*) '*ERROR in envtemp: at least one pressure'
        write(*,*) '       boundary condition must be given'
        call exit(201)
      endif
!     
!     check whether a specific gas constant was defined for all fluid
!     elements (except for purely thermal calculations)
!     
      if(network.gt.2) then
        do i=1,nflow
          nelem=ieg(i)
          if((lakon(nelem)(2:3).eq.'LI').or.
     &         (lakon(nelem)(2:3).eq.'LP').or.
     &         (lakon(nelem)(2:3).eq.'  ')) cycle
          imat=ielmat(1,nelem)
          r=shcon(3,1,imat)
          if(r.lt.1.d-10) then
            write(*,*)'*ERROR in envtemp: specific gas',
     &           'constant is close to zero'
            call exit(201)
          endif
        enddo
      endif
!     
!     check whether the temperature at each inlet or outlet node
!     is given
!     
      if(network.lt.4) then
        do i=1,nflow
          nelem=ieg(i)
          index=ipkon(nelem)
          node1=kon(index+1)
          node2=kon(index+3)
          if(node1.eq.0) then
            if(nactdog(0,node2).ne.0) then
              write(*,*) '*WARNING in envtemp: it is advised to'
              write(*,*) '         define the temperature at all'
              write(*,*) '         inlets and outlets by a boundary'
              write(*,*) '         condition. This is lacking for'
              write(*,*) '         node ',node2,' of element ',nelem
            endif
          endif
          if(node2.eq.0) then
            if(nactdog(0,node1).ne.0) then
              write(*,*) '*WARNING in envtemp: it is advised to'
              write(*,*) '         define the temperature at all'
              write(*,*) '         inlets and outlets by a boundary'
              write(*,*) '         condition. This is lacking for'
              write(*,*) '         node ',node1,' of element ',nelem
            endif
          endif
        enddo
      endif
!     
!     numbering the active equations      
!     
      nteq=0
      do i=1,ntg
        node=itg(i)
        do j=0,2
          if (nacteq(j,node).ne.0) then
            nteq=nteq+1
            nacteq(j,node)=nteq
          endif 
        enddo
!     write(30,*) 'unknowns ',node,(nactdog(j,node),j=0,3)
      enddo
!     do i=1,ntg
!     node=itg(i)
!     write(30,*) 'equations',node,(nacteq(j,node),j=0,2)
!     enddo
!     
!     taking network MPC's into account
!     
      if(networkmpcs.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:7).eq.'NETWORK') nteq=nteq+1
        enddo
      endif
!     
!     numbering the active degrees of freedom
!     
      ntq=0
      do i=1,ntg
        node=itg(i)
        do  j=0,3
          if (nactdog(j,node).ne.0) then
            ntq=ntq+1
            nactdog(j,node)=ntq
          endif 
        enddo
      enddo
c     
c      open(30,file='dummy',status='unknown')
c      write(30,*) 'nactdog'
c      do i=1,ntg
c        write(30,*) itg(i),(nactdog(j,itg(i)),j=0,3)
c      enddo
c
c      write(30,*) ''
c      write(30,*) 'nacteq'
c      do i=1,ntg
c        write(30,*) itg(i),(nacteq(j,itg(i)),j=0,3)
c      enddo
c      close(30)
!     
      if(ntq.ne.nteq) then
        write(*,*) '*ERROR in envtemp:'
        write(*,*) '*****number of network equations is not equal to'
        write(*,*) ' number of active degrees of freedom*****'
        write(*,*) ' # of network equations = ',nteq
        write(*,*) ' # of active degrees of freedom= ',ntq
        call exit(201)
      endif   
!     
!     for isothermal gas pipes the energy equation in the
!     topologically downstream node is replaced by an equation
!     expressing the equality of the static temperature at both
!     ends of the pipe. To this end these downstream nodes are
!     referring in nacteq(3,*) to the topologically upstream node
!     
!     if the temperature in the downstream node is a boundary
!     condition (i.e. the energy equation is not built), the
!     energy equation in the upstream node is replaced. In that 
!     case nacteq(3,upstreamnode) refers to the downstream node.
!     
!     if both nodes are boundary conditions, nothing is done
!     
      do i=1,nflow
        nelem=ieg(i)
        if((lakon(nelem)(1:4).eq."DGAP")
     &       .and.(lakon(nelem)(6:6).eq."I")) then
!     
          index=ipkon(nelem)
          node1=kon(index+1)
          node2=kon(index+3)
          if((node1.eq.0).or.(node2.eq.0)) cycle
!     
          if(nacteq(0,node2).ne.0) then
            nacteq(3,node2)=node1
          elseif(nacteq(0,node1).ne.0) then
            nacteq(3,node1)=node2
          endif
        endif
      enddo
!     
      return
      end
      
      
