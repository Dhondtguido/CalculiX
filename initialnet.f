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
!     calculate the initial conditions for the gas 
!         - the initial pressure
!         - identifying the chambers and gas pipe nodes
!           for gas networks
!         - the initial flow
!         - calculating the static temperature for gas networks
!     
      subroutine initialnet(itg,ieg,ntg,ac,bc,lakon,v,
     &     ipkon,kon,nflow,ikboun,nboun,prop,ielprop,
     &     nactdog,ndirboun,nodeboun,xbounact,
     &     ielmat,ntmat_,shcon,nshcon,physcon,ipiv,nteq,
     &     rhcon,nrhcon,ipobody,ibody,xbodyact,co,nbody,network,
     &     iin_abs,vold,set,istep,iit,mi,ineighe,ilboun,ichannel,
     &     iaxial,nmpc,labmpc,ipompc,nodempc,coefmpc,ttime,time,
     &     iponoel,inoel)
!     
      implicit none
     
      logical identity,gravity,gaspipe
!
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 set(*)
!           
      integer mi(*),ieg(*),nflow,i,j,ntg,ielmat(mi(3),*),ntmat_,
     &     node1,node2,ider,iaxial,nmpc,ipompc(*),nodempc(3,*),
     &     nelem,index,nshcon(*),ipkon(*),kon(*),ikboun(*),nboun,idof,
     &     nodem,idirf(8),nactdog(0:3,*),imat,ielprop(*),id1,id2,
     &     nodef(8),ndirboun(*),nodeboun(*),itg(*),node,kflag,ipiv(*),
     &     nrhs,info,idof1,idof2,nteq,nrhcon(*),ipobody(2,*),ibody(3,*),
     &     nbody,numf,network,iin_abs,icase,index2,index1,nelem1,nelem2,
     &     node11,node21,node12,node22,istep,iit,ineighe(*),iponoel(*),
     &     ilboun(*),idir,ichannel,inoel(2,*),indexe,iel,iplausi,id
!     
      real*8 ac(nteq,nteq), bc(nteq),prop(*),shcon(0:3,ntmat_,*),
     &     f,df(8),xflow,xbounact(*),v(0:mi(2),*),cp,r,tg1,
     &     tg2,gastemp,physcon(*),pressmin,dvi,rho,g(3),z1,z2,
     &     rhcon(0:1,ntmat_,*),co(3,*),xbodyact(7,*),kappa,
     &     A,Tt,pt,Ts,pressmax,constant,vold(0:mi(2),*),
     &     coefmpc(*),ttime,time,xflow360,A2,d,l,s
!
      kflag=1
      ider=0
      ichannel=0
      gaspipe=.false.
!
!     applying the boundary conditions
!
      do j=1,nboun
        v(ndirboun(j),nodeboun(j))=xbounact(j)
        vold(ndirboun(j),nodeboun(j))=xbounact(j)
      enddo
!     
!     turning a boundary pressure of 0.d0 into 1.235711d-10
!     (if the pressure is 0.d0 it is assumed in the present
!     routine to be unknown and not to be set by a boundary
!     condition; therefore it is changed here to a very small
!     value and reset at the end of the routine to zero; this
!     should only be relevant in liquids, where only pressure
!     differences matter)
!     
       do i=1,ntg
         node=itg(i)
         idof=8*(node-1)+2
         call nident(ikboun,idof,nboun,id)
         if(id.gt.0) then
           if(ikboun(id).eq.idof) then
             if(v(2,node).eq.0.d0) v(2,node)=1.235711d-10
           endif
         endif
       enddo
!
!     check for channel elements
!     
      do i=1,nflow
         nelem=ieg(i)
         if(lakon(nelem)(2:5).eq.'LICH') then
            ichannel=1
            return
         endif
      enddo
!     
!     initializing ac and bc
!     
      do i=1,nteq
         do j=1,nteq
            ac(i,j)=0.d0
         enddo
         bc(i)=0.d0
      enddo
!     
!     for all but purely thermal networks:
!     determining the initial pressure
!     and identifying the chamber and gas pipe nodes
!     
      if(network.gt.2) then
!
         call preinitialnet(ieg,lakon,v,ipkon,kon,nflow,prop,ielprop,
     &     ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,mi,iponoel,inoel,
     &        itg,ntg)
!     
!     determining whether pressure initial conditions 
!     are provided for all nodes
!     
         pressmin=-1.d0
         pressmax=0.d0
         constant=1.55d0
         
         do i=1,ntg
            node=itg(i)
            if(v(2,node).lt.1.d-10) then
               v(2,node)=0.d0
            else
               if(pressmin.lt.0.d0) then
                  pressmin=v(2,node)
               elseif(v(2,node).lt.pressmin) then
                  pressmin=v(2,node)
               endif
!     
               if(v(2,node).gt.pressmax)then
                  pressmax=v(2,node)
               endif
!     
            endif
         enddo
!     
         if(pressmin.lt.0.d0) then
            write(*,*) 
     &           '*ERROR in initialnet: minimum initial pressure'
            write(*,*) '       is smaller than zero'
            call exit(201)
         endif
!     
!     in nodes in which no initial pressure is given v(2,*)
!     is replaced by -n, where n is the number of elements the
!     node belongs to: allows to find boundary nodes of the 
!     network
!     
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            call nident(itg,node1,ntg,id1)
            call nident(itg,node2,ntg,id2)
!     
            if(iin_abs.eq.0) then
               if ((lakon(nelem)(2:5).eq.'GAPF').or.
     &             (lakon(nelem)(2:5).eq.'GAPR').or.
     &             ((lakon(nelem)(2:3).eq.'RE')
     &              .and.(lakon(nelem)(2:7).ne.'REWAOR')).or.
     &             (lakon(nelem)(2:3).eq.'UP')) then 
!     
!     In the case of a element of type GASPIPE or RESTRICTOR 
!     (except TYPE= RESTRICTOR WALL ORIFICE) or USER
!     the number of pipes connected to node 1 and 2
!     are computed and stored in ineighe(id1)
!     respectively ineighe(id2)
!     
                  gaspipe=.true.
                  if(node1.ne.0) then
                     if (ineighe(id1).ge.0) then
!     
                        if(node2.ne.0)then
                           ineighe(id1)=ineighe(id1)+1
                        endif
                     endif
                  endif
                  if(node2.ne.0) then
                     if (ineighe(id2).ge.0) then
                        if(node1.ne.0) then
                           ineighe(id2)=ineighe(id2)+1
                        endif
                     endif
                  endif
               else
!     
!     for all other elements (different from GASPIPE or 
!     RESTRICTOR), including RESTRICTOR WALL ORIFICE 
!     ineighe(idi)=-1
!     which means that they are connected to chambers
!     i.e. static and total values are equal
!     
                  if (node1.ne.0) then
                     ineighe(id1)=-1
                  endif
                  if(node2.ne.0) then
                     ineighe(id2)=-1
                  endif
               endif
            endif
!     
            if((node1.eq.0).or.(node2.eq.0)) cycle
            if(v(2,node1).lt.1.d-10) then
               v(2,node1)=v(2,node1)-1.d0
            endif
            if(v(2,node2).lt.1.d-10) then
               v(2,node2)=v(2,node2)-1.d0
            endif
         enddo
!
!              Guido 17.06.2022 Option 1  (cf. resultnet.f)
!
c!     
c!     for each end node i: if ineighe(i)<0: chamber
c!     else: ineighe(i)=number of pipe connections
c!
c!     check that isothermal gaspipe elements are not connected
c!     to chambers
c!     
c         do i=1,nflow
c           nelem=ieg(i)
c           index=ipkon(nelem)
c           node1=kon(index+1)
c           node2=kon(index+3)
c           call nident(itg,node1,ntg,id1)
c           call nident(itg,node2,ntg,id2)
c!     
c           if(iin_abs.eq.0) then
c             if (lakon(nelem)(2:6).eq.'GAPFI') then 
c               if((ineighe(id1).eq.-1).or.(ineighe(id2).eq.-1)) then
c                 write(*,*) '*ERROR in initialnet: element',nelem
c                 write(*,*) '       is an isothermal gas pipe'
c                 write(*,*) '       element but is connected to'
c                 write(*,*) '       at least one chamber'
c                 call exit(201)
c               endif
c             endif
c           endif
c         enddo
!     
!     assigning values to the boundary nodes of the network
!     (i.e. nodes belonging to only one element)
!     
         do i=1,ntg
            node=itg(i)
!     boundary condition nodes
            if(nactdog(2,node).eq.0) then 
               cycle 
            endif
!     initial condition nodes
            if((nactdog(2,node).ne.0)
     &           .and.(v(2,node).gt.0d0)) then
               cycle
            endif
!     nodes neither defined as *BOUNDARY nor as *INITIAL CONDITIONS
            if(abs(v(2,node)+1.d0).lt.1.d-10) then
!
!              if a nonzero flux is leaving the node increase
!              the maximum pressure, else decrease the minimum
!              pressure
!
!              determining the exit element (= element with one
!              zero node)
!
               index=iponoel(node)
               do
                  iel=inoel(1,index)
                  indexe=ipkon(iel)
                  node1=kon(indexe+1)
                  nodem=kon(indexe+2)
                  node2=kon(indexe+3)
                  if((node1.eq.0).or.(node2.eq.0)) exit
                  index=inoel(2,index)
               enddo
!
               if(((node1.eq.node).and.(v(1,nodem).ge.0.d0)).or.
     &            ((node2.eq.node).and.(v(1,nodem).le.0.d0))) then
                  v(2,node)=0.95d0*pressmin
                  pressmin=0.95d0*pressmin
               else
                  v(2,node)=1.05d0*pressmax
                  pressmax=1.05d0*pressmax
               endif
            endif              
         enddo
!
!        transforming the pressures using the tangent
!        function
!
         do i=1,ntg
            node=itg(i)
!     boundary condition nodes
            if(nactdog(2,node).eq.0) then 
               v(2,node)=dtan(v(2,node)
     &              *constant/pressmax) 
!     initial condition nodes
            elseif(v(2,node).gt.0d0) then
               v(2,node)=dtan(v(2,node)
     &              *constant/pressmax) 
            endif
         enddo
!     
!     mass flow and geometry dofs: 1 on the diagonal
!     
         do i=1,ntg
            node=itg(i)
            do j=1,1
               if(nactdog(j,node).ne.0) then
                  idof=nactdog(j,node)
                  ac(idof,idof)=1.d0
               endif
            enddo
            do j=3,3
               if(nactdog(j,node).ne.0) then
                  idof=nactdog(j,node)
                  ac(idof,idof)=1.d0
               endif
            enddo
         enddo
!     
!     pressure dofs
!     
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            if((node1.eq.0).or.(node2.eq.0)) cycle
            idof1=nactdog(2,node1)
            idof2=nactdog(2,node2)
            if(idof1.ne.0) then
               if(v(2,node1).gt.0.d0) then
!     initial pressure given in node1
                  ac(idof1,idof1)=1.d0
                  bc(idof1)=v(2,node1)
               else
                  ac(idof1,idof1)=ac(idof1,idof1)+1.d0
                  if(v(2,node2).gt.0.d0) then
!     initial pressure given in node2
                     bc(idof1)=bc(idof1)+v(2,node2)
                  else
                     ac(idof1,idof2)=ac(idof1,idof2)-1.d0
                  endif
               endif
            endif
            if(idof2.ne.0) then
               if(v(2,node2).gt.0.d0) then
!     initial pressure given in node2
                  ac(idof2,idof2)=1.d0
                  bc(idof2)=v(2,node2)
               else
                  ac(idof2,idof2)=ac(idof2,idof2)+1.d0
                  if(v(2,node1).gt.0.d0) then
!     initial pressure given in node1
                     bc(idof2)=bc(idof2)+v(2,node1)
                  else
                     ac(idof2,idof1)=ac(idof2,idof1)-1.d0
                  endif
               endif
            endif
         enddo
!     
!     checking for nodes not belonging to network elements
!     
         do i=1,ntg
            node=itg(i)
            if(nactdog(2,node).ne.0) then
               idof=nactdog(2,node)
               if(ac(idof,idof).eq.0.d0) then
                  ac(idof,idof)=1.d0
                  bc(idof)=v(2,node)
               endif
            endif
         enddo
      endif
!     
!     temperature conditions for all but purely hydrodynamic
!     networks
!     
      if(network.ne.4) then
!     
!     temperature dofs
!     
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            nodem=kon(index+2)
            node2=kon(index+3)
            if((node1.eq.0).or.(node2.eq.0)) cycle
            idof1=nactdog(0,node1)
            idof2=nactdog(0,node2)
            if(idof1.ne.0) then
               if(v(0,node1).gt.0.d0) then
!     initial temperature given in node1
                  ac(idof1,idof1)=1.d0
                  bc(idof1)=v(0,node1)
!
!                 temperature taken into account only for incoming
!                 flux (from the viewpoint of node1).
!
               elseif(v(1,nodem).lt.0.d0) then
                  ac(idof1,idof1)=ac(idof1,idof1)+1.d0
                  if(v(0,node2).gt.0.d0) then
!     initial temperature given in node2
                     bc(idof1)=bc(idof1)+v(0,node2)
                  else
                     ac(idof1,idof2)=ac(idof1,idof2)-1.d0
                  endif
               endif
            endif
            if(idof2.ne.0) then
               if(v(0,node2).gt.0.d0) then
!     initial temperature given in node2
                  ac(idof2,idof2)=1.d0
                  bc(idof2)=v(0,node2)
!
!                 temperature taken into account only for incoming
!                 flux (from the viewpoint of node2). Default (if
!                 no flux is predefined) is positive flux.
!
               elseif(v(1,nodem).ge.0.d0) then
                  ac(idof2,idof2)=ac(idof2,idof2)+1.d0
                  if(v(0,node1).gt.0.d0) then
!     initial temperature given in node1
                     bc(idof2)=bc(idof2)+v(0,node1)
                  else
                     ac(idof2,idof1)=ac(idof2,idof1)-1.d0
                  endif
               endif
            endif
         enddo
!     
!     checking for nodes not belonging to network elements
!     
         do i=1,ntg
            node=itg(i)
            if(nactdog(0,node).ne.0) then
               idof=nactdog(0,node)
               if(ac(idof,idof).eq.0.d0) then
                  ac(idof,idof)=1.d0
                  bc(idof)=v(0,node)
               endif
            endif
         enddo
      endif
!
!     additional multiple point constraints
!
      do j=1,nteq
         if(dabs(ac(j,j)).lt.1.d-20) ac(j,j)=1.d0
      enddo
!     
!     solving the system
!     
      if(nteq.gt.0) then
         nrhs=1
         call dgesv(nteq,nrhs,ac,nteq,ipiv,bc,nteq,info)
         if(info.ne.0) then
            write(*,*) '*ERROR in initialnet: singular matrix'
            call exit(201)
         endif
      endif
!     
!     storing the initial pressure in v
!     (not for purely thermal networks)
!     
      if(network.gt.2) then
         do i=1,ntg
            node=itg(i)
            if(nactdog(2,node).eq.0) then
               v(2,node)=pressmax
     &              *datan(v(2,node))/constant
               cycle
            endif
            v(2,node)=pressmax
     &           *datan(bc(nactdog(2,node)))/constant
         enddo
      endif
!     
!     storing the initial temperature in v
!     (not for purely hydrodynamic networks)
!     
      if(network.ne.4) then
         do i=1,ntg
            node=itg(i)
            if(nactdog(0,node).ne.0) v(0,node)=bc(nactdog(0,node))
         enddo
      endif
!     
!     for ELEMENT TYPE BRANCH
!     initial pressures in branch 1 and 2 are equal and set to the 
!     maximum of the two pressures
!     
      do i=1,nflow
         nelem=ieg(i)
         if(lakon(i)(2:5).ne.'REBR') cycle
         index=ielprop(nelem)
!     
         nelem1=nint(prop(index+2))
         index1=ipkon(nelem1)
         node11=kon(index1+1)
         node21=kon(index1+3)
!     
         nelem2=nint(prop(index+3))
         index2=ipkon(nelem2)
         node12=kon(index2+1)
         node22=kon(index2+3)
!     
         if(node11.eq.node12) then
            node1=node21
            node2=node22
         elseif(node11.eq.node22) then
            node1=node21
            node2=node12
         elseif(node21.eq.node12) then
            node1=node11
            node2=node22
         elseif(node21.eq.node22) then
            node1=node11
            node2=node12
         endif
!     
         if(lakon(i)(6:6).eq.'S') then
!
!           maximum
!
            if(v(2,node1).ge.v(2,node2)) then
               v(2,node2)=v(2,node1)
            else
               v(2,node1)=v(2,node2)
            endif
         elseif(lakon(i)(6:6).eq.'J') then
!
!           minimum
!
            if(v(2,node1).ge.v(2,node2)) then
               v(2,node1)=v(2,node2)
            else
               v(2,node2)=v(2,node1)
            endif
         endif
      enddo
!     
!     identifying the chamber nodes for purely thermal
!     gas networks (needed to determine the static
!     temperature which is used for the material properties)
!     
c      if(network.le.2) then
      if((network.le.2).and.(iin_abs.eq.0)) then
         do i=1,nflow
            nelem=ieg(i)
            if((lakon(nelem)(2:3).eq.'LP').or.
     &           (lakon(nelem)(2:3).eq.'LI')) cycle
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            if(node1.ne.0) then
               call nident(itg,node1,ntg,id1)
               ineighe(id1)=-1
            endif
            if(node2.ne.0) then
               call nident(itg,node2,ntg,id2)
               ineighe(id2)=-1
            endif
         enddo
      endif
!     
!     initialisation of bc
!     
      do i=1,nteq
         bc(i)=0.d0
      enddo
!     
!     determining the initial mass flow in those nodes for which no
!     flux boundary conditions are defined
!     
      if(network.gt.2)then
         do i=1,nflow
            nelem=ieg(i)
!     
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            if((node1.eq.0).or.(node2.eq.0)) cycle
            nodem=kon(index+2)
!     
!     cycle if both the geometry and the mass flow are known
!     
            if((nactdog(1,nodem).eq.0).and.(nactdog(3,nodem).eq.0)) 
     &           cycle
!     
!     for liquids: determine the gravity vector
!     
            if((lakon(nelem)(2:3).eq.'LI').or.
     &         (lakon(nelem)(2:3).eq.'LP')) then
               gravity=.false.
               do j=1,3
                  g(j)=0.d0
               enddo
               if(nbody.gt.0) then
                  index=nelem
                  do
                     j=ipobody(1,index)
                     if(j.eq.0) exit
                     if(ibody(1,j).eq.2) then
                        g(1)=g(1)+xbodyact(1,j)*xbodyact(2,j)
                        g(2)=g(2)+xbodyact(1,j)*xbodyact(3,j)
                        g(3)=g(3)+xbodyact(1,j)*xbodyact(4,j)
                        z1=-g(1)*co(1,node1)-g(2)*co(2,node1)
     &                       -g(3)*co(3,node1)
                        z2=-g(1)*co(1,node2)-g(2)*co(2,node2)
     &                       -g(3)*co(3,node2)
                        gravity=.true.
                     endif
                     index=ipobody(2,index)
                     if(index.eq.0) exit
                  enddo
               endif
               if(.not.gravity) then
                  write(*,*) 
     &                 '*WARNING in initialnet: no gravity vector'
                  write(*,*) 
     &                 '       was defined for liquid element',nelem
                  z1=0.d0
                  z2=0.d0
c                  call exit(201)
               endif
            endif
!     
            tg1=v(0,node1)
            tg2=v(0,node2)
!     
!     For liquid pipes elements the upstream temperature 
!     is used to determine the material properties
!     For all other elements the averaged value of the 
!     temperature between inlet and outlet is used
!     
            if((lakon(nelem)(2:3).ne.'LP').and.
     &           (lakon(nelem)(2:3).ne.'LI')) then
               gastemp=(tg1+tg2)/2.d0
            else
               xflow=v(1,nodem)
               if(xflow.gt.0) then
                  gastemp=tg1
               else
                  gastemp=tg2
               endif
            endif
!     
            imat=ielmat(1,nelem)
            call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,
     &           r,dvi,rhcon,nrhcon,rho)
!     
!     If inlet and outlet pressure are equal this leads to a massflow rate 
!     equal to 0 and in turn possibly to a singular matrix configuration
!     
!     gas
            if(((lakon(nelem)(2:3).ne.'LI').and.
     &           (v(2,node1).eq.v(2,node2))).or.
!     liquid
     &           ((lakon(nelem)(2:3).eq.'LI').and.
     &           (v(2,node1)+z1.eq.v(2,node2)+z2))) then
            
!     
!     if neither inlet nor outlet pressure are active D.O.F: error-> stop
!     
               if((nactdog(2,node1).eq.0)
     &           .and.(nactdog(2,node2).eq.0)) then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*ERROR: in subroutine initialnet.f'
                  write(*,*) '        no pressure gradient'
                  write(*,*) '        in element', nelem
                  write(*,*) '        Inlet and outlet pressures are '
                  write(*,*) '        boundary conditions '
                  write(*,*) '        node1',node1,' pressure',
     &                 v(2,node1)
                  write(*,*) '        node2',node2,' pressure',
     &                 v(2,node2)
                  call exit(201)
!     
!     if inlet pressure is an active degree of freedom
!     
               else if((nactdog(2,node1).ne.0)
     &                 .and.(nactdog(2,node2).eq.0))then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*WARNING: in subroutine initialnet.f'
                  write(*,*) '          no pressure gradient'
                  write(*,*) '          in element', nelem
                  write(*,*) 
     &               '          Inlet pressure initial condition '
                  write(*,*) '          is changed '
                  write(*,*) '          node1',node1,
     &                 ' given initial pressure',v(2,node1)
                  v(2,node1)=1.1d0*v(2,node1)
                  write(*,*) '          node1',node1,
     &                 ' new initial pressure',v(2,node1)
                  write(*,*) '          node2',node2,' pressure',
     &                 v(2,node2)
!     
!     if outlet pressure is an active D.O.F.
!     
               else if((nactdog(2,node1).eq.0)
     &                 .and.(nactdog(2,node2).ne.0))then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*WARNING: in subroutine initialnet.f'
                  write(*,*) '          no pressure gradient'
                  write(*,*) '          in element', nelem
                  write(*,*) 
     &                 '          Outlet pressure initial condition '
                  write(*,*) '          is changed '
                  write(*,*) '          node1',node1,' pressure'
     &                 ,v(2,node1)    
                  write(*,*) '          node2',node2,
     &                 'given intial pressure',
     &                 v(2,node2)
                  v(2,node2)=0.9d0*v(2,node2)
                  write(*,*) '          node2',node2,
     &                 ' new initial pressure',v(2,node2)
!     
!     if both inlet and outlet pressures are active D.O.F.
!     
               else if((nactdog(2,node1).ne.0)
     &                 .and.(nactdog(2,node2).ne.0))then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*WARNING: in subroutine initialnet.f'
                  write(*,*) '          no pressure gradient'
                  write(*,*) '          in element', nelem
                  write(*,*) '          Inlet and outlet pressure '
                  write(*,*) '          initial condition are changed '
                  write(*,*) '          node1',node1,
     &                 ' given initial pressure',v(2,node1)  
                  v(2,node1)=1.05d0*v(2,node2)
                  write(*,*) '          node1',node1,
     &                 ' new intial pressure',v(2,node1)
                  write(*,*) '          node2',node2,
     &                 ' given initial pressure',v(2,node2)
                  v(2,node2)=0.95d0*v(2,node2)
                  write(*,*) '          node2',node2,
     &                 ' new intial pressure',v(2,node2)
               endif
            endif
!     
!            calculating flux if the flux is an unknown; any initial
!            flux defined by the user is checked on plausibility            
!
            if(nactdog(1,nodem).ne.0) then
               call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &           nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &           nodef,idirf,df,cp,r,rho,physcon,g,co,dvi,numf,
     &           vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &           ttime,time,iaxial,iplausi)
               v(1,nodem)=xflow
            endif
!     
            if((lakon(nelem)(2:4).ne.'LIP').and.
     &         (lakon(nelem)(2:3).ne.'VO').and.
     &         (lakon(nelem)(2:5).ne.'GAPR').and.
     &         (lakon(nelem)(2:4).ne.'ATR').and.
     &         (lakon(nelem)(2:4).ne.'RTA')) then
               if(v(1,nodem).eq.0d0) then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*ERROR: in subroutine initialnet.f'
                  write(*,*) '        in element', nelem,
     &                 lakon(nelem)(1:6)
                  write(*,*) '        mass flow rate value = 0 !'
                  write(*,*) '        node1',node1,' pressure',
     &                 v(2,node1)
                  write(*,*) '        node2',node2,' pressure',
     &                 v(2,node2)
                  call exit(201)
               endif
               if (v(1,nodem)*(v(2,node1)-v(2,node2)).lt.0) then
                  WRITE(*,*) '**************************************'
                  write(*,*) '*WARNING: in subroutine initialnet.f'
                  write(*,*) '          in element', nelem
                  write(*,*) 
     &                 '          initial mass flow rate value does not'
                  write(*,*) '          correspond to pressure gradient'
                  write(*,*) '          node1',node1,'pressure',
     &                 v(2,node1)
                  write(*,*) '          node2',node2,'pressure',
     &                 v(2,node2)
                  write(*,*) '          nodem',nodem,'mass flow',
     &                 v(1,nodem)
                  write(*,*) 
     &               '          the sign of the initial mass flow'
                  write(*,*) '          rate will be changed'
               endif
            endif
         enddo
      endif
!
      call postinitialnet(ieg,lakon,v,ipkon,kon,nflow,prop,ielprop,
     &     ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,mi,iponoel,inoel,
     &     itg,ntg,nactdog)
!     
!     calculating the static temperature for nodes belonging to gas pipes
!     and restrictors (except RESTRICTOR WALL ORIFICE)
!     
      if (gaspipe.and.(iin_abs.eq.0)) then
!     
!     ineighe(i) is set to -1 (= chamber node) if more than 2 pipes
!     are connected to the node
!     
         do i=1,ntg
            if(ineighe(i).gt.2) then
               ineighe(i)=-1
               write(*,*) '*WARNING :in subroutine initialnet.f'
               write(*,*) '          more than 2 elements GASPIPE'
               write(*,*) '          or RESTRICTOR are connected '
               write(*,*) '          to node',itg(i),'. The common'
               write(*,*) 
     &              '          node is converted into a chamber.'
               write(*,*) '          Total and static parameters are'
               write(*,*) '          equal'
            endif
         enddo
!     
         do i=1,nflow
            nelem=ieg(i)
            index=ipkon(nelem)
            node1=kon(index+1)
            node2=kon(index+3)
            call nident(itg,node1,ntg,id1)
            call nident(itg,node2,ntg,id2)
!     
!           for each end node i: 
!             if ineighe(i)=-1: chamber
!                          = 0: middle node
!                          > 0: number of pipe connections (max. 2 allowed)
!     
!           the number of pipe connections for ineighe(i)>0 is 
!           replaced by the number of an element the node belongs to:
!     
            if(node1.gt.0) then
               if((ineighe(id1).eq.1).or.
     &              (ineighe(id1).eq.2)) then
                  if(node2.ne.0)then
                     ineighe(id1)=nelem
                  endif
               endif
            endif
!     
            if(node2.gt.0) then
               if((ineighe(id2).eq.1).or.
     &              (ineighe(id2).eq.2)) then
                  if (node1.ne.0) then
                     ineighe(id2)=nelem
                  endif
               endif
            endif
         enddo
!     
!     The static temperature is calculated and stored in v(3,node)
!     total temperatures are supposed equal (adiabatic pipe)
!     
         do i=1,ntg
            node=itg(i)
            if(ineighe(i).gt.0) then 
!     
               nelem=ineighe(i)
               index=ielprop(nelem)
               nodem=kon(ipkon(nelem)+2)
!     
               imat=ielmat(1,nelem)
               call materialdata_tg(imat,ntmat_,v(0,node),
     &              shcon,nshcon,cp,r,dvi,rhcon,nrhcon,rho)
               kappa=cp/(cp-R)
               xflow=v(1,nodem)
               Tt=v(0,node)
               pt=v(2,node)
!     
               if(lakon(nelem)(2:5).eq.'GAPF') then
!
!                 distinguish between standard pipes and flexible
!                 pipes
! 
                  if((lakon(nelem)(6:8).eq.'AFR').or.
     &               (lakon(nelem)(6:8).eq.'ARL').or.
     &               (lakon(nelem)(6:8).eq.'ARG').or.
     &               (lakon(nelem)(6:8).eq.'IFR').or.
     &               (lakon(nelem)(6:8).eq.'IRL')) then
                     call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &                   ttime,time,ielprop,mi,A,A2,d,l,s)
                  else
                     A=prop(index+1)
                  endif
!
                  if(lakon(nelem)(2:6).eq.'GAPFA') then
                     icase=0
                  elseif(lakon(nelem)(2:6).eq.'GAPFI') then
                     icase=1
                  endif     
!
               elseif(lakon(nelem)(2:5).eq.'GAPR') then
!
!                 rotating adiabatic gas pipe with variable cross section:
!                 check whether section 1 or 2                  
!
                  if(kon(ipkon(nelem)+1).eq.node) then
                     A=prop(index+1)
                  else
                     A=prop(index+2)
                  endif
                  icase=0
               elseif(lakon(nelem)(2:3).eq.'RE') then
                  index2=ipkon(nelem)
                  node1=kon(index2+1)
                  node2=kon(index2+3)
                  if(lakon(nelem)(4:5).eq.'EX') then
                     if(lakon(nint(prop(index+4)))(2:6).eq.'GAPFA')
     &                    then
                        icase=0
                     elseif
     &                   (lakon(nint(prop(index+4)))(2:6).eq.'GAPFI') 
     &                       then
                        icase=1
                     endif
                  else
                     icase=0
                  endif
!     
                  if(lakon(nelem)(4:5).eq.'BE') then
                     A=prop(index+1)
!     
                  elseif(lakon(nelem)(4:5).eq.'BR') then
                     if(lakon(nelem)(4:6).eq.'BRJ') then
                        if(nelem.eq.nint(prop(index+2)))then
                           A=prop(index+5)
                        elseif(nelem.eq.nint(prop(index+3))) then
                           A=prop(index+6)
                        endif
                     elseif(lakon(nelem)(4:6).eq.'BRS') then
                        if(nelem.eq.nint(prop(index+2)))then
                           A=prop(index+5)
                        elseif(nelem.eq.nint(prop(index+3))) then
                           A=prop(index+6)
                        endif
                     endif
!     
                  else
!
                     if((lakon(nelem)(2:5).eq.'RBSF')) then
                        call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &                       ttime,time,ielprop,mi,A,A2,d,l,s)
                        if(node.eq.node2) then
                           A=A2
                        endif
                     else
                        if(node.eq.node1) then
                           A=prop(index+1)
                        elseif(node.eq.node2) then
                           A=prop(index+2)
                        endif
                     endif
                  endif
               elseif(lakon(nelem)(2:3).eq.'UP') then
!
!                 user elements whose names start with UP are assumed
!                 to be Pipe-like (static and total temperatures differ)
!
                  call calcgeomelemnet(vold,co,prop,lakon,nelem,
     &                       ttime,time,ielprop,mi,A,A2,d,l,s)
                  icase=0
               endif
!     
               if(v(3,node).eq.0.d0) then
                  xflow360=xflow*iaxial
                  call ts_calc(xflow360,Tt,pt,kappa,r,A,Ts,icase)
                  v(3,node)=Ts
               endif
!     
!     if the element is not of gaspipe or branch type,
!     total and static temperatures are equal for all endnodes
!     
            endif
         enddo
      endif
!     
!     for chambers the static temperature equals the total
!     temperature
!     
      do i=1,ntg
         if(ineighe(i).eq.-1) v(3,itg(i))=v(0,itg(i))
      enddo
!
       do i=1,ntg
         node=itg(i)
         if(v(2,node).eq.1.235711d-10) v(2,node)=0.d0
       enddo
c      write(*,*) 'initialnet '
c      do i=1,ntg
c         write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c      enddo
!     
      return
      end
      
      
