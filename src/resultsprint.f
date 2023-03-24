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
      subroutine resultsprint(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
     &     stx,ielorien,norien,orab,t1,ithermal,filab,een,iperturb,fn,
     &     nactdof,iout,vold,nodeboun,ndirboun,nboun,nmethod,ttime,
     &     xstate,
     &     epn,mi,nstate_,ener,enern,xstaten,eei,set,nset,istartset,
     &     iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
     &     nelemload,nload,ikin,ielmat,thicke,eme,emn,rhcon,nrhcon,
     &     shcon,
     &     nshcon,cocon,ncocon,ntmat_,sideload,icfd,inomat,pslavsurf,
     &     islavact,cdn,mortar,islavnode,nslavnode,ntie,islavsurf,time,
     &     ielprop,prop,veold,ne0,nmpc,ipompc,nodempc,labmpc,energyini,
     &     energy,orname,xload,itiefac,pmastsurf,springarea,tieset,
     &     ipobody,ibody,xbody,nbody,iinc)
!     
!     - stores the results in the .dat file, if requested
!     - nodal quantities at the nodes
!     - element quantities at the integration points
!     - calculates the extrapolation of element quantities to
!     the nodes (if requested for .frd output)
!     - calculates 1d/2d results for 1d/2d elements by
!     interpolation
!     
      implicit none
!     
      logical force,rfprint
!     
      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*),labmpc(*)
      character*80 orname(*)
      character*81 set(*),prset(*),tieset(3,*)
      character*87 filab(*)
!     
      integer kon(*),inum(*),iperm(20),mi(*),ielorien(mi(3),*),
     &     ipkon(*),nactdof(0:mi(2),*),nodeboun(*),compressible,
     &     nelemload(2,*),ndirboun(*),ielmat(mi(3),*),nrhcon(*),
     &     inotr(2,*),iorienloc,iflag,nload,mt,nk,ne,ithermal(*),i,
     &     norien,iperturb(*),iout,nboun,nmethod,node,nshcon(*),iinc,
     &     nfield,ndim,nstate_,nset,istartset(*),iendset(*),ialset(*),
     &     nprint,ntrans,ikin,ncocon(2,*),ntmat_,icfd,inomat(*),mortar,
     &     islavact(*),islavnode(*),nslavnode(*),ntie,islavsurf(2,*),
     &     ielprop(*),ne0,index,nmpc,ipompc(*),nodempc(3,*),nactdoh,
     &     iextrapolate,itiefac(2,*),ipobody(2,*),ibody(3,*),nbody 
!     
      real*8 co(3,*),v(0:mi(2),*),stx(6,mi(1),*),stn(6,*),cdn(6,*),
     &     qfx(3,mi(1),*),qfn(3,*),orab(7,*),fn(0:mi(2),*),
     &     t1(*),een(6,*),vold(0:mi(2),*),epn(*),thicke(mi(3),*),time,
     &     ener(2,mi(1),*),enern(*),eei(6,mi(1),*),rhcon(0:1,ntmat_,*),
     &     ttime,xstate(nstate_,mi(1),*),trab(7,*),xstaten(nstate_,*),
     &     eme(6,mi(1),*),emn(6,*),shcon(0:3,ntmat_,*),
     &     prop(*),veold(0:mi(2),*),energy(*),energyini(*),xload(2,*),
     &     pmastsurf,springarea(2,*),xbody(7,*),pslavsurf(3,*),
     &     cocon(0:6,ntmat_,*)
!     
      data iflag /3/
      data iperm /5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/
!     
      mt=mi(2)+1
      iextrapolate=0
!     
!     no print requests
!     
      if(iout.le.0) then
!     
!     2d basic dof results (displacements, temperature) are
!     calculated in each iteration, so that they are available
!     in the user subroutines
!     
        if(filab(1)(5:5).ne.' ') then
          nfield=mt
          call map3dto1d2d_v(v,ipkon,inum,kon,lakon,nfield,nk,
     &         ne,nactdof)
        endif
!     
!     the total energy should not be calculated:
!     - for non-dynamical calculations (nmethod!=4)
!     - for modal dynamics (iperturb(1)<=1)
!     - for thermal and thermomechanical calculations (ithermal(1)>1)
!     - for electromagnetic calculations (mi(2)=5)
!     
        return
      endif
!     
!     printing the increment number in the .dat file; only if
!     data output was requested AND the procedure is static/visco or
!     dynamic
!
      if((nprint.gt.0).and.((abs(nmethod).eq.1).or.(nmethod.eq.4))) then
        write(5,*)
        write(5,100) iinc
 100    format(32x,'INCREMENT',1x,i5)
        write(5,*)
      endif
!
!     output in dat file (with *NODE PRINT or *EL PRINT)
!
        call printout(set,nset,istartset,iendset,ialset,nprint,
     &     prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &     mi(1),nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &     orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,
     &     thicke,eme,islavsurf,mortar,time,ielprop,prop,veold,orname,
     &     nelemload,nload,sideload,xload,rhcon,nrhcon,ntmat_,ipobody,
     &     ibody,xbody,nbody,nmethod)
!     
!     for facial information (*section print): if forces and/or
!     moments in sections are requested, the stresses have to be
!     extrapolated from the integration points to the nodes first
!     
      do i=1,nprint
        if(prlab(i)(1:3).eq.'SOF') then
          nfield=6
          ndim=6
          if((norien.gt.0).and.(filab(3)(6:6).eq.'L')) then
            iorienloc=1
          else
            iorienloc=0
          endif
          cflag=filab(3)(5:5)
          force=.false.
!     
          call extrapolate(stx,stn,ipkon,inum,kon,lakon,nfield,nk,
     &         ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &         vold,force,ielmat,thicke,ielprop,prop)
          iextrapolate=1
          exit
        endif
      enddo
!     
      compressible=0
      call printoutface(co,ntmat_,v,
     &     cocon,ncocon,istartset,iendset,ipkon,lakon,kon,
     &     ialset,prset,ttime,nset,set,nprint,prlab,ielmat,mi,
     &     time,stn,iperturb)
!     
      call printoutcontact(co,vold,lakon,ne0,ne,pslavsurf,stx,
     &     prset,ttime,nprint,prlab,mi,ipkon,kon,springarea,
     &     time,tieset,itiefac,ntie,pmastsurf)
!     
!     interpolation in the original nodes of 1d and 2d elements
!     this operation has to be performed in any case since
!     the interpolated values may be needed as boundary conditions
!     in the next step (e.g. the temperature in a heat transfer
!     calculation as boundary condition in a subsequent static
!     step)
!     
      if(filab(1)(5:5).ne.' ') then
        nfield=mt
        cflag=filab(1)(5:5)
        force=.false.
        call map3dto1d2d(v,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,cflag,co,vold,force,mi,ielprop,prop)
      endif
!     
      if((filab(2)(1:4).eq.'NT  ').and.(ithermal(1).le.1)) then
        if(filab(2)(5:5).eq.'I') then
          nfield=1
          cflag=filab(2)(5:5)
          force=.false.
          call map3dto1d2d(t1,ipkon,inum,kon,lakon,nfield,nk,
     &         ne,cflag,co,vold,force,mi,ielprop,prop)
        endif
      endif
!     
!     check whether forces are requested in the frd-file. If so, but
!     none are requested in the .dat file, and output=2d, 
!     map3dto1d2d has to be called
!     
      if(filab(5)(1:2).eq.'RF') then
        if(filab(5)(5:5).eq.'I') then
          rfprint=.false.
          do i=1,nprint
            if(prlab(i)(1:2).eq.'RF') then
              rfprint=.true.
              exit
            endif
          enddo
          if(.not.rfprint) then
            nfield=mt
            cflag=' '
            force=.true.
            call map3dto1d2d(fn,ipkon,inum,kon,lakon,nfield,nk,
     &           ne,cflag,co,vold,force,mi,ielprop,prop)
          endif
        endif
      endif
!     
!     for composites:
!     interpolation of the displacements and temperatures
!     from the expanded nodes to the layer nodes
!     
      if(mi(3).gt.1) then
        if((filab(1)(1:3).eq.'U  ').or.
     &       ((filab(2)(1:4).eq.'NT  ').and.(ithermal(1).gt.1))) then
          nfield=mt
          call map3dtolayer(v,ipkon,kon,lakon,nfield,
     &         ne,co,ielmat,mi)
        endif
        if((filab(2)(1:4).eq.'NT  ').and.(ithermal(1).le.1)) then
          nfield=1
          call map3dtolayer(t1,ipkon,kon,lakon,nfield,
     &         ne,co,ielmat,mi)
        endif
      endif
!     
!     determining the contact differential displacements and stresses
!     in the contact nodes for output in frd format (only for face-
!     to-face penalty; for node-to-face penalty these quantities are
!     determined in the slave nodes and no extrapolation is necessary)
!     
!     This block must precede all calls to extrapolate, since the
!     field inum from extrapolatecontact.f is not correct; by a
!     subsequent call to extrapolate inum is corrected.
!     
      if((filab(26)(1:4).eq.'CONT').or.(filab(46)(1:4).eq.'PCON')) then
        if(mortar.eq.1) then
          nfield=6
          ndim=6
          cflag=filab(3)(5:5)
          force=.false.
          call extrapolatecontact(stx,cdn,ipkon,inum,kon,lakon,nfield,
     &         nk,ne,mi(1),ndim,co,cflag,vold,force,pslavsurf,
     &         islavact,islavnode,nslavnode,ntie,islavsurf,ielprop,prop,
     &         ielmat,ne0)
        endif
      endif
!     
!     determining the stresses in the nodes for output in frd format
!     
      if((filab(3)(1:4).eq.'S   ').or.(filab(18)(1:4).eq.'PHS ').or.
     &     (filab(20)(1:4).eq.'MAXS').or.
     &     (((filab(44)(1:4).eq.'EMFE').or.(filab(45)(1:4).eq.'EMFB'))
     &     .and.(ithermal(1).ne.2))) then
        nfield=6
        ndim=6
        if((norien.gt.0).and.(filab(3)(6:6).eq.'L')) then
          iorienloc=1
        else
          iorienloc=0
        endif
        cflag=filab(3)(5:5)
        force=.false.
!     
        call extrapolate(stx,stn,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
!     
      endif
!     
!     determining the total strains in the nodes for output in frd format
!     
      if((filab(4)(1:4).eq.'E   ').or.(filab(30)(1:4).eq.'MAXE'))
     &     then
        nfield=6
        ndim=6
        if((norien.gt.0).and.(filab(4)(6:6).eq.'L')) then
          iorienloc=1
        else
          iorienloc=0
        endif
        cflag=filab(4)(5:5)
        force=.false.
        call extrapolate(eei,een,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     determining the mechanical strains in the nodes for output in 
!     frd format
!     
      if(filab(32)(1:4).eq.'ME  ') then
        nfield=6
        ndim=6
        if((norien.gt.0).and.(filab(4)(6:6).eq.'L')) then
          iorienloc=1
        else
          iorienloc=0
        endif
        cflag=filab(4)(5:5)
        force=.false.
        call extrapolate(eme,emn,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     determining the plastic equivalent strain in the nodes 
!     for output in frd format
!     
      if(filab(6)(1:4).eq.'PEEQ') then
        nfield=1
        ndim=nstate_
        iorienloc=0
        cflag=filab(6)(5:5)
        force=.false.
        call extrapolate(xstate,epn,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     determining the internal energy in the nodes 
!     for output in frd format
!     
      if(filab(7)(1:4).eq.'ENER') then
        nfield=1
        ndim=2
        iorienloc=0
        cflag=filab(7)(5:5)
        force=.false.
        call extrapolate(ener,enern,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     determining the internal state variables in the nodes 
!     for output in frd format
!     
      if(filab(8)(1:4).eq.'SDV ') then
        nfield=nstate_
        ndim=nstate_
        if((norien.gt.0).and.(filab(9)(6:6).eq.'L')) then
          write(*,*) '*WARNING in results: SDV variables cannot'
          write(*,*) '         be stored in a local frame;'
          write(*,*) '         the global frame will be used'
        endif
        iorienloc=0
        cflag=filab(8)(5:5)
        force=.false.
        call extrapolate(xstate,xstaten,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     determining the heat flux in the nodes for output in frd format
!     
      if(((filab(9)(1:4).eq.'HFL ').and.(ithermal(1).gt.1)).or.
     &     ((filab(42)(1:3).eq.'ECD').and.(ithermal(1).eq.2))) then
        nfield=3
        ndim=3
        if((norien.gt.0).and.(filab(9)(6:6).eq.'L')) then
          iorienloc=1
        else
          iorienloc=0
        endif
        cflag=filab(9)(5:5)
        force=.false.
        call extrapolate(qfx,qfn,ipkon,inum,kon,lakon,nfield,nk,
     &       ne,mi(1),ndim,orab,ielorien,co,iorienloc,cflag,
     &       vold,force,ielmat,thicke,ielprop,prop)
        iextrapolate=1
      endif
!     
!     if no element quantities requested in the nodes: calculate
!     inum if nodal quantities are requested: used in subroutine frd
!     to determine which nodes are active in the model 
!     
      if((iextrapolate.eq.0).and.
     &     ((nmethod.ne.4).or.(iperturb(1).ge.2))) then
!     
        nfield=0
        ndim=0
        iorienloc=0
        cflag=filab(1)(5:5)
        call createinum(ipkon,inum,kon,lakon,nk,ne,cflag,nelemload,
     &       nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
     &       ielprop,prop)
      endif
!     
      if(ithermal(2).gt.1) then
!     
!     next section is executed if at least one step is thermal
!     or thermomechanical
!     
!     extrapolation for the network
!     -interpolation for the total pressure and temperature
!     in the middle nodes
!     -extrapolation for the mass flow in the end nodes
!     
        call networkextrapolate(v,ipkon,inum,kon,lakon,ne,mi)
!     
!     printing values for environmental film and
!     pressure nodes (these nodes are considered to be network
!     nodes)
!     
        do i=1,nload
          if((sideload(i)(3:4).ne.'FC').and.
     &         (sideload(i)(3:4).ne.'NP')) cycle
          node=nelemload(2,i)
          if(icfd.ne.0) then
            if(node.gt.0) then
              if(inomat(node).ne.0) cycle
            endif
          endif
          if((node.gt.0).and.(sideload(i)(1:1).ne.' ')) then
            if(inum(node).lt.0) cycle
            inum(node)=-1
          endif
        enddo
!     
!     printing values radiation 
!     (these nodes are considered to be network nodes, unless
!     they were already assigned to the structure)
!     
        do i=1,nload
          if((sideload(i)(3:4).ne.'CR')) cycle
          node=nelemload(2,i)
          if(icfd.ne.0) then
            if(node.gt.0) then
              if(inomat(node).ne.0) cycle
            endif
          endif
          if((node.gt.0).and.(sideload(i)(1:1).ne.' ')) then
            if(inum(node).ne.0) cycle
            inum(node)=-1
          endif
        enddo
!     
!     printing values for nodes belonging to network MPC's
!     (these nodes are considered to be network nodes)
!     
        do i=1,nmpc
          if(labmpc(i)(1:7).eq.'NETWORK') then
            index=ipompc(i)
            do
              node=nodempc(1,index)
              if(inum(node).ge.0) inum(node)=-1
              index=nodempc(3,index)
              if(index.eq.0) exit
            enddo
          endif
        enddo
!     
!     printing values of prescribed boundary conditions (these
!     nodes are considered to be network nodes)
!     
        do i=1,nboun
          node=nodeboun(i)
          if(inum(node).ne.0) cycle
          if(icfd.ne.0) then
            if(inomat(node).ne.0) cycle
          endif
          if((cflag.ne.' ').and.(ndirboun(i).eq.3)) cycle
          inum(node)=-1
        enddo
      endif
!     
      return
      end
