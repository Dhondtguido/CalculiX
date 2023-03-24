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
!     Solve the Bresse equation for the turbulent stationary flow
!     in channels with a non-erosive bottom
!     
      subroutine initialchannel(itg,ieg,ntg,lakon,v,ipkon,kon,nflow,
     &     ikboun,nboun,prop,ielprop,ndirboun,nodeboun,xbounact,
     &     ielmat,ntmat_,shcon,nshcon,physcon,rhcon,nrhcon,ipobody,
     &     ibody,xbodyact,co,nbody,network,vold,set,istep,iit,mi,
     &     ineighe,ilboun,ttime,time,itreated,iponoel,inoel,istack,
     &     sfr,hfr,sba,hba,ndata,jumpup,jumpdo,istackb,nelemload,
     &     ixnode,iyload,nload,sideload,xloadact,cocon,ncocon,iinc,
     &     nforc,ikforc,ilforc,xforcact)
!
      implicit none
!
      character*1 mode
      character*8 lakon(*),lakonl
      character*20 sideload(*)
      character*81 set(*)
!      
      integer mi(*),itg(*),ieg(*),ntg,nflow,ipkon(*),kon(*),ikboun(*),
     &     nboun,ielprop(*),ndirboun(*),nodeboun(*),ielmat(mi(3),*),
     &     ntmat_,nshcon(*),nrhcon(*),ipobody(2,*),ibody(3,*),nbody,
     &     network,istep,iit,ineighe(*),ilboun(*),i,j,nelem,indexe,
     &     node1,node2,id,itreated(*),id1,id2,nup,index,iponoel(*),
     &     inoel(2,*),nmid,ndo,inv,nelemio,nelup,node,imat,neldo,
     &     istack(2,*),nstack,nel,ndata,jumpup(*),jumpdo(*),
     &     istackb(2,*),nstackb,nel1,nup1,nentry,newel,nelemload(2,*),
     &     ixnode(*),iyload(*),nload,kflag,ii,nelemwall,ncocon(2,*),
     &     nfield,iflag,ig,jltyp,iinc,k,m,mint2d,nope,nopes,konl(20),
     &     ifaceq(8,6),ifacet(6,4),ifacew(8,5),iloop,idof,nforc,
     &     ikforc(*),ilforc(*),iemchange
!
      real*8 v(0:mi(2),*),prop(*),xbounact(*),shcon(0:3,ntmat_,*),
     &     physcon(*),rhcon(0:1,ntmat_,*),xbodyact(7,*),co(3,*),
     &     vold(0:mi(2),*),ttime,time,xflow,g(3),dg,temp,cp,dvi,r,
     &     rho,sfr(*),hfr(*),sba(*),hba(*),epsilon,heatflux,temp1,
     &     xflow1,xflowact,xlin,xconst,areaj,cocon(0:6,ntmat_,*),
     &     dxsj2,xi,et,field(1),heatnod,heatfac,h(2),sinktemp,tvar(2),
     &     weight,xs2(3,7),xsj2(3),tl2(8),xl2(3,8),coords(3),shp2(7,8),
     &     xloadact(2,*),walltemp,xforcact(*),xquart,eps,timeend(2),
     &     dtemp,x,dx
!     
      include "gauss.f"
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
      data iflag /2/
!
      if(network.le.2) then
        write(*,*) '*ERROR: a network channel canot be used for'
        write(*,*) '        temperature calculations only'
        write(*,*)
        call exit(201)
      endif
!
      epsilon=1.d-5
!
!     identify the number of "real" channel neighbors (i.e. without
!     CHANNEL INOUT elements) per node: stored in field ineighe(*)
!
      do i=1,ntg
        ineighe(i)=0
      enddo
!
      do i=1,nflow
        nelem=ieg(i)
        if(lakon(nelem)(2:5).ne.'LICH') cycle
        if(lakon(nelem)(6:7).eq.'IO') cycle
        indexe=ipkon(nelem)
!
        node1=kon(indexe+1)
        call nident(itg,node1,ntg,id)
        ineighe(id)=ineighe(id)+1
!
        node2=kon(indexe+3)
        call nident(itg,node2,ntg,id)
        ineighe(id)=ineighe(id)+1
!
      enddo
!
!     setting all values to zero
!
      do i=1,ntg
        node=itg(i)
        v(0,i)=physcon(1)-1.d0
        do j=1,mi(2)
          v(j,i)=0.d0
        enddo
      enddo
!
!     applying the boundary conditions
!
      do j=1,nboun
        v(ndirboun(j),nodeboun(j))=xbounact(j)
      enddo
!     
!     determine the gravity vector
!     
      do j=1,3
        g(j)=0.d0
      enddo
      if(nbody.gt.0) then
        loop3: do i=1,nflow
          nelem=ieg(i)
          if(lakon(nelem)(2:5).ne.'LICH') cycle
          if(lakon(nelem)(6:7).eq.'IO') cycle
!          
          index=nelem
          do
            j=ipobody(1,index)
            if(j.eq.0) exit
            if(ibody(1,j).eq.2) then
              g(1)=g(1)+xbodyact(1,j)*xbodyact(2,j)
              g(2)=g(2)+xbodyact(1,j)*xbodyact(3,j)
              g(3)=g(3)+xbodyact(1,j)*xbodyact(4,j)
              exit loop3
            endif
            index=ipobody(2,index)
            if(index.eq.0) exit
          enddo
        enddo loop3
      endif
      dg=dsqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3))
!
      nstack=0
      nstackb=0
!
!     mechanical computations
!
!     major loop: looking for SLUICE GATE and WEAR elements
!
      loop1: do i=1,nflow
        nelem=ieg(i)
        if(((lakon(nelem)(6:7).ne.'SG').and.
     &     (lakon(nelem)(6:7).ne.'WE')).or.
     &       (itreated(i).eq.1)) then
          if(i.lt.nflow) then
            cycle
!
!         if no more branches starting with a SLUICE GATE or WEAR
!         element: look for non-treated backwater curves starting
!         upstream of joints
!
          elseif(nstackb.gt.0) then
            neldo=istackb(1,nstackb)
            ndo=istackb(2,nstackb)
            nstackb=nstackb-1
            nelem=0
            xflow=dabs(v(1,kon(ipkon(neldo)+2)))
            mode='B'
          else
            exit
          endif
        else
!
!       untreated SLUICE GATE or WEAR element found
!
          indexe=ipkon(nelem)
          node1=kon(indexe+1)
          call nident(itg,node1,ntg,id1)
          node2=kon(indexe+3)
          call nident(itg,node2,ntg,id2)
!
!         the SLUICE GATE or WEAR element should be connected on one side to
!         a CHANNEL INOUT element (as only element)
!
          if((ineighe(id1).gt.1).and.(ineighe(id2).gt.1)) cycle
!
!         new branch found
!
!         determine the upstream node nup of the element
!
          if(ineighe(id1).eq.1) then
            nup=node1
            inv=1
          else
            nup=node2
            inv=-1
          endif
!
!         determine the upstream element nelup
!
          index=iponoel(nup)
          do
            if(index.eq.0) then
              write(*,*) '*ERROR: node',nup
              write(*,*) '        is only connected to one element'
              write(*,*) '        with number',nelem
              write(*,*)
              call exit(201)
            endif
            if(inoel(1,index).ne.nelem) then
              nelup=inoel(1,index)
              exit
            else
              index=inoel(2,index)
            endif
          enddo
!
!         define mode to be "forward"
!
          mode='F'
!
!         mass flow is taken from the IO element upstream of the sluice gate
!
          xflow=dabs(v(1,kon(ipkon(nelup)+2)))
        endif
!
!       loop over all elements in present branch
!
        loop2: do
!
!         if F (forward): nelup and nup known
!
          if(mode.eq.'F') then
            if(nelem.eq.0) then
              call nident(itg,nup,ntg,id)
!
!             end of branch
!
              if(ineighe(id).eq.1) then
                write(*,*) '*INFO: branch finished'
                write(*,*)
!     
!     IO-element: determine the mass flow
!     
                index=iponoel(nup)
                do
                  if(inoel(1,index).ne.nelup) then
                    nelemio=inoel(1,index)
                    if(nup.eq.kon(ipkon(nelemio)+1)) then
                      v(1,kon(ipkon(nelemio)+2))=xflow
                    else
                      v(1,kon(ipkon(nelemio)+2))=-xflow
                    endif
                    cycle loop1
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) then
                    write(*,*) '*ERROR in initialchannel: no IO'
                    write(*,*) '       element at end of branch'
                    write(*,*) '       most downstream node: ',nup
                    call exit(201)
                  endif
                enddo
!
!               branch continues 
!
              elseif(ineighe(id).eq.2) then
!     
!     one "true" element connected downstream
!     loop over all elements connected to nup
!     
                index=iponoel(nup)
                do
                  if(inoel(1,index).ne.nelup) then
                    if(lakon(inoel(1,index))(6:7).ne.'IO') then
                      nelem=inoel(1,index)
                    else
!     
!     add flow
!     
                      nelemio=inoel(1,index)
                      if((lakon(nelup)(6:7).eq.'SG').or.
     &                     (lakon(nelup)(6:7).eq.'WE')) then
                        write(*,*)
     &                       '*ERROR in initialchannel: no IO element'
                        write(*,*)
     &                       '       allowed immediately downstream'
                        write(*,*) '       of a SLUICE GATE or WEAR'
                        write(*,*) '       element; faulty element:',
     &                       nelemio
                        write(*,*) 
                        call exit(201)
                      endif
                      if(nup.eq.kon(ipkon(nelemio)+3)) then
                        xflow=xflow+v(1,kon(ipkon(nelemio)+2))
                      else
                        xflow=xflow-v(1,kon(ipkon(nelemio)+2))
                      endif
                    endif
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) exit
                enddo
!
!               joint of three channels
!
              elseif(ineighe(id).eq.3) then
!     
!               taking the temperature of the upstream node for the
!               material properties
!     
                temp=v(0,nup)
                imat=ielmat(1,nelup)
!     
                call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &               dvi,rhcon,nrhcon,rho)
!     
                call channeljointfront(nelem,nelup,nup,iponoel,inoel,
     &               ielprop,prop,ipkon,kon,mi,v,g,dg,nstackb,istackb,
     &               rho,xflow,co,lakon)
!
!               if nelem=0 more than one branch has no mass flux
!
                if(nelem.eq.0) then
                  cycle loop1
                endif
!
!               joint of more than three channels
!
              elseif(ineighe(id).gt.3) then
                write(*,*) '*ERROR in initialchannel: branch joint'
                write(*,*) '       of more than 3 channels'
                write(*,*)
                call exit(201)
              endif
!
            endif
!     
!     actual element = nelem
!     determining the middle and downstream node nmid and ndo
!     
!     if the actual flow is from kon(indexe+1) to kon(indexe+3) then
!     inv=1
!     if the actual flow is from kon(indexe+3) to kon(indexe+1) then
!     inv=-1
!     
            indexe=ipkon(nelem)
            nmid=kon(indexe+2)
            if(kon(indexe+1).eq.nup) then
              ndo=kon(indexe+3)
              inv=1
            else
              ndo=kon(indexe+1)
              inv=-1
            endif
          else
!     
!     B (backward): either nelem and ndo known or
!                   neldo and ndo known
!
            if(nelem.eq.0) then
              call nident(itg,ndo,ntg,id)
              if(ineighe(id).eq.1) then
!
!              IO-element: determine the mass flow
!     
                index=iponoel(ndo)
                do
                  if(inoel(1,index).ne.neldo) then
                    nelemio=inoel(1,index)
                    if(ndo.eq.kon(ipkon(nelemio)+3)) then
                      v(1,kon(ipkon(nelemio)+2))=xflow
                    else
                      v(1,kon(ipkon(nelemio)+2))=-xflow
                    endif
                    exit
                  endif
                  index=inoel(2,index)
                  if(index.eq.0) exit
                enddo
!
                if(nstack.gt.0) then
                  mode='F'
                  nelup=istack(1,nstack)
                  nup=istack(2,nstack)
                  nelem=0
                  xflow=dabs(v(1,kon(ipkon(nelup)+2)))
                  nstack=nstack-1
                  cycle loop2
                else
                  write(*,*) '*INFO: branch finished'
                  write(*,*)
                  cycle loop1
                endif
!
!               joint of three channels
!
              elseif(ineighe(id).eq.3) then
!     
!               taking the temperature of the downstream node for the
!               material properties
!     
                temp=v(0,ndo)
                imat=ielmat(1,neldo)
!     
                call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &               dvi,rhcon,nrhcon,rho)
!     
                call channeljointback(neldo,ndo,iponoel,inoel,ipkon,kon,
     &               mi,v,istackb,nstackb,co,ielprop,prop,g,dg,xflow,
     &               rho)
!
!               joint of more than three channels
!
              elseif(ineighe(id).gt.3) then
                write(*,*) '*ERROR in initialchannel: branch joint'
                write(*,*) '       of more than 3 channels'
                write(*,*)
                call exit(201)
              endif
!     
!     if nelem is zero, neldo is known and nelem has to be
!     determined
!     
              index=iponoel(ndo)
              do
                if(inoel(1,index).ne.neldo) then
                  if(lakon(inoel(1,index))(6:7).ne.'IO') then
                    nelem=inoel(1,index)
                  else
!     
!     add flow
!     
                    nelemio=inoel(1,index)
                    if(lakon(nelup)(6:7).eq.'RE') then
                      write(*,*)
     &                     '*ERROR in initialchannel: no IO element'
                      write(*,*) '       allowed immediately upstream'
                      write(*,*) '       of a RESERVOIR element.'
                      write(*,*) '       faulty element:',nelemio
                      write(*,*) 
                      call exit(201)
                    endif
                    if(kon(ipkon(nelemio)+3).eq.ndo) then
                      xflow=xflow-v(1,kon(ipkon(nelemio)+2))
                    else
                      xflow=xflow+v(1,kon(ipkon(nelemio)+2))
                    endif
                  endif
                endif
                index=inoel(2,index)
                if(index.eq.0) exit
              enddo
            endif
!
!           determine nmid and nup
!
            indexe=ipkon(nelem)
            nmid=kon(indexe+2)
            if(kon(indexe+3).eq.ndo) then
              nup=kon(indexe+1)
              inv=1
            else
              nup=kon(indexe+3)
              inv=-1
            endif
!     
!           determining the upstream element (needed in istack
!           in certain cases
!     
            nelup=0
            index=iponoel(nup)
            do
              if(index.eq.0) exit
              if(inoel(1,index).ne.nelem) then
                if(lakon(inoel(1,index))(6:7).ne.'IO') then
                  nelup=inoel(1,index)
                  exit
                endif
              endif
              index=inoel(2,index)
            enddo
          endif
!
!         taking the temperature of the upstream node for the
!         material properties
!
          temp=v(0,nup)
          imat=ielmat(1,nelem)
!     
          call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &         dvi,rhcon,nrhcon,rho)
!     
!     treating the specific element type
!     
          call nident(ieg,nelem,nflow,nel)
          if((lakon(nelem)(6:7).eq.'SG').or.
     &       (lakon(nelem)(6:7).eq.'WE')) then
            call sluicegate(nelem,ielprop,prop,nup,nmid,ndo,co,g,dg,
     &           mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,
     &           nboun,mi,v,ipkon,kon,inv,epsilon,lakon)
          elseif(lakon(nelem)(6:7).eq.'  ') then
            call straightchannel(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &           dg,mode,xflow,rho,dvi,nelup,neldo,istack,nstack,ikboun,
     &           nboun,mi,v,ipkon,kon,ndata,nel,sfr((nel-1)*ndata+1),
     &           hfr((nel-1)*ndata+1),sba((nel-1)*ndata+1),
     &           hba((nel-1)*ndata+1),jumpup,jumpdo,inv,epsilon)
          elseif(lakon(nelem)(6:7).eq.'RE') then
            call reservoir(nelem,ielprop,prop,nup,nmid,ndo,co,g,
     &           dg,mode,xflow,rho,dvi,nelup,mi,v,inv,epsilon,istack,
     &           nstack)
          elseif((lakon(nelem)(6:7).eq.'CO').or.
     &           (lakon(nelem)(6:7).eq.'EL').or.
     &           (lakon(nelem)(6:7).eq.'ST')) then
c     &           (lakon(nelem)(6:7).eq.'DR')) then
            call contraction(nelem,ielprop,prop,nup,nmid,ndo,dg,
     &           mode,xflow,rho,nelup,neldo,istack,nstack,
     &           mi,v,inv,epsilon,co)
          else
            write(*,*) '*ERROR in initialchannel:'
            write(*,*) '       element of type ',lakon(nelem)
            write(*,*) '       is not known'
            write(*,*)
            call exit(201)
          endif
          itreated(nel)=1
        enddo loop2
      enddo loop1
!
c      i=1
c      if(i.eq.1) return
!     
!     thermal computations
!     
      timeend(1)=time
      timeend(2)=ttime+time
!
!     storing the distributed facial loads and corresponding
!     reference nodes
!
      do i=1,nload
        ixnode(i)=nelemload(2,i)
        iyload(i)=i
      enddo
!
!     sorting ixnode and iyload along
!
      kflag=2
      call isortii(ixnode,iyload,nload,kflag)
!
!     major loop: looking for SLUICE GATE and WEAR elements
!
      do i=1,nflow
        itreated(i)=0
      enddo
!
      loop4: do i=1,nflow
        nelem=ieg(i)
        if(((lakon(nelem)(6:7).ne.'SG').and.
     &       (lakon(nelem)(6:7).ne.'WE')).or.
     &       (itreated(i).eq.1)) cycle
!
!       untreated SLUICE GATE or WEAR element found
!
        indexe=ipkon(nelem)
        node1=kon(indexe+1)
        call nident(itg,node1,ntg,id1)
        node2=kon(indexe+3)
        call nident(itg,node2,ntg,id2)
!     
!       the SLUICE GATE or WEAR element should be connected on one side to
!       a CHANNEL INOUT element (as only element)
!     
        if((ineighe(id1).gt.1).and.(ineighe(id2).gt.1)) cycle
!     
!       new branch found
!       the temperature of the upstream node of the SLUICE GATE
!       or WEAR element should be known (boundary condition)
!     
!       determine the upstream node nup of the element downstream
!       of the SLUICE GATE or WEAR element
!     
        if(ineighe(id1).eq.1) then
          nup=node2
        else
          nup=node1
        endif
!     
!     determine the upstream element nelup
!     
        nelup=nelem
!
!       loop over all elements in present branch
!
        loop5: do
!
!         nelup and nup known
!     
!         mass flow in nelup
!
          xflow=dabs(v(1,kon(ipkon(nelup)+2)))
!     
!         heat flux flow in nelup
!
          if(nup.eq.kon(ipkon(nelup)+3)) then
            nentry=kon(ipkon(nelup)+1)
          else
            nentry=kon(ipkon(nelup)+3)
          endif
          temp=v(0,nentry)
          if(temp-physcon(1).lt.0.d0) then
            write(*,*) '*WARNING in initialchannel: no thermal'
            write(*,*) '         calculation can be performed'
            write(*,*) '         since the temperature in upstream'
            write(*,*) '         node ',nentry,' is lacking'
            write(*,*)
            exit loop4
          endif
          imat=ielmat(1,nelup)
          call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &         dvi,rhcon,nrhcon,rho)
          heatflux=cp*(temp-physcon(1))*xflow
!
!         starting value for the temperature in nup
!
          v(0,nup)=temp
!     
          call nident(itg,nup,ntg,id)
!
!         end of branch
!
          if(ineighe(id).eq.1) then
            write(*,*) '*INFO: branch finished'
            write(*,*)
            v(0,nup)=heatflux/(cp*xflow)+physcon(1)
            exit loop5
!
!           branch continues 
!
          elseif(ineighe(id).eq.2) then
!     
!           one "true" element connected downstream
!           loop over all elements connected to nup
!     
            index=iponoel(nup)
            do
              if(inoel(1,index).ne.nelup) then
                if(lakon(inoel(1,index))(6:7).ne.'IO') then
!     
!                 actual element
!     
                  nelem=inoel(1,index)
                else
!     
!                 IO element
!     
                  nelemio=inoel(1,index)
!
!                 mass flow
!
                  if(nup.eq.kon(ipkon(nelemio)+3)) then
                    xflowact=v(1,kon(ipkon(nelemio)+2))
                  else
                    xflowact=-v(1,kon(ipkon(nelemio)+2))
                  endif
                  xflow=xflow+xflowact
!
!                 heat flux
!
!                 the temperature is exceptionally stored in
!                 the middle node of the IO-element, since the
!                 upstream node has node number zero
!
                  nentry=kon(ipkon(nelemio)+2)
                  temp=v(0,nentry)
                  if(temp-physcon(1).lt.0.d0) then
                    write(*,*) '*WARNING in initialchannel: no thermal'
                    write(*,*) '         calculation can be performed'
                    write(*,*) '         since the temperature at '
                    write(*,*) '         entrance node ',nentry,
     &                   ' is lacking'
                    exit loop4
                  endif
                  imat=ielmat(1,nelemio)
                  call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,
     &                 r,dvi,rhcon,nrhcon,rho)
                  heatflux=heatflux+cp*(temp-physcon(1))*xflowact
                endif
              endif
              index=inoel(2,index)
              if(index.eq.0) exit
            enddo
!     
!           joint of three channels
!
          elseif(ineighe(id).eq.3) then
!     
            nel1=0
            index=iponoel(nup)
            do
              if(inoel(1,index).eq.nelup) then
                index=inoel(2,index)
                if(index.eq.0) exit
                cycle
              endif
!
!             element not equal to nelup
!
              newel=inoel(1,index)
              indexe=ipkon(newel)
!
!             temperature at the end node of element newel which
!             is not equal to nup
!
              if(kon(indexe+1).eq.nup) then
                temp1=v(0,kon(indexe+3))
              else
                temp1=v(0,kon(indexe+1))
              endif
!
!             if absolute temperature is negative: element not
!             treated yet
!
              if(temp1-physcon(1).lt.0.d0) then
                nelem=newel
                index=inoel(2,index)
                if(index.eq.0) exit
                cycle
              endif
!
!             temperature is not negative in new branch: must be
!             branch 1
!
              nel1=newel
              imat=ielmat(1,nel1)
              call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,
     &             r,dvi,rhcon,nrhcon,rho)
              if(kon(indexe+1).eq.nup) then
                nup1=kon(indexe+3)
                xflow1=-v(1,kon(indexe+2))
              else
                nup1=kon(indexe+1)
                xflow1=v(1,kon(indexe+2))
              endif
              xflow=xflow+xflow1
              heatflux=heatflux+cp*(temp1-physcon(1))*xflow1
!     
              index=inoel(2,index)
              if(index.eq.0) exit
            enddo
!     
!           if nel1=0 more than one branch has no mass flux
!
            if(nel1.eq.0) then
              nelem=0
              cycle loop4
            endif
!     
!           joint of more than three channels
!
          elseif(ineighe(id).gt.3) then
            write(*,*) '*ERROR in initialchannel: branch joint'
            write(*,*) '       of more than 3 channels'
            write(*,*)
            call exit(201)
          endif
!
!         iteration loop for the temperature in nup;
!         the linearity is very mild (material properties)
!
          iloop=0
          do
            iloop=iloop+1
            xquart=0.d0
!     
!           taking the temperature of the upstream node for the
!           material properties
!
            temp=v(0,nup)
            imat=ielmat(1,nelem)
!     
            call materialdata_tg(imat,ntmat_,temp,shcon,nshcon,cp,r,
     &           dvi,rhcon,nrhcon,rho)
!     
            xconst=heatflux
            xlin=cp*xflow
!     
!           check for convection with walls
!     
            call nident(ixnode,nup,nload,id)
!     
            if(id.gt.0) then
              do
                if(ixnode(id).eq.nup) then
                  ii=iyload(id)
!     
!                 forced convection or radiation  
!     
                  if((sideload(ii)(3:4).eq.'FC').or.
     &               (sideload(ii)(1:1).eq.'R')) then
                    nelemwall=nelemload(1,ii)
                    index=ipkon(nelemwall)
                    if(index.lt.0) cycle
                    lakonl=lakon(nelemwall)
!     
!                   calculate the area
!     
                    read(sideload(ii)(2:2),'(i1)') ig
!     
!                   number of nodes and integration points in the face
!     
                    if(lakonl(4:4).eq.'2') then
                      nope=20
                      nopes=8
                    elseif(lakonl(4:4).eq.'8') then
                      nope=8
                      nopes=4
                    elseif(lakonl(4:5).eq.'10') then
                      nope=10
                      nopes=6
                    elseif(lakonl(4:4).eq.'4') then
                      nope=4
                      nopes=3
                    elseif(lakonl(4:5).eq.'15') then
                      nope=15
                    else
                      nope=6
                    endif
!     
                    if(lakonl(4:5).eq.'8R') then
                      mint2d=1
                    elseif((lakonl(4:4).eq.'8').or.
     &                     (lakonl(4:6).eq.'20R')) then
                      if(lakonl(7:7).eq.'A') then
                        mint2d=2
                      else
                        mint2d=4
                      endif
                    elseif(lakonl(4:4).eq.'2') then
                      mint2d=9
                    elseif(lakonl(4:5).eq.'10') then
                      mint2d=3
                    elseif(lakonl(4:4).eq.'4') then
                      mint2d=1
                    endif
!     
                    if(lakonl(4:4).eq.'6') then
                      mint2d=1
                      if(ig.le.2) then
                        nopes=3
                      else
                        nopes=4
                      endif
                    endif
                    if(lakonl(4:5).eq.'15') then
                      if(ig.le.2) then
                        mint2d=3
                        nopes=6
                      else
                        mint2d=4
                        nopes=8
                      endif
                    endif
!     
!                   connectivity of the element
!     
                    do k=1,nope
                      konl(k)=kon(index+k)
                    enddo
!     
!                   coordinates of the nodes belonging to the face
!     
                    if((nope.eq.20).or.(nope.eq.8)) then
                      do k=1,nopes
                        tl2(k)=v(0,konl(ifaceq(k,ig)))
                        do j=1,3
                          xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &                         v(j,konl(ifaceq(k,ig)))
                        enddo
                      enddo
                    elseif((nope.eq.10).or.(nope.eq.4)) then
                      do k=1,nopes
                        tl2(k)=v(0,konl(ifacet(k,ig)))
                        do j=1,3
                          xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &                         v(j,konl(ifacet(k,ig)))
                        enddo
                      enddo
                    else
                      do k=1,nopes
                        tl2(k)=v(0,konl(ifacew(k,ig)))
                        do j=1,3
                          xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &                         v(j,konl(ifacew(k,ig)))
                        enddo
                      enddo
                    endif
!     
!                   integration to obtain the area and the mean
!                   temperature
!     
                    do m=1,mint2d
                      if((lakonl(4:5).eq.'8R').or.
     &                     ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                        xi=gauss2d1(1,m)
                        et=gauss2d1(2,m)
                        weight=weight2d1(m)
                      elseif((lakonl(4:4).eq.'8').or.
     &                       (lakonl(4:6).eq.'20R').or.
     &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                        xi=gauss2d2(1,m)
                        et=gauss2d2(2,m)
                        weight=weight2d2(m)
                      elseif(lakonl(4:4).eq.'2') then
                        xi=gauss2d3(1,m)
                        et=gauss2d3(2,m)
                        weight=weight2d3(m)
                      elseif((lakonl(4:5).eq.'10').or.
     &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                        xi=gauss2d5(1,m)
                        et=gauss2d5(2,m)
                        weight=weight2d5(m)
                      elseif((lakonl(4:4).eq.'4').or.
     &                     ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
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
                      else
                        call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                      endif
!     
                      dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                     xsj2(3)*xsj2(3))
                      areaj=dxsj2*weight
!     
                      walltemp=0.d0
                      do k=1,3
                        coords(k)=0.d0
                      enddo
                      do j=1,nopes
                        walltemp=walltemp+tl2(j)*shp2(4,j)
                        do k=1,3
                          coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                        enddo
                      enddo
!
                      if(sideload(ii)(1:1).eq.'F') then
                        if(sideload(ii)(5:6).ne.'NU') then
                          h(1)=xloadact(1,ii)
                        else
                          read(sideload(ii)(2:2),'(i1)') jltyp
                          jltyp=jltyp+10
                          call film(h,temp,walltemp,istep,iinc,tvar,
     &                         nelemwall,m,coords,jltyp,field,nfield,
     &                         sideload(ii),node,areaj,v,mi,ipkon,
     &                         kon,lakon,iponoel,inoel,ielprop,prop,
     &                         ielmat,shcon,nshcon,rhcon,nrhcon,ntmat_,
     &                         cocon,ncocon,ipobody,xbodyact,ibody,
     &                         heatnod,heatfac)
                        endif
                        xconst=xconst+h(1)*(walltemp-physcon(1))*areaj
                        xlin=xlin+h(1)*areaj
                      else
                        if(sideload(ii)(5:6).ne.'NU') then
                          eps=xloadact(1,ii)
                        else
                          read(sideload(ii)(2:2),'(i1)') jltyp
                          jltyp=jltyp+10
                          call radiate(eps,temp,walltemp,
     &                         istep,iinc,timeend,nelem,i,coords,jltyp,
     &                         field,nfield,sideload(id),node,areaj,
     &                         vold,mi,iemchange)
                        endif
                        xconst=xconst+physcon(2)*eps*
     &                       (walltemp-physcon(1))**4*areaj
                        xquart=xquart+physcon(2)*eps*areaj
c                        write(*,*) 'xquart',xquart,areaj
                       endif
                    enddo
                  endif
                else
                  exit
                endif
                id=id-1
                if(id.eq.0) exit
              enddo
            endif
!     
!           check for heat sources in node "nup" (energy/time)
!     
            idof=8*(nup-1)
            call nident(ikforc,idof,nforc,id)
            if(id.gt.0) then
              if(ikforc(id).eq.idof) then
                xconst=xconst+xforcact(ilforc(id))
              endif
            endif
!
!           determining the temperature
!
            if(xquart.eq.0.d0) then
!
!             no radiation: solving linear equation
!             xlin*x=xconst
!             temperatures in the equation are absolute temperatures
!
              v(0,nup)=xconst/xlin+physcon(1)
              if((dabs(xconst/xlin-(temp-physcon(1))).lt.1.d-3).or.
     &             (dabs(xconst/xlin-(temp-physcon(1))).lt.
     &             1.d-3*(temp-physcon(1))))
     &             then
                exit
              endif
            else
!
!             radiation: solving quartic equation
!             xquart*x**4+xlin*x=xconst
!             temperatures in the equation are absolute temperatures
!
              x=temp-physcon(1)
              do
                dx=-(xquart*x**4+xlin*x-xconst)/
     &               (4.d0*xquart*x**3+xlin)
c                write(*,*) iloop,x,dx
c                v(0,nup)=x+dx+physcon(1)
                if((dabs(dx).lt.1.d-3).or.
     &               (dabs(dx).lt.1.d-3*x)) then
                  exit
                endif
                x=x+dx
              enddo
              dtemp=x-(temp-physcon(1))
              v(0,nup)=x+physcon(1)
              if((dabs(dtemp).lt.1.d-3).or.
     &             (dabs(dtemp).lt.1.d-3*(temp-physcon(1))))
     &             then
                exit
              endif
c              write(*,*) 'temp ',iloop,temp,dtemp,xconst,xlin,xquart
            endif
!
            if(iloop.gt.100) then
              write(*,*) '*ERROR in initialchannel:'
              write(*,*) '       temperature does not converge.'
              call exit(201)
            endif
!     
          enddo
!     
!         determining new nelup and nup
!     
          indexe=ipkon(nelem)
          if(kon(indexe+1).eq.nup) then
            nup=kon(indexe+3)
          else
            nup=kon(indexe+1)
          endif
          nelup=nelem
!     
          call nident(ieg,nelem,nflow,nel)
          itreated(nel)=1
        enddo loop5
      enddo loop4
!
      return
      end
      
