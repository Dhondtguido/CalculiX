!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine initialcfdfem(yy,nk,co,ne,ipkon,kon,lakon,x,y,z,xo,
     &     yo,zo,nx,ny,nz,isolidsurf,neighsolidsurf,xsolidsurf,dh,
     &     nshcon,shcon,nrhcon,rhcon,vold,vcon,ntmat_,iponoel,inoel,
     &     nsolidsurf,iturbulent,physcon,compressible,
     &     matname,inomat,mi,ithermal,dhel,jyy,ifreesurface,
     &     nbody,ipobody,ibody,xbody,depth,nodfreesurf,dgravity,xg)
!     
!     initial calculations for cfd applicatons:
!     - determine the distance from the nearest solid surface
!     (stored in yy and corresponding wall node in jyy)
!     - determing the distance from the nearest in-flow node
!     for solid surface nodes (stored in xsolidsurf)
!     - determine the adjacent element height for each node 
!     (stored in field dh)
!     - calculate the value of the conservative variables (vcon)
!     - calculate initial values for the turbulence parameters
!       (if not defined by the user in the input deck with
!        *INITIAL CONDITIONS,TYPE=TURBULENCE)
!     
      implicit none
!
      logical iniconbyuser
!     
      character*8 lakon(*)
      character*80 matname(*)
!     
      integer ne,ipkon(*),kon(*),indexe,ifaceq(8,6),ifacet(7,4),
     &     ifacew(8,5),kflag,isolidsurf(*),nsolidsurf,nope,node1,node2,
     &     nk,node,i,j,k,iponoel(*),inoel(2,*),nx(*),ny(*),index,nelem,
     &     nz(*),neighsolidsurf(*),kneigh,jyy(*),ifreesurface,
     &     nshcon(*),nrhcon(*),ntmat_,neigh,mi(*),
     &     imat,inomat(*),ithermal(*),iturbulent,
     &     iflag,konl(26),nopes,nfaces,ig,compressible,ipobody(2,*),
     &     nbody,ibody(3,*),m,neinode(3),nodedepth,nodfreesurf(*),
     &     neigh4(3,4),neigh6(3,6),neigh8(3,8)
!     
      real*8 x(*),y(*),z(*),xo(*),yo(*),zo(*),xsolidsurf(*),xsj2(3),
     &     yy(*),co(3,*),dh(*),r,cp,rho,shcon(0:3,ntmat_,*),
     &     rhcon(0:1,ntmat_,*),vold(0:mi(2),*),vcon(nk,0:mi(2)),
     &     temp,physcon(*),xtu,xkin,dvi,dhel(*),px,py,pz,
     &     xl2(3,9),xl(3,26),xs2(3,7),xi,et,ze,xsj,weight,volume,
     &     shp(4,26),shp2(7,9),factor,area,hmin,dgravity,
     &     xg(3),xbody(7,*),depth(*),cosangle,cosanglemax,dp(3),dd
!     
!     nodes belonging to the element faces
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &     1,2,4,5,9,8,12,
     &     2,3,4,6,10,9,13,
     &     1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
!
!     neighbors of a given node within the element (local numbers)
!
      data neigh4 /2,3,4,1,3,4,1,2,4,1,2,3/
      data neigh6 /2,3,4,1,3,5,1,2,6,1,5,6,2,4,6,3,4,5/
      data neigh8 /2,4,5,1,3,6,2,4,7,1,3,8,1,6,8,2,5,7,3,6,8,4,5,7/
!     
      data iflag /2/
      data kflag /2/
!     
      if(iturbulent.ne.0) then
!     
!     determining the nearest solid boundary node
!     
        kneigh=1
!     
        do i=1,nk
          yy(i)=1.d30
        enddo
!
        if(nsolidsurf.ne.0) then
          do i=1,nsolidsurf
            node=isolidsurf(i)
            x(i)=co(1,node)
            y(i)=co(2,node)
            z(i)=co(3,node)
            xo(i)=x(i)
            yo(i)=y(i)
            zo(i)=z(i)
            nx(i)=i
            ny(i)=i
            nz(i)=i
          enddo
          call dsort(x,nx,nsolidsurf,kflag)
          call dsort(y,ny,nsolidsurf,kflag)
          call dsort(z,nz,nsolidsurf,kflag)
!     
          do node=1,nk
            index=iponoel(node)
            if(index.le.0) cycle
            px=co(1,node)
            py=co(2,node)
            pz=co(3,node)
!     
!     determining the neighboring solid surface node
!     
            call near3d(xo,yo,zo,x,y,z,nx,ny,nz,px,py,pz,
     &           nsolidsurf,neigh,kneigh)
!     
            neigh=isolidsurf(neigh)
            jyy(node)=neigh
            yy(node)=dsqrt((co(1,node)-co(1,neigh))**2+
     &           (co(2,node)-co(2,neigh))**2+
     &           (co(3,node)-co(3,neigh))**2)
          enddo
!     
!     determining the distance to the nearest in-flow node for
!     solid surface nodes (middle nodes do not count as valid
!     in-flow nodes)
!     
          do i=1,nsolidsurf
            node1=isolidsurf(i)
            node2=neighsolidsurf(i)
            xsolidsurf(i)=dsqrt((co(1,node1)-co(1,node2))**2+
     &           (co(2,node1)-co(2,node2))**2+
     &           (co(3,node1)-co(3,node2))**2)
          enddo
!     
        endif
      endif
!     
!     determining the smallest height in each element
!     (volume divided by area times factor)
!     cf. calcstabletimeincvol.f       
!     
c     write(*,*) 'initialcfdfem element height'
!     
      do nelem=1,ne
        if((ipkon(nelem).lt.0).or.(lakon(nelem)(1:1).ne.'F')) cycle
        dhel(nelem)=1.d30
        indexe=ipkon(nelem)
!     
!     volume=area*height*factor
!     
        if(lakon(nelem)(4:4).eq.'8') then
          nope=8
          nopes=4
          nfaces=6
          factor=1.d0
        elseif(lakon(nelem)(4:4).eq.'4') then
          nope=4
          nopes=3
          nfaces=4
          factor=1.d0/3.d0
        elseif(lakon(nelem)(4:4).eq.'6') then
          nope=6
          nfaces=5
        endif
!     
!     computation of the coordinates of the local nodes
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
!     
        if(lakon(nelem)(4:4).eq.'8') then
          xi=0.d0
          et=0.d0
          ze=0.d0
          weight=8.d0
          call shape8h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(nelem)(4:4).eq.'4') then
          xi=0.25d0
          et=0.25d0
          ze=0.25d0
          weight=1.d0/6.d0
          call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(nelem)(4:4).eq.'6') then
          xi=1.d0/3.d0
          et=1.d0/3.d0
          ze=0.d0
          weight=1.d0
          call shape6w(xi,et,ze,xl,xsj,shp,iflag)
        endif
!     
        volume=xsj*weight
!     
!     loop over the faces
!     
        do ig=1,nfaces
          if(lakon(nelem)(4:4).eq.'6')then
            if(ig.le.2)then
              nopes=3
              factor=1.d0
            else
              nopes=4
              factor=0.5d0
            endif
          endif
!     
          if(nope.eq.8)then
            do i=1,nopes
              do j=1,3
                xl2(j,i)=co(j,konl(ifaceq(i,ig)))
              enddo
            enddo
          elseif(nope.eq.4)then
            do i=1,nopes
              do j=1,3
                xl2(j,i)=co(j,konl(ifacet(i,ig)))
              enddo
            enddo
          else
            do i=1,nopes
              do j=1,3
                xl2(j,i)=co(j,konl(ifacew(i,ig)))
              enddo
            enddo
          endif
!     
          if(nopes.eq.4)then
            xi=0.d0
            et=0.d0
            weight=4.d0
          else
            xi=1.d0/3.d0
            et=1.d0/3.d0
            weight=0.5d0
          endif
!     
          if(nopes.eq.4) then
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
          else
            call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
          endif
!     
          area=weight*dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &         xsj2(3)*xsj2(3))

          dhel(nelem)=min(dhel(nelem),(volume/(area*factor)))
!     
        enddo
!     ENDDO over sides
c     write(*,*) nelem,dhel(nelem)
      enddo
!      
      hmin=1.d30
      loop:do i=1,nk
        index=iponoel(i)
        if(index.le.0) cycle
        dh(i)=1.d30
!     
!     loop over all elements belonging to the edge node:
!     determining the minimum height
!     
        do
          if(index.le.0) exit
          nelem=inoel(1,index)
          dh(i)=min(dh(i),dhel(nelem))
          index=inoel(2,index)
        enddo
        hmin=min(hmin,dh(i))
!     
      enddo loop
!
!     shallow water equations: determine the depth at all fluid nodes
!     (element length in the direction of the gravity vector)      
!     and label the nodes which are on the free surface
!
      if(ifreesurface.eq.1) then
        dgravity=-1.d0
        do nelem=1,ne
          if((ipkon(nelem).lt.0).or.(lakon(nelem)(1:1).ne.'F')) cycle
!
!     determine the gravity vector: is assumed to be constant for
!     the complete model
!
          if(dgravity.lt.0.d0) then
            if(nbody.le.0) then
              write(*,*) '*ERROR in initialcfdfem: no gravity vector'
              write(*,*) '       defined in a shallow water calculation'
              call exit(201)
            endif
            index=nelem
            do
              j=ipobody(1,index)
              if(j.eq.0) exit
              if(ibody(1,j).eq.2) then
                dgravity=xbody(1,j)
                xg(1)=xbody(2,j)
                xg(2)=xbody(3,j)
                xg(3)=xbody(4,j)
              endif
              index=ipobody(2,index)
              if(index.eq.0) exit
            enddo
            if(dgravity.lt.0.d0) then
              write(*,*) '*ERROR in initialcfdfem: no gravity vector'
              write(*,*) '       defined in a shallow water calculation'
              call exit(201)
            endif
          endif
!
          nope=ichar(lakon(nelem)(4:4))-48
          indexe=ipkon(nelem)
!
!     loop over all nodes belonging to the element
!
          do j=1,nope
            node=kon(indexe+j)
            if(depth(node).gt.0.d0) cycle
!
!           neighboring nodes of node
!
            if(nope.eq.8) then
              do k=1,3
                neinode(k)=kon(indexe+neigh8(k,j))
              enddo
            elseif(nope.eq.6) then
              do k=1,3
                neinode(k)=kon(indexe+neigh6(k,j))
              enddo
            else
              do k=1,3
                neinode(k)=kon(indexe+neigh4(k,j))
              enddo
            endif
            cosanglemax=0.d0
!     
!     loop over all neighbors of the node within the element
!     determine the neighbor most parallel to the direction of the
!     gravity vector
!     
            do k=1,3
              do m=1,3
                dp(m)=co(m,neinode(k))-co(m,node)
              enddo
              dd=dsqrt(dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3))
              cosangle=(xg(1)*dp(1)+xg(2)*dp(2)+xg(3)*dp(3))/dd
              if(dabs(cosangle).gt.cosanglemax) then
                cosanglemax=dabs(cosangle)
                nodedepth=neinode(k)
                depth(node)=dd*cosangle/dabs(cosangle)
              endif
            enddo
!
            if(depth(node).gt.0.d0) then
              nodfreesurf(node)=1
              depth(nodedepth)=depth(node)
            else
              depth(node)=-depth(node)
              depth(nodedepth)=depth(node)
              nodfreesurf(nodedepth)=1
            endif
          enddo
        enddo
      endif
!
!     calculate conservative fields
!     
      do node=1,nk
        if(inomat(node).eq.0) cycle
        imat=inomat(node)
        temp=vold(0,node)
        call materialdata_cp_sec(imat,ntmat_,temp,shcon,nshcon,cp,
     &       physcon)
!     
!     different treatment for gases and liquids
!     
        if(compressible.eq.1) then
          if(ifreesurface.eq.0) then
            r=shcon(3,1,imat)
            if(r.lt.1.d-10) then
              write(*,*) '*ERROR in initialcfd: specific gas '
              write(*,*) 'constant for material ',matname(imat)
              write(*,*) 'is close to zero; maybe it has'
              write(*,*) 'not been defined'
              stop
            endif
            if(vold(0,node)-physcon(1).le.1.d-10) then
              write(*,*) '*ERROR in initialcfd: absolute temperature '
              write(*,*) '       is nearly zero; maybe absolute zero '
              write(*,*) '       was wrongly defined or not defined'
              write(*,*) '       at all (*PHYSICAL CONSTANTS card)'
              stop
            endif
            if(vold(4,node).le.0.d0) then
              write(*,*) '*ERROR in inicialcfd: initial pressure'
              write(*,*) '       must be strictly positive'
              stop
            endif
            rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
            vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)-vold(4,node)
          else
!
!           shallow water equations: H=0.1 fixed
!
            rho=2.d0*vold(4,node)/dgravity+depth(node)**2
            if(rho.le.0.d0) then
              write(*,*) '*ERROR in initialcfd: height  '
              write(*,*) '       is les or equal to zero; initial '
              write(*,*) '       "pressure" too low'
              stop
            endif
            rho=dsqrt(rho)
            if(ithermal(1).gt.1) then
              vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &             (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &             /2.d0)
            endif
          endif
        else
!
!         incompressible calculations
!
          call materialdata_rho(rhcon,nrhcon,imat,rho,
     &         temp,ntmat_,ithermal)
          if(ithermal(1).gt.1) then
            vcon(node,0)=rho*(cp*(temp-physcon(1))+
     &           (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &           /2.d0)
          endif
        endif
        vcon(node,4)=rho
        do k=1,3
          vcon(node,k)=rho*vold(k,node)
        enddo
      enddo
!     
!     initial conditions for turbulence parameters:
!     freestream conditions
!     
      if(iturbulent.ne.0) then
!
!     check whether initial turbulence conditions were defined
!     by user
!
        iniconbyuser=.false.
        do node=1,nk
          imat=inomat(node)
          if(imat.eq.0) cycle
          if((vold(5,node).gt.0.d0).or.(vold(6,node).gt.0.d0))
     &         iniconbyuser=.true.
          exit
        enddo
!     
        if(.not.iniconbyuser) then
          if(dabs(physcon(8)).lt.1.d-20) then
            write(*,*) '*ERROR in initialcfd: typical length of the'
            write(*,*) '       computational domain is too small'
            write(*,*) '       use *VALUES AT INFINITY to define a'
            write(*,*) '       finite value'
            stop
          endif
          xtu=10.d0*physcon(5)/physcon(8)
          xkin=10.d0**(-3.5d0)*xtu
          do node=1,nk
            imat=inomat(node)
            if(imat.eq.0) cycle
            temp=vold(0,node)
            call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
            rho=vcon(node,4)
!     
            vcon(node,5)=xkin*dvi
            vcon(node,6)=xtu*rho
!     
            vold(5,node)=vcon(node,5)/rho
            vold(6,node)=xtu
          enddo
        else
          do node=1,nk
            imat=inomat(node)
            if(imat.eq.0) cycle
!     
            rho=vcon(node,4)
!     
            vcon(node,5)=vold(5,node)*rho
            vcon(node,6)=vold(6,node)*rho
          enddo
        endif
      endif
!     
      return
      end
