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
      subroutine cyclicsymmetrymodels(inpc,textpart,set,
     &     istartset,iendset,
     &     ialset,nset,tieset,tietol,co,nk,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,nr,nz,
     &     rcs0,zcs0,ncs_,cs,labmpc,istep,istat,n,iline,ipol,inl,
     &     ipoinp,inp,ntie,mcs,lprev,ithermal,rcscg,rcs0cg,zcscg,
     &     zcs0cg,nrcg,nzcg,jcs,kontri,straight,ne,ipkon,kon,
     &     lakon,lcs,ifacetet,inodface,ipoinpc,maxsectors,
     &     trab,ntrans,ntrans_,jobnamec,vold,nef,mi,iaxial,ier)
!     
!     reading the input deck: *CYCLIC SYMMETRY MODEL
!     
!     several cyclic symmetry parts can be defined for one and the
!     same model; for each part there must be a *CYCLIC SYMMETRY MODEL
!     card
!     
!     cs(1,mcs): # segments in 360 degrees
!     cs(2,mcs): minimum node diameter
!     cs(3,mcs): maximum node diameter
!     cs(4,mcs): # nodes on the independent side (for fluids: upper bound)
!     cs(5,mcs): # sectors to be plotted
!     cs(6..8,mcs): first point on cyclic symmetry axis
!     cs(9..11,mcs): second point on cylic symmetry axis; turning
!     the slave surface clockwise about the cyclic symmetry axis
!     while looking from the first point to the second point one
!     arrives at the master surface without leaving the body
!     cs(12,mcs): -1 (denotes a cylindrical coordinate system)
!     cs(13,mcs): number of the element set
!     cs(14,mcs): sum of previous independent nodes
!     cs(15,mcs): cos(angle); angle = 2*pi/cs(1,mcs)
!     cs(16,mcs): sin(angle)
!     cs(17,mcs): number of tie constraint
!     
!     notice that in this routine ics, zcs and rcs start for 1 for
!     each *cyclic symmetry model card (look at the pointer in
!     the cyclicsymmetrymodels call in calinput.f)
!     
      implicit none
!     
      logical triangulation,calcangle,nodesonaxis,check,exist
!     
      character*1 inpc(*),depkind,indepkind
      character*8 lakon(*)
      character*20 labmpc(*)
      character*80 tie
      character*81 set(*),depset,indepset,tieset(3,*),elset
      character*132 textpart(16),jobnamec(*)
      character*256 fn
!     
      integer istartset(*),iendset(*),ialset(*),ipompc(*),ifaces,
     &     nodempc(3,*),itiecyc,ntiecyc,iaxial,nopes,nelems,m,indexe,
     &     nset,istep,istat,n,key,i,j,k,nk,nmpc,nmpc_,mpcfree,ics(*),
     &     nr(*),nz(*),jdep,jindep,l,noded,ikmpc(*),ilmpc(*),lcs(*),
     &     kflag,node,ncsnodes,ncs_,iline,ipol,inl,ipoinp(2,*),nneigh,
     &     inp(3,*),itie,iset,ipos,mcs,lprev,ntie,ithermal(*),
     &     nrcg(*),nzcg(*),jcs(*),kontri(3,*),ne,ipkon(*),kon(*),nodei,
     &     ifacetet(*),inodface(*),ipoinpc(0:*),maxsectors,id,jfaces,
     &     noden(2),ntrans,ntrans_,nef,mi(*),ifaceq(8,6),ifacet(6,4),
     &     ifacew1(4,5),ifacew2(8,5),idof,ier,icount,nodeaxd,nodeaxi,
     &     ilen
!     
      real*8 tolloc,co(3,*),coefmpc(*),rcs(*),zcs(*),rcs0(*),zcs0(*),
     &     csab(7),xn,yn,zn,dd,xap,yap,zap,tietol(4,*),cs(17,*),
     &     gsectors,x3,y3,z3,phi,rcscg(*),rcs0cg(*),zcscg(*),zcs0cg(*),
     &     straight(9,*),x1,y1,z1,x2,y2,z2,zp,rp,dist,trab(7,*),rpd,zpd,
     &     vold(0:mi(2),*),calculated_angle,user_angle,xsectors,
     &     axdistmin
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
      if(istep.gt.0) then
        write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '       *CYCLIC SYMMETRY MODEL should'
        write(*,*) '       be placed before all step definitions'
        write(*,*)
        ier=1
        return
      endif
!     
      check=.true.
      gsectors=1
      elset='
     &'
      tie='
     &'
      do i=2,n
        if(textpart(i)(1:2).eq.'N=') then
          read(textpart(i)(3:22),'(f20.0)',iostat=istat) xsectors
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC SYMMETRY MODEL%",ier)
            return
          endif
        elseif(textpart(i)(1:8).eq.'CHECK=NO') then
          check=.false.
        elseif(textpart(i)(1:7).eq.'NGRAPH=') then
          read(textpart(i)(8:27),'(f20.0)',iostat=istat) gsectors
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC SYMMETRY MODEL%",ier)
            return
          endif
        elseif(textpart(i)(1:4).eq.'TIE=') then
          read(textpart(i)(5:84),'(a80)',iostat=istat) tie
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC SYMMETRY MODEL%",ier)
            return
          endif
        elseif(textpart(i)(1:6).eq.'ELSET=') then
          read(textpart(i)(7:86),'(a80)',iostat=istat) elset
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CYCLIC SYMMETRY MODEL%",ier)
            return
          endif
          elset(81:81)=' '
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
        else
          write(*,*) 
     &         '*WARNING reading *CYCLIC SYMMETRY MODEL:'
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*CYCLIC SYMMETRY MODEL%")
        endif
      enddo
!     
      mcs=mcs+1
      cs(2,mcs)=-0.5d0
      cs(3,mcs)=-0.5d0
      cs(14,mcs)=lprev+0.5d0
!     
!     determining the tie constraint
!     
      itie=0
      ntiecyc=0
      do i=1,ntie
        if((tieset(1,i)(1:80).eq.tie).and.
     &         (tieset(1,i)(81:81).ne.'C').and.
     &         (tieset(1,i)(81:81).ne.'T').and.
     &         (tieset(1,i)(81:81).ne.'M').and.
     &         (tieset(1,i)(81:81).ne.'S').and.
     &         (tieset(1,i)(81:81).ne.'D')) then
          itie=i
          exit
        elseif((tieset(1,i)(81:81).ne.'C').and.
     &         (tieset(1,i)(81:81).ne.'T').and.
     &         (tieset(1,i)(81:81).ne.'M').and.
     &         (tieset(1,i)(81:81).ne.'S').and.
     &         (tieset(1,i)(81:81).ne.'D')) then
          ntiecyc=ntiecyc+1
          itiecyc=i
        endif
      enddo
      if(itie.eq.0) then
        if(ntiecyc.eq.1) then
          itie=itiecyc
        else
          write(*,*)
     &         '*ERROR reading *CYCLIC SYMMETRY MODEL:'
          write(*,*) '       tie constraint is nonexistent'
          call inputerror(inpc,ipoinpc,iline,
     &         "*CYCLIC SYMMETRY MODEL%",ier)
          return
        endif
      endif
!     
      if(xsectors.le.0) then
        write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '       the required parameter N'
        write(*,*) '       is lacking on the *CYCLIC SYMMETRY MODEL'
        write(*,*) '       keyword card or has a value <=0'
        write(*,*)
        ier=1
        return
      endif
      if(gsectors.lt.1) then
        write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '         cannot plot less than'
        write(*,*) '         one sector: one sector will be plotted'
        write(*,*)
        gsectors=1
      endif
      if(gsectors.gt.xsectors) then
        write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '         cannot plot more than'
        write(*,*) '         ',xsectors,'sectors;',
     &       xsectors,' sectors will'
        write(*,*) '       be plotted'
        write(*,*)
        gsectors=xsectors
      endif
!     
      maxsectors=max(maxsectors,int(xsectors+0.5d0))
!     
      cs(1,mcs)=xsectors
      cs(5,mcs)=gsectors+0.5d0
      cs(17,mcs)=itie+0.5d0
      depset=tieset(2,itie)
      indepset=tieset(3,itie)
      tolloc=tietol(1,itie)
!     
!     determining the element set
!     
      iset=0
      if(elset.eq.'                     ') then
        write(*,*) '*INFO reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '      no element set given'
        call inputinfo(inpc,ipoinpc,iline,
     &       "*CYCLIC SYMMETRY MODEL%")
      else
        call cident81(set,elset,nset,id)
        if(id.gt.0) then
          if(elset.eq.set(id)) then
            iset=id
          endif
        endif
        if(iset.eq.0) then
          write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
          write(*,*) '       element set does not'
          write(*,*) '       exist; '
          call inputerror(inpc,ipoinpc,iline,
     &         "*CYCLIC SYMMETRY MODEL%",ier)
          return
        endif
      endif
      cs(13,mcs)=iset+0.5d0
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      if((istat.lt.0).or.(key.eq.1)) then
        write(*,*)'*ERROR reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '      definition of the cyclic'
        write(*,*) '      symmetry model is not complete'
        write(*,*)
        ier=1
        return
      endif
!     
      ntrans=ntrans+1
      if(ntrans.gt.ntrans_) then
        write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '       increase ntrans_'
        write(*,*)
        ier=1
        return
      endif
!     
      do i=1,6
        read(textpart(i)(1:20),'(f20.0)',iostat=istat) csab(i)
        trab(i,ntrans)=csab(i)
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*CYCLIC SYMMETRY MODEL%",ier)
          return
        endif
      enddo
!     
!     cyclic coordinate system
!     
      csab(7)=-1.d0
!     
!     marker for cyclic symmetry axis
!     
      trab(7,ntrans)=2
!     
!     check whether depset and indepset exist
!     determine the kind of set (nodal or facial)
!     
      ipos=index(depset,' ')
      depkind='S'
      depset(ipos:ipos)=depkind
      call cident81(set,depset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(depset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
        depkind='T'
        depset(ipos:ipos)=depkind
        call cident81(set,depset,nset,id)
        i=nset+1
        if(id.gt.0) then
          if(depset.eq.set(id)) then
            i=id
          endif
        endif
        if(i.gt.nset) then
          write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
          write(*,*) '       surface ',depset(1:ipos-1)
          write(*,*) '       has not yet been defined.' 
          write(*,*)
          ier=1
          return
        endif
      endif
      jdep=i
!     
      ipos=index(indepset,' ')
      indepkind='S'
      indepset(ipos:ipos)=indepkind
      call cident81(set,indepset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(indepset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
        indepkind='T'
        indepset(ipos:ipos)=indepkind
        call cident81(set,indepset,nset,id)
        i=nset+1
        if(id.gt.0) then
          if(indepset.eq.set(id)) then
            i=id
          endif
        endif
        if(i.gt.nset) then
          write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
          write(*,*) '       surface ',indepset(1:ipos-1)
          write(*,*) '       has not yet been defined.' 
          write(*,*)
          ier=1
          return
        endif
      endif
      jindep=i
!     
!     unit vector along the rotation axis (xn,yn,zn)
!     
      xn=csab(4)-csab(1)
      yn=csab(5)-csab(2)
      zn=csab(6)-csab(3)
      dd=dsqrt(xn*xn+yn*yn+zn*zn)
      xn=xn/dd
      yn=yn/dd
      zn=zn/dd
!     
!     defining the indepset as a 2-D data field (axes: r=radial
!     coordinate, z=axial coordinate): needed to allocate a node
!     of the depset to a node of the indepset for the cyclic
!     symmetry equations
!     
      l=0
      do j=istartset(jindep),iendset(jindep)
        if(ialset(j).gt.0) then
          if(indepkind.eq.'T') then
!     
!     facial independent surface
!     
            ifaces=ialset(j)
            nelems=int(ifaces/10)
            jfaces=ifaces - nelems*10
            indexe=ipkon(nelems)
!     
            if(lakon(nelems)(4:5).eq.'20') then
              nopes=8
            elseif(lakon(nelems)(4:4).eq.'2') then
              nopes=9
            elseif(lakon(nelems)(4:4).eq.'8') then
              nopes=4
            elseif(lakon(nelems)(4:5).eq.'10') then
              nopes=6
            elseif(lakon(nelems)(4:4).eq.'4') then
              nopes=3
            endif
!     
            if(lakon(nelems)(4:4).eq.'6') then
              if(jfaces.le.2) then
                nopes=3
              else
                nopes=4
              endif
            endif
            if(lakon(nelems)(4:5).eq.'15') then
              if(jfaces.le.2) then
                nopes=6
              else
                nopes=8
              endif
            endif   
          else
            nopes=1
          endif
!     
          do m=1,nopes
            if(indepkind.eq.'T') then
              if((lakon(nelems)(4:4).eq.'2').or.
     &             (lakon(nelems)(4:4).eq.'8')) then
                node=kon(indexe+ifaceq(m,jfaces))
              elseif((lakon(nelems)(4:4).eq.'4').or.
     &               (lakon(nelems)(4:5).eq.'10')) then
                node=kon(indexe+ifacet(m,jfaces))
              elseif(lakon(nelems)(4:4).eq.'6') then
                node=kon(indexe+ifacew1(m,jfaces))
              elseif(lakon(nelems)(4:5).eq.'15') then
                node=kon(indexe+ifacew2(m,jfaces))
              endif
              call nident(ics,node,l,id)
              exist=.FALSE.
              if(id.gt.0) then
                if(ics(id).eq.node) then
                  exist=.TRUE.
                endif
              endif
              if(exist) cycle
              l=l+1
              if(lprev+l.gt.ncs_) then
                write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
                write(*,*) '       increase ncs_'
                write(*,*)
                ier=1
                return
              endif
              do k=l,id+2,-1
                ics(k)=ics(k-1)
                zcs(k)=zcs(k-1)
                rcs(k)=rcs(k-1)
              enddo
!     
              xap=co(1,node)-csab(1)
              yap=co(2,node)-csab(2)
              zap=co(3,node)-csab(3)
!     
              ics(id+1)=node
              zcs(id+1)=xap*xn+yap*yn+zap*zn
              rcs(id+1)=dsqrt((xap-zcs(id+1)*xn)**2+
     &             (yap-zcs(id+1)*yn)**2+
     &             (zap-zcs(id+1)*zn)**2)
            else
              l=l+1
              if(lprev+l.gt.ncs_) then
                write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
                write(*,*) '       increase ncs_'
                write(*,*)
                ier=1
                return
              endif
              node =ialset(j)
!     
              xap=co(1,node)-csab(1)
              yap=co(2,node)-csab(2)
              zap=co(3,node)-csab(3)
!     
              ics(l)=node
              zcs(l)=xap*xn+yap*yn+zap*zn
              rcs(l)=dsqrt((xap-zcs(l)*xn)**2+
     &             (yap-zcs(l)*yn)**2+
     &             (zap-zcs(l)*zn)**2)
            endif
          enddo
        else
          k=ialset(j-2)
          do
            k=k-ialset(j)
            if(k.ge.ialset(j-1)) exit
            l=l+1
            if(l.gt.ncs_) then
              write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL:'
              write(*,*) '       increase ncs_'
              write(*,*)
              ier=1
              return
            endif
            node=k
!     
            xap=co(1,node)-csab(1)
            yap=co(2,node)-csab(2)
            zap=co(3,node)-csab(3)
!     
            ics(l)=node
            zcs(l)=xap*xn+yap*yn+zap*zn
            rcs(l)=dsqrt((xap-zcs(l)*xn)**2+
     &           (yap-zcs(l)*yn)**2+
     &           (zap-zcs(l)*zn)**2)
          enddo
        endif
      enddo
!     
      ncsnodes=l
!     
!     initialization of near2d
!     
      do i=1,ncsnodes
        nr(i)=i
        nz(i)=i
        rcs0(i)=rcs(i)
        zcs0(i)=zcs(i)
      enddo
      kflag=2
      call dsort(rcs,nr,ncsnodes,kflag)
      call dsort(zcs,nz,ncsnodes,kflag)
!     
!     check whether a tolerance was defined. If not, a tolerance
!     is calculated as 1.e-10 times the mean of the distance of every
!     independent node to its nearest neighbour
!     
      if(tolloc.lt.0.d0) then
        nneigh=2
        dist=0.d0
        do i=1,ncsnodes
          nodei=ics(i)
!     
          xap=co(1,nodei)-csab(1)
          yap=co(2,nodei)-csab(2)
          zap=co(3,nodei)-csab(3)
!     
          zp=xap*xn+yap*yn+zap*zn
          rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!     
          call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rp,zp,ncsnodes,noden,
     &         nneigh)
!     
          dist=dist+dsqrt((co(1,nodei)-co(1,noden(2)))**2+
     &         (co(2,nodei)-co(2,noden(2)))**2+
     &         (co(3,nodei)-co(3,noden(2)))**2)
        enddo
        tolloc=1.d-10*dist/ncsnodes
        write(*,*) '*INFO reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '      no tolerance was defined'
        write(*,*) '      in the *TIE option; a tolerance of ',
     &       tolloc
        write(*,*) '      will be used'
        write(*,*)
      endif
!     
!     calculating the angle between dependent and independent
!     side and check for nodes on the axis
!     
!     this angle may be different from 2*pi/xsectors: in that way
!     the user can simulate fractional nodal diameters
!     
!     (x2,y2,z2): unit vector on the dependent side and orthogonal
!     to the rotation axis
!     (x3,y3,z3): unit vector on the independent side and orthogonal
!     to the rotation axis
!     (x1,y1,z1)=(x2,y2,z2)x(x3,y3,z3)
!     points in the same direction of xn if the independent
!     side is on the clockwise side of the dependent side if
!     looking in the direction of xn
!     
      calcangle=.false.
      nodesonaxis=.false.
      phi=0.d0
      axdistmin=1.d30
!     
      nneigh=1
      loop1: do i=istartset(jdep),iendset(jdep)
      if(ialset(i).gt.0) then
!     
!     check whether dependent side is node based or
!     face based
!     
        if(depkind.eq.'T') then
          ifaces=ialset(i)
          nelems=int(ifaces/10)
          jfaces=ifaces - nelems*10
          indexe=ipkon(nelems)
!     
          if(lakon(nelems)(4:5).eq.'20') then
            nopes=8
          elseif(lakon(nelems)(4:4).eq.'2') then
            nopes=9
          elseif(lakon(nelems)(4:4).eq.'8') then
            nopes=4
          elseif(lakon(nelems)(4:5).eq.'10') then
            nopes=6
          elseif(lakon(nelems)(4:4).eq.'4') then
            nopes=3
          endif
!     
          if(lakon(nelems)(4:4).eq.'6') then
            if(jfaces.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
          if(lakon(nelems)(4:5).eq.'15') then
            if(jfaces.le.2) then
              nopes=6
            else
              nopes=8
            endif
          endif 
        else
          nopes=1
        endif
!     
        do m=1,nopes
          if(depkind.eq.'T') then
            if((lakon(nelems)(4:4).eq.'2').or.
     &           (lakon(nelems)(4:4).eq.'8')) then
              noded=kon(indexe+ifaceq(m,jfaces))
            elseif((lakon(nelems)(4:4).eq.'4').or.
     &             (lakon(nelems)(4:5).eq.'10')) then
              noded=kon(indexe+ifacet(m,jfaces))
            elseif(lakon(nelems)(4:4).eq.'6') then
              noded=kon(indexe+ifacew1(m,jfaces))
            elseif(lakon(nelems)(4:5).eq.'15') then
              noded=kon(indexe+ifacew2(m,jfaces))
            endif
          else
            if(i.gt.istartset(jdep)) then
              if(ialset(i).eq.ialset(i-1)) cycle loop1
            endif
            noded=ialset(i)
          endif
!     
          xap=co(1,noded)-csab(1)
          yap=co(2,noded)-csab(2)
          zap=co(3,noded)-csab(3)
!     
          zpd=xap*xn+yap*yn+zap*zn
          rpd=dsqrt((xap-zpd*xn)**2+(yap-zpd*yn)**2+
     &         (zap-zpd*zn)**2)
!     
          if((.not.calcangle).and.(rpd.gt.1.d-10)) then
            x2=(xap-zpd*xn)/rpd
            y2=(yap-zpd*yn)/rpd
            z2=(zap-zpd*zn)/rpd
          endif
!     
          call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rpd,zpd,ncsnodes,
     &         noden,nneigh)
          node=noden(1)
!     
          nodei=ics(node)
          if(nodei.lt.0) cycle
          if(nodei.eq.noded) then
            ics(node)=-nodei
            nodesonaxis=.true.
            cycle
          endif
!     
          xap=co(1,nodei)-csab(1)
          yap=co(2,nodei)-csab(2)
          zap=co(3,nodei)-csab(3)
!     
          zp=xap*xn+yap*yn+zap*zn
          rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
!         store the minimum axial distance between dependent and
!         independent nodes
!
          if(dabs(zp-zpd).lt.axdistmin) then
            axdistmin=dabs(zp-zpd)
            nodeaxd=noded
            nodeaxi=nodei
          endif
!     
!     in order for the angle to be correct the axial position
!     of the dependent and independent node must be the same
!     (important for non-coincident meshes)
!     
          if((.not.calcangle).and.(rp.gt.1.d-10).and.
     &         (dabs(zp-zpd).lt.1.d-10)) then
            x3=(xap-zp*xn)/rp
            y3=(yap-zp*yn)/rp
            z3=(zap-zp*zn)/rp
!     
            x1=y2*z3-y3*z2
            y1=x3*z2-x2*z3
            z1=x2*y3-x3*y2
!     
            phi=(x1*xn+y1*yn+z1*zn)/dabs(x1*xn+y1*yn+z1*zn)*
     &           dacos(x2*x3+y2*y3+z2*z3)
            if(check) then
              calculated_angle=dacos(x2*x3+y2*y3+z2*z3)
              user_angle=6.28318531d0/cs(1,mcs)
              if(dabs(calculated_angle-user_angle)/
     &             calculated_angle.gt.0.01d0) then
                write(*,*) 
     &               '*ERROR reading *CYCLIC SYMMETRY MODEL'
                write(*,*) '       number of segments does not'
                write(*,*) '       agree with the geometry'
                write(*,*) '       angle based on N:',
     &               user_angle*57.29577951d0
                write(*,*)'       angle based on the geometry:',
     &               calculated_angle*57.29577951d0
                write(*,*)
                ier=1
                return
              endif
            else
              write(*,*) '*INFO in cyclicsymmetrymodels: angle'
              write(*,*)'      check is deactivated by the user;'
              write(*,*) '      the real geometry is used for'
              write(*,*) '      the calculation of the segment'
              write(*,*) '      angle'
              write(*,*)
            endif
            calcangle=.true.
          endif
        enddo
!     
      else
        k=ialset(i-2)
        do
          k=k-ialset(i)
          if(k.ge.ialset(i-1)) exit
          noded=k
!     
          xap=co(1,noded)-csab(1)
          yap=co(2,noded)-csab(2)
          zap=co(3,noded)-csab(3)
!     
          zpd=xap*xn+yap*yn+zap*zn
          rpd=dsqrt((xap-zpd*xn)**2+(yap-zpd*yn)**2+
     &         (zap-zpd*zn)**2)
!     
          if((.not.calcangle).and.(rpd.gt.1.d-10)) then
            x2=(xap-zpd*xn)/rpd
            y2=(yap-zpd*yn)/rpd
            z2=(zap-zpd*zn)/rpd
          endif
!     
          call near2d(rcs0,zcs0,rcs,zcs,nr,nz,rpd,zpd,ncsnodes,
     &         noden,nneigh)
          node=noden(1)
!     
          nodei=ics(node)
          if(nodei.lt.0) cycle
          if(nodei.eq.noded) then
            ics(node)=-nodei
            nodesonaxis=.true.
            cycle
          endif
!     
          xap=co(1,nodei)-csab(1)
          yap=co(2,nodei)-csab(2)
          zap=co(3,nodei)-csab(3)
!     
          zp=xap*xn+yap*yn+zap*zn
          rp=dsqrt((xap-zp*xn)**2+(yap-zp*yn)**2+(zap-zp*zn)**2)
!
!         store the minimum axial distance between dependent and
!         independent nodes
!
          if(dabs(zp-zpd).lt.axdistmin) then
            axdistmin=dabs(zp-zpd)
            nodeaxd=noded
            nodeaxi=nodei
          endif
!     
!     in order for the angle to be correct the axial position
!     of the dependent and independent node must be the same
!     (important for non-coincident meshes)
!     
          if((.not.calcangle).and.(rp.gt.1.d-10).and.
     &         (dabs(zp-zpd).lt.1.d-10)) then
            x3=(xap-zp*xn)/rp
            y3=(yap-zp*yn)/rp
            z3=(zap-zp*zn)/rp
!     
            x1=y2*z3-y3*z2
            y1=x3*z2-x2*z3
            z1=x2*y3-x3*y2
!     
            phi=(x1*xn+y1*yn+z1*zn)/dabs(x1*xn+y1*yn+z1*zn)*
     &           dacos(x2*x3+y2*y3+z2*z3)
            if(check) then
              calculated_angle=dacos(x2*x3+y2*y3+z2*z3)
              user_angle=6.28318531d0/cs(1,mcs)
              if(dabs(calculated_angle-user_angle)
     &             /calculated_angle.gt.0.01d0) then
                write(*,*) 
     &               '*ERROR reading *CYCLIC SYMMETRY MODEL'
                write(*,*) '       number of segments does not'
                write(*,*) '       agree with the geometry'
                write(*,*) '       angle based on N:',
     &               user_angle*57.29577951d0
                write(*,*) '       angle based on the geometry:'
     &               ,calculated_angle*57.29577951d0
                write(*,*)
                ier=1
                return
              endif
            endif
            calcangle=.true.
          endif
!     
        enddo
      endif
!     
      enddo loop1
!     
      if(phi.eq.0.d0) then
        write(*,*) '*ERROR reading *CYCLIC SYMMETRY MODEL'
        write(*,*) '       sector angle cannot be determined:'
        write(*,*) '       there exists no dependent node'
        write(*,*) '       with the same axial position as'
        write(*,*) '       an independent node (a tolerance'
        write(*,*) '       of 1.e-10 length units is applied).'
        write(*,*) '       The smallest axial distance of ',
     &       axdistmin
        write(*,*) '       exists between dependent node ',nodeaxd
        write(*,*) '       and independent node ',nodeaxi
        write(*,*)
        ier=1
        return
      endif
!     
!     allocating a node of the depset to each node of the indepset 
!     
      triangulation=.false.
!     
!     opening a file to store the nodes which are not connected
!
      do ilen=1,132
        if(ichar(jobnamec(1)(ilen:ilen)).eq.0) exit
      enddo
      ilen=ilen-1
      fn=jobnamec(1)(1:ilen)//'_WarnNodeMissCyclicSymmetry.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeCyclicSymmetry'
      icount=0
!     
!     generating the thermal MPC's; the generated MPC's are for nodal
!     diameter 0. BETTER: based on ithermal(2), cf. gen3dfrom2d.f
!     
!     about next info write statement:
!     in tempload cyclic symmetry is enforced for field t1, but not
!     for field t0. This may lead to stresses if t1 is not cyclic
!     symmetric. If there is a *initial conditions,type=temperature
!     card in the input deck but no *temperature card t1 is copied from
!     t0 before ensuring the cyclic symmetry for t1. So also in this
!     case a non-cyclic symmetric field t0 can lead to stresses.
!     
      if(ithermal(1).eq.1) then
        write(*,*) '*INFO reading *CYCLIC SYMMETRY MODEL'
        write(*,*) '      cyclic symmetry equations are generated'
        write(*,*) '      for the temperature; if the initial'
        write(*,*) '      temperatures are not cyclic symmetric'
        write(*,*) '      and/or the applied temperature is not'
        write(*,*) '      cyclic symmetric this may lead to'
        write(*,*) '      additional stresses'
        write(*,*)
      endif
!     
      loop2: do i=istartset(jdep),iendset(jdep)
      if(ialset(i).gt.0) then
!     
!     check whether dependent side is node based or
!     face based
!     
        if(depkind.eq.'T') then
          ifaces=ialset(i)
          nelems=int(ifaces/10)
          jfaces=ifaces - nelems*10
          indexe=ipkon(nelems)
!     
          if(lakon(nelems)(4:5).eq.'20') then
            nopes=8
          elseif(lakon(nelems)(4:4).eq.'2') then
            nopes=9
          elseif(lakon(nelems)(4:4).eq.'8') then
            nopes=4
          elseif(lakon(nelems)(4:5).eq.'10') then
            nopes=6
          elseif(lakon(nelems)(4:4).eq.'4') then
            nopes=3
          endif
!     
          if(lakon(nelems)(4:4).eq.'6') then
            if(jfaces.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
          if(lakon(nelems)(4:5).eq.'15') then
            if(jfaces.le.2) then
              nopes=6
            else
              nopes=8
            endif
          endif 
        else
          nopes=1
        endif
!     
        do m=1,nopes
          if(depkind.eq.'T') then
            if((lakon(nelems)(4:4).eq.'2').or.
     &           (lakon(nelems)(4:4).eq.'8')) then
              noded=kon(indexe+ifaceq(m,jfaces))
            elseif((lakon(nelems)(4:4).eq.'4').or.
     &             (lakon(nelems)(4:5).eq.'10')) then
              noded=kon(indexe+ifacet(m,jfaces))
            elseif(lakon(nelems)(4:4).eq.'6') then
              noded=kon(indexe+ifacew1(m,jfaces))
            elseif(lakon(nelems)(4:5).eq.'15') then
              noded=kon(indexe+ifacew2(m,jfaces))
            endif
          else
            if(i.gt.istartset(jdep)) then
              if(ialset(i).eq.ialset(i-1)) cycle loop2
            endif
            noded=ialset(i)
          endif
!     
!     check whether cyclic MPC's have already been
!     generated (e.g. for nodes belonging to several
!     faces for face based dependent surfaces)
!     
          idof=8*(noded-1)+1
          call nident(ikmpc,idof,nmpc,id)
          if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
              if(labmpc(ilmpc(id))(1:6).eq.'CYCLIC') cycle
            endif
          endif
!     
          call generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &         coefmpc,nmpc,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,
     &         nr,nz,rcs0,zcs0,labmpc,
     &         mcs,triangulation,csab,xn,yn,zn,phi,noded,
     &         ncsnodes,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,
     &         nzcg,jcs,lcs,kontri,straight,ne,ipkon,kon,lakon,
     &         ifacetet,inodface,vold,nef,mi,
     &         indepset,ithermal,icount)
        enddo
!     
      else
        k=ialset(i-2)
        do
          k=k-ialset(i)
          if(k.ge.ialset(i-1)) exit
          noded=k
!     
!     check whether cyclic MPC's have already been
!     generated (e.g. for nodes belonging to several
!     faces for face based dependent surfaces)
!     
          idof=8*(noded-1)+1
          call nident(ikmpc,idof,nmpc,id)
          if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
              if(labmpc(ilmpc(id))(1:6).eq.'CYCLIC') cycle
            endif
          endif
!     
          call generatecycmpcs(tolloc,co,nk,ipompc,nodempc,
     &         coefmpc,nmpc,ikmpc,ilmpc,mpcfree,rcs,zcs,ics,
     &         nr,nz,rcs0,zcs0,labmpc,
     &         mcs,triangulation,csab,xn,yn,zn,phi,noded,
     &         ncsnodes,rcscg,rcs0cg,zcscg,zcs0cg,nrcg,
     &         nzcg,jcs,lcs,kontri,straight,ne,ipkon,kon,lakon,
     &         ifacetet,inodface,vold,nef,mi,
     &         indepset,ithermal,icount)
        enddo
      endif
!     
      enddo loop2
!
      if(icount.gt.0) then
        write(*,*) '*WARNING reading *CYCLIC SYMMETRY MODEL:'
        write(*,*) '        for at least one dependent'
        write(*,*) '        node in a cyclic symmetry definition no '
        write(*,*) '        independent counterpart was found.'
        write(*,*) '        Failed nodes are stored in file '
        write(*,*) '        ',fn(1:ilen+31)
        write(*,*) '        This file can be loaded into'
        write(*,*) '        an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+31),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!     
!     sorting ics
!     ics contains the master (independent) nodes
!     
      kflag=1
      call isortii(ics,nr,ncsnodes,kflag)
      cs(4,mcs)=ncsnodes+0.5d0
      lprev=lprev+ncsnodes
!     
!     check orientation of (xn,yn,zn) (important for copying of base
!     sector in arpackcs)
!     
      if(phi.lt.0.d0) then
        csab(4)=2.d0*csab(1)-csab(4)
        csab(5)=2.d0*csab(2)-csab(5)
        csab(6)=2.d0*csab(3)-csab(6)
      endif
!     
      do i=1,7
        cs(5+i,mcs)=csab(i)
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

