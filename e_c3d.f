!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine e_c3d(co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &     ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &     alzero,ielmat,ielorien,norien,orab,ntmat_,
     &     t0,t1,ithermal,vold,iperturb,nelemload,
     &     sideload,xload,nload,idist,sti,stx,iexpl,plicon,
     &     nplicon,plkcon,nplkcon,xstiff,npmat_,dtime,
     &     matname,mi,ncmat_,mass,stiffness,buckling,rhsi,intscheme,
     &     ttime,time,istep,iinc,coriolis,xloadold,reltime,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,veold,springarea,
     &     nstate_,xstateini,xstate,ne0,ipkon,thicke,
     &     integerglob,doubleglob,tieset,istartset,iendset,ialset,ntie,
     &     nasym,pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,
     &     kscale,smscalel,mscalmethod,set,nset,islavelinv,autloc,
     &     irowtloc,jqtloc,mortartrafoflag)
!     
!     computation of the element matrix and rhs for the element with
!     the topology in konl
!     
!     ff: rhs without temperature and eigenstress contribution
!     
!     nmethod=0: check for positive Jacobian
!     nmethod=1: stiffness matrix + right hand side
!     nmethod=2: stiffness matrix + mass matrix
!     nmethod=3: static stiffness + buckling stiffness
!     nmethod=4: stiffness matrix + mass matrix
!     
!     intscheme=0: use the integration scheme corresponding to the element     
!                  type, i.e. reduced integration for reduced integration    
!                  elements etc.
!     intscheme=1: use for C3D8R and C3D20R elements the integration point
!                  schemes for C3D8 and C3D20, respectively, for all other    
!                  elements the the integration scheme corresponding to the
!                  element type     
!     
      implicit none
!     
      integer mass,stiffness,buckling,rhsi,coriolis
!     
      character*1 entity
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*),amat
      character*81 tieset(3,*),set(*)
!     
      integer konl(26),ifaceq(8,6),nelemload(2,*),nbody,nelem,
     &     mi(*),jfaces,igauss,mortar,kon(*),ielprop(*),null,
     &     mattyp,ithermal(*),iperturb(*),nload,idist,i,j,k,l,i1,i2,j1,
     &     nmethod,k1,l1,ii,jj,ii1,jj1,id,ipointer,ig,m1,m2,m3,m4,kk,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),six,
     &     ielorien(mi(3),*),ilayer,nlayer,ki,kl,ipkon(*),indexe,
     &     ntmat_,nope,nopes,norien,ihyper,iexpl,kode,imat,mint2d,
     &     mint3d,ifacet(6,4),nopev,iorien,istiff,ncmat_,iface,
     &     ifacew(8,5),intscheme,n,ipointeri,ipointerj,istep,iinc,
     &     layer,kspt,jltyp,iflag,iperm(60),m,ipompc(*),nodempc(3,*),
     &     nmpc,ikmpc(*),ilmpc(*),iscale,nstate_,ne0,iselect(6),kscale,
     &     istartset(*),iendset(*),ialset(*),ntie,integerglob(*),nasym,
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,nopered,
     &     mscalmethod,nset,islavelinv(*),jqtloc(*),irowtloc(*),
     &     node1,node2,irowtloc1(16),jqtloc1(9),j2,mortartrafoflag
!     
      real*8 co(3,*),xl(3,26),shp(4,26),xs2(3,7),veold(0:mi(2),*),
     &     s(60,60),w(3,3),p1(3),p2(3),bodyf(3),bodyfx(3),ff(60),
     &     bf(3),q(3),shpj(4,26),elcon(0:ncmat_,ntmat_,*),t(3),
     &     rhcon(0:1,ntmat_,*),xkl(3,3),eknlsign,reltime,prop(*),
     &     alcon(0:6,ntmat_,*),alzero(*),orab(7,*),t0(*),t1(*),
     &     anisox(3,3,3,3),voldl(0:mi(2),26),vo(3,3),xloadold(2,*),
     &     xl2(3,9),xsj2(3),shp2(7,9),vold(0:mi(2),*),xload(2,*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     v(3,3,3,3),springarea(2,*),thickness,tlayer(4),dlayer(4),
     &     om,omx,e,un,al,um,xi,et,ze,tt,const,xsj,xsjj,sm(60,60),
     &     sti(6,mi(1),*),stx(6,mi(1),*),s11,s22,s33,s12,s13,s23,s11b,
     &     s22b,s33b,s12b,s13b,s23b,t0l,t1l,coefmpc(*),xlayer(mi(3),4),
     &     senergy,senergyb,rho,elas(21),summass,summ,thicke(mi(3),*),
     &     sume,factorm,factore,alp,elconloc(ncmat_),eth(6),
     &     weight,coords(3),dmass,xl1(3,9),term,clearini(3,9,*),
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),plconloc(802),dtime,ttime,time,tvar(2),
     &     sax(60,60),ffax(60),gs(8,4),a,stress(6),stre(3,3),
     &     pslavsurf(3,*),pmastsurf(6,*),xmass,xsjmass,shpmass(4,26),
     &     shpjmass(4,26),smscalel,smfactor,shptil(4,26),autloc(*),
     &     shptil2(7,9),autloc1(16),doubleglob(*)
!     
      include "gauss.f"
!     
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/),(/8,5/))
      iflag=3
      null=0
      iperm=(/13,14,-15,16,17,-18,19,20,-21,22,23,-24,
     &     1,2,-3,4,5,-6,7,8,-9,10,11,-12,
     &     37,38,-39,40,41,-42,43,44,-45,46,47,-48,
     &     25,26,-27,28,29,-30,31,32,-33,34,35,-36,
     &     49,50,-51,52,53,-54,55,56,-57,58,59,-60/)
!     
      tvar(1)=time
      tvar(2)=ttime+time
!     
      summass=0.d0
!     
      indexe=ipkon(nelem)
c     Bernhardi start
      if(lakonl(1:5).eq.'C3D8I') then
        nope=11
        nopev=8
        nopes=4
      elseif(lakonl(4:5).eq.'20') then
c     Bernhardi end
        nope=20
        nopev=8
        nopes=8
      elseif(lakonl(4:4).eq.'8') then
        nope=8
        nopev=8
        nopes=4
      elseif(lakonl(4:5).eq.'10') then
        nope=10
        nopev=4
        nopes=6
      elseif(lakonl(4:4).eq.'4') then
        nope=4
        nopev=4
        nopes=3
      elseif(lakonl(4:5).eq.'15') then
        nope=15
        nopev=6
      elseif(lakonl(4:4).eq.'6') then
        nope=6
        nopev=6
      elseif(lakonl(1:2).eq.'ES') then
        if(lakonl(7:7).eq.'C') then
          if(mortar.eq.0) then
            nope=ichar(lakonl(8:8))-47
            konl(nope+1)=kon(indexe+nope+1)
          elseif(mortar.eq.1) then
            nope=kon(indexe)
          endif
        else
          nope=ichar(lakonl(8:8))-47
        endif
      elseif(lakonl(1:4).eq.'MASS') then
        nope=1
      endif
!     
!     material and orientation
!     
      if(lakonl(7:8).ne.'LC') then
!     
        imat=ielmat(1,nelem)
        amat=matname(imat)
        if(norien.gt.0) then
          iorien=max(0,ielorien(1,nelem))
        else
          iorien=0
        endif
!     
        if(nelcon(1,imat).lt.0) then
          ihyper=1
        else
          ihyper=0
        endif
      elseif(lakonl(4:5).eq.'20') then
!     
!     composite materials
!     
!     determining the number of layers
!     
        nlayer=0
        do k=1,mi(3)
          if(ielmat(k,nelem).ne.0) then
            nlayer=nlayer+1
          endif
        enddo
        mint2d=4
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
        iflag=1
        do kk=1,mint2d
          xi=gauss3d2(1,kk)
          et=gauss3d2(2,kk)
          call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
          tlayer(kk)=0.d0
          do i=1,nlayer
            thickness=0.d0
            do j=1,nopes
              thickness=thickness+thicke(i,indexe+j)*shp2(4,j)
            enddo
            tlayer(kk)=tlayer(kk)+thickness
            xlayer(i,kk)=thickness
          enddo
        enddo
        iflag=3
!     
        ilayer=0
        do i=1,4
          dlayer(i)=0.d0
        enddo
!     
      elseif(lakonl(4:5).eq.'15') then
!     
!     composite materials
!     
!     determining the number of layers
!     
        nlayer=0
        do k=1,mi(3)
          if(ielmat(k,nelem).ne.0) then
            nlayer=nlayer+1
          endif
        enddo
        mint2d=3
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
        iflag=1
        do kk=1,mint2d
          xi=gauss3d10(1,kk)
          et=gauss3d10(2,kk)
          call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
          tlayer(kk)=0.d0
          do i=1,nlayer
            thickness=0.d0
            do j=1,6
              thickness=thickness+thicke(i,indexe+j)*shp2(4,j)
            enddo
            tlayer(kk)=tlayer(kk)+thickness
            xlayer(i,kk)=thickness
          enddo
        enddo
        iflag=3
!     
        ilayer=0
        do i=1,3
          dlayer(i)=0.d0
        enddo
!     
      endif
!     
        if(lakonl(4:5).eq.'8R') then
          mint2d=1
          if(intscheme.eq.0) then
            mint3d=1
          else
            mint3d=8
          endif
        elseif(lakonl(4:7).eq.'20RB') then
          if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint2d=4
            mint3d=50
          else
            mint2d=4
            call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &           null,xi,et,ze,weight)
          endif
        elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
          if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &         (lakonl(7:7).eq.'E')) then
            mint2d=2
            mint3d=4
          else
            mint2d=4
            if(lakonl(7:8).eq.'LC') then
              mint3d=8*nlayer
            else
              if((intscheme.eq.0).or.(lakonl(4:8).ne.'20R  ')) then
                mint3d=8
              else
                mint3d=27
              endif
            endif
          endif
        elseif(lakonl(4:4).eq.'2') then
          mint2d=9
          mint3d=27
        elseif(lakonl(4:5).eq.'10') then
          mint2d=3
          mint3d=4
        elseif(lakonl(4:4).eq.'4') then
          mint2d=1
          mint3d=1
        elseif(lakonl(4:5).eq.'15') then
          if(lakonl(7:8).eq.'LC') then
            mint3d=6*nlayer
          else
            mint3d=9
          endif
        elseif(lakonl(4:4).eq.'6') then
          mint3d=2
        else
          mint3d=0
        endif
!     
!     computation of the coordinates of the local nodes
!     
      do i=1,nope
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!     
!     initialisation for distributed forces
!     
      if(rhsi.eq.1) then
        if(idist.ne.0) then
          do i=1,3*nope
            ff(i)=0.d0
          enddo
        endif
      endif
!     
!     displacements for 2nd order static and modal theory
!     
      if(((iperturb(1).eq.1).or.(iperturb(2).eq.1)).and.
     &     (stiffness.eq.1).and.(buckling.eq.0)) then
        do i1=1,nope
          do i2=1,3
            voldl(i2,i1)=vold(i2,konl(i1))
          enddo
        enddo
      endif
!     
!     initialisation of sm
!     
      if((mass.eq.1).or.(buckling.eq.1).or.(coriolis.eq.1)) then
        do i=1,3*nope
          do j=1,3*nope
            sm(i,j)=0.d0
          enddo
        enddo
      endif
!     
!     initialisation of s
!     
      do i=1,3*nope
        do j=1,3*nope
          s(i,j)=0.d0
        enddo
      enddo
!     
!     calculating the stiffness matrix for the (contact) spring elements
!     and the mass matrix for the mass elements
!     
      if(mint3d.eq.0) then
!     
!     for a
!     first step in a displacement-driven geometrically nonlinear
!     calculation or a nonlinear calculation with linear strains
!     the next block has to be evaluated
!     
        if(iperturb(2).eq.0) then
          do i1=1,nope
            do i2=1,3
              voldl(i2,i1)=vold(i2,konl(i1))
            enddo
          enddo
        endif
!     
        kode=nelcon(1,imat)
        if(lakonl(1:2).eq.'ES') then
          if(lakonl(7:7).ne.'C') then
            t0l=0.d0
            t1l=0.d0
            if(ithermal(1).eq.1) then
              t0l=(t0(konl(1))+t0(konl(2)))/2.d0
              t1l=(t1(konl(1))+t1(konl(2)))/2.d0
            elseif(ithermal(1).ge.2) then
              t0l=(t0(konl(1))+t0(konl(2)))/2.d0
              t1l=(vold(0,konl(1))+vold(0,konl(2)))/2.d0
            endif
          endif
!     
          if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'1').or.
     &         (lakonl(7:7).eq.'2').or.(mortar.eq.0)) then
            call springstiff_n2f(xl,elas,konl,voldl,s,imat,elcon,
     &           nelcon,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &           plicon,nplicon,npmat_,iperturb,
     &           springarea(1,konl(nope+1)),nmethod,mi,ne0,nstate_,
     &           xstateini,xstate,reltime,nasym,ielorien,orab,norien,
     &           nelem)
          elseif(mortar.eq.1) then
            jfaces=kon(indexe+nope+2)
            igauss=kon(indexe+nope+1)
c     write(*,*) 'e_c3d ',nelem
            call springstiff_f2f(xl,elas,voldl,s,imat,elcon,nelcon,
     &           ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,plicon,
     &           nplicon,npmat_,iperturb,springarea(1,igauss),
     &           nmethod,mi,ne0,nstate_,xstateini,xstate,reltime,
     &           nasym,jfaces,igauss,pslavsurf,
     &           pmastsurf,clearini,kscale)
          endif
        elseif(lakonl(1:4).eq.'MASS') then
!     
!     mass matrix contribution
!     
          if(mass.eq.1) then
            do i1=1,3
              sm(i1,i1)=rhcon(1,1,imat)
            enddo
          endif
!     
!     contribution to the rhs
!     
          if(rhsi.eq.1) then
            if(nbody.ne.0) then
!     
              xmass=rhcon(1,1,imat)
!     
              do i1=1,3
                ff(i1)=bodyfx(i1)*xmass
              enddo
!     
              if(omx.gt.0.d0) then
!     
!     omega**2 * mass
!     
                om=omx*xmass
                do i1=1,3
!     
!     computation of the global coordinates of the
!     mass node
!     
                  if((iperturb(1).ne.1).and.(iperturb(2).ne.1)) 
     &                 then
                    q(i1)=xl(i1,1)
                  else
                    q(i1)=xl(i1,1)+voldl(i1,1)
                  endif
!     
                  q(i1)=q(i1)-p1(i1)
                enddo
                const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!     
!     inclusion of the centrifugal force into the body force
!     
                do i1=1,3
                  ff(i1)=ff(i1)+(q(i1)-const*p2(i1))*om
                enddo
              endif
            endif
          endif
!     
        endif
        return
      endif
!     
!     computation of the matrix: loop over the Gauss points
!     
      do kk=1,mint3d
          if(lakonl(4:5).eq.'8R') then
            if(intscheme.eq.0) then
              xi=gauss3d1(1,kk)
              et=gauss3d1(2,kk)
              ze=gauss3d1(3,kk)
              weight=weight3d1(kk)
            else
!
!     initial acceleration calculation for dynamic calculations:
!     full integration used
!
              xi=gauss3d2(1,kk)
              et=gauss3d2(2,kk)
              ze=gauss3d2(3,kk)
              weight=weight3d2(kk)
            endif
          elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
              xi=gauss3d13(1,kk)
              et=gauss3d13(2,kk)
              ze=gauss3d13(3,kk)
              weight=weight3d13(kk)
            else
              call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &             kk,xi,et,ze,weight)
            endif
          elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &           then
            if(lakonl(7:8).ne.'LC') then
              if((intscheme.eq.0).or.(lakonl(4:8).ne.'20R  ')) then
                xi=gauss3d2(1,kk)
                et=gauss3d2(2,kk)
                ze=gauss3d2(3,kk)
                weight=weight3d2(kk)
              else
!
!     initial acceleration calculation for dynamic calculations:
!     full integration used
!
                xi=gauss3d3(1,kk)
                et=gauss3d3(2,kk)
                ze=gauss3d3(3,kk)
                weight=weight3d3(kk)
              endif
            else
              kl=mod(kk,8)
              if(kl.eq.0) kl=8
!     
              xi=gauss3d2(1,kl)
              et=gauss3d2(2,kl)
              ze=gauss3d2(3,kl)
              weight=weight3d2(kl)
!     
              ki=mod(kk,4)
              if(ki.eq.0) ki=4
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do i=1,4
                    dlayer(i)=dlayer(i)+xlayer(ilayer-1,i)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,nelem)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,nelem))
              else
                iorien=0
              endif
!     
              if(nelcon(1,imat).lt.0) then
                ihyper=1
              else
                ihyper=0
              endif
            endif
          elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,kk)
            et=gauss3d3(2,kk)
            ze=gauss3d3(3,kk)
            weight=weight3d3(kk)
          elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,kk)
            et=gauss3d5(2,kk)
            ze=gauss3d5(3,kk)
            weight=weight3d5(kk)
          elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,kk)
            et=gauss3d4(2,kk)
            ze=gauss3d4(3,kk)
            weight=weight3d4(kk)
          elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d8(1,kk)
              et=gauss3d8(2,kk)
              ze=gauss3d8(3,kk)
              weight=weight3d8(kk)
            else
              kl=mod(kk,6)
              if(kl.eq.0) kl=6
!     
              xi=gauss3d10(1,kl)
              et=gauss3d10(2,kl)
              ze=gauss3d10(3,kl)
              weight=weight3d10(kl)
!     
              ki=mod(kk,3)
              if(ki.eq.0) ki=3
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do i=1,3
                    dlayer(i)=dlayer(i)+xlayer(ilayer-1,i)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,nelem)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,nelem))
              else
                iorien=0
              endif
!     
              if(nelcon(1,imat).lt.0) then
                ihyper=1
              else
                ihyper=0
              endif
            endif
          elseif(lakonl(4:4).eq.'6') then
            xi=gauss3d7(1,kk)
            et=gauss3d7(2,kk)
            ze=gauss3d7(3,kk)
            weight=weight3d7(kk)
          endif
!     
!     calculation of the shape functions and their derivatives
!     in the gauss point
!     
c     Bernhardi start
        if((lakonl(1:5).eq.'C3D8R').and.(intscheme.eq.0)) then
          call shape8hr(xl,xsj,shp,gs,a)
        elseif(lakonl(1:5).eq.'C3D8I') then
          call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
          call shape8humass(xi,et,ze,xl,xsjmass,shpmass,iflag)
        elseif(nope.eq.20) then
c     Bernhardi end
          if(lakonl(7:7).eq.'A') then
            call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
          elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
            call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
          else
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
          endif
        elseif(nope.eq.8) then
          call shape8h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.10) then
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.4) then
          call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(nope.eq.15) then
          call shape15w(xi,et,ze,xl,xsj,shp,iflag)
        else
          call shape6w(xi,et,ze,xl,xsj,shp,iflag)
        endif
!     
c     mortar start
        if(mortartrafoflag.eq.1) then
          do i1=1,nope
            shptil(1,i1)=shp(1,i1)
            shptil(2,i1)=shp(2,i1)
            shptil(3,i1)=shp(3,i1)
            shptil(4,i1)=shp(4,i1)
          enddo
          if(islavelinv(nelem).gt.0) then
            if((nope.eq.20).or.(nope.eq.10).or.(nope.eq.15)) then
              do i1=1,nope
                if(jqtloc(i1+1)-jqtloc(i1).gt.0) then
                  shptil(1,i1)=0.0
                  shptil(2,i1)=0.0
                  shptil(3,i1)=0.0
                  shptil(4,i1)=0.0
                endif
                do j1=jqtloc(i1),jqtloc(i1+1)-1
                  j2=irowtloc(j1)
                  shptil(1,i1)=shptil(1,i1)+autloc(j1)*shp(1,j2)
                  shptil(2,i1)=shptil(2,i1)+autloc(j1)*shp(2,j2)
                  shptil(3,i1)=shptil(3,i1)+autloc(j1)*shp(3,j2)
                  shptil(4,i1)=shptil(4,i1)+autloc(j1)*shp(4,j2)
                enddo
              enddo
            endif
          endif
        endif
c     mortar end
!     
!     check the jacobian determinant
!     
        if(xsj.lt.1.d-20) then
          write(*,*) '*ERROR in e_c3d: nonpositive jacobian'
          write(*,*) '       determinant in element',nelem
          write(*,*)
          xsj=dabs(xsj)
          nmethod=0
        endif
!     
        if(((iperturb(1).eq.1).or.(iperturb(2).eq.1)).and.
     &       (stiffness.eq.1).and.(buckling.eq.0)) then
!     
!     stresses for 2nd order static and modal theory
!     
          s11=sti(1,kk,nelem)
          s22=sti(2,kk,nelem)
          s33=sti(3,kk,nelem)
          s12=sti(4,kk,nelem)
          s13=sti(5,kk,nelem)
          s23=sti(6,kk,nelem)
        endif
!     
!     calculating the temperature in the integration
!     point
!     
        t0l=0.d0
        t1l=0.d0
        if(ithermal(1).eq.1) then
          if(lakonl(4:5).eq.'8 ') then
            do i1=1,nope
              t0l=t0l+t0(konl(i1))/8.d0
              t1l=t1l+t1(konl(i1))/8.d0
            enddo
          elseif(lakonl(4:6).eq.'20 ') then
            nopered=20
            call lintemp(t0,konl,nopered,kk,t0l)
            call lintemp(t1,konl,nopered,kk,t1l)
          elseif(lakonl(4:6).eq.'10T') then
            call linscal10(t0,konl,t0l,null,shp)
            call linscal10(t1,konl,t1l,null,shp)
          else
            do i1=1,nope
              t0l=t0l+shp(4,i1)*t0(konl(i1))
              t1l=t1l+shp(4,i1)*t1(konl(i1))
            enddo
          endif
        elseif(ithermal(1).ge.2) then
          if(lakonl(4:5).eq.'8 ') then
            do i1=1,nope
              t0l=t0l+t0(konl(i1))/8.d0
              t1l=t1l+vold(0,konl(i1))/8.d0
            enddo
          elseif(lakonl(4:6).eq.'20 ') then
            nopered=20
            call lintemp_th0(t0,konl,nopered,kk,t0l,mi)
            call lintemp_th1(vold,konl,nopered,kk,t1l,mi)
          elseif(lakonl(4:6).eq.'10T') then
            call linscal10(t0,konl,t0l,null,shp)
            call linscal10(vold,konl,t1l,mi(2),shp)
          else
            do i1=1,nope
              t0l=t0l+shp(4,i1)*t0(konl(i1))
              t1l=t1l+shp(4,i1)*vold(0,konl(i1))
            enddo
          endif
        endif
        tt=t1l-t0l
!     
!     calculating the coordinates of the integration point
!     for material orientation purposes (for cylindrical
!     coordinate systems)
!     
        if(iorien.gt.0) then
          do j=1,3
            coords(j)=0.d0
            do i1=1,nope
              coords(j)=coords(j)+shp(4,i1)*co(j,konl(i1))
            enddo
          enddo
        endif
!     
!     for deformation plasticity: calculating the Jacobian
!     and the inverse of the deformation gradient
!     needed to convert the stiffness matrix in the spatial
!     frame of reference to the material frame
!     
        kode=nelcon(1,imat)
!     
!     material data and local stiffness matrix
!     
        istiff=1
        call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
     &       imat,amat,iorien,coords,orab,ntmat_,elas,rho,
     &       nelem,ithermal,alzero,mattyp,t0l,t1l,
     &       ihyper,istiff,elconloc,eth,kode,plicon,
     &       nplicon,plkcon,nplkcon,npmat_,
     &       plconloc,mi(1),dtime,kk,
     &       xstiff,ncmat_)
!     
        if(mattyp.eq.1) then
c     write(*,*) 'elastic co', elas(1),elas(2)
          e=elas(1)
          un=elas(2)
          um=e/(1.d0+un)
          al=un*um/(1.d0-2.d0*un)
          um=um/2.d0
        elseif(mattyp.eq.2) then
c     call orthotropic(elas,anisox)
        else
          call anisotropic(elas,anisox)
        endif
!     
!     initialisation for the body forces
!     
        om=omx*rho
        if(rhsi.eq.1) then
          if(nbody.ne.0) then
            do ii=1,3
              bodyf(ii)=bodyfx(ii)*rho
            enddo
          endif
        endif
!     
        if(buckling.eq.1) then
!     
!     buckling stresses
!     
          s11b=stx(1,kk,nelem)
          s22b=stx(2,kk,nelem)
          s33b=stx(3,kk,nelem)
          s12b=stx(4,kk,nelem)
          s13b=stx(5,kk,nelem)
          s23b=stx(6,kk,nelem)
!     
        endif
!     
!     incorporating the jacobian determinant in the shape
!     functions
!     
        xsjj=dsqrt(xsj)
!        
c     mortar start
        if(mortartrafoflag.eq.1) then
          do i1=1,nope
            shpj(1,i1)=shptil(1,i1)*xsjj
            shpj(2,i1)=shptil(2,i1)*xsjj
            shpj(3,i1)=shptil(3,i1)*xsjj
            shpj(4,i1)=shptil(4,i1)*xsj
          enddo
c     mortar end
        else
          do i1=1,nope
            shpj(1,i1)=shp(1,i1)*xsjj
            shpj(2,i1)=shp(2,i1)*xsjj
            shpj(3,i1)=shp(3,i1)*xsjj
            shpj(4,i1)=shp(4,i1)*xsj
          enddo
        endif
!        
        if(lakonl(4:5).eq.'8I') then
          do i1=1,nope
            shpjmass(4,i1)=shpmass(4,i1)*xsjmass
          enddo
        endif
!     
!     determination of the stiffness, and/or mass and/or
!     buckling matrix
!     
        if((stiffness.eq.1).or.(mass.eq.1).or.(buckling.eq.1).or.
     &       (coriolis.eq.1)) then
!     
          if(((iperturb(1).ne.1).and.(iperturb(2).ne.1)).or.
     &         (buckling.eq.1)) then
!     
!     calculating the total mass of the element for
!     lumping purposes: only for explicit
!     dynamic calculations
!     
            if((mass.eq.1).and.(iexpl.gt.1)) then
              summass=summass+rho*xsj*weight
            endif
!     
            jj1=1
            do jj=1,nope
!     
              ii1=1
              do ii=1,jj
!     
!     all products of the shape functions for a given ii
!     and jj
!     
                do i1=1,3
                  do j1=1,3
                    w(i1,j1)=shpj(i1,ii)*shpj(j1,jj)
                  enddo
                enddo
!     
!     the following section calculates the static
!     part of the stiffness matrix which, for buckling 
!     calculations, is done in a preliminary static
!     call
!     
                if(buckling.eq.0) then
!     
                  if(mattyp.eq.1) then
!     
                    s(ii1,jj1)=s(ii1,jj1)+(al*w(1,1)+
     &                   um*(2.d0*w(1,1)+w(2,2)+w(3,3)))*weight
!     
                    s(ii1,jj1+1)=s(ii1,jj1+1)+(al*w(1,2)+
     &                   um*w(2,1))*weight
                    s(ii1,jj1+2)=s(ii1,jj1+2)+(al*w(1,3)+
     &                   um*w(3,1))*weight
                    s(ii1+1,jj1)=s(ii1+1,jj1)+(al*w(2,1)+
     &                   um*w(1,2))*weight
                    s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+(al*w(2,2)+
     &                   um*(2.d0*w(2,2)+w(1,1)+w(3,3)))*weight
                    s(ii1+1,jj1+2)=s(ii1+1,jj1+2)+(al*w(2,3)+
     &                   um*w(3,2))*weight
                    s(ii1+2,jj1)=s(ii1+2,jj1)+(al*w(3,1)+
     &                   um*w(1,3))*weight
                    s(ii1+2,jj1+1)=s(ii1+2,jj1+1)+(al*w(3,2)+
     &                   um*w(2,3))*weight
                    s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+(al*w(3,3)+
     &                   um*(2.d0*w(3,3)+w(2,2)+w(1,1)))*weight
!     
                  elseif(mattyp.eq.2) then
!     
                    s(ii1,jj1)=s(ii1,jj1)+(elas(1)*w(1,1)+
     &                   elas(7)*w(2,2)+elas(8)*w(3,3))*weight
                    s(ii1,jj1+1)=s(ii1,jj1+1)+(elas(2)*w(1,2)+
     &                   elas(7)*w(2,1))*weight
                    s(ii1,jj1+2)=s(ii1,jj1+2)+(elas(4)*w(1,3)+
     &                   elas(8)*w(3,1))*weight
                    s(ii1+1,jj1)=s(ii1+1,jj1)+(elas(7)*w(1,2)+
     &                   elas(2)*w(2,1))*weight
                    s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+
     &                   (elas(7)*w(1,1)+
     &                   elas(3)*w(2,2)+elas(9)*w(3,3))*weight
                    s(ii1+1,jj1+2)=s(ii1+1,jj1+2)+
     &                   (elas(5)*w(2,3)+
     &                   elas(9)*w(3,2))*weight
                    s(ii1+2,jj1)=s(ii1+2,jj1)+
     &                   (elas(8)*w(1,3)+
     &                   elas(4)*w(3,1))*weight
                    s(ii1+2,jj1+1)=s(ii1+2,jj1+1)+
     &                   (elas(9)*w(2,3)+
     &                   elas(5)*w(3,2))*weight
                    s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+
     &                   (elas(8)*w(1,1)+
     &                   elas(9)*w(2,2)+elas(6)*w(3,3))*weight
!     
                  else
!     
                    do i1=1,3
                      do j1=1,3
                        do k1=1,3
                          do l1=1,3
                            s(ii1+i1-1,jj1+j1-1)=
     &                           s(ii1+i1-1,jj1+j1-1)
     &                           +anisox(i1,k1,j1,l1)
     &                           *w(k1,l1)*weight
                          enddo
                        enddo
                      enddo
                    enddo
!     
                  endif
!     
!     mass matrix
!     
                  if(mass.eq.1) then
                    if(lakonl(4:5).ne.'8I') then
c     mortar start
                      if(mortartrafoflag.eq.1) then
                        sm(ii1,jj1)=sm(ii1,jj1)
     &                       +rho*shpj(4,ii)*shptil(4,jj)*weight
                      else
c     mortar end                      
                        sm(ii1,jj1)=sm(ii1,jj1)
     &                       +rho*shpj(4,ii)*shp(4,jj)*weight
                      endif
                    else
                      sm(ii1,jj1)=sm(ii1,jj1)
     &                     +rho*shpjmass(4,ii)*shpmass(4,jj)*weight
                    endif
                    sm(ii1+1,jj1+1)=sm(ii1,jj1)
                    sm(ii1+2,jj1+2)=sm(ii1,jj1)
                  endif
!     
!     Coriolis matrix
!     
                  if(coriolis.eq.1) then
c     mortar start
                    if(mortartrafoflag.eq.1) then
                      dmass=2.d0*
     &                     rho*shpj(4,ii)*shptil(4,jj)*weight
     &                     *dsqrt(omx)
                    else
c     mortar end                      
                      dmass=2.d0*
     &                     rho*shpj(4,ii)*shp(4,jj)*weight*dsqrt(omx)
                    endif
                    sm(ii1,jj1+1)=sm(ii1,jj1+1)-p2(3)*dmass
                    sm(ii1,jj1+2)=sm(ii1,jj1+2)+p2(2)*dmass
                    sm(ii1+1,jj1)=sm(ii1+1,jj1)+p2(3)*dmass
                    sm(ii1+1,jj1+2)=sm(ii1+1,jj1+2)-p2(1)*dmass
                    sm(ii1+2,jj1)=sm(ii1+2,jj1)-p2(2)*dmass
                    sm(ii1+2,jj1+1)=sm(ii1+2,jj1+1)+p2(1)*dmass
                  endif
!     
                else
!     
!     buckling matrix  
!     
                  senergyb=
     &                 (s11b*w(1,1)+s12b*(w(1,2)+w(2,1))
     &                 +s13b*(w(1,3)+w(3,1))+s22b*w(2,2)
     &                 +s23b*(w(2,3)+w(3,2))+s33b*w(3,3))*weight
                  sm(ii1,jj1)=sm(ii1,jj1)-senergyb
                  sm(ii1+1,jj1+1)=sm(ii1+1,jj1+1)-senergyb
                  sm(ii1+2,jj1+2)=sm(ii1+2,jj1+2)-senergyb
!     
                endif
!     
                ii1=ii1+3
              enddo
              jj1=jj1+3
            enddo
          else
!     
!     stiffness matrix for static and modal
!     2nd order calculations
!     
!     large displacement stiffness
!     
            do i1=1,3
              do j1=1,3
                vo(i1,j1)=0.d0
                do k1=1,nope
                  vo(i1,j1)=vo(i1,j1)+shp(j1,k1)*voldl(i1,k1)
                enddo
              enddo
            enddo
!     
            if(mattyp.eq.1) then
              call wcoef(v,vo,al,um)
            endif
!     
!     calculating the total mass of the element for
!     lumping purposes: only for explicit
!     dynamic calculations
!     
            if((mass.eq.1).and.(iexpl.gt.1)) then
              summass=summass+rho*xsj*weight
            endif
!     
            jj1=1
            do jj=1,nope
!     
              ii1=1
              do ii=1,jj
!     
!     all products of the shape functions for a given ii
!     and jj
!     
                do i1=1,3
                  do j1=1,3
                    w(i1,j1)=shpj(i1,ii)*shpj(j1,jj)
                  enddo
                enddo
!     
                if(mattyp.eq.1) then
!     
                  do m1=1,3
                    do m2=1,3
                      do m3=1,3
                        do m4=1,3
                          s(ii1+m2-1,jj1+m1-1)=
     &                         s(ii1+m2-1,jj1+m1-1)
     &                         +v(m4,m3,m2,m1)*w(m4,m3)*weight
                        enddo
                      enddo
                    enddo
                  enddo
!     
                elseif(mattyp.eq.2) then
!     
                  call orthonl(w,vo,elas,s,ii1,jj1,weight)
!     
                else
!     
                  call anisonl(w,vo,elas,s,ii1,jj1,weight)
!     
                endif
!     
!     stress stiffness
!     
                senergy=
     &               (s11*w(1,1)+s12*(w(1,2)+w(2,1))
     &               +s13*(w(1,3)+w(3,1))+s22*w(2,2)
     &               +s23*(w(2,3)+w(3,2))+s33*w(3,3))*weight
                s(ii1,jj1)=s(ii1,jj1)+senergy
                s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+senergy
                s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+senergy
!     
!     stiffness contribution of centrifugal forces
!     
                if((mass.eq.1).and.(om.gt.0.d0)) then
c     mortar start                  
                  if(mortartrafoflag.eq.1) then
                    dmass=shpj(4,ii)*shptil(4,jj)*weight*om
                  else
c     mortar end                  
                    dmass=shpj(4,ii)*shp(4,jj)*weight*om
                  endif
                  do m1=1,3
                    s(ii1+m1-1,jj1+m1-1)=s(ii1+m1-1,jj1+m1-1)-
     &                   dmass
                    do m2=1,3
                      s(ii1+m1-1,jj1+m2-1)=s(ii1+m1-1,jj1+m2-1)+
     &                     dmass*p2(m1)*p2(m2)
                    enddo
                  enddo
                endif
!     
!     mass matrix
!     
                if(mass.eq.1) then
                  if(lakonl(4:5).ne.'8I') then
c     mortar start                    
                    if(mortartrafoflag.eq.1) then
                      sm(ii1,jj1)=sm(ii1,jj1)
     &                     +rho*shpj(4,ii)*shptil(4,jj)*weight
                    else
c     mortar end                    
                      sm(ii1,jj1)=sm(ii1,jj1)
     &                     +rho*shpj(4,ii)*shp(4,jj)*weight
                    endif
                  else
                    sm(ii1,jj1)=sm(ii1,jj1)
     &                   +rho*shpjmass(4,ii)*shpmass(4,jj)*weight
                  endif
                  sm(ii1+1,jj1+1)=sm(ii1,jj1)
                  sm(ii1+2,jj1+2)=sm(ii1,jj1)
                endif
!     
!     Coriolis matrix
!     
                if(coriolis.eq.1) then
c     mortar start                    
                  if(mortartrafoflag.eq.1) then
                    dmass=2.d0*
     &                   rho*shpj(4,ii)*shptil(4,jj)*
     &                   weight*dsqrt(omx)
                  else
c     mortar end                    
                    dmass=2.d0*
     &                   rho*shpj(4,ii)*shp(4,jj)*weight*dsqrt(omx)
                  endif
                  sm(ii1,jj1+1)=sm(ii1,jj1+1)-p2(3)*dmass
                  sm(ii1,jj1+2)=sm(ii1,jj1+2)+p2(2)*dmass
                  sm(ii1+1,jj1)=sm(ii1+1,jj1)+p2(3)*dmass
                  sm(ii1+1,jj1+2)=sm(ii1+1,jj1+2)-p2(1)*dmass
                  sm(ii1+2,jj1)=sm(ii1+2,jj1)-p2(2)*dmass
                  sm(ii1+2,jj1+1)=sm(ii1+2,jj1+1)+p2(1)*dmass
                endif
!     
                ii1=ii1+3
              enddo
              jj1=jj1+3
            enddo
          endif
!     
        endif
!     
!     add hourglass control stiffnesses: C3D8R only. 
        if(lakonl(1:5).eq.'C3D8R') then 
          call hgstiffness(s,elas,a,gs)
        endif
!     
!     computation of the right hand side
!     
        if(rhsi.eq.1) then
!     
!     distributed body flux
!     
          if(nload.gt.0) then
            call nident2(nelemload,nelem,nload,id)
            do
              if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
              if((sideload(id)(1:2).ne.'BX').and.
     &             (sideload(id)(1:2).ne.'BY').and.
     &             (sideload(id)(1:2).ne.'BZ')) then
                id=id-1
                cycle
              endif
              if(sideload(id)(3:4).eq.'NU') then
                do j=1,3
                  coords(j)=0.d0
                  do i1=1,nope
                    coords(j)=coords(j)+
     &                   shp(4,i1)*xl(j,i1)
                  enddo
                enddo
                if(sideload(id)(1:2).eq.'BX') then
                  jltyp=1
                elseif(sideload(id)(1:2).eq.'BY') then
                  jltyp=2
                elseif(sideload(id)(1:2).eq.'BZ') then
                  jltyp=3
                endif
                iscale=1
                call dload(xload(1,id),istep,iinc,tvar,nelem,i,
     &               layer,kspt,coords,jltyp,sideload(id),vold,co,
     &               lakonl,konl,ipompc,nodempc,coefmpc,nmpc,ikmpc,
     &               ilmpc,iscale,veold,rho,amat,mi)
                if((nmethod.eq.1).and.(iscale.ne.0))
     &               xload(1,id)=xloadold(1,id)+
     &               (xload(1,id)-xloadold(1,id))*reltime
              endif
              jj1=1
              do jj=1,nope
                if(sideload(id)(1:2).eq.'BX')
     &               ff(jj1)=ff(jj1)+xload(1,id)*shpj(4,jj)*weight
                if(sideload(id)(1:2).eq.'BY')
     &               ff(jj1+1)=ff(jj1+1)+xload(1,id)*shpj(4,jj)*weight
                if(sideload(id)(1:2).eq.'BZ')
     &               ff(jj1+2)=ff(jj1+2)+xload(1,id)*shpj(4,jj)*weight
                jj1=jj1+3
              enddo
              id=id-1
            enddo
          endif
!     
!     body forces
!     
          if(nbody.ne.0) then
            if(om.gt.0.d0) then
              do i1=1,3
!     
!     computation of the global coordinates of the gauss
!     point
!     
                q(i1)=0.d0
                if((iperturb(1).ne.1).and.(iperturb(2).ne.1)) then
                  do j1=1,nope
                    q(i1)=q(i1)+shp(4,j1)*xl(i1,j1)
                  enddo
                else
                  do j1=1,nope
                    q(i1)=q(i1)+shp(4,j1)*
     &                   (xl(i1,j1)+voldl(i1,j1))
                  enddo
                endif
!     
                q(i1)=q(i1)-p1(i1)
              enddo
              const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!     
!     inclusion of the centrifugal force into the body force
!     
              do i1=1,3
                bf(i1)=bodyf(i1)+(q(i1)-const*p2(i1))*om
              enddo
            else
              do i1=1,3
                bf(i1)=bodyf(i1)
              enddo
            endif
            jj1=1
            do jj=1,nope
              ff(jj1)=ff(jj1)+bf(1)*shpj(4,jj)*weight
              ff(jj1+1)=ff(jj1+1)+bf(2)*shpj(4,jj)*weight
              ff(jj1+2)=ff(jj1+2)+bf(3)*shpj(4,jj)*weight
              jj1=jj1+3
            enddo
          endif
!     
        endif
!     
      enddo
!     
c     write(*,*) nelem
c     write(*,'(6(1x,e11.4))') ((s(i1,j1),i1=1,j1),j1=1,60)
c     write(*,*)
c     
      if((buckling.eq.0).and.(nload.ne.0)) then
!     
!     distributed loads
!     
        call nident2(nelemload,nelem,nload,id)
        do
          if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
          if(sideload(id)(1:1).ne.'P') then
            id=id-1
            cycle
          endif
c     read(sideload(id)(2:2),'(i1)') ig
          ig=ichar(sideload(id)(2:2))-48
!     
!     
!     treatment of wedge faces
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
c     Bernhardi start
          if((nope.eq.20).or.(nope.eq.8).or.
     &         (nope.eq.11)) then
c     Bernhardi end
            if((iperturb(1).ne.1).and.(iperturb(2).ne.1)) then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                enddo
              enddo
            else
              if(mass.eq.1) then
                do i=1,nopes
                  do j=1,3
                    xl1(j,i)=co(j,konl(ifaceq(i,ig)))
                  enddo
                enddo
              endif
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))+
     &                 vold(j,konl(ifaceq(i,ig)))
                enddo
              enddo
            endif
          elseif((nope.eq.10).or.(nope.eq.4)) then
            if((iperturb(1).ne.1).and.(iperturb(2).ne.1)) then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))
                enddo
              enddo
            else
              if(mass.eq.1) then
                do i=1,nopes
                  do j=1,3
                    xl1(j,i)=co(j,konl(ifacet(i,ig)))
                  enddo
                enddo
              endif
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))+
     &                 vold(j,konl(ifacet(i,ig)))
                enddo
              enddo
            endif
          else
            if((iperturb(1).ne.1).and.(iperturb(2).ne.1)) then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))
                enddo
              enddo
            else
              if(mass.eq.1) then
                do i=1,nopes
                  do j=1,3
                    xl1(j,i)=co(j,konl(ifacew(i,ig)))
                  enddo
                enddo
              endif
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))+
     &                 vold(j,konl(ifacew(i,ig)))
                enddo
              enddo
            endif
          endif
!     
c     mortar start
          if(mortartrafoflag.eq.1) then
!
!     generate autloc1
!
            if(islavelinv(nelem).gt.0) then
              if((nope.eq.20).or.(nope.eq.10).or.(nope.eq.15)) then
                jqtloc1(1)=1
                ii=1
                do i1=1,nopes
                  if(nope.eq.20) then
                    ipointer=ifaceq(i1,ig)
                  elseif(nope.eq.10) then
                    ipointer=ifacet(i1,ig)
                  else
                    ipointer=ifacew(i1,ig)
                  endif
                  node1=ipointer
                  do j1=jqtloc(node1),jqtloc(node1+1)-1
                    node2=irowtloc(j1)
                    do j2=1,nopes
                      if(nope.eq.20) then
                        ipointer=ifaceq(j2,ig)
                      elseif(nope.eq.10) then
                        ipointer=ifacet(j2,ig)
                      else
                        ipointer=ifacew(j2,ig)
                      endif
                      if(ipointer.eq.node2) then
                        autloc1(ii)=autloc(j1)
                        irowtloc1(ii)=j2
                        ii=ii+1
                      endif
                    enddo
                  enddo
                  jqtloc1(i1+1)=ii
                enddo
              endif
            endif
          endif
c     mortar end
!          
          do i=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
              xi=gauss2d1(1,i)
              et=gauss2d1(2,i)
              weight=weight2d1(i)
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
              xi=gauss2d2(1,i)
              et=gauss2d2(2,i)
              weight=weight2d2(i)
            elseif(lakonl(4:4).eq.'2') then
              xi=gauss2d3(1,i)
              et=gauss2d3(2,i)
              weight=weight2d3(i)
            elseif((lakonl(4:5).eq.'10').or.
     &             ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
              xi=gauss2d5(1,i)
              et=gauss2d5(2,i)
              weight=weight2d5(i)
            elseif((lakonl(4:4).eq.'4').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
              xi=gauss2d4(1,i)
              et=gauss2d4(2,i)
              weight=weight2d4(i)
            endif
!     
            if(rhsi.eq.1) then
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
c     mortar start
              if(mortartrafoflag.eq.1) then
                do i1=1,nopes
                  shptil2(1,i1)=shp2(1,i1)
                  shptil2(2,i1)=shp2(2,i1)
                  shptil2(3,i1)=shp2(3,i1)
                  shptil2(4,i1)=shp2(4,i1)
                enddo
                if(islavelinv(nelem).gt.0) then
                  if((nopes.eq.8).or.(nopes.eq.6)) then
                    do i1=1,nopes
                      if(jqtloc1(i1+1)-jqtloc1(i1).gt.0) then
                        shptil2(1,i1)=0.0
                        shptil2(2,i1)=0.0
                        shptil2(3,i1)=0.0
                        shptil2(4,i1)=0.0
                      endif
                      do j1=jqtloc1(i1),jqtloc1(i1+1)-1
                        j2=irowtloc1(j1)
                        shptil2(1,i1)=shptil2(1,i1)
     &                       +autloc1(j1)*shp2(1,j2)
                        shptil2(2,i1)=shptil2(2,i1)
     &                       +autloc1(j1)*shp2(2,j2)
                        shptil2(3,i1)=shptil2(3,i1)
     &                       +autloc1(j1)*shp2(3,j2)
                        shptil2(4,i1)=shptil2(4,i1)
     &                       +autloc1(j1)*shp2(4,j2)
                      enddo
                    enddo
                  endif
                endif
              endif
c     mortar end
!     
!     for nonuniform load: determine the coordinates of the
!     point (transferred into the user subroutine)
!     
              if(sideload(id)(3:4).eq.'NU') then
                do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                    coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
                enddo
c     read(sideload(id)(2:2),'(i1)') jltyp
                jltyp=ichar(sideload(id)(2:2))-48
                jltyp=jltyp+20
                iscale=1
                call dload(xload(1,id),istep,iinc,tvar,nelem,i,layer,
     &               kspt,coords,jltyp,sideload(id),vold,co,lakonl,
     &               konl,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &               iscale,veold,rho,amat,mi)
                if((nmethod.eq.1).and.(iscale.ne.0))
     &               xload(1,id)=xloadold(1,id)+
     &               (xload(1,id)-xloadold(1,id))*reltime
              elseif(sideload(id)(3:4).eq.'SM') then
!     
!     submodel boundary: interpolation from the
!     global model
!     
                do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                    coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
                enddo
c     read(sideload(id)(2:2),'(i1)') jltyp
                jltyp=ichar(sideload(id)(2:2))-48
!     
                entity='T'
                six=6
                do k=1,6
                  iselect(k)=k+4
                enddo
                iface=10*nelem+jltyp
                call interpolsubmodel(integerglob,doubleglob,stress,
     &               coords,iselect,six,iface,tieset,istartset,
     &               iendset,ialset,ntie,entity,set,nset)
c     write(*,*) 'e_c3d ',(stress(k),k=1,6)
!     
!     cave: stress order of cgx: xx,yy,zz,xy,yz,xz
!     
                t(1)=stress(1)*xsj2(1)+stress(4)*xsj2(2)+
     &               stress(6)*xsj2(3)
                t(2)=stress(4)*xsj2(1)+stress(2)*xsj2(2)+
     &               stress(5)*xsj2(3)
                t(3)=stress(6)*xsj2(1)+stress(5)*xsj2(2)+
     &               stress(3)*xsj2(3)
!     
                xload(1,id)=-1.d0
                do k=1,3
                  xsj2(k)=t(k)
                enddo
              endif
!     
              do k=1,nopes
c     Bernhardi start
                if((nope.eq.20).or.(nope.eq.8).or.
     &               (nope.eq.11)) then
c     Bernhardi end
                  ipointer=(ifaceq(k,ig)-1)*3
                elseif((nope.eq.10).or.(nope.eq.4)) then
                  ipointer=(ifacet(k,ig)-1)*3
                else
                  ipointer=(ifacew(k,ig)-1)*3
                endif
c     mortar start                
                if(mortartrafoflag.eq.1) then
                  ff(ipointer+1)=ff(ipointer+1)
     &                 -shptil2(4,k)*xload(1,id)
     &                 *xsj2(1)*weight
                  ff(ipointer+2)=ff(ipointer+2)
     &                 -shptil2(4,k)*xload(1,id)
     &                 *xsj2(2)*weight
                  ff(ipointer+3)=ff(ipointer+3)
     &                 -shptil2(4,k)*xload(1,id)
     &                 *xsj2(3)*weight
                else
c     mortar end                  
                  ff(ipointer+1)=ff(ipointer+1)-shp2(4,k)*xload(1,id)
     &                 *xsj2(1)*weight
                  ff(ipointer+2)=ff(ipointer+2)-shp2(4,k)*xload(1,id)
     &                 *xsj2(2)*weight
                  ff(ipointer+3)=ff(ipointer+3)-shp2(4,k)*xload(1,id)
     &                 *xsj2(3)*weight
                endif
              enddo
!     
!     stiffness< contribution of the distributed load 
!     reference: Dhondt G., The Finite Element Method for
!     three-dimensional thermomechanical Applications,
!     Wiley, 2004, p 153, eqn. (3.54).
!     
            elseif((mass.eq.1).and.
     &             ((iperturb(1).eq.1).or.(iperturb(2).eq.1))) then
              if(nopes.eq.8) then
                call shape8q(xi,et,xl1,xsj2,xs2,shp2,iflag)
              elseif(nopes.eq.4) then
                call shape4q(xi,et,xl1,xsj2,xs2,shp2,iflag)
              elseif(nopes.eq.6) then
                call shape6tri(xi,et,xl1,xsj2,xs2,shp2,iflag)
              else
                call shape3tri(xi,et,xl1,xsj2,xs2,shp2,iflag)
              endif
c     mortar start
              if(mortartrafoflag.eq.1) then
                do i1=1,nopes
                  shptil2(1,i1)=shp2(1,i1)
                  shptil2(2,i1)=shp2(2,i1)
                  shptil2(3,i1)=shp2(3,i1)
                  shptil2(4,i1)=shp2(4,i1)
                enddo
                if(islavelinv(nelem).gt.0) then
                  if((nopes.eq.8).or.(nopes.eq.6)) then
                    do i1=1,nopes
                      if(jqtloc1(i1+1)-jqtloc1(i1).gt.0) then
                        shptil2(1,i1)=0.0
                        shptil2(2,i1)=0.0
                        shptil2(3,i1)=0.0
                        shptil2(4,i1)=0.0
                      endif
                      do j1=jqtloc1(i1),jqtloc1(i1+1)-1
                        j2=irowtloc1(j1)
                        shptil2(1,i1)=shptil2(1,i1)
     &                       +autloc1(j1)*shp2(1,j2)
                        shptil2(2,i1)=shptil2(2,i1)
     &                       +autloc1(j1)*shp2(2,j2)
                        shptil2(3,i1)=shptil2(3,i1)
     &                       +autloc1(j1)*shp2(3,j2)
                        shptil2(4,i1)=shptil2(4,i1)
     &                       +autloc1(j1)*shp2(4,j2)
                      enddo
                    enddo
                  endif
                endif
              endif
c     mortar end
!     
!     for nonuniform load: determine the coordinates of the
!     point (transferred into the user subroutine)
!     
              if(sideload(id)(3:4).eq.'NU') then
                do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                    coords(k)=coords(k)+xl1(k,j)*shp2(4,j)
                  enddo
                enddo
c     read(sideload(id)(2:2),'(i1)') jltyp
                jltyp=ichar(sideload(id)(2:2))-48
                jltyp=jltyp+20
                iscale=1
                call dload(xload(1,id),istep,iinc,tvar,nelem,i,layer,
     &               kspt,coords,jltyp,sideload(id),vold,co,lakonl,
     &               konl,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &               iscale,veold,rho,amat,mi)
                if((nmethod.eq.1).and.(iscale.ne.0))
     &               xload(1,id)=xloadold(1,id)+
     &               (xload(1,id)-xloadold(1,id))*reltime
              elseif(sideload(id)(3:4).eq.'SM') then
!     
!     submodel boundary: interpolation from the
!     global model
!     
                do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                    coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
                enddo
c     read(sideload(id)(2:2),'(i1)') jltyp
                jltyp=ichar(sideload(id)(2:2))-48
!     
                entity='T'
                six=6
                do k=1,6
                  iselect(k)=k+4
                enddo
                iface=10*nelem+jltyp
                call interpolsubmodel(integerglob,doubleglob,stress,
     &               coords,iselect,six,iface,tieset,istartset,
     &               iendset,ialset,ntie,entity,set,nset)
c     write(*,*) 'e_c3d ',(stress(k),k=1,6)
!     
                stre(1,1)=stress(1)
                stre(1,2)=stress(4)
                stre(1,3)=stress(6)
                stre(2,1)=stress(4)
                stre(2,2)=stress(2)
                stre(2,3)=stress(5)
                stre(3,1)=stress(6)
                stre(3,2)=stress(5)
                stre(3,3)=stress(3)
              endif
!     
!     calculation of the deformation gradient
!     
              do k=1,3
                do l=1,3
                  xkl(k,l)=0.d0
                  do ii=1,nopes
                    xkl(k,l)=xkl(k,l)+shp2(l,ii)*xl2(k,ii)
                  enddo
                enddo
              enddo
!     
              do ii=1,nopes
c     Bernhardi start
                if((nope.eq.20).or.(nope.eq.8).or.
     &               (nope.eq.11)) then
c     Bernhardi end
                  ipointeri=(ifaceq(ii,ig)-1)*3
                elseif((nope.eq.10).or.(nope.eq.4)) then
                  ipointeri=(ifacet(ii,ig)-1)*3
                else
                  ipointeri=(ifacew(ii,ig)-1)*3
                endif
                do jj=1,nopes
c     Bernhardi start
                  if((nope.eq.20).or.(nope.eq.8)
     &                 .or.(nope.eq.11)) then
c     Bernhardi end
                    ipointerj=(ifaceq(jj,ig)-1)*3
                  elseif((nope.eq.10).or.(nope.eq.4)) then
                    ipointerj=(ifacet(jj,ig)-1)*3
                  else
                    ipointerj=(ifacew(jj,ig)-1)*3
                  endif
!     
!     if no submodel: only pressure
!     else: complete stress vector
!     
                  if(sideload(id)(3:4).ne.'SM') then
                    do k=1,3
                      do l=1,3
                        if(k.eq.l) cycle
                        eknlsign=1.d0
                        if(k*l.eq.2) then
                          n=3
                          if(k.lt.l) eknlsign=-1.d0
                        elseif(k*l.eq.3) then
                          n=2
                          if(k.gt.l) eknlsign=-1.d0
                        else
                          n=1
                          if(k.lt.l) eknlsign=-1.d0
                        endif
c     mortar start                        
                        if(mortartrafoflag.eq.1) then
                          term=weight*xload(1,id)*shptil2(4,ii)*
     &                         eknlsign*(xsj2(1)*
     &                         (xkl(n,2)*shptil2(3,jj)-xkl(n,3)*
     &                         shptil2(2,jj))+xsj2(2)*
     &                         (xkl(n,3)*shptil2(1,jj)-xkl(n,1)*
     &                         shptil2(3,jj))+xsj2(3)*
     &                         (xkl(n,1)*shptil2(2,jj)-xkl(n,2)*
     &                         shptil2(1,jj)))
                        else
c     mortar end                                  
                          term=weight*xload(1,id)*shp2(4,ii)*
     &                         eknlsign*(xsj2(1)*
     &                         (xkl(n,2)*shp2(3,jj)-xkl(n,3)*
     &                         shp2(2,jj))+xsj2(2)*
     &                         (xkl(n,3)*shp2(1,jj)-xkl(n,1)*
     &                         shp2(3,jj))+xsj2(3)*
     &                         (xkl(n,1)*shp2(2,jj)-xkl(n,2)*
     &                         shp2(1,jj)))
                        endif
                        s(ipointeri+k,ipointerj+l)=
     &                       s(ipointeri+k,ipointerj+l)+term/2.d0
                        s(ipointerj+l,ipointeri+k)=
     &                       s(ipointerj+l,ipointeri+k)+term/2.d0
                      enddo
                    enddo
                  else
                    do kk=1,3
                      do k=1,3
                        do l=1,3
                          if(k.eq.l) cycle
                          eknlsign=1.d0
                          if(k*l.eq.2) then
                            n=3
                            if(k.lt.l) eknlsign=-1.d0
                          elseif(k*l.eq.3) then
                            n=2
                            if(k.gt.l) eknlsign=-1.d0
                          else
                            n=1
                            if(k.lt.l) eknlsign=-1.d0
                          endif
c     mortar start                          
                          if(mortartrafoflag.eq.1) then
                            term=-weight*stre(kk,k)*
     &                           shptil2(4,ii)*
     &                           eknlsign*(xsj2(1)*
     &                           (xkl(n,2)*shptil2(3,jj)-xkl(n,3)*
     &                           shptil2(2,jj))+xsj2(2)*
     &                           (xkl(n,3)*shptil2(1,jj)-xkl(n,1)*
     &                           shptil2(3,jj))+xsj2(3)*
     &                           (xkl(n,1)*shptil2(2,jj)-xkl(n,2)*
     &                           shptil2(1,jj)))
                          else
c     mortar end                                    
                            term=-weight*stre(kk,k)*shp2(4,ii)*
     &                           eknlsign*(xsj2(1)*
     &                           (xkl(n,2)*shp2(3,jj)-xkl(n,3)*
     &                           shp2(2,jj))+xsj2(2)*
     &                           (xkl(n,3)*shp2(1,jj)-xkl(n,1)*
     &                           shp2(3,jj))+xsj2(3)*
     &                           (xkl(n,1)*shp2(2,jj)-xkl(n,2)*
     &                           shp2(1,jj)))
                          endif
                          s(ipointeri+kk,ipointerj+l)=
     &                         s(ipointeri+kk,ipointerj+l)+term/2.d0
                          s(ipointerj+l,ipointeri+kk)=
     &                         s(ipointerj+l,ipointeri+kk)+term/2.d0
                        enddo
                      enddo
                    enddo
                  endif
                enddo
              enddo
!     
            endif
          enddo
!     
          id=id-1
        enddo
      endif
!     
!     for axially symmetric and plane stress/strain elements: 
!     complete s and sm
!     
      if(intscheme.eq.0) then
        if(((lakonl(4:5).eq.'8 ').or.
     &       ((lakonl(4:6).eq.'20R').and.(lakonl(7:8).ne.'BR'))).and.
     &       ((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &       (lakonl(7:7).eq.'E'))) then
          do i=1,60
            do j=i,60
              k=abs(iperm(i))
              l=abs(iperm(j))
              if(k.gt.l) then
                m=k
                k=l
                l=m
              endif
              sax(i,j)=s(k,l)*iperm(i)*iperm(j)/(k*l)
            enddo
          enddo
          do i=1,60
            do j=i,60
              s(i,j)=s(i,j)+sax(i,j)
            enddo
          enddo
!     
          if((nload.ne.0).or.(nbody.ne.0)) then
            do i=1,60
              k=abs(iperm(i))
              ffax(i)=ff(k)*iperm(i)/k
            enddo
            do i=1,60
              ff(i)=ff(i)+ffax(i)
            enddo
          endif
!     
          if(mass.eq.1) then
            summass=2.d0*summass
            do i=1,60
              do j=i,60
                k=abs(iperm(i))
                l=abs(iperm(j))
                if(k.gt.l) then
                  m=k
                  k=l
                  l=m
                endif
                sax(i,j)=sm(k,l)*iperm(i)*iperm(j)/(k*l)
              enddo
            enddo
            do i=1,60
              do j=i,60
                sm(i,j)=sm(i,j)+sax(i,j)
              enddo
            enddo
          endif
        endif
      endif
!
!     lumping the mass matrix for explicit dynamics
!     (not for the massless method)
!
      if((mass.eq.1).and.(iexpl.gt.1).and.(mortar.ne.-1)) then
!     
!     scaling the diagonal terms of the mass matrix such that the total mass
!     is right (LUMPING; for explicit dynamic calculations)
!     
        sume=0.d0
        summ=0.d0
        do i=1,3*nopev,3
          sume=sume+sm(i,i)
        enddo
        do i=3*nopev+1,3*nope,3
          summ=summ+sm(i,i)
        enddo
!     
        if(nope.eq.20) then
          alp=.2917d0
        elseif(nope.eq.10) then
          alp=0.1203d0
        elseif(nope.eq.15) then
          alp=0.2141d0
        elseif(nope.eq.11) then
          alp=0.1852d0
        endif
!     
        if((nope.eq.20).or.(nope.eq.10).or.
     &       (nope.eq.15).or.(nope.eq.11)) then
          factore=summass*alp/(1.d0+alp)/sume
          factorm=summass/(1.d0+alp)/summ
        else
          factore=summass/sume
        endif
!     
        do i=1,3*nopev,3
          sm(i,i)=sm(i,i)*factore
          sm(i+1,i+1)=sm(i,i)
          sm(i+2,i+2)=sm(i,i)
        enddo
        do i=3*nopev+1,3*nope,3
          sm(i,i)=sm(i,i)*factorm
          sm(i+1,i+1)=sm(i,i)
          sm(i+2,i+2)=sm(i,i)
        enddo
!     
c        if((mscalmethod.ne.1).and.(mscalmethod.ne.3)) then
!     
!     setting all off-diagonal terms to zero     
!     
          do i=1,3*nope
            do j=1,3*nope
              if(i.eq.j) cycle
              sm(i,j)=0.d0
            enddo
          enddo
c        endif
      endif
!
!     mscalmethod = 1 or 3: selective mass scaling SMS
!     
      if((mass.eq.1).and.(iexpl.gt.1)) then
        if((mscalmethod.eq.1).or.(mscalmethod.eq.3)) then
!     
!     beta = smscalel
!     
          smfactor=smscalel*summass/((nope-1)*nope)
!     
          do i=1,3*nope
c            do j=1,3*nope
c              if(i.ne.j) then
c!     set non diagonals to zero           
c                sm(i,j)=0.d0
c!     diagonal terms of M for SMS
c              else
c                sm(i,j)=sm(i,j)+(nope-1)*smfactor
                sm(i,i)=sm(i,i)+(nope-1)*smfactor
c              endif
c            enddo
          enddo
!     
!     nondiagonal terms of M for SMS
!     
!         no massless method: regular mass matrix was lumped     
!     
c          if(mortar.ne.-1) then
c            i=0
c            do j=0,nope*3-1,3
c              i=i+1
c              do k=1,(nope-i)
c                do l=1,3
c                  sm(j+l+k*3,j+l)=-smfactor
c                  sm(j+l,j+l+k*3)=-smfactor
c                enddo
c              enddo
c            enddo
c          else
c!     
c!          massless method: regular mass matrix was not lumped     
c!     
            i=0
            do j=0,nope*3-1,3
              i=i+1
              do k=1,(nope-i)
                do l=1,3
                  sm(j+l+k*3,j+l)=sm(j+l+k*3,j+l)-smfactor
                  sm(j+l,j+l+k*3)=sm(j+l,j+l+k*3)-smfactor
                enddo
              enddo
            enddo
c          endif
!     to calculate additional energy in resultsmech.f:
          smscalel=smfactor
        endif
!     
      endif
!     
      return
      end
