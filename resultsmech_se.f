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
      subroutine resultsmech_se(co,kon,ipkon,lakon,ne,v,
     &     stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &     ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,eme,iperturb,fn,iout,vold,nmethod,
     &     veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &     xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &     ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &     springarea,reltime,calcul_fn,calcul_cauchy,nener,
     &     ikin,ne0,thicke,emeini,pslavsurf,
     &     pmastsurf,mortar,clearini,nea,neb,ielprop,prop,dfn,
     &     idesvar,nodedesi,fn0,sti,icoordinate,dxstiff,
     &     ialdesi,xdesi)
!     
!     calculates stresses and the material tangent at the integration
!     points and the internal forces at the nodes
!     
      implicit none
!     
      integer cauchy
!     
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!     
      integer kon(*),konl(26),nea,neb,mi(*),mint2d,nopes,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &     ielorien(mi(3),*),ntmat_,ipkon(*),ne0,iflag,null,
     &     istep,iinc,mt,ne,mattyp,ithermal(*),iprestr,i,j,k,m1,m2,jj,
     &     i1,m3,m4,kk,nener,indexe,nope,norien,iperturb(*),iout,
     &     icmd,ihyper,nmethod,kode,imat,mint3d,iorien,ielas,
     &     istiff,ncmat_,nstate_,ikin,ilayer,nlayer,ki,kl,ielprop(*),
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &     calcul_cauchy,nopered,mortar,jfaces,igauss,
     &     idesvar,node,nodedesi(*),kscale,idir,nlgeom_undo,
     &     iactive,icoordinate,ialdesi(*),ii,nr,mscalmethod,
     &     node1,node2,ifaceqexp(2,20),ifacewexp(2,15)
!     
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),stiini(6,mi(1),*),
     &     stx(6,mi(1),*),xl(3,26),vl(0:mi(2),26),stre(6),prop(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),xs2(3,7),
     &     alcon(0:6,ntmat_,*),vini(0:mi(2),*),thickness,
     &     alzero(*),orab(7,*),elas(21),rho,fn(0:mi(2),*),
     &     fnl(3,10),skl(3,3),beta(6),xl2(3,8),qa(4),
     &     vkl(0:3,3),t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &     ckl(3,3),vold(0:mi(2),*),eloc(9),veold(0:mi(2),*),
     &     springarea(2,*),elconloc(ncmat_),eth(6),xkl(3,3),
     &     voldl(0:mi(2),26),xikl(3,3),ener(2,mi(1),*),emec(6),
     &     eei(6,mi(1),*),enerini(2,mi(1),*),senergy,venergy,
     &     emec0(6),veoldl(0:mi(2),26),xsj2(3),shp2(7,8),
     &     e,un,al,um,am1,xi,et,ze,tt,exx,eyy,ezz,exy,exz,eyz,
     &     xsj,vj,t0l,t1l,dtime,weight,pgauss(3),vij,time,ttime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &     vokl(3,3),xstateini(nstate_,mi(1),*),vikl(3,3),
     &     gs(8,4),a,reltime,tlayer(4),dlayer(4),xlayer(mi(3),4),
     &     thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &     pslavsurf(3,*),pmastsurf(6,*),sti(6,mi(1),*),
     &     fn0(0:mi(2),*),dfn(0:mi(2),*),hglf(3,4),ahr,
     &     dxstiff(27,mi(1),ne,*),xdesi(3,*),smscale(1)
!     
      include "gauss.f"
!     
      iflag=3
      null=0
      qa(3)=-1.d0
      qa(4)=0.d0
      mscalmethod=0
      kscale=1
!     
      mt=mi(2)+1
!     
!     nodes in expansion direction of hex element 
!     
      data ifaceqexp /5,17,
     &     6,18,
     &     7,19,
     &     8,20,
     &     1,17,
     &     2,18,
     &     3,19,
     &     4,20,
     &     13,0,
     &     14,0,
     &     15,0,
     &     16,0,
     &     9,0,
     &     10,0,
     &     11,0,
     &     12,0,
     &     0,0,
     &     0,0,
     &     0,0,
     &     0,0/
!     
!     nodes in expansion direction of wedge element 
!     
      data ifacewexp /4,13,
     &     5,14,
     &     6,16,
     &     1,13,
     &     2,14,
     &     3,15,
     &     10,0,
     &     11,0,
     &     12,0,
     &     7,0,
     &     8,0,
     &     9,0,
     &     0,0,
     &     0,0,
     &     0,0/
!     
!     -------------------------------------------------------------
!     Initialisation of the loop for the variation of 
!     the internal forces
!     -------------------------------------------------------------
!     
!     
!     Loop over all elements in thread
      do ii=nea,neb
        if(idesvar.gt.0) then
          i=ialdesi(ii)
        else
          i=ii
        endif
!     
        lakonl=lakon(i)
!     
        if(ipkon(i).lt.0) cycle
!     
!     no 3D-fluid elements
!     
        if(lakonl(1:1).eq.'F') cycle
        if(lakonl(1:7).eq.'DCOUP3D') cycle
!     
        if(lakonl(7:8).ne.'LC') then
!     
          imat=ielmat(1,i)
          amat=matname(imat)
          if(norien.gt.0) then
            iorien=max(0,ielorien(1,i))
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
          mint2d=4
          nopes=8
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
          iflag=1
          indexe=ipkon(i)
          do kk=1,mint2d
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
              thickness=0.d0
              do j=1,nopes
                thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
              enddo
              tlayer(kk)=tlayer(kk)+thickness
              xlayer(k,kk)=thickness
            enddo
          enddo
          iflag=3
!     
          ilayer=0
          do k=1,4
            dlayer(k)=0.d0
          enddo
        elseif(lakonl(4:5).eq.'15') then
!     
!     composite materials
!     
          mint2d=3
          nopes=6
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
          iflag=1
          indexe=ipkon(i)
          do kk=1,mint2d
            xi=gauss3d10(1,kk)
            et=gauss3d10(2,kk)
            call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
              thickness=0.d0
              do j=1,nopes
                thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
              enddo
              tlayer(kk)=tlayer(kk)+thickness
              xlayer(k,kk)=thickness
            enddo
          enddo
          iflag=3
!     
          ilayer=0
          do k=1,3
            dlayer(k)=0.d0
          enddo
!     
        endif
!     
        indexe=ipkon(i)
c     Bernhardi start
        if(lakonl(1:5).eq.'C3D8I') then
          nope=11
        elseif(lakonl(4:5).eq.'20') then
c     Bernhardi end
          nope=20
        elseif(lakonl(4:4).eq.'8') then
          nope=8
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:4).eq.'6') then
          nope=6
        elseif((lakonl(1:1).eq.'E').and.(lakonl(7:7).ne.'F')) then
!     
!     spring elements, contact spring elements and
!     dashpot elements
!     
          if(lakonl(7:7).eq.'C') then
!     
!     contact spring elements
!     
            if(mortar.eq.1) then
!     
!     face-to-face penalty
!     
              nope=kon(ipkon(i))
            elseif(mortar.eq.0) then
!     
!     node-to-face penalty
!     
              nope=ichar(lakonl(8:8))-47
              konl(nope+1)=kon(indexe+nope+1)
            endif
          else
!     
!     genuine spring elements and dashpot elements
!     
            nope=ichar(lakonl(8:8))-47
          endif
        else
          cycle
        endif
!     
!     computation of the coordinates of the local nodes w/ variation
!     if the designnode belongs to the considered element
!     
        if(idesvar.gt.0) then 
!     
          if(icoordinate.eq.1) then
!     
!     the coordinates of the nodes are the design variables
!     
            do j=1,nope
              node=kon(indexe+j)
              if(node.eq.nodedesi(idesvar)) then
                iactive=j
                exit
              endif
            enddo
          endif
!     
        endif
!     
        if(lakonl(4:5).eq.'8R') then
          mint3d=1
        elseif(lakonl(4:7).eq.'20RB') then
          if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint3d=50
          else
            call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &           null,xi,et,ze,weight)
          endif
        elseif((lakonl(4:4).eq.'8').or.
     &         (lakonl(4:6).eq.'20R')) then
          if(lakonl(7:8).eq.'LC') then
            mint3d=8*nlayer
          else
            mint3d=8
          endif
        elseif(lakonl(4:4).eq.'2') then
          mint3d=27
        elseif(lakonl(4:5).eq.'10') then
          mint3d=4
        elseif(lakonl(4:4).eq.'4') then
          mint3d=1
        elseif(lakonl(4:5).eq.'15') then
          if(lakonl(7:8).eq.'LC') then
            mint3d=6*nlayer
          else
            mint3d=9
          endif
        elseif(lakonl(4:4).eq.'6') then
          mint3d=2
        elseif(lakonl(1:1).eq.'E') then
          mint3d=0
        endif

        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
            vl(k,j)=v(k,konl(j))
            voldl(k,j)=vold(k,konl(j))
          enddo
        enddo
!     
!     computation of the coordinates of the local nodes w/ variation
!     if the designnode belongs to the considered element
!     
        if((idesvar.gt.0).and.(icoordinate.eq.1)) then 
          do j=1,3
            xl(j,iactive)=xl(j,iactive)+xdesi(j,idesvar)
          enddo
          if(lakonl(1:5).eq.'C3D20') then
            if((lakonl(7:7).eq.'A').or.
     &           (lakonl(7:7).eq.'L').or.
     &           (lakonl(7:7).eq.'S').or.
     &           (lakonl(7:7).eq.'E')) then
              node1=ifaceqexp(1,iactive)
              node2=ifaceqexp(2,iactive)
              do j=1,3
                if(iactive.le.8) then
                  xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)
                  xl(j,node2)=xl(j,node2)+xdesi(j,idesvar)
                else
                  xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)
                endif      
              enddo      
            endif  
          elseif(lakonl(1:5).eq.'C3D15') then
            if((lakonl(7:7).eq.'A').or.
     &           (lakonl(7:7).eq.'L').or.
     &           (lakonl(7:7).eq.'S').or.
     &           (lakonl(7:7).eq.'E')) then
              node1=ifacewexp(1,iactive)
              node2=ifacewexp(2,iactive)
              do j=1,3
                if(iactive.le.6) then
                  xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)
                  xl(j,node2)=xl(j,node2)+xdesi(j,idesvar)
                else
                  xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)
                endif      
              enddo      
            endif  
          elseif(lakonl(1:4).eq.'C3D8') then
            if((lakonl(7:7).eq.'A').or.
     &           (lakonl(7:7).eq.'L').or.
     &           (lakonl(7:7).eq.'S').or.
     &           (lakonl(7:7).eq.'E')) then
              node1=ifaceqexp(1,iactive)
              do j=1,3
                xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)      
              enddo      
            endif  
          elseif(lakonl(1:4).eq.'C3D6') then
            if((lakonl(7:7).eq.'A').or.
     &           (lakonl(7:7).eq.'L').or.
     &           (lakonl(7:7).eq.'S').or.
     &           (lakonl(7:7).eq.'E')) then
              node1=ifaceqexp(1,iactive)
              do j=1,3
                xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)       
              enddo  
            endif     
          endif
        endif
!     
!     calculating the forces for the contact elements
!     
        if(mint3d.eq.0) then
!     
!     "normal" spring and dashpot elements
!     
          kode=nelcon(1,imat)
          if(lakonl(7:7).eq.'A') then
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
!     spring elements (including contact springs)
!     
          if(lakonl(2:2).eq.'S') then
            if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'1').or.
     &           (lakonl(7:7).eq.'2').or.((mortar.eq.0).and.
     &         ((nmethod.ne.1).or.(iperturb(1).ge.2).or.(iout.ne.-1))))
     &           then
!
!     nr is the element number in the assumption that all
!     slave nodes lead to contact, i.e. that in each slave node a
!     spring element was generated in gencontelem_n2f.f
!
              if(nener.eq.1) then
                if(lakonl(7:7).ne.'C') then
                  nr=i
                else
                  nr=ne0+konl(nope+1)
                endif
                venergy=enerini(2,1,nr)
              endif
              call springforc_n2f(xl,konl,vl,imat,elcon,nelcon,elas,
     &             fnl,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &             plicon,nplicon,npmat_,senergy,nener,
     &             stx(1,1,i),mi,springarea(1,konl(nope+1)),nmethod,
     &             ne0,nstate_,xstateini,xstate,reltime,
     &             ielas,venergy,ielorien,orab,norien,i,
     &             smscale,mscalmethod)
              if(nener.eq.1) then
                ener(1,1,nr)=senergy
                ener(2,1,nr)=venergy
              endif
            elseif((mortar.eq.1).and.
     &           ((nmethod.ne.1).or.(iperturb(1).ge.2).or.(iout.ne.-1)))
     &             then
              jfaces=kon(indexe+nope+2)
              igauss=kon(indexe+nope+1)
              if(nener.eq.1) venergy=enerini(2,1,ne0+igauss)
              call springforc_f2f(xl,vl,imat,elcon,nelcon,elas,
     &             fnl,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &             plicon,nplicon,npmat_,senergy,nener,
     &             stx(1,1,i),mi,springarea(1,igauss),nmethod,
     &             ne0,nstate_,xstateini,xstate,reltime,
     &             ielas,jfaces,igauss,pslavsurf,pmastsurf,
     &             clearini,venergy,kscale,konl,iout,i,
     &             smscale,mscalmethod)
            endif
!     
!     next lines are not executed in linstatic.c before the
!     setup of the stiffness matrix (i.e. nmethod=1 and
!     iperturb(1)<1 and iout=-1).
!     
            if((lakonl(7:7).eq.'A').or.
     &           ((nmethod.ne.1).or.(iperturb(1).ge.2).or.(iout.ne.-1)))
     &           then
              if(idesvar.eq.0) then
                do j=1,nope
                  do k=1,3
                    fn0(k,indexe+j)=fn0(k,indexe+j)+fnl(k,j)
                  enddo
                enddo
              else
                do j=1,nope
                  do k=1,3
                    dfn(k,konl(j))=dfn(k,konl(j))+fnl(k,j)
                  enddo
                enddo
              endif
            endif
!     
!     dashpot elements (including contact dashpots)
!     
          elseif((nmethod.eq.4).or.(nmethod.eq.5).or.
     &           ((abs(nmethod).eq.1).and.(iperturb(1).ge.2))) then
          endif
        elseif(ikin.eq.1) then
          do j=1,nope
            do k=1,3
              veoldl(k,j)=veold(k,konl(j))
            enddo
          enddo            
        endif
!     
        do jj=1,mint3d
          if(lakonl(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
          elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
              xi=gauss3d13(1,jj)
              et=gauss3d13(2,jj)
              ze=gauss3d13(3,jj)
              weight=weight3d13(jj)
            else
              call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &             jj,xi,et,ze,weight)
            endif
          elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d2(1,jj)
              et=gauss3d2(2,jj)
              ze=gauss3d2(3,jj)
              weight=weight3d2(jj)
            else
              kl=mod(jj,8)
              if(kl.eq.0) kl=8
!     
              xi=gauss3d2(1,kl)
              et=gauss3d2(2,kl)
              ze=gauss3d2(3,kl)
              weight=weight3d2(kl)
!     
              ki=mod(jj,4)
              if(ki.eq.0) ki=4
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do k=1,4
                    dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,i)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,i))
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
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
          elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
          elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
          elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d8(1,jj)
              et=gauss3d8(2,jj)
              ze=gauss3d8(3,jj)
              weight=weight3d8(jj)
            else
              kl=mod(jj,6)
              if(kl.eq.0) kl=6
!     
              xi=gauss3d10(1,kl)
              et=gauss3d10(2,kl)
              ze=gauss3d10(3,kl)
              weight=weight3d10(kl)
!     
              ki=mod(jj,3)
              if(ki.eq.0) ki=3
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do k=1,3
                    dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,i)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,i))
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
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
          endif
!     
c     Bernhardi start
          if(lakonl(1:5).eq.'C3D8R') then
            call shape8hr(xl,xsj,shp,gs,a)
          elseif(lakonl(1:5).eq.'C3D8I') then
            call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.20) then
c     Bernhardi end
            if(lakonl(7:7).eq.'A') then
              call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
            elseif((lakonl(7:7).eq.'E').or.
     &             (lakonl(7:7).eq.'S')) then
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
!     vkl(m2,m3) contains the derivative of the m2-
!     component of the displacement with respect to
!     direction m3
!     
          do m2=1,3
            do m3=1,3
              vkl(m2,m3)=0.d0
            enddo
          enddo
!     
          do m1=1,nope
            do m2=1,3
              do m3=1,3
                vkl(m2,m3)=vkl(m2,m3)+shp(m3,m1)*vl(m2,m1)
              enddo
c     write(*,*) 'vnoeie',i,konl(m1),(vkl(m2,k),k=1,3)
            enddo
          enddo
!     
!     for frequency analysis or buckling with preload the
!     strains are calculated with respect to the deformed
!     configuration
!     for a linear iteration within a nonlinear increment:
!     the tangent matrix is calculated at strain at the end
!     of the previous increment
!     
          if((iperturb(1).eq.1).or.(iperturb(1).eq.-1))then
            do m2=1,3
              do m3=1,3
                vokl(m2,m3)=0.d0
              enddo
            enddo
!     
            do m1=1,nope
              do m2=1,3
                do m3=1,3
                  vokl(m2,m3)=vokl(m2,m3)+
     &                 shp(m3,m1)*voldl(m2,m1)
                enddo
              enddo
            enddo
          endif
!     
          kode=nelcon(1,imat)
!     
!     calculating the strain
!     
!     attention! exy,exz and eyz are engineering strains!
!     
          exx=vkl(1,1)
          eyy=vkl(2,2)
          ezz=vkl(3,3)
          exy=vkl(1,2)+vkl(2,1)
          exz=vkl(1,3)+vkl(3,1)
          eyz=vkl(2,3)+vkl(3,2)
!     
          if(iperturb(2).eq.1) then
!     
!     Lagrangian strain
!     
            exx=exx+(vkl(1,1)**2+vkl(2,1)**2+vkl(3,1)**2)/2.d0
            eyy=eyy+(vkl(1,2)**2+vkl(2,2)**2+vkl(3,2)**2)/2.d0
            ezz=ezz+(vkl(1,3)**2+vkl(2,3)**2+vkl(3,3)**2)/2.d0
            exy=exy+vkl(1,1)*vkl(1,2)+vkl(2,1)*vkl(2,2)+
     &           vkl(3,1)*vkl(3,2)
            exz=exz+vkl(1,1)*vkl(1,3)+vkl(2,1)*vkl(2,3)+
     &           vkl(3,1)*vkl(3,3)
            eyz=eyz+vkl(1,2)*vkl(1,3)+vkl(2,2)*vkl(2,3)+
     &           vkl(3,2)*vkl(3,3)
!     
!     for frequency analysis or buckling with preload the
!     strains are calculated with respect to the deformed
!     configuration
!     
          elseif(iperturb(1).eq.1) then
            exx=exx+vokl(1,1)*vkl(1,1)+vokl(2,1)*vkl(2,1)+
     &           vokl(3,1)*vkl(3,1)
            eyy=eyy+vokl(1,2)*vkl(1,2)+vokl(2,2)*vkl(2,2)+
     &           vokl(3,2)*vkl(3,2)
            ezz=ezz+vokl(1,3)*vkl(1,3)+vokl(2,3)*vkl(2,3)+
     &           vokl(3,3)*vkl(3,3)
            exy=exy+vokl(1,1)*vkl(1,2)+vokl(1,2)*vkl(1,1)+
     &           vokl(2,1)*vkl(2,2)+vokl(2,2)*vkl(2,1)+
     &           vokl(3,1)*vkl(3,2)+vokl(3,2)*vkl(3,1)
            exz=exz+vokl(1,1)*vkl(1,3)+vokl(1,3)*vkl(1,1)+
     &           vokl(2,1)*vkl(2,3)+vokl(2,3)*vkl(2,1)+
     &           vokl(3,1)*vkl(3,3)+vokl(3,3)*vkl(3,1)
            eyz=eyz+vokl(1,2)*vkl(1,3)+vokl(1,3)*vkl(1,2)+
     &           vokl(2,2)*vkl(2,3)+vokl(2,3)*vkl(2,2)+
     &           vokl(3,2)*vkl(3,3)+vokl(3,3)*vkl(3,2)
          endif
!     
!     storing the local strains
!     
          if(iperturb(1).ne.-1) then
            eloc(1)=exx
            eloc(2)=eyy
            eloc(3)=ezz
            eloc(4)=exy/2.d0
            eloc(5)=exz/2.d0
            eloc(6)=eyz/2.d0
          else
!     
!     linear iteration within a nonlinear increment:
!     
            eloc(1)=vokl(1,1)+
     &           (vokl(1,1)**2+vokl(2,1)**2+vokl(3,1)**2)/2.d0
            eloc(2)=vokl(2,2)+
     &           (vokl(1,2)**2+vokl(2,2)**2+vokl(3,2)**2)/2.d0
            eloc(3)=vokl(3,3)+
     &           (vokl(1,3)**2+vokl(2,3)**2+vokl(3,3)**2)/2.d0
            eloc(4)=(vokl(1,2)+vokl(2,1)+vokl(1,1)*vokl(1,2)+
     &           vokl(2,1)*vokl(2,2)+vokl(3,1)*vokl(3,2))/2.d0
            eloc(5)=(vokl(1,3)+vokl(3,1)+vokl(1,1)*vokl(1,3)+
     &           vokl(2,1)*vokl(2,3)+vokl(3,1)*vokl(3,3))/2.d0
            eloc(6)=(vokl(2,3)+vokl(3,2)+vokl(1,2)*vokl(1,3)+
     &           vokl(2,2)*vokl(2,3)+vokl(3,2)*vokl(3,3))/2.d0
          endif
!     
!     calculating the deformation gradient (needed to
!     convert the element stiffness matrix from spatial
!     coordinates to material coordinates
!     deformation plasticity)
!     
          if((kode.eq.-50).or.(kode.le.-100)) then
!     
!     calculating the deformation gradient
!     
c     Bernhardi start
            xkl(1,1)=vkl(1,1)+1.0d0
            xkl(2,2)=vkl(2,2)+1.0d0
            xkl(3,3)=vkl(3,3)+1.0d0
c     Bernhardi end
            xkl(1,2)=vkl(1,2)
            xkl(1,3)=vkl(1,3)
            xkl(2,3)=vkl(2,3)
            xkl(2,1)=vkl(2,1)
            xkl(3,1)=vkl(3,1)
            xkl(3,2)=vkl(3,2)
!     
!     calculating the Jacobian
!     
            vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &           -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &           +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
c!     
c!     inversion of the deformation gradient (only for
c!     deformation plasticity)
c!     
c            if(kode.eq.-50) then
c!     
c              ckl(1,1)=(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))/vj
c              ckl(2,2)=(xkl(1,1)*xkl(3,3)-xkl(1,3)*xkl(3,1))/vj
c              ckl(3,3)=(xkl(1,1)*xkl(2,2)-xkl(1,2)*xkl(2,1))/vj
c              ckl(1,2)=(xkl(1,3)*xkl(3,2)-xkl(1,2)*xkl(3,3))/vj
c              ckl(1,3)=(xkl(1,2)*xkl(2,3)-xkl(2,2)*xkl(1,3))/vj
c              ckl(2,3)=(xkl(2,1)*xkl(1,3)-xkl(1,1)*xkl(2,3))/vj
c              ckl(2,1)=(xkl(3,1)*xkl(2,3)-xkl(2,1)*xkl(3,3))/vj
c              ckl(3,1)=(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))/vj
c              ckl(3,2)=(xkl(3,1)*xkl(1,2)-xkl(1,1)*xkl(3,2))/vj
c!     
c!     converting the Lagrangian strain into Eulerian
c!     strain (only for deformation plasticity)
c!     
c              cauchy=0
c              call str2mat(eloc,ckl,vj,cauchy)
c            endif
!     
          endif
!     
!     calculating fields for incremental plasticity
!     
          if(kode.le.-100) then
!     
!     calculating the deformation gradient at the
!     start of the increment
!     
!     calculating the displacement gradient at the
!     start of the increment
!     
            do m2=1,3
              do m3=1,3
                vikl(m2,m3)=0.d0
              enddo
            enddo
!     
            do m1=1,nope
              do m2=1,3
                do m3=1,3
                  vikl(m2,m3)=vikl(m2,m3)
     &                 +shp(m3,m1)*vini(m2,konl(m1))
                enddo
              enddo
            enddo
!     
!     calculating the deformation gradient of the old
!     fields
!     
            xikl(1,1)=vikl(1,1)+1.d0
            xikl(2,2)=vikl(2,2)+1.d0
            xikl(3,3)=vikl(3,3)+1.d0
            xikl(1,2)=vikl(1,2)
            xikl(1,3)=vikl(1,3)
            xikl(2,3)=vikl(2,3)
            xikl(2,1)=vikl(2,1)
            xikl(3,1)=vikl(3,1)
            xikl(3,2)=vikl(3,2)
!     
!     calculating the Jacobian
!     
            vij=xikl(1,1)*(xikl(2,2)*xikl(3,3)
     &           -xikl(2,3)*xikl(3,2))
     &           -xikl(1,2)*(xikl(2,1)*xikl(3,3)
     &           -xikl(2,3)*xikl(3,1))
     &           +xikl(1,3)*(xikl(2,1)*xikl(3,2)
     &           -xikl(2,2)*xikl(3,1))
!     
!     stresses at the start of the increment
!     
            do m1=1,6
              stre(m1)=sti(m1,jj,i)
            enddo
!     
          endif
!     
!     prestress values
!     
          if(iprestr.eq.1) then
            do kk=1,6
              beta(kk)=-prestr(kk,jj,i)
            enddo
          else
            do kk=1,6
              beta(kk)=0.d0
            enddo
          endif
!     
          if(ithermal(1).ge.1) then
!     
!     calculating the temperature difference in
!     the integration point
!     
            t0l=0.d0
            t1l=0.d0
            if(ithermal(1).eq.1) then
              if((lakonl(4:5).eq.'8 ').or.
     &             (lakonl(4:5).eq.'8I')) then
                do i1=1,8
                  t0l=t0l+t0(konl(i1))/8.d0
                  t1l=t1l+t1(konl(i1))/8.d0
                enddo
              elseif(lakonl(4:6).eq.'20 ') then
                nopered=20
                call lintemp(t0,konl,nopered,jj,t0l)
                call lintemp(t1,konl,nopered,jj,t1l)
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
              if((lakonl(4:5).eq.'8 ').or.
     &             (lakonl(4:5).eq.'8I')) then
                do i1=1,8
                  t0l=t0l+t0(konl(i1))/8.d0
                  t1l=t1l+vold(0,konl(i1))/8.d0
                enddo
              elseif(lakonl(4:6).eq.'20 ') then
                nopered=20
                call lintemp_th0(t0,konl,nopered,jj,t0l,mi)
                call lintemp_th1(vold,konl,nopered,jj,t1l,mi)
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
          endif
!     
!     calculating the coordinates of the integration point
!     for material orientation purposes (for cylindrical
!     coordinate systems)
!     
          if((iorien.gt.0).or.(kode.le.-100)) then
            do j=1,3
              pgauss(j)=0.d0
              do i1=1,nope
                pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
              enddo
            enddo
          endif
!     
!     material data; for linear elastic materials
!     this includes the calculation of the stiffness
!     matrix
!     
          istiff=0
!     
          call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &         nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &         elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &         istiff,elconloc,eth,kode,plicon,nplicon,
     &         plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jj,
     &         xstiff,ncmat_)
!     
!     determining the mechanical strain
!     
          if(ithermal(1).ne.0) then
            do m1=1,6
              emec(m1)=eloc(m1)-eth(m1)
c     emec0(m1)=emeini(m1,jj,i)
            enddo
          else
            do m1=1,6
              emec(m1)=eloc(m1)
c     emec0(m1)=emeini(m1,jj,i)
            enddo
          endif
          if(kode.le.-100) then
            do m1=1,6
              emec0(m1)=emeini(m1,jj,i)
            enddo
          endif
!     
!     subtracting the plastic initial strains
!     
          if(iprestr.eq.2) then
            do m1=1,6
              emec(m1)=emec(m1)-prestr(m1,jj,i)
            enddo
          endif
!     
!     calculating the local stiffness and stress
!     
          nlgeom_undo=0
          call mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &         icmd,beta,stre,xkl,ckl,vj,xikl,vij,
     &         plconloc,xstate,xstateini,ielas,
     &         amat,t1l,dtime,time,ttime,i,jj,nstate_,mi(1),
     &         iorien,pgauss,orab,eloc,mattyp,qa(3),istep,iinc,
     &         ipkon,nmethod,iperturb,qa(4),nlgeom_undo)
!     
          if(((nmethod.ne.4).or.(iperturb(1).ne.0)).and.
     &         (nmethod.ne.5)) then
            if(idesvar.eq.0) then
              do m1=1,21
                xstiff(m1,jj,i)=elas(m1)
              enddo
            elseif(icoordinate.ne.1) then
!     
!     for orientation design variables: store the modified 
!     stiffness matrix for use in mafillsmse
!     
              idir=idesvar-3*((idesvar-1)/3)
              do m1=1,21
                dxstiff(m1,jj,i,idir)=elas(m1)
              enddo
            endif
          endif
!     
          if((iperturb(1).eq.-1).and.(nlgeom_undo.eq.0)) then
!     
!     if the forced displacements were changed at
!     the start of a nonlinear step, the nodal
!     forces due do this displacements are 
!     calculated in a purely linear way, and
!     the first iteration is purely linear in order
!     to allow the displacements to redistribute
!     in a quasi-static way (only applies to
!     quasi-static analyses (*STATIC))
!     
            eloc(1)=exx-vokl(1,1)
            eloc(2)=eyy-vokl(2,2)
            eloc(3)=ezz-vokl(3,3)
            eloc(4)=exy-(vokl(1,2)+vokl(2,1))
            eloc(5)=exz-(vokl(1,3)+vokl(3,1))
            eloc(6)=eyz-(vokl(2,3)+vokl(3,2))
!     
            if(mattyp.eq.1) then
              e=elas(1)
              un=elas(2)
              um=e/(1.d0+un)
              al=un*um/(1.d0-2.d0*un)
              um=um/2.d0
              am1=al*(eloc(1)+eloc(2)+eloc(3))
              stre(1)=am1+2.d0*um*eloc(1)
              stre(2)=am1+2.d0*um*eloc(2)
              stre(3)=am1+2.d0*um*eloc(3)
              stre(4)=um*eloc(4)
              stre(5)=um*eloc(5)
              stre(6)=um*eloc(6)
            elseif(mattyp.eq.2) then
              stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)
     &             +eloc(3)*elas(4)
              stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)
     &             +eloc(3)*elas(5)
              stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)
     &             +eloc(3)*elas(6)
              stre(4)=eloc(4)*elas(7)
              stre(5)=eloc(5)*elas(8)
              stre(6)=eloc(6)*elas(9)
            elseif(mattyp.eq.3) then
              stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)+
     &             eloc(3)*elas(4)+eloc(4)*elas(7)+
     &             eloc(5)*elas(11)+eloc(6)*elas(16)
              stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)+
     &             eloc(3)*elas(5)+eloc(4)*elas(8)+
     &             eloc(5)*elas(12)+eloc(6)*elas(17)
              stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)+
     &             eloc(3)*elas(6)+eloc(4)*elas(9)+
     &             eloc(5)*elas(13)+eloc(6)*elas(18)
              stre(4)=eloc(1)*elas(7)+eloc(2)*elas(8)+
     &             eloc(3)*elas(9)+eloc(4)*elas(10)+
     &             eloc(5)*elas(14)+eloc(6)*elas(19)
              stre(5)=eloc(1)*elas(11)+eloc(2)*elas(12)+
     &             eloc(3)*elas(13)+eloc(4)*elas(14)+
     &             eloc(5)*elas(15)+eloc(6)*elas(20)
              stre(6)=eloc(1)*elas(16)+eloc(2)*elas(17)+
     &             eloc(3)*elas(18)+eloc(4)*elas(19)+
     &             eloc(5)*elas(20)+eloc(6)*elas(21)
            endif
          endif
!     
          skl(1,1)=stre(1)
          skl(2,2)=stre(2)
          skl(3,3)=stre(3)
          skl(2,1)=stre(4)
          skl(3,1)=stre(5)
          skl(3,2)=stre(6)
!     
          skl(1,2)=skl(2,1)
          skl(1,3)=skl(3,1)
          skl(2,3)=skl(3,2)
!     
!     calculation of the nodal forces
!     
          if(calcul_fn.eq.1)then
            if(idesvar.eq.0) then
!     
!     unperturbed state
!     
              do m1=1,nope
                do m2=1,3
!     
!     linear elastic part
!     
                  do m3=1,3
                    fn0(m2,indexe+m1)=fn0(m2,indexe+m1)+
     &                   xsj*skl(m2,m3)*shp(m3,m1)*weight
                  enddo
!     
!     nonlinear geometric part
!     
                  if((iperturb(2).eq.1).and.(nlgeom_undo.eq.0)) 
     &                 then
                    do m3=1,3
                      do m4=1,3
                        fn0(m2,indexe+m1)=fn0(m2,indexe+m1)+
     &                       xsj*skl(m4,m3)*weight*
     &                       (vkl(m2,m4)*shp(m3,m1)+
     &                       vkl(m2,m3)*shp(m4,m1))/2.d0
                      enddo
                    enddo
                  endif
!     
                enddo
              enddo
c     Bernhardi start
              if(lakonl(1:5).eq.'C3D8R') then
                ahr=elas(1)*a
!     
                do i1=1,3
                  do k=1,4    
                    hglf(i1,k)=0.0d0
                    do j=1,8
                      hglf(i1,k)=hglf(i1,k)+gs(j,k)*vl(i1,j)
                    enddo
                    hglf(i1,k)=hglf(i1,k)*ahr
                  enddo
                enddo
                do i1=1,3
                  do j=1,8
                    do k=1,4
                      fn0(i1,indexe+j)=fn0(i1,indexe+j)+
     &                     hglf(i1,k)*gs(j,k)
                    enddo
                  enddo
                enddo
              endif
c     Bernhardi end
            else
!     
!     perturbed state     
!     
              do m1=1,nope
                do m2=1,3
!     
!     linear elastic part
!     
                  do m3=1,3
                    dfn(m2,konl(m1))=dfn(m2,konl(m1))+
     &                   xsj*skl(m2,m3)*shp(m3,m1)*weight
                  enddo
!     
!     nonlinear geometric part
!     
                  if((iperturb(2).eq.1).and.(nlgeom_undo.eq.0)) 
     &                 then
                    do m3=1,3
                      do m4=1,3
                        dfn(m2,konl(m1))=dfn(m2,konl(m1))+
     &                       xsj*skl(m4,m3)*weight*
     &                       (vkl(m2,m4)*shp(m3,m1)+
     &                       vkl(m2,m3)*shp(m4,m1))/2.d0
                      enddo
                    enddo
                  endif
!     
                enddo
              enddo
c     Bernhardi start
              if(lakonl(1:5).eq.'C3D8R') then
!     
                ahr=elas(1)*a
!     
                do i1=1,3
                  do k=1,4    
                    hglf(i1,k)=0.0d0
                    do j=1,8
                      hglf(i1,k)=hglf(i1,k)+gs(j,k)*vl(i1,j)
                    enddo
                    hglf(i1,k)=hglf(i1,k)*ahr
                  enddo
                enddo
                do i1=1,3
                  do j=1,8
                    do k=1,4
                      dfn(i1,konl(j))=dfn(i1,konl(j))
     &                     +hglf(i1,k)*gs(j,k)
                    enddo
                  enddo
                enddo
              endif
c     Bernhardi end
            endif
          endif
!     
        enddo                   ! enddo loop over the integration points
!     
!     subtracting the unperturbed state of the element at stake
!     
        if(idesvar.gt.0) then
          do m1=1,nope
            do m2=1,3
              dfn(m2,konl(m1))=dfn(m2,konl(m1))
     &             -fn0(m2,indexe+m1)
            enddo
          enddo
        endif
!     
!     end of loop over all elements in thread
      enddo
!     
      return
      end
