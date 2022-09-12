!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine objective_mass_dx(co,kon,ipkon,lakon,nelcon,rhcon,
     &  ielmat,ielorien,norien,ntmat_,matname,mi,
     &  thicke,mortar,nea,neb,ielprop,prop,distmin,
     &  ndesi,nodedesi,nobject,g0,dgdx,iobject,xmass,
     &  istartdesi,ialdesi,xdesi,idesvar)
!
!     calculates the mass and its derivative w.r.t. the coordinates
!     of the mesh
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),nea,neb,mi(*),mint2d,nopes,nelcon(2,*),
     &  ielmat(mi(3),*),ielorien(mi(3),*),ntmat_,ipkon(*),iflag,null,
     &  mt,i,ii,j,k,jj,kk,indexe,nope,norien,ihyper,iactive,
     &  imat,mint3d,iorien,ilayer,nlayer,ki,kl,istartdesi(*),
     &  ielprop(*),mortar,idesvar,node,ndesi,nodedesi(*),
     &  nobject,iobject,ialdesi(*),ij,node1,node2,
     &  ifaceqexp(2,20),ifacewexp(2,15) 
!
      real*8 co(3,*),shp(4,26),xl(3,26),vl(0:mi(2),26),
     &  prop(*),rhcon(0:1,ntmat_,*),xs2(3,7),thickness,rho,xl2(3,8),
     &  xsj2(3),shp2(7,8),xi,et,ze,
     &  xsj,weight,gs(8,4),a,tlayer(4),dlayer(4),xlayer(mi(3),4),
     &  thicke(mi(3),*),distmin,xmassel,g0(*),xmass(*),xdesi(3,*),
     &  dgdx(ndesi,nobject)
!
!
!
      include "gauss.f"
!
      iflag=3
      null=0
!
      mt=mi(2)+1
!
!     nodes in expansion direction of hex element 
!
      data ifaceqexp /5,17,
     &                6,18,
     &                7,19,
     &                8,20,
     &                1,17,
     &                2,18,
     &                3,19,
     &                4,20,
     &                13,0,
     &                14,0,
     &                15,0,
     &                16,0,
     &                9,0,
     &                10,0,
     &                11,0,
     &                12,0,
     &                0,0,
     &                0,0,
     &                0,0,
     &                0,0/
!
!     nodes in expansion direction of wedge element 
!
      data ifacewexp /4,13,
     &                5,14,
     &                6,16,
     &                1,13,
     &                2,14,
     &                3,15,
     &                10,0,
     &                11,0,
     &                12,0,
     &                7,0,
     &                8,0,
     &                9,0,
     &                0,0,
     &                0,0,
     &                0,0/
!
!     -------------------------------------------------------------
!     Initialisation of the loop for the variation of 
!     the internal forces
!     -------------------------------------------------------------
!
!     
!     Loop over all elements in thread
         do ij=nea,neb
            if(idesvar.gt.0) then
               i=ialdesi(ij)
            else
               i=ij
            endif
     
!     initialisation of mass
!     
            xmassel=0.d0
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
            else
!     
!     composite materials
!
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
               if(lakonl(4:5).eq.'20') then
!
                  mint2d=4
                  nopes=8
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
                           thickness=thickness+thicke(k,indexe+j)*
     &                        shp2(4,j)
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
                  mint2d=3
                  nopes=6
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
                           thickness=thickness+thicke(k,indexe+j)*
     &                        shp2(4,j)
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
               endif
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
            elseif(lakonl(4:4).eq.'2') then
               nope=26
            elseif(lakonl(4:4).eq.'8') then
               nope=8
            elseif(lakonl(4:5).eq.'10') then
               nope=10
            elseif(lakonl(4:5).eq.'14') then
               nope=14
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
!     check whether a design variable belongs to the element    
!     
            if(idesvar.gt.0) then
               do ii=1,nope
                  node=kon(indexe+ii)
                  if(node.eq.nodedesi(idesvar)) then
                     iactive=ii
                     exit
                  endif
               enddo
            endif
!     
            if(lakonl(4:5).eq.'8R') then
               mint3d=1
            elseif(lakonl(4:7).eq.'20RB') then
               if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
                  mint3d=50
               else
                  call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &                 null,xi,et,ze,weight)
               endif
            elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'26R').or.
     &              (lakonl(4:6).eq.'20R')) then
               if(lakonl(7:8).eq.'LC') then
                  mint3d=8*nlayer
               else
                  mint3d=8
               endif
            elseif(lakonl(4:4).eq.'2') then
               mint3d=27
            elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
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
               enddo
            enddo
!     
!     computation of the objective function and the derivate 
!     if the designnode belongs to the considered element
!     if not, next element in loop will be taken 
!     
            if(idesvar.gt.0) then
               do j=1,3
                  xl(j,iactive)=xl(j,iactive)+xdesi(j,idesvar)
               enddo
               if(lakonl(1:5).eq.'C3D20') then
                  if((lakonl(7:7).eq.'A').or.
     &               (lakonl(7:7).eq.'L').or.
     &               (lakonl(7:7).eq.'S').or.
     &               (lakonl(7:7).eq.'E')) then
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
     &               (lakonl(7:7).eq.'L').or.
     &               (lakonl(7:7).eq.'S').or.
     &               (lakonl(7:7).eq.'E')) then
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
     &               (lakonl(7:7).eq.'L').or.
     &               (lakonl(7:7).eq.'S').or.
     &               (lakonl(7:7).eq.'E')) then
                     node1=ifaceqexp(1,iactive)
                     do j=1,3
                        xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)       
                     enddo  
                  endif     
               elseif(lakonl(1:4).eq.'C3D6') then
                  if((lakonl(7:7).eq.'A').or.
     &               (lakonl(7:7).eq.'L').or.
     &               (lakonl(7:7).eq.'S').or.
     &               (lakonl(7:7).eq.'E')) then
                     node1=ifaceqexp(1,iactive)
                     do j=1,3
                        xl(j,node1)=xl(j,node1)+xdesi(j,idesvar)       
                     enddo  
                  endif     
               endif
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
     &                    jj,xi,et,ze,weight)
                  endif
               elseif((lakonl(4:4).eq.'8').or.
     &                 (lakonl(4:6).eq.'20R').or.(lakonl(4:6).eq.'26R'))
     &                 then
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
                     ze=2.d0*(dlayer(ki)+(ze+1.d0)/
     &                    2.d0*xlayer(ilayer,ki))/tlayer(ki)-1.d0
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
               elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14'))
     &               then
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
                     ze=2.d0*(dlayer(ki)+(ze+1.d0)/
     &                    2.d0*xlayer(ilayer,ki))/tlayer(ki)-1.d0
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
     &                    (lakonl(7:7).eq.'S')) then
                     call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
                  else
                     call shape20h(xi,et,ze,xl,xsj,shp,iflag)
                  endif
c               elseif(nope.eq.26) then
c                  call shape26h(xi,et,ze,xl,xsj,shp,iflag,konl)
               elseif(nope.eq.8) then
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.10) then
                  call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
c               elseif(nope.eq.14) then
c                  call shape14tet(xi,et,ze,xl,xsj,shp,iflag,konl)
               elseif(nope.eq.4) then
                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.15) then
                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
               endif
!     
!     calculation of the objective function
!     
               rho=rhcon(1,1,imat)
               xmassel=xmassel+weight*xsj*rho
!
!               write(5,'(i5,3x,i5,3x,i5,3x,e14.8,3x,e14.7,3x,e14.8,
!     &                    3x,e14.8,3x,e14.8,3x,e20.14)') 
!     &                    idesvar,i,jj,weight,xsj,xi,et,ze,
!     &                    xmassel
!     
!     end of loop over all integration points
            enddo
!     
!     summation of the objective function over all elements
!     
            if(idesvar.eq.0) then
               xmass(i)=xmassel
               g0(iobject)=g0(iobject)+xmassel
            else
               dgdx(idesvar,iobject)=dgdx(idesvar,iobject)+
     &              (xmassel-xmass(i))/distmin
            endif
!     
!     end of loop over all elements in thread
!     
         enddo
!     
!     
      return
      end
      
