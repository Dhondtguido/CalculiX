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
      subroutine resultsem(co,kon,ipkon,lakon,v,elcon,nelcon,ielmat,
     &  ntmat_,vini,dtime,matname,mi,ncmat_,nea,neb,sti,alcon,
     &  nalcon,h0,istartset,iendset,ialset,iactive,fn,eei,iout,nmethod)
!
!     calculates the heat flux and the material tangent at the integration
!     points and the internal concentrated flux at the nodes
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),mi(*),nelcon(2,*),ielmat(mi(3),*),
     &  ntmat_,ipkon(*),null,three,iflag,mt,i,j,k,m1,kk,i1,m3,indexe,
     &  nope,imat,mint3d,ncmat_,nea,neb,nalcon(2,*),mm,l,istart,iset,
     &  isurf,ilength,istartset(*),iendset(*),ialset(*),nopes,m,
     &  m2,ig,id,iactive(3),konl2(9),ifaceq(8,6),ifacet(6,4),iel,
     &  ifacew(8,5),mint2d,nfaces,one,iout,nmethod
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),xl(3,26),vl(0:mi(2),26),
     &  elcon(0:ncmat_,ntmat_,*),vkl(0:mi(2),3),vini(0:mi(2),*),c1,
     &  elconloc(ncmat_),xi,et,ze,xsj,t1l,dtime,weight,xsj2(3),
     &  vinikl(0:mi(2),3),alpha(6),vl2(0:mi(2),9),xl2(1:3,9),
     &  h0l(3),al(3),ainil(3),sti(6,mi(1),*),xs2(3,7),phi,
     &  um,alcon(0:6,ntmat_,*),h0(3,*),vinil(0:mi(2),26),
     &  fn(0:mi(2),*),eei(6,mi(1),*),shp2(7,9)
!
      include "gauss.f"
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      iflag=3
      null=0
      one=1
      three=3
!
      mt=mi(2)+1
!
!     calculation of temperatures and thermal flux
!
      do i=nea,neb
!
         lakonl(1:8)=lakon(i)(1:8)
!
         if(ipkon(i).lt.0) cycle
!
         imat=ielmat(1,i)
         amat=matname(imat)
!
         indexe=ipkon(i)
         if(lakonl(4:5).eq.'20') then
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
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
         if(lakonl(4:5).eq.'8R') then
            mint3d=1
            mint2d=1
         elseif((lakonl(4:4).eq.'8').or.
     &          (lakonl(4:6).eq.'20R')) then
            mint3d=8
            mint2d=4
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
            mint2d=9
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
            mint2d=3
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
            mint2d=1
         elseif(lakonl(4:5).eq.'15') then
            mint3d=9
         elseif(lakonl(4:4).eq.'6') then
            mint3d=2
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
            enddo
            do k=1,5
               vl(k,j)=v(k,konl(j))
            enddo
            vinil(4,j)=vini(4,konl(j))
         enddo
!
         do kk=1,mint3d
            if(lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,kk)
               et=gauss3d1(2,kk)
               ze=gauss3d1(3,kk)
               weight=weight3d1(kk)
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,kk)
               et=gauss3d2(2,kk)
               ze=gauss3d2(3,kk)
               weight=weight3d2(kk)
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
               xi=gauss3d8(1,kk)
               et=gauss3d8(2,kk)
               ze=gauss3d8(3,kk)
               weight=weight3d8(kk)
            elseif(lakonl(4:4).eq.'6') then
               xi=gauss3d7(1,kk)
               et=gauss3d7(2,kk)
               ze=gauss3d7(3,kk)
               weight=weight3d7(kk)
            endif
!
            if(nope.eq.20) then
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
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
            c1=xsj*weight
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!     
            do k=1,5
               do m3=1,3
                  vkl(k,m3)=0.d0
               enddo
!     
               do m1=1,nope
                  do m3=1,3
                     vkl(k,m3)=vkl(k,m3)+shp(m3,m1)*vl(k,m1)
                  enddo
               enddo
            enddo
!     
            do m3=1,3
               vinikl(4,m3)=0.d0
            enddo
            do m1=1,nope
               do m3=1,3
                  vinikl(4,m3)=vinikl(4,m3)+shp(m3,m1)*vinil(4,m1)
               enddo
            enddo
!
!              calculating the temperature difference in
!              the integration point
!
            t1l=0.d0
            do j=1,3
               h0l(j)=0.d0
               al(j)=0.d0
               ainil(j)=0.d0
            enddo
            if(lakonl(4:5).eq.'8 ') then
               do i1=1,8
                  do j=1,3
                     h0l(j)=h0l(j)+h0(j,konl(i1))/8.d0
                     al(j)=al(j)+v(j,konl(i1))/8.d0
                     ainil(j)=ainil(j)+vini(j,konl(i1))/8.d0
                  enddo
                  t1l=t1l+v(0,konl(i1))/8.d0
               enddo
            elseif(lakonl(4:6).eq.'20 ') then
               call linscal(v,konl,nope,kk,t1l,mi(2))
               call linvec(h0,konl,nope,kk,h0l,one,three)
               call linvec(v,konl,nope,kk,al,null,mi(2))
               call linvec(vini,konl,nope,kk,ainil,null,mi(2))
            elseif(lakonl(4:6).eq.'10T') then
               call linscal10(v,konl,t1l,mi(2),shp)
               call linvec10(h0,konl,h0l,one,three,shp)
               call linvec10(v,konl,al,null,mi(2),shp)
               call linvec10(vini,konl,ainil,null,mi(2),shp)
            else
               do i1=1,nope
                  t1l=t1l+shp(4,i1)*v(0,konl(i1))
                  do j=1,3
                     h0l(j)=h0l(j)+shp(4,i1)*h0(j,konl(i1))
                     al(j)=al(j)+shp(4,i1)*v(j,konl(i1))
                     ainil(j)=ainil(j)+shp(4,i1)*vini(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data (permeability)
!
            call materialdata_em(elcon,nelcon,alcon,nalcon,
     &           imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!
            um=elconloc(1)
!
            if(int(elconloc(2)).eq.1) then
!
!              magnetic field B in phi-domain
!
               do k=1,3
                  sti(k+3,kk,i)=um*(h0l(k)-vkl(5,k))
               enddo
!
!     calculating the electromagnetic force K_phiphi
!     
               do m1=1,nope
                  do m3=1,3
                     fn(5,konl(m1))=fn(5,konl(m1))-c1*
     &                  um*shp(m3,m1)*vkl(5,m3)
                  enddo
               enddo
            else
!
!              magnetic field B in A and A-V domain
!
               sti(4,kk,i)=vkl(3,2)-vkl(2,3)
               sti(5,kk,i)=vkl(1,3)-vkl(3,1)
               sti(6,kk,i)=vkl(2,1)-vkl(1,2)
!
!              electric field E in A-V domain
!
               if(int(elconloc(2)).eq.2) then
                 if((nmethod.eq.2).and.(iout.eq.1)) then
                   do k=1,3
                     sti(k,kk,i)=-al(k)-vkl(4,k)
                   enddo
                 else
                   do k=1,3
                     sti(k,kk,i)=(ainil(k)-al(k)+
     &                    vinikl(4,k)-vkl(4,k))/dtime
                   enddo
                 endif
               endif
!
!     calculating the electromagnetic force K_AA
!     
               do m1=1,nope
                  do m2=1,3
                     do m3=1,3
                        fn(m2,konl(m1))=fn(m2,konl(m1))+c1*
     &                       (shp(m3,m1)*(vkl(m2,m3)-vkl(m3,m2))+
     &                       shp(m2,m1)*vkl(m3,m3))/um
                     enddo
                  enddo
               enddo
!
            endif
!
         enddo
!     
!        surface integrals
!     
!        determining the number of faces per element
!
         if((lakonl(4:4).eq.'8').or.(lakonl(4:4).eq.'2')) then
            nfaces=6
         elseif((lakonl(4:4).eq.'6').or.(lakonl(4:5).eq.'15')) then
            nfaces=5
         elseif((lakonl(4:4).eq.'4').or.(lakonl(4:5).eq.'10')) then
            nfaces=4
         endif
!
         m=int(elconloc(2))
            if(iactive(m).eq.0) cycle
            iset=iactive(m)
            istart=istartset(iset)
            ilength=iendset(iset)-istart+1
!     
            isurf=10*i+nfaces
            call nident(ialset(istart),isurf,ilength,id)
!     
            do
               if(id.eq.0) exit
               isurf=ialset(istart+id-1)
               iel=int(isurf/10.d0)
               if(iel.ne.i) exit
               ig=isurf-10*iel
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
!              face connectivity is stored in konl2
!  
               if((nope.eq.20).or.(nope.eq.8)) then
                  do j=1,nopes
                     konl2(j)=konl(ifaceq(j,ig))
                  enddo
               elseif((nope.eq.10).or.(nope.eq.4)) then
                  do j=1,nopes
                     konl2(j)=konl(ifacet(j,ig))
                  enddo
               else
                  do j=1,nopes
                     konl2(j)=konl(ifacew(j,ig))
                  enddo
               endif
!
!              face coordinates and solution fields
!
               do j=1,nopes
                  do k=1,3
                     xl2(k,j)=co(k,konl2(j))
                     vl2(k,j)=v(k,konl2(j))
                  enddo
                  vl2(5,j)=v(5,konl2(j))
               enddo
!     
               do kk=1,mint2d
!     
                  if((lakonl(4:5).eq.'8R').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                     xi=gauss2d1(1,kk)
                     et=gauss2d1(2,kk)
                     weight=weight2d1(kk)
                  elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')
     &                .or.((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                     xi=gauss2d2(1,kk)
                     et=gauss2d2(2,kk)
                     weight=weight2d2(kk)
                  elseif(lakonl(4:4).eq.'2') then
                     xi=gauss2d3(1,kk)
                     et=gauss2d3(2,kk)
                     weight=weight2d3(kk)
                  elseif((lakonl(4:5).eq.'10').or.
     &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                     xi=gauss2d5(1,kk)
                     et=gauss2d5(2,kk)
                     weight=weight2d5(kk)
                  elseif((lakonl(4:4).eq.'4').or.
     &                    ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                     xi=gauss2d4(1,kk)
                     et=gauss2d4(2,kk)
                     weight=weight2d4(kk)
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
!                 determining the electromagnetic force
!
                  if(m.gt.1) then
!
!                    determining the values of phi
!
                     phi=0.d0
                     do m1=1,nopes
                        phi=phi+shp2(4,m1)*vl2(5,m1)
                     enddo
!
                     do m1=1,nopes
                        do k=1,3
                           l=k+1
                           if(l.gt.3) l=1
                           mm=l+1
                           if(mm.gt.3) mm=1
                           fn(k,konl2(m1))=fn(k,konl2(m1))-
     &                        (shp2(l,m1)*xsj2(mm)-shp2(mm,m1)*xsj2(l))
     &                          *phi*weight
                        enddo
                     enddo
                  elseif(m.eq.1) then
!     
!     determining the derivative of A w.r.t. x, y and z
!     
                     do k=1,3
                        do m3=1,3
                           vkl(k,m3)=0.d0
                        enddo
                        do m1=1,nopes
                           do m3=1,3
                              vkl(k,m3)=vkl(k,m3)+shp2(m3,m1)*vl2(k,m1)
                           enddo
                        enddo
                     enddo
!
                     do m1=1,nopes
                        do k=1,3
                           l=k+1
                           if(l.gt.3) l=1
                           mm=l+1
                           if(mm.gt.3) mm=1
                           fn(5,konl2(m1))=fn(5,konl2(m1))+
     &                        shp2(4,m1)*(vkl(mm,l)-vkl(l,mm))*xsj2(k)
     &                        *weight
                        enddo
                     enddo
                  endif
!                        
               enddo
               id=id-1
            enddo
      enddo
!
      return
      end
