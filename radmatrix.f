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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine radmatrix(ntr,adrad,aurad,bcr,sideload,nelemload,
     &     xloadact,lakon,vold,ipkon,kon,co,nloadtr,tarea,tenv,physcon,
     &     erad,adview,auview,ithermal,iinc,iit,fenv,istep,dtime,ttime,
     &     time,iviewfile,xloadold,reltime,nmethod,mi,iemchange,nam,
     &     iamload,jqrad,irowrad,nzsrad)
!     
      implicit none
!
      character*8 lakonl,lakon(*)
      character*20 sideload(*)
!     
      integer ntr,nelemload(2,*),nope,nopes,mint2d,i,j,k,l,
     &     node,ifaceq(8,6),ifacet(6,4),iviewfile,mi(*),
     &     ifacew(8,5),nelem,ig,index,konl(20),iflag,
     &     ipkon(*),kon(*),nloadtr(*),i1,ithermal(*),iinc,iit,
     &     istep,jltyp,nfield,nmethod,iemchange,nam,
     &     iamload(2,*),irowrad(*),jqrad(*),nzsrad
!     
      real*8 adrad(*),aurad(*),bcr(ntr,1),xloadact(2,*),h(2),
     &     xl2(3,8),coords(3),dxsj2,temp,xi,et,weight,xsj2(3),
     &     vold(0:mi(2),*),co(3,*),shp2(7,8),xs2(3,7),
     &     areamean,tarea(*),tenv(*),erad(*),fenv(*),physcon(*),
     &     adview(*),auview(*),tl2(8),tvar(2),field,
     &     dtime,ttime,time,areaj,xloadold(2,*),reltime,
     &     fform
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
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
      external fform
!     
      include "gauss.f"
!     
      tvar(1)=time
      tvar(2)=ttime+time
!     
!     filling acr and bcr
!     
      do i1=1,ntr
         i=nloadtr(i1)
         nelem=nelemload(1,i)
         lakonl=lakon(nelem)
!     
!     calculate the mean temperature of the face
!     
         read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
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
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &           then
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
!     connectivity of the element
!     
         index=ipkon(nelem)
         if(index.lt.0) then
            write(*,*) '*ERROR in radmatrix: element ',nelem
            write(*,*) '       is not defined'
            call exit(201)
         endif
         do k=1,nope
            konl(k)=kon(index+k)
         enddo
!     
!     coordinates of the nodes belonging to the face
!     
         if((nope.eq.20).or.(nope.eq.8)) then
            do k=1,nopes
               tl2(k)=vold(0,konl(ifaceq(k,ig)))
!     
               do j=1,3
                  xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &                 vold(j,konl(ifaceq(k,ig)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do k=1,nopes
               tl2(k)=vold(0,konl(ifacet(k,ig)))
               do j=1,3
                  xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &                 vold(j,konl(ifacet(k,ig)))
               enddo
            enddo
         else
            do k=1,nopes
               tl2(k)=vold(0,konl(ifacew(k,ig)))
               do j=1,3
                  xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &                 vold(j,konl(ifacew(k,ig)))
               enddo
            enddo
         endif
!     
!     integration to obtain the center of gravity and the mean
!     temperature; radiation coefficient
!     
         areamean=0.d0
         tarea(i1)=0.d0
!     
         read(sideload(i)(2:2),'(i1)') jltyp
         jltyp=jltyp+10
         if(sideload(i)(5:6).ne.'NU') then
            erad(i1)=xloadact(1,i)
!
!           if an amplitude was defined for the emissivity it is
!           assumed that the emissivity changes with the step, so
!           acr has to be calculated anew in every iteration
!
            if(nam.gt.0) then
               if(iamload(1,i).ne.0) then
                  iemchange=1
               endif
            endif
         else
            erad(i1)=0.d0
         endif
!     
         do l=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
               xi=gauss2d1(1,l)
               et=gauss2d1(2,l)
               weight=weight2d1(l)
            elseif((lakonl(4:4).eq.'8').or.
     &              (lakonl(4:6).eq.'20R').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
               xi=gauss2d2(1,l)
               et=gauss2d2(2,l)
               weight=weight2d2(l)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss2d3(1,l)
               et=gauss2d3(2,l)
               weight=weight2d3(l)
            elseif((lakonl(4:5).eq.'10').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
               xi=gauss2d5(1,l)
               et=gauss2d5(2,l)
               weight=weight2d5(l)
            elseif((lakonl(4:4).eq.'4').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
               xi=gauss2d4(1,l)
               et=gauss2d4(2,l)
               weight=weight2d4(l)
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
     &           xsj2(3)*xsj2(3))
!     
            temp=0.d0
            do j=1,nopes
               temp=temp+tl2(j)*shp2(4,j)
            enddo
!     
            tarea(i1)=tarea(i1)+temp*dxsj2*weight
            areamean=areamean+dxsj2*weight
!     
            if(sideload(i)(5:6).eq.'NU') then
               areaj=dxsj2*weight
               do k=1,3
                  coords(k)=0.d0
               enddo
               do j=1,nopes
                  do k=1,3
                     coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
               enddo
               call radiate(h(1),tenv(i1),temp,istep,
     &              iinc,tvar,nelem,l,coords,jltyp,field,nfield,
     &              sideload(i),node,areaj,vold,mi,iemchange)
               if(nmethod.eq.1) h(1)=xloadold(1,i)+
     &              (h(1)-xloadold(1,i))*reltime
               erad(i1)=erad(i1)+h(1)
            endif
!     
         enddo
         tarea(i1)=tarea(i1)/areamean-physcon(1)
         if(sideload(i)(5:6).eq.'NU') then
            erad(i1)=erad(i1)/mint2d
         endif
!     
!        updating the right hand side
!     
         bcr(i1,1)=physcon(2)*(erad(i1)*tarea(i1)**4+
     &        (1.d0-erad(i1))*fenv(i1)*tenv(i1)**4)
c         write(*,*) 'radmatrix ',bcr(i1,1)
!     
      enddo
!
!     adrad and aurad is recalculated only if the viewfactors
!     or the emissivity changed
!
      if(((ithermal(1).eq.3).and.(iviewfile.ge.0)).or.
     &      ((iit.eq.-1).and.(iviewfile.ne.-2)).or.(iemchange.eq.1).or.
     &      ((iit.eq.0).and.(abs(nmethod).eq.1))) then
!
         do i=1,ntr
            erad(i)=1.d0-erad(i)
         enddo
!
!     diagonal entries
!     
         do i=1,ntr
            adrad(i)=1.d0-erad(i)*adview(i)
         enddo
!     
!     lower triangular entries
!     
         do i=1,nzsrad
            aurad(i)=-erad(irowrad(i))*auview(i)
         enddo
!     
!     upper triangular entries
!     
         do i=1,ntr
            do j=nzsrad+jqrad(i),nzsrad+jqrad(i+1)-1
               aurad(j)=-erad(i)*auview(j)
            enddo
         enddo
!
         do i=1,ntr
            erad(i)=1.d0-erad(i)
         enddo
      endif
!     
      return
      end
