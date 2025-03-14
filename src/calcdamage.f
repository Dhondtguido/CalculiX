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
      subroutine calcdamage(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold)
!     
!     calculates the damage due to plastic strain
!     
      implicit none
!     
      character*8 lakon(*),lakonl
!     
      integer ipkon(*),kon(*),mi(*),nope,indexe,i,j,k,ii,
     &     konl(20),mint3d,jj,iflag,ki,kl,ilayer,nlayer,kk,
     &     nopes,ielmat(mi(3),*),mint2d,null,ielprop(*),ne0,
     &     ndmcon(2,*),ndmat_,ntmat_,i1,nopered,ithermal(*),
     &     nstate_,id,imat,ndmconst
!     
      real*8 co(3,*),prop(*),xl(3,20),xi,et,ze,xsj,shp(4,20),weight,
     &     a,gs(8,4),dlayer(4),tlayer(4),thickness,skl(3,3),s(3,3),
     &     thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &     xl2(3,8),eps0RT,xlimit,d1,d2,d3,d4,d5,Tmelt,Ttrans,eps0p,shy,
     &     svm,triax,t1l,dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),
     &     dpeq,dpeqdt,dtime,ef,sti(6,mi(1),*),t1(*),vold(0:mi(2),*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     dmconloc(ndmat_)
!     
      include "gauss.f"
!     
      do i=1,ne0
!     
!     element must exist and be a volume element
!     
        if((ipkon(i).lt.0).or.(lakon(i)(1:1).ne.'C')) cycle
!     
        lakonl=lakon(i)
        indexe=ipkon(i)
!     
        if(lakonl(1:5).eq.'C3D8I') then
          nope=11
        elseif(lakonl(4:4).eq.'2') then
          nope=20
        elseif(lakonl(4:4).eq.'8') then
          nope=8
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:5).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
!     material
!     
        if(lakonl(7:8).ne.'LC') then
          imat=ielmat(1,i)
          if(ndmcon(2,imat).eq.0) cycle
!     
!     determining the model for this element
!     
          if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
            eps0RT=dmcon(2,1,imat)
            xlimit=dmcon(3,1,imat)
          elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
            d1=dmcon(2,1,imat)
            d2=dmcon(3,1,imat)
            d3=dmcon(4,1,imat)
            d4=dmcon(5,1,imat)
            d5=dmcon(6,1,imat)
            Tmelt=dmcon(7,1,imat)
            Ttrans=dmcon(8,1,imat)
            eps0p=dmcon(9,1,imat)
            xlimit=dmcon(10,1,imat)
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
!     the thickness of the composite layers is only needed for
!     models requiring the temperature at the integration points
!     (so far only for the Johnson-Cook model)
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            if(lakonl(4:4).eq.'2') then
              mint2d=4
              nopes=8
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
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,4
                dlayer(ii)=0.d0
              enddo
            elseif(lakonl(4:5).eq.'15') then
              mint2d=3
              nopes=6
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
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,3
                dlayer(ii)=0.d0
              enddo
            endif
          endif
!     
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
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
        elseif(lakonl(4:5).eq.'6') then
          mint3d=2
        else
          cycle
        endif
!     
        do jj=1,mint3d
!     
!     the calculation of the integration point coordinates and the
!     shape functions is only needed for models with temperature-dependent
!     data (Johnson-Cook).
!     
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
     &             kk,xi,et,ze,weight)
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
                  do ii=1,4
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
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
                  do ii=1,3
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
            endif
          else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
          endif
!     
!     determining the material model for this layer
!     (for composites only)
!     
          if(lakonl(7:8).ne.'LC') then
            if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
              eps0RT=dmcon(2,1,imat)
              xlimit=dmcon(3,1,imat)
            elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
              d1=dmcon(2,1,imat)
              d2=dmcon(3,1,imat)
              d3=dmcon(4,1,imat)
              d4=dmcon(5,1,imat)
              d5=dmcon(6,1,imat)
              Tmelt=dmcon(7,1,imat)
              Ttrans=dmcon(8,1,imat)
              eps0p=dmcon(9,1,imat)
              xlimit=dmcon(10,1,imat)
            endif
          endif
!     
!     shape functions need only be determined if the     
!     temperature is needed, i.e. for the Johnson-Cook model
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            if(lakonl(1:5).eq.'C3D8R') then
              call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
              call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
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
          endif
!     
!     stress at the integration point
!     
          skl(1,1)=sti(1,jj,i)
          skl(2,2)=sti(2,jj,i)
          skl(3,3)=sti(3,jj,i)
          skl(1,2)=sti(4,jj,i)
          skl(1,3)=sti(5,jj,i)
          skl(2,3)=sti(6,jj,i)
!     
!     hydrostatic stress
!     
          shy=(skl(1,1)+skl(2,2)+skl(3,3))/3.d0
!     
!     deviatoric stress tensor
!     
          s(1,1)=skl(1,1)-shy
          s(2,2)=skl(2,2)-shy
          s(3,3)=skl(3,3)-shy
          s(1,2)=skl(1,2)
          s(1,3)=skl(1,3)
          s(2,3)=skl(2,3)
!     
!     von Mises stress
!     
          svm=dsqrt(3.d0/2.d0*(
     &         s(1,1)*s(1,1)+s(2,2)*s(2,2)+s(3,3)*s(3,3)+
     &         2.d0*(s(1,2)*s(1,2)+s(1,3)*s(1,3)+s(2,3)*s(2,3))))
!     
!     triaxiality
!     
          triax=shy/svm
!     
!     calculate the temperature
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            if(ithermal(1).ge.1) then
              t1l=0.d0
              if(ithermal(1).eq.1) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+t1(konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp(t1,konl,nopered,jj,t1l)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(t1,konl,t1l,null,shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*t1(konl(i1))
                  enddo
                endif
              elseif(ithermal(1).ge.2) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+vold(0,konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp_th1(vold,konl,nopered,jj,t1l,mi)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(vold,konl,t1l,mi(2),shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                  enddo
                endif
              endif
            endif
          endif
!     
!     interpolating the material data
!     (for models with temperature dependent parameters;
!     no such model is implemented so far; a model
!     number of 3 or higher is assumed)          
!     
          if(int(dmcon(1,1,imat)).gt.2) then
!     
!           number of constants in this model    
!     
            ndmconst=ndmcon(1,imat)
            if(ithermal(1).eq.0) then
              do k=1,ndmconst
                dmconloc(k)=dmcon(k,1,imat)
              enddo
            else
              call ident2(dmcon(0,1,imat),t1l,ndmcon(2,imat),ndmat_+1,
     &             id)
              if(ndmcon(2,imat).eq.1) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.0) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.ndmcon(2,imat)) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)
                enddo
              else
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)+
     &                 (dmcon(k,id+1,imat)-dmcon(k,id,imat))*
     &                 (t1l-dmcon(0,id,imat))/
     &                 (dmcon(0,id+1,imat)-dmcon(0,id,imat))
                enddo
              endif
            endif
          endif          
!     
!     change in equivalent plastic strain
!     
          dpeq=xstate(1,j,i)-xstateini(1,j,i)
          dpeqdt=dpeq/dtime
!     
          if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
            ef=1.65d0*eps0RT*dexp(-3.d0*triax/2.d0)
            xlimit=dmcon(3,2,imat)
!     
          elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
            ef=(d1+d2*dexp(d3*triax))*(1.d0+d4*dlog(dpeqdt/eps0p))*
     &           (1.d0+d5*(t1l-Ttrans)/(Tmelt-Ttrans))
          endif
!     
!     damage
!     
          dam(j,i)=dam(j,i)+dpeq/ef
!     
          if(dam(j,i).gt.xlimit) then
            ipkon(i)=-ipkon(i)-2
            exit
          endif
        enddo
      enddo
!     
      return
      end
      
