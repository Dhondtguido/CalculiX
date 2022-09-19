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
      subroutine e_c3d_em(co,konl,lakonl,s,sm,
     &     ff,nelem,nmethod,ielmat,
     &     ntmat_,t1,ithermal,vold,idist,
     &     matname,mi,mass,rhsi,
     &     ncmat_,elcon,nelcon,h0,iactive,
     &     alcon,nalcon,istartset,iendset,ialset)
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
      implicit none
!     
      logical mass,rhsi
!     
      character*8 lakonl
      character*80 matname(*),amat
!     
      integer konl(26),ifaceq(8,6),nelem,nmethod,iactive(3),mi(*),
     &     ithermal(*),idist,i,j,k,i1,m,one,ii,jj,id,ipointer,ig,kk,
     &     ielmat(mi(3),*),ntmat_,nope,nopes,imat,mint2d,mint3d,
     &     ifacet(6,4),nopev,ifacew(8,5),ipointeri,ipointerj,iflag,
     &     nelcon(2,*),ncmat_,nalcon(2,*),iel,ii1,ilength,istart,iset,
     &     isurf,jj1,istartset(*),iendset(*),ialset(*),three,nfaces,
     &     null
!     
      real*8 co(3,*),xl(3,26),shp(4,26),s(100,100),ff(100),xs2(3,7),
     &     t1(*),h0(3,*),xl2(3,9),xsj2(3),shp2(7,9),vold(0:mi(2),*),
     &     xi,et,ze,xsj,sm(100,100),t1l,weight,elcon(0:ncmat_,ntmat_,*),
     &     elconloc(ncmat_),sigma,h0l(3),h0l2(3,8),alcon(0:6,ntmat_,*),
     &     diag,um,uminv,alpha(6),d(3,3)
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
      data d /1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
!     
      null=0
      one=1
      three=3
      iflag=3
!     
      imat=ielmat(1,nelem)
      amat=matname(imat)
!     
      if(lakonl(4:5).eq.'20') then
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
      endif
!     
      if(lakonl(4:5).eq.'8R') then
        mint2d=1
        mint3d=1
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
        mint2d=4
        mint3d=8
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
        mint3d=9
      elseif(lakonl(4:4).eq.'6') then
        mint3d=2
      else
        mint3d=0
      endif
!     
!     computation of the coordinates of the local nodes
!     
      do i=1,nope
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!     
!     initialisation for distributed forces
!     
      if(rhsi) then
        if(idist.ne.0) then
          do i=1,5*nope
            ff(i)=0.d0
          enddo
        endif
      endif
!     
!     initialisation of sm
!     
      if(mass) then
        do i=1,5*nope
          do j=i,5*nope
            sm(i,j)=0.d0
          enddo
        enddo
      endif
!     
!     initialisation of s
!     
      do i=1,5*nope
        do j=i,5*nope
          s(i,j)=0.d0
        enddo
      enddo
!     
!     computation of the matrix: loop over the Gauss points
!     
      do kk=1,mint3d
        if(lakonl(4:5).eq.'8R') then
          xi=gauss3d1(1,kk)
          et=gauss3d1(2,kk)
          ze=gauss3d1(3,kk)
          weight=weight3d1(kk)
        elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &         then
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
        else
          xi=gauss3d7(1,kk)
          et=gauss3d7(2,kk)
          ze=gauss3d7(3,kk)
          weight=weight3d7(kk)
        endif
!     
!     calculation of the shape functions and their derivatives
!     in the gauss point
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
!     check the jacobian determinant
!     
        if(xsj.lt.1.d-20) then
          write(*,*) '*ERROR in e_c3d_th: nonpositive jacobian'
          write(*,*) '       determinant in element',nelem
          write(*,*)
          xsj=dabs(xsj)
          nmethod=0
        endif
!     
!     calculating the temperature and the magnetic intensity
!     in the integration point (one order lower than the
!     interpolation of the primary unknowns)
!     
        do j=1,3
          h0l(j)=0.d0
        enddo
        t1l=0.d0
!     
!     calculating the magnetic intensity
!     
        if(lakonl(4:5).eq.'8 ') then
          do i1=1,nope
            do j=1,3
              h0l(j)=h0l(j)+h0(j,konl(i1))/8.d0
            enddo
          enddo
        elseif(lakonl(4:6).eq.'20 ') then
          call linvec(h0,konl,nope,kk,h0l,one,three)
        elseif(lakonl(4:6).eq.'10T') then
          call linvec10(h0,konl,h0l,one,three,shp)
        else
          do i1=1,nope
            do j=1,3
              h0l(j)=h0l(j)+shp(4,i1)*h0(j,konl(i1))
            enddo
          enddo
        endif
!     
!     calculating the temperature
!     
        if(ithermal(1).eq.1) then
          if(lakonl(4:5).eq.'8 ') then
            do i1=1,nope
              t1l=t1l+t1(konl(i1))/8.d0
            enddo
          elseif(lakonl(4:6).eq.'20 ')then
            call linscal(t1,konl,nope,kk,t1l,one)
          elseif(lakonl(4:6).eq.'10T') then
            call linscal10(t1,konl,t1l,null,shp)
          else
            do i1=1,nope
              t1l=t1l+shp(4,i1)*t1(konl(i1))
            enddo
          endif
        elseif(ithermal(1).ge.2) then
          if(lakonl(4:5).eq.'8 ') then
            do i1=1,nope
              t1l=t1l+vold(0,konl(i1))/8.d0
            enddo
          elseif(lakonl(4:6).eq.'20 ')then
            call linscal(vold,konl,nope,kk,t1l,mi(2))
          elseif(lakonl(4:6).eq.'10T') then
            call linscal10(vold,konl,t1l,mi(2),shp)
          else
            do i1=1,nope
              t1l=t1l+shp(4,i1)*vold(0,konl(i1))
            enddo
          endif
        endif
!     
!     material data (electric conductivity and
!     magnetic permeability)
!     
        call materialdata_em(elcon,nelcon,alcon,nalcon,
     &       imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!     
        weight=weight*xsj
!     
        uminv=weight/elconloc(1)
        um=weight*elconloc(1)
        sigma=weight*alpha(1)
!     
        jj1=0
        do jj=1,nope
          ii1=0
          do ii=1,jj
!     
            if((int(elconloc(2)).eq.2).or.(int(elconloc(2)).eq.3))
     &           then
!     
!     K_AA matrix
!     
              diag=shp(1,ii)*shp(1,jj)+shp(2,ii)*shp(2,jj)+
     &             shp(3,ii)*shp(3,jj)
              do k=1,3
                do m=1,3
                  s(ii1+k,jj1+m)=s(ii1+k,jj1+m)+
     &                 (diag*d(k,m)-shp(m,ii)*shp(k,jj)
     &                 +shp(k,ii)*shp(m,jj))*uminv
                enddo
              enddo
!     
!     M_AA matrix
!     
              if(mass) then
                do k=1,3
                  sm(ii1+k,jj1+k)=sm(ii1+k,jj1+k)+
     &                 sigma*shp(4,ii)*shp(4,jj)
                enddo
              endif
              if((int(elconloc(2)).eq.2).and.mass) then
!     
!     M_AV and M_VA matrix
!     
                do k=1,3
                  sm(ii1+k,jj1+4)=sm(ii1+k,jj1+4)+
     &                 sigma*shp(4,ii)*shp(k,jj)
                  sm(ii1+4,jj1+k)=sm(ii1+4,jj1+k)+
     &                 sigma*shp(k,ii)*shp(4,jj)
                enddo
!     
!     M_VV matrix
!     
                sm(ii1+4,jj1+4)=sm(ii1+4,jj1+4)+
     &               sigma*(shp(1,ii)*shp(1,jj)+
     &               shp(2,ii)*shp(2,jj)+
     &               shp(3,ii)*shp(3,jj))
              endif
            else if(int(elconloc(2)).eq.1) then
!     
!     K_phiphi matrix
!     
              s(ii1+5,jj1+5)=s(ii1+5,jj1+5)-
     &             um*(shp(1,ii)*shp(1,jj)+
     &             shp(2,ii)*shp(2,jj)+
     &             shp(3,ii)*shp(3,jj))
            endif
!     
            ii1=ii1+5
          enddo
!     
!     F_phi matrix
!     
          if((rhsi).and.(int(elconloc(2)).eq.1)) then
            ff(jj1+5)=ff(jj1+5)-um*(shp(1,jj)*h0l(1)+
     &           shp(2,jj)*h0l(2)+
     &           shp(3,jj)*h0l(3))
          endif
!     
          jj1=jj1+5
        enddo
      enddo
!     
!     surface integrals
!     
!     determining the number of faces per element
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
      if(iactive(m).eq.0) return
      iset=iactive(m)
      istart=istartset(iset)
      ilength=iendset(iset)-istart+1
!     
      isurf=10*nelem+nfaces
      call nident(ialset(istart),isurf,ilength,id)
!     
      do
        if(id.eq.0) exit
        isurf=ialset(istart+id-1)
        iel=int(isurf/10.d0)
        if(iel.ne.nelem) exit
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
        if((nope.eq.20).or.(nope.eq.8)) then
          do i=1,nopes
            do j=1,3
              h0l2(j,i)=h0(j,konl(ifaceq(i,ig)))
              xl2(j,i)=co(j,konl(ifaceq(i,ig)))
            enddo
          enddo
        elseif((nope.eq.10).or.(nope.eq.4)) then
          do i=1,nopes
            do j=1,3
              h0l2(j,i)=h0(j,konl(ifacet(i,ig)))
              xl2(j,i)=co(j,konl(ifacet(i,ig)))
            enddo
          enddo
        else
          do i=1,nopes
            do j=1,3
              h0l2(j,i)=h0(j,konl(ifacew(i,ig)))
              xl2(j,i)=co(j,konl(ifacew(i,ig)))
            enddo
          enddo
        endif
!     
        do i=1,mint2d
!     
          if((lakonl(4:5).eq.'8R').or.
     &         ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
            xi=gauss2d1(1,i)
            et=gauss2d1(2,i)
            weight=weight2d1(i)
          elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R').or.
     &           ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
            xi=gauss2d2(1,i)
            et=gauss2d2(2,i)
            weight=weight2d2(i)
          elseif(lakonl(4:4).eq.'2') then
            xi=gauss2d3(1,i)
            et=gauss2d3(2,i)
            weight=weight2d3(i)
          elseif((lakonl(4:5).eq.'10').or.
     &           ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
            xi=gauss2d5(1,i)
            et=gauss2d5(2,i)
            weight=weight2d5(i)
          elseif((lakonl(4:4).eq.'4').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
            xi=gauss2d4(1,i)
            et=gauss2d4(2,i)
            weight=weight2d4(i)
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
          do k=1,3
            h0l(k)=0.d0
            do j=1,nopes
              h0l(k)=h0l(k)+h0l2(k,j)*shp2(4,j)
            enddo
          enddo
!     
          if((rhsi).and.(m.gt.1)) then
            do k=1,nopes
              if((nope.eq.20).or.(nope.eq.8)) then
                ipointer=5*(ifaceq(k,ig)-1)
              elseif((nope.eq.10).or.(nope.eq.4)) then
                ipointer=5*(ifacet(k,ig)-1)
              else
                ipointer=5*(ifacew(k,ig)-1)
              endif
!     
!     F_A vector
!     
              ff(ipointer+1)=ff(ipointer+1)+
     &             shp2(4,k)*(h0l(2)*xsj2(3)-h0l(3)*xsj2(2))*weight
              ff(ipointer+2)=ff(ipointer+2)+
     &             shp2(4,k)*(h0l(3)*xsj2(1)-h0l(1)*xsj2(3))*weight
              ff(ipointer+3)=ff(ipointer+3)+
     &             shp2(4,k)*(h0l(1)*xsj2(2)-h0l(2)*xsj2(1))*weight
!     
            enddo
          endif
!     
          do ii=1,nopes
            if((nope.eq.20).or.(nope.eq.8)) then
              ipointeri=5*(ifaceq(ii,ig)-1)
            elseif((nope.eq.10).or.(nope.eq.4)) then
              ipointeri=5*(ifacet(ii,ig)-1)
            else
              ipointeri=5*(ifacew(ii,ig)-1)
            endif
            do jj=1,nopes
              if((nope.eq.20).or.(nope.eq.8)) then
                ipointerj=5*(ifaceq(jj,ig)-1)
              elseif((nope.eq.10).or.(nope.eq.4)) then
                ipointerj=5*(ifacet(jj,ig)-1)
              else
                ipointerj=5*(ifacew(jj,ig)-1)
              endif
!     
              if(m.gt.1) then
!     
!     K_Aphi matrix
!     
                if((ipointeri+1).gt.ipointerj+5) cycle
                s(ipointeri+1,ipointerj+5)=
     &               s(ipointeri+1,ipointerj+5)
     &               +shp2(4,jj)*(shp2(3,ii)*xsj2(2)
     &               -shp2(2,ii)*xsj2(3))
     &               *weight
!     
                if((ipointeri+2).gt.ipointerj+5) cycle
                s(ipointeri+2,ipointerj+5)=
     &               s(ipointeri+2,ipointerj+5)
     &               +shp2(4,jj)*(shp2(1,ii)*xsj2(3)
     &               -shp2(3,ii)*xsj2(1))
     &               *weight
!     
                if((ipointeri+3).gt.ipointerj+5) cycle
                s(ipointeri+3,ipointerj+5)=
     &               s(ipointeri+3,ipointerj+5)
     &               +shp2(4,jj)*(shp2(2,ii)*xsj2(1)
     &               -shp2(1,ii)*xsj2(2))
     &               *weight
!     
!     K_phiA matrix
!     
                if(ipointeri+5.gt.(ipointerj+1)) cycle
                s(ipointeri+5,ipointerj+1)=
     &               s(ipointeri+5,ipointerj+1)
     &               +shp2(4,ii)*(shp2(3,jj)*xsj2(2)
     &               -shp2(2,jj)*xsj2(3))
     &               *weight
!     
                if(ipointeri+5.gt.(ipointerj+2)) cycle
                s(ipointeri+5,ipointerj+2)=
     &               s(ipointeri+5,ipointerj+2)
     &               +shp2(4,ii)*(shp2(1,jj)*xsj2(3)
     &               -shp2(3,jj)*xsj2(1))
     &               *weight
!     
                if(ipointeri+5.gt.(ipointerj+3)) cycle
                s(ipointeri+5,ipointerj+3)=
     &               s(ipointeri+5,ipointerj+3)
     &               +shp2(4,ii)*(shp2(2,jj)*xsj2(1)
     &               -shp2(1,jj)*xsj2(2))
     &               *weight
              endif
            enddo
          enddo
!     
        enddo
!     
        id=id-1
      enddo
!     
      return
      end
