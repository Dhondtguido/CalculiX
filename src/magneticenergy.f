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
      subroutine magneticenergy(ipkon,lakon,kon,co,elcon,nelcon,
     &     mi,ne,sti,ielmat,ncmat_,ntmat_,
     &     alcon,nalcon,ithermal,vold,t1,nmethod,prlab,nprint,
     &     ttime,time,prset,istartset,iendset,ialset,set,nset)
!     
!     determines the magnetic energy
!     
      implicit none
!     
      logical magnetic
!     
      character*6 prlab(*)
      character*8 lakon(*)
      character*81 set(*),elset,prset(*)
!     
      integer ipkon(*),nelem,kon(*),mi(*),nope,indexe,j,k,null,
     &     mint3d,jj,iflag,ne,ielmat(mi(3),*),konl(20),
     &     one,nelcon(2,*),nalcon(2,*),ithermal(*),i1,ncmat_,ntmat_,
     &     imat,nmethod,ii,nprint,ipos,kk,istartset(*),iendset(*),
     &     ialset(*),id,iset,nstart,nend,ninc,nset
!     
      real*8 co(3,*),xl(3,20),xi,et,ze,xsj,shp(4,20),weight,
     &     sti(6,mi(1),*),alpha(6),elcon(0:ncmat_,ntmat_,*),
     &     elconloc(ncmat_),t1l,alcon(0:6,ntmat_,*),vold(0:mi(2),*),
     &     t1(*),emagn,emagntot,time,ttime
!     
      include "gauss.f"
!     
      data iflag /2/
!     
      magnetic=.false.
      do ii=1,nprint
        if(prlab(ii)(1:4).eq.'ELME') then
          magnetic=.true.
          ipos=index(prset(ii),' ')
          elset=' '
          elset(1:ipos-1)=prset(ii)(1:ipos-1)
          write(5,*)
          write(5,110) elset(1:ipos-2),ttime+time
 110      format(' magnetic energy (elem, energy) for set '
     &         ,A,' and time ',e14.7)
          write(5,*)
          exit
        endif
      enddo
      if(.not.magnetic) return
!     
      emagntot=0.d0
      null=0
      one=1
!     
      call cident81(set,prset(ii),nset,id)
      iset=nset+1
      if(id.gt.0) then
        if(prset(ii).eq.set(id)) then
          iset=id
        endif
      endif
      do kk=istartset(iset),iendset(iset)
        if(ialset(kk).lt.0) cycle
        if(kk.eq.iendset(iset)) then
          nstart=ialset(kk)
          nend=ialset(kk)
          ninc=1
        elseif(ialset(kk+1).gt.0) then
          nstart=ialset(kk)
          nend=ialset(kk)
          ninc=1
        else
          nstart=ialset(kk-1)-ialset(kk+1)
          nend=ialset(kk)
          ninc=-ialset(kk+1)
        endif
!     
        do nelem=nstart,nend,ninc
          if(ipkon(nelem).lt.0) cycle
!     
c     if(int(elcon(2,1,ielmat(1,nelem))).ne.2) cycle
!     
          imat=ielmat(1,nelem)
          indexe=ipkon(nelem)
!     
          if(lakon(nelem)(1:5).eq.'C3D8I') then
            nope=11
          elseif(lakon(nelem)(4:4).eq.'2') then
            nope=20
          elseif(lakon(nelem)(4:4).eq.'8') then
            nope=8
          elseif(lakon(nelem)(4:5).eq.'10') then
            nope=10
          elseif(lakon(nelem)(4:4).eq.'4') then
            nope=4
          elseif(lakon(nelem)(4:5).eq.'15') then
            nope=15
          elseif(lakon(nelem)(4:5).eq.'6') then
            nope=6
          endif
!     
          do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
              xl(k,j)=co(k,konl(j))
            enddo
          enddo
!     
          emagn=0.d0
!     
          if(lakon(nelem)(4:5).eq.'8R') then
            mint3d=1
          elseif((lakon(nelem)(4:4).eq.'8').or.
     &           (lakon(nelem)(4:6).eq.'20R')) then
            mint3d=8
          elseif(lakon(nelem)(4:4).eq.'2') then
            mint3d=27
          elseif(lakon(nelem)(4:5).eq.'10') then
            mint3d=4
          elseif(lakon(nelem)(4:4).eq.'4') then
            mint3d=1
          elseif(lakon(nelem)(4:5).eq.'15') then
            mint3d=9
          elseif(lakon(nelem)(4:5).eq.'6') then
            mint3d=2
          endif
!     
!     loop over the integration points
!     
          do jj=1,mint3d
            if(lakon(nelem)(4:5).eq.'8R') then
              xi=gauss3d1(1,jj)
              et=gauss3d1(2,jj)
              ze=gauss3d1(3,jj)
              weight=weight3d1(jj)
            elseif((lakon(nelem)(4:4).eq.'8').or.
     &             (lakon(nelem)(4:6).eq.'20R'))
     &             then
              xi=gauss3d2(1,jj)
              et=gauss3d2(2,jj)
              ze=gauss3d2(3,jj)
              weight=weight3d2(jj)
            elseif(lakon(nelem)(4:4).eq.'2') then
              xi=gauss3d3(1,jj)
              et=gauss3d3(2,jj)
              ze=gauss3d3(3,jj)
              weight=weight3d3(jj)
            elseif(lakon(nelem)(4:5).eq.'10') then
              xi=gauss3d5(1,jj)
              et=gauss3d5(2,jj)
              ze=gauss3d5(3,jj)
              weight=weight3d5(jj)
            elseif(lakon(nelem)(4:4).eq.'4') then
              xi=gauss3d4(1,jj)
              et=gauss3d4(2,jj)
              ze=gauss3d4(3,jj)
              weight=weight3d4(jj)
            elseif(lakon(nelem)(4:5).eq.'15') then
              xi=gauss3d8(1,jj)
              et=gauss3d8(2,jj)
              ze=gauss3d8(3,jj)
              weight=weight3d8(jj)
            else
              xi=gauss3d7(1,jj)
              et=gauss3d7(2,jj)
              ze=gauss3d7(3,jj)
              weight=weight3d7(jj)
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
!     calculating the temperature
!     
            t1l=0.d0
            if(ithermal(1).eq.1) then
              if(lakon(nelem)(4:5).eq.'8 ') then
                do i1=1,nope
                  t1l=t1l+t1(konl(i1))/8.d0
                enddo
              elseif(lakon(nelem)(4:6).eq.'20 ')then
                call linscal(t1,konl,nope,jj,t1l,one)
              elseif(lakon(nelem)(4:6).eq.'10T') then
                call linscal10(t1,konl,t1l,null,shp)
              else
                do i1=1,nope
                  t1l=t1l+shp(4,i1)*t1(konl(i1))
                enddo
              endif
            elseif(ithermal(1).ge.2) then
              if(lakon(nelem)(4:5).eq.'8 ') then
                do i1=1,nope
                  t1l=t1l+vold(0,konl(i1))/8.d0
                enddo
              elseif(lakon(nelem)(4:6).eq.'20 ')then
                call linscal(vold,konl,nope,jj,t1l,mi(2))
              elseif(lakon(nelem)(4:6).eq.'10T') then
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
     &           imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!     
            if(nmethod.ne.2) then
!     
!     time dependent current: 
!     magneticenergy=B**2/(2*mu)
!     
              emagn=emagn+weight*xsj*
     &             (sti(4,jj,nelem)**2+
     &             sti(5,jj,nelem)**2+
     &             sti(6,jj,nelem)**2)/elconloc(1)
            else
!     
!     alternating current: 
!     magneticenergy=(B_real**2+B_imaginary**2)/(2*mu)
!     
              emagn=emagn+weight*xsj*
     &             (sti(4,jj,nelem)**2+
     &             sti(5,jj,nelem)**2+
     &             sti(6,jj,nelem)**2+
     &             sti(4,jj,nelem+ne)**2+
     &             sti(5,jj,nelem+ne)**2+
     &             sti(6,jj,nelem+ne)**2)/elconloc(1)
            endif
          enddo
        enddo
        emagn=emagn/2.d0
        write(5,'(i10,1p,1x,e13.6)') nelem,emagn
        emagntot=emagntot+emagn
      enddo
!     
      write(5,*)
      write(5,132) ttime+time
 132  format(' total magnetic energy for time ',
     &     e14.7)
      write(5,*)
      write(5,'(6x,1p,1x,e13.6)') emagntot
!     
      return
      end
