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
      subroutine initialconditionss(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,t0,t1,prestr,iprestr,ithermal,veold,inoelfree,
     &     nk_,
     &     mi,istep,istat,n,iline,ipol,inl,ipoinp,inp,lakon,kon,co,ne,
     &     ipkon,vold,ipoinpc,xstate,nstate_,nk,t0g,t1g,iaxial,ielprop,
     &     prop,ier,nuel_)
!     
!     reading the input deck: *INITIAL CONDITIONS
!     
      implicit none
!     
      logical user
!     
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 rebarn
      character*81 set(*),noset
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),nset,iprestr,
     &     ithermal(*),m,id,nuel_,
     &     istep,istat,n,i,j,k,l,ii,key,idir,ipos,inoelfree,nk_,mi(*),
     &     iline,ipol,inl,ipoinp(2,*),inp(3,*),ij,jj,ntens,ncrds,layer,
     &     kspt,lrebar,iflag,i1,mint3d,nope,kon(*),konl(20),indexe,
     &     ipkon(*),ne,ipoinpc(0:*),nstate_,nk,jmax,ntot,numberoflines,
     &     iaxial,null,ielprop(*),ier
!     
      real*8 t0(*),t1(*),beta(8),prestr(6,mi(1),*),veold(0:mi(2),*),
     &     temperature,velocity,tempgrad1,tempgrad2,pgauss(3),
     &     shp(4,20),xsj,xl(3,20),xi,et,ze,weight,co(3,*),pressure,
     &     vold(0:mi(2),*),xstate(nstate_,mi(1),*),dispvelo,totpres,
     &     xmassflow,t0g(2,*),t1g(2,*),prop(*),turbini(0:mi(2))
!     
      include "gauss.f"
!     
      null=0
!     
      if(istep.gt.0) then
        write(*,*) 
     &       '*ERROR reading *INITIAL CONDITIONS: *INITIAL CONDITIONS'
        write(*,*) '  should be placed before all step definitions'
        ier=1
        return
      endif
!     
      do ij=2,n
        if(textpart(ij)(1:16).eq.'TYPE=TEMPERATURE') then
!     
          ithermal(1)=1
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &           temperature
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            temperature=1.d-6*int(1.d6*temperature+0.5d0)
!     
            if((inoelfree.ne.0).or.(nuel_.gt.0)) then
              tempgrad1=0.d0
              tempgrad2=0.d0
              if(n.gt.2) then
                read(textpart(3)(1:20),'(f20.0)',iostat=istat) 
     &               tempgrad1
                if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*INITIAL CONDITIONS%",ier)
                  return
                endif
              endif
              if(n.gt.3) then
                read(textpart(4)(1:20),'(f20.0)',iostat=istat) 
     &               tempgrad2
                if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*INITIAL CONDITIONS%",ier)
                  return
                endif
              endif
            endif
!     
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              t0(l)=temperature
              t1(l)=temperature
              vold(0,l)=temperature
              if((inoelfree.ne.0).or.(nuel_.gt.0)) then
                t0g(1,l)=tempgrad1
                t0g(2,l)=tempgrad2
                t1g(1,l)=tempgrad1
                t1g(2,l)=tempgrad2
              endif
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  t0(ialset(j))=temperature
                  t1(ialset(j))=temperature
                  vold(0,ialset(j))=temperature
                  if((inoelfree.ne.0).or.(nuel_.gt.0)) then
                    t0g(1,ialset(j))=tempgrad1
                    t0g(2,ialset(j))=tempgrad2
                    t1g(1,ialset(j))=tempgrad1
                    t1g(2,ialset(j))=tempgrad2
                  endif
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    t0(k)=temperature
                    t1(k)=temperature
                    vold(0,k)=temperature
                    if((inoelfree.ne.0).or.(nuel_.gt.0)) then
                      t0g(1,k)=tempgrad1
                      t0g(2,k)=tempgrad2
                      t1g(1,k)=tempgrad1
                      t1g(2,k)=tempgrad2
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:11).eq.'TYPE=STRESS') then
!     
          iprestr=1
          do jj=1,n
            if(textpart(jj)(1:4).eq.'USER') then
!     
!     residual stresses are defined by user subroutine
!     sigini
!     
              iflag=1
              ntens=6
              ncrds=3
              lrebar=0
              do i=1,ne
                indexe=ipkon(i)
                if(lakon(i)(4:4).eq.'2') then
                  nope=20
                elseif(lakon(i)(4:4).eq.'8') then
                  nope=8
                elseif(lakon(i)(4:5).eq.'10') then
                  nope=10
                elseif(lakon(i)(4:4).eq.'4') then
                  nope=4
                elseif(lakon(i)(4:5).eq.'15') then
                  nope=15
                elseif(lakon(i)(4:4).eq.'6') then
                  nope=6
                else
                  cycle
                endif
!     
                if(lakon(i)(4:5).eq.'8R') then
                  mint3d=1
                elseif(lakon(i)(4:7).eq.'20RB') then
                  if((lakon(i)(8:8).eq.'R').or.
     &                 (lakon(i)(8:8).eq.'C')) then
                    mint3d=50
                  else
                    call beamintscheme(lakon(i),mint3d,
     &                   ielprop(i),prop,
     &                   null,xi,et,ze,weight)
                  endif
                elseif((lakon(i)(4:4).eq.'8').or.
     &                 (lakon(i)(4:6).eq.'20R')) then
                  mint3d=8
                elseif(lakon(i)(4:4).eq.'2') then
                  mint3d=27
                elseif(lakon(i)(4:5).eq.'10') then
                  mint3d=4
                elseif(lakon(i)(4:4).eq.'4') then
                  mint3d=1
                elseif(lakon(i)(4:5).eq.'15') then
                  mint3d=9
                elseif(lakon(i)(4:4).eq.'6') then
                  mint3d=2
                endif
!     
                do j=1,nope
                  konl(j)=kon(indexe+j)
                  do k=1,3
                    xl(k,j)=co(k,konl(j))
                  enddo
                enddo
!     
                do j=1,mint3d
                  if(lakon(i)(4:5).eq.'8R') then
                    xi=gauss3d1(1,j)
                    et=gauss3d1(2,j)
                    ze=gauss3d1(3,j)
                    weight=weight3d1(j)
                  elseif(lakon(i)(4:7).eq.'20RB') then
                    if((lakon(i)(8:8).eq.'R').or.
     &                   (lakon(i)(8:8).eq.'C')) then
                      xi=gauss3d13(1,j)
                      et=gauss3d13(2,j)
                      ze=gauss3d13(3,j)
                      weight=weight3d13(j)
                    else
                      call beamintscheme(lakon(i),mint3d,
     &                     ielprop(i),prop,
     &                     j,xi,et,ze,weight)
                    endif
                  elseif((lakon(i)(4:4).eq.'8').or.
     &                   (lakon(i)(4:6).eq.'20R'))
     &                   then
                    xi=gauss3d2(1,j)
                    et=gauss3d2(2,j)
                    ze=gauss3d2(3,j)
                    weight=weight3d2(j)
                  elseif(lakon(i)(4:4).eq.'2') then
                    xi=gauss3d3(1,j)
                    et=gauss3d3(2,j)
                    ze=gauss3d3(3,j)
                    weight=weight3d3(j)
                  elseif(lakon(i)(4:5).eq.'10') then
                    xi=gauss3d5(1,j)
                    et=gauss3d5(2,j)
                    ze=gauss3d5(3,j)
                    weight=weight3d5(j)
                  elseif(lakon(i)(4:4).eq.'4') then
                    xi=gauss3d4(1,j)
                    et=gauss3d4(2,j)
                    ze=gauss3d4(3,j)
                    weight=weight3d4(j)
                  elseif(lakon(i)(4:5).eq.'15') then
                    xi=gauss3d8(1,j)
                    et=gauss3d8(2,j)
                    ze=gauss3d8(3,j)
                    weight=weight3d8(j)
                  elseif(lakon(i)(4:4).eq.'6') then
                    xi=gauss3d7(1,j)
                    et=gauss3d7(2,j)
                    ze=gauss3d7(3,j)
                    weight=weight3d7(j)
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
                  do k=1,3
                    pgauss(k)=0.d0
                    do i1=1,nope
                      pgauss(k)=pgauss(k)+
     &                     shp(4,i1)*co(k,konl(i1))
                    enddo
                  enddo
!     
                  call sigini(prestr(1,j,i),pgauss,ntens,ncrds,
     &                 i,j,layer,kspt,lrebar,rebarn)
!     
                enddo
              enddo
              call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &             inl,ipoinp,inp,ipoinpc)
              return
            endif
          enddo
!     
!     residual stresses are written explicitly in the input deck
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            do j=1,6
              read(textpart(j+2)(1:20),'(f20.0)',iostat=istat) 
     &             beta(j)
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
            enddo
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.ne.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            if(l.gt.ne) then
              write(*,*) 
     &             '*WARNING reading *INITIAL CONDITIONS: element ',l
              write(*,*)'          exceeds the largest defined ',
     &             'element number'
              cycle
            endif
            read(textpart(2)(1:10),'(i10)',iostat=istat) k
            if(istat.eq.0) then
              do j=1,6
                prestr(j,k,l)=beta(j)
              enddo
            else
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
          enddo
          return
        elseif(textpart(ij)(1:18).eq.'TYPE=PLASTICSTRAIN') then
!     
          iprestr=2
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            do j=1,6
              read(textpart(j+2)(1:20),'(f20.0)',iostat=istat) 
     &             beta(j)
              if(istat.gt.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
            enddo
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.ne.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            if(l.gt.ne) then
              write(*,*) 
     &             '*WARNING reading *INITIAL CONDITIONS: element ',l
              write(*,*)'          exceeds the largest defined ',
     &             'element number'
              cycle
            endif
            read(textpart(2)(1:10),'(i10)',iostat=istat) k
            if(istat.eq.0) then
              do j=1,6
                prestr(j,k,l)=beta(j)
              enddo
            else
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
          enddo
          return
        elseif((textpart(ij)(1:17).eq.'TYPE=DISPLACEMENT').or.
     &         (textpart(ij)(1:18).eq.'TYPE=FLUIDVELOCITY')) then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:10),'(i10)',iostat=istat) idir
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) dispvelo
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              vold(idir,l)=dispvelo
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  vold(idir,ialset(j))=dispvelo
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    vold(idir,k)=dispvelo
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:13).eq.'TYPE=VELOCITY') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:10),'(i10)',iostat=istat) idir
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) velocity
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              veold(idir,l)=velocity
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  veold(idir,ialset(j))=velocity
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    veold(idir,k)=velocity
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:13).eq.'TYPE=PRESSURE') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) pressure
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              vold(4,l)=pressure
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  vold(4,ialset(j))=pressure
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    vold(4,k)=pressure
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:15).eq.'TYPE=TURBULENCE') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            do m=5,mi(2)
              read(textpart(m-3)(1:20),'(f20.0)',iostat=istat)
     &             turbini(m)
            enddo
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              do m=5,mi(2)
                vold(m,l)=turbini(m)
              enddo
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  do m=5,mi(2)
                    vold(m,ialset(j))=turbini(m)
                  enddo
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    do m=5,mi(2)
                      vold(m,k)=turbini(m)
                    enddo
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:18).eq.'TYPE=TOTALPRESSURE') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) totpres
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              vold(2,l)=totpres
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  vold(2,ialset(j))=totpres
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    vold(2,k)=totpres
                  enddo
                endif
              enddo
            endif
          enddo
          return
        elseif(textpart(ij)(1:13).eq.'TYPE=MASSFLOW') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) xmassflow
            if(iaxial.eq.180) xmassflow=xmassflow/iaxial
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*INITIAL CONDITIONS%",ier)
              return
            endif
            read(textpart(1)(1:10),'(i10)',iostat=istat) l
            if(istat.eq.0) then
              if(l.gt.nk) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: node ',l
                write(*,*)'          exceeds the largest defined ',
     &               'node number'
                cycle
              endif
              vold(1,l)=xmassflow
            else
              read(textpart(1)(1:80),'(a80)',iostat=istat) noset
              noset(81:81)=' '
              ipos=index(noset,' ')
              noset(ipos:ipos)='N'
c              do ii=1,nset
c                if(set(ii).eq.noset) exit
c              enddo
              call cident81(set,noset,nset,id)
              ii=nset+1
              if(id.gt.0) then
                if(noset.eq.set(id)) then
                  ii=id
                endif
              endif
              if(ii.gt.nset) then
                noset(ipos:ipos)=' '
                write(*,*) 
     &               '*ERROR reading *INITIAL CONDITIONS: node set '
     &               ,noset
                write(*,*)'  has not yet been defined. '
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              do j=istartset(ii),iendset(ii)
                if(ialset(j).gt.0) then
                  vold(1,ialset(j))=xmassflow
                else
                  k=ialset(j-2)
                  do
                    k=k-ialset(j)
                    if(k.ge.ialset(j-1)) exit
                    vold(1,k)=xmassflow
                  enddo
                endif
              enddo
            endif
          enddo
          return
!     
        elseif(textpart(ij)(1:13).eq.'TYPE=SOLUTION') then
          user=.false.
          do j=2,n
            if(textpart(j)(1:4).eq.'USER') user=.true.
          enddo
c     if(.not.user) then
c     write(*,*) 
c     &            '*ERROR reading *INITIAL CONDITIONS: TYPE=SOLUTION'
c     write(*,*) '       can only be used in combination with'
c     write(*,*) '       USER'
c     ier=1
c     return
c     endif
          if(user) then
!     
!     internal state variables are read in file sdvini.f
!     
            iflag=1
            ncrds=3
            do i=1,ne
              indexe=ipkon(i)
              if(lakon(i)(4:4).eq.'2') then
                nope=20
              elseif(lakon(i)(4:4).eq.'8') then
                nope=8
              elseif(lakon(i)(4:5).eq.'10') then
                nope=10
              elseif(lakon(i)(4:4).eq.'4') then
                nope=4
              elseif(lakon(i)(4:5).eq.'15') then
                nope=15
              elseif(lakon(i)(4:4).eq.'6') then
                nope=6
              else
                cycle
              endif
!     
              if(lakon(i)(4:5).eq.'8R') then
                mint3d=1
              elseif((lakon(i)(4:4).eq.'8').or.
     &               (lakon(i)(4:6).eq.'20R')) then
                mint3d=8
              elseif(lakon(i)(4:4).eq.'2') then
                mint3d=27
              elseif(lakon(i)(4:5).eq.'10') then
                mint3d=4
              elseif(lakon(i)(4:4).eq.'4') then
                mint3d=1
              elseif(lakon(i)(4:5).eq.'15') then
                mint3d=9
              elseif(lakon(i)(4:4).eq.'6') then
                mint3d=2
              endif
!     
              do j=1,nope
                konl(j)=kon(indexe+j)
                do k=1,3
                  xl(k,j)=co(k,konl(j))
                enddo
              enddo
!     
              do j=1,mint3d
                if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,j)
                  et=gauss3d1(2,j)
                  ze=gauss3d1(3,j)
                  weight=weight3d1(j)
                elseif((lakon(i)(4:4).eq.'8').or.
     &                 (lakon(i)(4:6).eq.'20R'))
     &                 then
                  xi=gauss3d2(1,j)
                  et=gauss3d2(2,j)
                  ze=gauss3d2(3,j)
                  weight=weight3d2(j)
                elseif(lakon(i)(4:4).eq.'2') then
                  xi=gauss3d3(1,j)
                  et=gauss3d3(2,j)
                  ze=gauss3d3(3,j)
                  weight=weight3d3(j)
                elseif(lakon(i)(4:5).eq.'10') then
                  xi=gauss3d5(1,j)
                  et=gauss3d5(2,j)
                  ze=gauss3d5(3,j)
                  weight=weight3d5(j)
                elseif(lakon(i)(4:4).eq.'4') then
                  xi=gauss3d4(1,j)
                  et=gauss3d4(2,j)
                  ze=gauss3d4(3,j)
                  weight=weight3d4(j)
                elseif(lakon(i)(4:5).eq.'15') then
                  xi=gauss3d8(1,j)
                  et=gauss3d8(2,j)
                  ze=gauss3d8(3,j)
                  weight=weight3d8(j)
                elseif(lakon(i)(4:4).eq.'6') then
                  xi=gauss3d7(1,j)
                  et=gauss3d7(2,j)
                  ze=gauss3d7(3,j)
                  weight=weight3d7(j)
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
                do k=1,3
                  pgauss(k)=0.d0
                  do i1=1,nope
                    pgauss(k)=pgauss(k)+
     &                   shp(4,i1)*co(k,konl(i1))
                  enddo
                enddo
!     
                call sdvini(xstate(1,j,i),pgauss,nstate_,ncrds,
     &               i,j,layer,kspt)
!     
              enddo
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            return
          else
!     
!     internal variables are written explicitly in the input deck
!     
            do
              call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &             inl,ipoinp,inp,ipoinpc)
              if((istat.lt.0).or.(key.eq.1)) return
!     
              if(nstate_.lt.6) then
                ntot=nstate_
              else
                ntot=6
              endif
!     
              do j=1,ntot
                read(textpart(j+2)(1:20),'(f20.0)',iostat=istat) 
     &               beta(j)
                if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*INITIAL CONDITIONS%",ier)
                  return
                endif
              enddo
              read(textpart(1)(1:10),'(i10)',iostat=istat) l
              if(istat.ne.0) then
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
              if(l.gt.ne) then
                write(*,*) 
     &               '*WARNING reading *INITIAL CONDITIONS: element ',l
                write(*,*)'          exceeds the largest defined ',
     &               'element number'
                cycle
              endif
              read(textpart(2)(1:10),'(i10)',iostat=istat) k
              if(istat.eq.0) then
                do j=1,ntot
                  xstate(j,k,l)=beta(j)
                enddo
              else
                call inputerror(inpc,ipoinpc,iline,
     &               "*INITIAL CONDITIONS%",ier)
                return
              endif
!     
              if(nstate_.gt.6) then
                numberoflines=(nstate_-7)/8+1
                do ii=1,numberoflines
                  if(ii.lt.numberoflines) then
                    jmax=8
                  else
                    jmax=nstate_-ntot
                  endif
                  call getnewline(inpc,textpart,istat,n,key,iline,
     &                 ipol,inl,ipoinp,inp,ipoinpc)
                  if((istat.lt.0).or.(key.eq.1)) return
                  do j=1,jmax
                    read(textpart(j+2)(1:20),'(f20.0)',
     &                   iostat=istat) beta(j)
                    if(istat.gt.0) 
     &                   call inputerror(inpc,ipoinpc,iline,
     &                   "*INITIAL CONDITIONS%",ier)
                    return
                    xstate(ntot+j,k,l)=beta(j)
                  enddo
                  ntot=ntot+jmax
                enddo
              endif
!     
            enddo
            return
          endif
!     
        else
          write(*,*) 
     &'*WARNING reading *INITIAL CONDITIONS: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(ij)(1:index(textpart(ij),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*INITIAL CONDITIONS%")
        endif
      enddo
!     
      write(*,*) '*ERROR reading *INITIAL CONDITIONS: unknown type'
      write(*,*) '  '
      call inputerror(inpc,ipoinpc,iline,
     &     "*INITIAL CONDITIONS%",ier)
!     
      return
      end

