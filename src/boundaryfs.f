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
      subroutine boundaryfs(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nodeboun,ndirboun,xboun,nboun,nboun_,nk,
     &  iamboun,amname,nam,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,trab,ntrans,ikboun,ilboun,ikmpc,ilmpc,nk_,
     &  co,labmpc,typeboun,istat,n,iline,ipol,
     &  inl,ipoinp,inp,nam_,namtot_,namta,amta,nmethod,iperturb,
     &  ipoinpc,vold,mi,xload,sideload,nload,nelemload,lakon,kon,
     &  ipkon,ne,iamplitudedefault,namtot,ier)
!
!     reading the input deck: *BOUNDARYF
!
!     (boundary conditions for CFD-calculations)
!
      implicit none
!
      logical user,fixed,surface
!
      character*1 typeboun(*),type,inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*80 amname(*),amplitude
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nodeboun(*),
     &  ndirboun(*),iface,nload,nelemload(2,*),kon(*),ipkon(*),
     &  nset,nboun,nboun_,istat,n,i,j,k,l,ibounstart,ibounend,
     &  key,nk,iamboun(*),nam,iamplitude,ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,mpcfree,ikboun(*),ilboun(*),ikmpc(*),
     &  ilmpc(*),ntrans,nk_,ipos,m,ne,id,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),nam_,namtot,namtot_,
     &  namta(3,*),idelay,nmethod,iperturb(*),ipoinpc(0:*),
     &  mi(*),iamplitudedefault,ier
!
      real*8 xboun(*),bounval,coefmpc(*),trab(7,*),co(3,*),amta(2,*),
     &  vold(0:mi(2),*),xload(2,*)
!
      type='F'
      iamplitude=iamplitudedefault
      idelay=0
      user=.false.
      fixed=.false.
      surface=.false.
!
      do i=2,n
         if(textpart(i)(1:10).eq.'AMPLITUDE=') then
            read(textpart(i)(11:90),'(a80)') amplitude
            do j=nam,1,-1
               if(amname(j).eq.amplitude) then
                  iamplitude=j
                  exit
               endif
            enddo
            if(j.eq.0) then
               write(*,*)
     &           '*ERROR reading *BOUNDARYF: nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
            iamplitude=j
         elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
            if(idelay.ne.0) then
               write(*,*)'*ERROR reading *BOUNDARYF: the parameter TIME'
               write(*,*) '       DELAY is used twice in the same'
               write(*,*) '       keyword; '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            else
               idelay=1
            endif
            nam=nam+1
            if(nam.gt.nam_) then
               write(*,*) '*ERROR reading *BOUNDARYF: increase nam_'
               ier=1
               return
            endif
            amname(nam)='
     &                                 '
            if(iamplitude.eq.0) then
               write(*,*)'*ERROR reading *BOUNDARYF: time delay must be'
               write(*,*) '       preceded by the amplitude parameter'
               ier=1
               return
            endif
            namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
            iamplitude=nam
            namtot=namtot+1
            if(namtot.gt.namtot_) then
               write(*,*) '*ERROR boundaries: increase namtot_'
               ier=1
               return
            endif
            namta(1,nam)=namtot
            namta(2,nam)=namtot
            read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &           amta(1,namtot)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
         elseif(textpart(i)(1:4).eq.'USER') then
            user=.true.
         else
            write(*,*) 
     &        '*WARNING reading *BOUNDARYF: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BOUNDARYF%")
         endif
      enddo
!
      if(user.and.(iamplitude.ne.0)) then
         write(*,*) '*WARNING: no amplitude definition is allowed'
         write(*,*) '          for temperatures defined by a'
         write(*,*) '          user routine'
         iamplitude=0
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(3)(1:10),'(i10)',iostat=istat) ibounstart
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*BOUNDARYF%",ier)
            return
         endif
!     
         if(textpart(4)(1:1).eq.' ') then
            ibounend=ibounstart
         else
            read(textpart(4)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
         endif
!     
         if(textpart(5)(1:1).eq.' ') then
            bounval=0.d0
         else
            read(textpart(5)(1:20),'(f20.0)',iostat=istat) bounval
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
         endif
!     
!        dummy boundary condition consisting of the first primes
!
         if(user) bounval=1.2357111317d0
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if((l.gt.ne).or.(l.le.0)) then
               write(*,*) '*ERROR reading *BOUNDARYF:'
               write(*,*) '       element ',l,' is not defined'
               ier=1
               return
            endif
            read(textpart(2)(2:2),'(i1)',iostat=istat) iface
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
            l=10*l+iface
            call bounaddf(l,ibounstart,ibounend,bounval,
     &        nodeboun,ndirboun,xboun,nboun,nboun_,
     &        iamboun,iamplitude,nam,ipompc,nodempc,
     &        coefmpc,nmpc,nmpc_,mpcfree,trab,
     &        ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &        type,typeboun,nmethod,iperturb,vold,mi,
     &        nelemload,sideload,xload,nload,lakon,ipkon,kon)
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
c            do i=1,nset
c               if(set(i).eq.elset) exit
c            enddo
            call cident81(set,elset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(elset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
!
!              check for facial surface
!
               surface=.true.
               elset(ipos:ipos)='T'
c               do i=1,nset
c                  if(set(i).eq.elset) exit
c               enddo
               call cident81(set,elset,nset,id)
               i=nset+1
               if(id.gt.0) then
                 if(elset.eq.set(id)) then
                   i=id
                 endif
               endif
               if(i.gt.nset) then
                  elset(ipos:ipos)=' '
                  write(*,*) '*ERROR reading *BOUNDARYF: surface ',elset
                  write(*,*) '       has not yet been defined. '
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*BOUNDARYF%",ier)
                  return
               endif
            endif
            read(textpart(2)(2:2),'(i1)',iostat=istat) iface
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BOUNDARYF%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  k=ialset(j)
                  if(.not.surface) k=10*k+iface
                  call bounaddf(k,ibounstart,ibounend,bounval,
     &               nodeboun,ndirboun,xboun,nboun,nboun_,
     &               iamboun,iamplitude,nam,ipompc,nodempc,
     &               coefmpc,nmpc,nmpc_,mpcfree,trab,
     &               ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &               type,typeboun,nmethod,iperturb,vold,mi,
     &               nelemload,sideload,xload,nload,lakon,ipkon,kon)
               else
                  m=ialset(j-2)
                  do
                     m=m-ialset(j)
                     if(m.ge.ialset(j-1)) exit
                     k=10*m+iface
                     call bounaddf(k,ibounstart,ibounend,bounval,
     &                 nodeboun,ndirboun,xboun,nboun,nboun_,
     &                 iamboun,iamplitude,nam,ipompc,nodempc,
     &                 coefmpc,nmpc,nmpc_,mpcfree,trab,
     &                 ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &                 labmpc,type,typeboun,nmethod,iperturb,
     &                 vold,mi,
     &                 nelemload,sideload,xload,nload,lakon,ipkon,kon)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

