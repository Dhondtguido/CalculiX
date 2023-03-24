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
      subroutine cfluxs(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &  amname,nam,ntrans,trab,inotr,co,ikforc,ilforc,nk,
     &  cflux_flag,istep,istat,n,iline,ipol,inl,ipoinp,inp,nam_,
     &  namtot_,namta,amta,iaxial,ipoinpc,idefforc,ipompc,nodempc,
     &  nmpc,ikmpc,ilmpc,labmpc,iamplitudedefault,namtot,ier)
!
!     reading the input deck: *CFLUX
!
      implicit none
!
      logical cflux_flag,user,add
!
      character*1 inpc(*)
      character*20 labmpc(*)
      character*80 amplitude,amname(*)
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nodeforc(2,*),
     &  nset,nforc,nforc_,istep,istat,n,i,j,k,l,iforcdir,key,
     &  iamforc(*),nam,iamplitude,ntrans,inotr(2,*),ipos,ikforc(*),
     &  ilforc(*),nk,iline,ipol,inl,ipoinp(2,*),inp(3,*),nam_,namtot,
     &  namtot_,namta(3,*),idelay,ndirforc(*),isector,iaxial,
     &  ipoinpc(0:*),idefforc(*),ipompc(*),id,
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),iamplitudedefault,ier
!
      real*8 xforc(*),forcval,co(3,*),trab(7,*),amta(2,*)
!
      iamplitude=iamplitudedefault
      idelay=0
      user=.false.
      add=.false.
      isector=0
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *CFLUX: *CFLUX should only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.cflux_flag)) then
            do j=1,nforc
               if(ndirforc(j).eq.0) xforc(j)=0.d0
            enddo
         elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
            read(textpart(i)(11:90),'(a80)') amplitude
            do j=nam,1,-1
               if(amname(j).eq.amplitude) then
                  iamplitude=j
                  exit
               endif
            enddo
            if(j.eq.0) then
               write(*,*)'*ERROR reading *CFLUX: nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*CFLUX%",ier)
               return
            endif
            iamplitude=j
         elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
            if(idelay.ne.0) then
               write(*,*) 
     &           '*ERROR reading *CFLUX: the parameter TIME DELAY'
               write(*,*) '       is used twice in the same keyword'
               write(*,*) '       '
               call inputerror(inpc,ipoinpc,iline,
     &              "*CFLUX%",ier)
               return
            else
               idelay=1
            endif
            nam=nam+1
            if(nam.gt.nam_) then
               write(*,*) '*ERROR reading *CFLUX: increase nam_'
               ier=1
               return
            endif
            amname(nam)='
     &                                 '
            if(iamplitude.eq.0) then
               write(*,*) '*ERROR reading *CFLUX: time delay must be'
               write(*,*) '       preceded by the amplitude parameter'
               ier=1
               return
            endif
            namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
            iamplitude=nam
c            if(nam.eq.1) then
c               namtot=0
c            else
c               namtot=namta(2,nam-1)
c            endif
            namtot=namtot+1
            if(namtot.gt.namtot_) then
               write(*,*) '*ERROR cfluxes: increase namtot_'
               ier=1
               return
            endif
            namta(1,nam)=namtot
            namta(2,nam)=namtot
c            call reorderampl(amname,namta,nam)
            read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &           amta(1,namtot)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CFLUX%",ier)
               return
            endif
         elseif(textpart(i)(1:4).eq.'USER') then
            user=.true.
         elseif(textpart(i)(1:3).eq.'ADD') then
            add=.true.
         else
            write(*,*) 
     &        '*WARNING reading *CFLUX: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CFLUX%")
         endif
      enddo
!
      if(user.and.(iamplitude.ne.0)) then
         write(*,*) '*WARNING: no amplitude definition is allowed'
         write(*,*) '          for heat fluxes defined by a'
         write(*,*) '          user routine'
         iamplitude=0
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:10),'(i10)',iostat=istat) iforcdir
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*CFLUX%",ier)
            return
         endif
         if((iforcdir.ne.0).and.(iforcdir.ne.11)) then
            write(*,*) '*ERROR reading *CFLUX: nonexistent degree of '
            write(*,*) '       freedom. '
            call inputerror(inpc,ipoinpc,iline,
     &           "*CFLUX%",ier)
            return
         endif
         iforcdir=0
!
         if(textpart(3)(1:1).eq.' ') then
            forcval=0.d0
         else
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) forcval
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CFLUX%",ier)
               return
            endif
            if(iaxial.eq.180) forcval=forcval/iaxial
         endif
!
!        dummy flux consisting of the first primes
!
         if(user) forcval=1.2357111317d0
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if(l.gt.nk) then
               write(*,*) '*ERROR reading *CFLUX: node ',l
               write(*,*) '       is not defined'
               ier=1
               return
            endif
            call forcadd(l,iforcdir,forcval,
     &        nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &        iamplitude,nam,ntrans,trab,inotr,co,ikforc,ilforc,
     &        isector,add,user,idefforc,ipompc,nodempc,
     &        nmpc,ikmpc,ilmpc,labmpc)
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c            do i=1,nset
c               if(set(i).eq.noset) exit
c            enddo
            call cident81(set,noset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *CFLUX: node set ',noset
               write(*,*) '  has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*CFLUX%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
               call forcadd(ialset(j),iforcdir,forcval,
     &           nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &           iamplitude,nam,ntrans,trab,inotr,co,ikforc,ilforc,
     &           isector,add,user,idefforc,ipompc,nodempc,
     &           nmpc,ikmpc,ilmpc,labmpc)
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     call forcadd(k,iforcdir,forcval,
     &                 nodeforc,ndirforc,xforc,nforc,nforc_,
     &                 iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &                 ikforc,ilforc,isector,add,user,idefforc,
     &                 ipompc,nodempc,nmpc,ikmpc,ilmpc,labmpc)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

