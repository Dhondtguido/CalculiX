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
      subroutine dfluxs(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,nelemload,sideload,xload,nload,nload_,
     &     ielmat,ntmat_,iamload,
     &     amname,nam,lakon,ne,dflux_flag,istep,istat,n,iline,ipol,inl,
     &     ipoinp,inp,nam_,namtot_,namta,amta,ipoinpc,mi,idefload,
     &     iamplitudedefault,namtot,ier)
!     
!     reading the input deck: *DFLUX
!     
      implicit none
!     
      logical dflux_flag,surface
!     
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 sideload(*),label
      character*80 amname(*),amplitude
      character*81 set(*),elset
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),nelemload(2,*),mi(*),
     &     ielmat(mi(3),*),nset,nload,nload_,ntmat_,istep,istat,n,i,j,l,
     &     key,idefload(*),ier,iamload(2,*),nam,iamplitude,ipos,ne,
     &     iline,ipol,inl,ipoinp(2,*),id,
     &     inp(3,*),nam_,namtot,namtot_,namta(3,*),idelay,isector,
     &     ipoinpc(0:*),iamplitudedefault
!     
      real*8 xload(2,*),xmagnitude,amta(2,*)
!     
      iamplitude=iamplitudedefault
      idelay=0
      isector=0
      surface=.false.
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *DFLUX: *DFLUX should only be used'
        write(*,*) '  within a STEP'
        ier=1
        return
      endif
!     
      do i=2,n
        if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.dflux_flag)) then
          do j=1,nload
            if((sideload(j)(1:1).eq.'S').or.
     &           (sideload(j)(1:2).eq.'BF')) then
              sideload(j)(3:4)='  '
              xload(1,j)=0.d0
            endif
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
            write(*,*)'*ERROR reading *DFLUX: nonexistent amplitude'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DFLUX%",ier)
            return
          endif
          iamplitude=j
        elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
          if(idelay.ne.0) then
            write(*,*) 
     &           '*ERROR reading *DFLUX: the parameter TIME DELAY'
            write(*,*) '       is used twice in the same keyword'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DFLUX%",ier)
            return
          else
            idelay=1
          endif
          nam=nam+1
          if(nam.gt.nam_) then
            write(*,*) '*ERROR reading *DFLUX: increase nam_'
            ier=1
            return
          endif
          amname(nam)='
     &'
          if(iamplitude.eq.0) then
            write(*,*) '*ERROR reading *DFLUX: time delay must be'
            write(*,*) '       preceded by the amplitude parameter'
            ier=1
            return
          endif
          namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
          iamplitude=nam
          namtot=namtot+1
          if(namtot.gt.namtot_) then
            write(*,*) '*ERROR dfluxes: increase namtot_'
            ier=1
            return
          endif
          namta(1,nam)=namtot
          namta(2,nam)=namtot
          read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &         amta(1,namtot)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DFLUX%",ier)
            return
          endif
        else
          write(*,*) 
     &         '*WARNING reading *DFLUX: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DFLUX%")
        endif
      enddo
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
!     
        read(textpart(2)(1:20),'(a20)',iostat=istat) label
!     
!     compatibility with ABAQUS for shells
!     
        if(label(2:4).eq.'NEG') label(2:4)='1  '
        if(label(2:4).eq.'POS') label(2:4)='2  '
!     
!     for plane stress elements: 'N' and 'P' are converted
!     into '5' and '6' and further down in '1' and '2'
!     
        if(label(2:2).eq.'N') label(2:2)='5'
        if(label(2:2).eq.'P') label(2:2)='6'
!     
        read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmagnitude
!     
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DFLUX%",ier)
          return
        endif
        if(((label(1:2).ne.'S1').and.(label(1:2).ne.'S2').and.
     &       (label(1:2).ne.'S0').and.
     &       (label(1:2).ne.'S3').and.(label(1:2).ne.'S4').and.
     &       (label(1:2).ne.'S5').and.(label(1:2).ne.'S6').and.
     &       (label(1:2).ne.'BF').and.(label(1:2).ne.'S ')).or.
     &       ((label(3:4).ne.'  ').and.(label(3:4).ne.'NU'))) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DFLUX%",ier)
          return
        endif
!     
        read(textpart(1)(1:10),'(i10)',iostat=istat) l
        if(istat.eq.0) then
          if(l.gt.ne) then
            write(*,*) '*ERROR reading *DFLUX: element ',l
            write(*,*) '       is not defined'
            ier=1
            return
          endif
!     
          if((lakon(l)(1:2).eq.'CP').or.
     &         (lakon(l)(2:2).eq.'A').or.
     &         (lakon(l)(7:7).eq.'E').or.
     &         (lakon(l)(7:7).eq.'S').or.
     &         (lakon(l)(7:7).eq.'A')) then
            if(label(1:2).eq.'S1') then
              label(1:2)='S3'
            elseif(label(1:2).eq.'S2') then
              label(1:2)='S4'
            elseif(label(1:2).eq.'S3') then
              label(1:2)='S5'
            elseif(label(1:2).eq.'S4') then
              label(1:2)='S6'
            elseif(label(1:2).eq.'S5') then
              label(1:2)='S1'
            elseif(label(1:2).eq.'S6') then
              label(1:2)='S2'
            endif
c          elseif((lakon(l)(1:1).eq.'B').or.
c     &           (lakon(l)(7:7).eq.'B')) then
c          elseif((lakon(l)(1:1).eq.'S').or.
c     &           (lakon(l)(7:7).eq.'L')) then
          endif
          call loadadd(l,label,xmagnitude,nelemload,sideload,
     &         xload,nload,nload_,iamload,iamplitude,
     &         nam,isector,idefload)
        else
          read(textpart(1)(1:80),'(a80)',iostat=istat) elset
          elset(81:81)=' '
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
c          do i=1,nset
c            if(set(i).eq.elset) exit
c          enddo
          call cident81(set,elset,nset,id)
          i=nset+1
          if(id.gt.0) then
            if(elset.eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
!     
!     check for facial surface
!     
            surface=.true.
            elset(ipos:ipos)='T'
c            do i=1,nset
c              if(set(i).eq.elset) exit
c            enddo
            call cident81(set,elset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(elset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
              elset(ipos:ipos)=' '
              write(*,*) '*ERROR reading *DFLUX: element set '
              write(*,*) '       or facial surface ',elset
              write(*,*) '       has not yet been defined. '
              call inputerror(inpc,ipoinpc,iline,
     &             "*DFLUX%",ier)
              return
            endif
          endif
!     
          l=ialset(istartset(i))
          if(surface) then
            write(label(2:2),'(i1)') l-10*(l/10)
            l=l/10
          endif
!     
          if((lakon(l)(1:2).eq.'CP').or.
     &         (lakon(l)(2:2).eq.'A').or.
     &         (lakon(l)(7:7).eq.'E').or.
     &         (lakon(l)(7:7).eq.'S').or.
     &         (lakon(l)(7:7).eq.'A')) then
            if(label(1:2).eq.'S1') then
              label(1:2)='S3'
            elseif(label(1:2).eq.'S2') then
              label(1:2)='S4'
            elseif(label(1:2).eq.'S3') then
              label(1:2)='S5'
            elseif(label(1:2).eq.'S4') then
              label(1:2)='S6'
            elseif(label(1:2).eq.'S5') then
              label(1:2)='S1'
            elseif(label(1:2).eq.'S6') then
              label(1:2)='S2'
            endif
c          elseif((lakon(l)(1:1).eq.'B').or.
c     &           (lakon(l)(7:7).eq.'B')) then
c          elseif((lakon(l)(1:1).eq.'S').or.
c     &           (lakon(l)(7:7).eq.'L')) then
          endif
          do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
              l=ialset(j)
              if(surface) then
                write(label(2:2),'(i1)') l-10*(l/10)
                l=l/10
              endif
              call loadadd(l,label,xmagnitude,nelemload,sideload,
     &             xload,nload,nload_,iamload,iamplitude,
     &             nam,isector,idefload)
            else
              l=ialset(j-2)
              do
                l=l-ialset(j)
                if(l.ge.ialset(j-1)) exit
                call loadadd(l,label,xmagnitude,nelemload,
     &               sideload,xload,nload,nload_,
     &               iamload,iamplitude,nam,isector,idefload)
              enddo
            endif
          enddo
        endif
      enddo
!     
      return
      end
