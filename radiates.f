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
      subroutine radiates(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,nelemload,sideload,xload,nload,nload_,
     &     ielmat,ntmat_,iamload,amname,nam,lakon,ne,radiate_flag,
     &     istep,istat,n,iline,ipol,inl,ipoinp,inp,physcon,nam_,namtot_,
     &     namta,amta,ipoinpc,mi,iamplitudedefault,namtot,ier)
!     
!     reading the input deck: *RADIATE
!     
      implicit none
!     
      logical radiate_flag,environmentnode,surface
!     
      character*1 inpc(*)
      character*3 cavlabel
      character*8 lakon(*)
      character*20 sideload(*),label
      character*80 amname(*),amplitude
      character*81 set(*),elset
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),nelemload(2,*),mi(*),
     &     ielmat(mi(3),*),nset,nload,nload_,ntmat_,istep,istat,n,id,
     &     i,j,l,key,iload,
     &     iamload(2,*),nam,iamptemp,ipos,ne,node,iampradi,iline,ipol,
     &     inl,ipoinp(2,*),inp(3,*),nam_,namtot,namtot_,namta(3,*),
     &     idelay1,idelay2,ipoinpc(0:*),iamplitudedefault,ier
!     
      real*8 xload(2,*),xmagradi,xmagtemp,physcon(*),amta(2,*)
!     
      iamptemp=iamplitudedefault
      iampradi=0
      idelay1=0
      idelay2=0
      cavlabel='   '
!     
      environmentnode=.false.
      surface=.false.
!     
      if(istep.lt.1) then
        write(*,*) 
     &       '*ERROR reading *RADIATE: *RADIATE should only be used'
        write(*,*) '  within a STEP'
        ier=1
        return
      endif
!     
      if(physcon(2).le.0.d0) then
        write(*,*) 
     &       '*ERROR reading *RADIATE: *RADIATE card was selected'
        write(*,*) '       but no *PHYSICAL CONSTANTS card encountered'
        write(*,*) '       or incomplete (Stefan Boltzmann constant'
        write(*,*) '       is lacking)'
        write(*,*)
        ier=1
        return
      endif
!     
      do i=2,n
        if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.radiate_flag)) then
          do j=1,nload
            if(sideload(j)(1:1).eq.'R') then
              sideload(j)(3:4)='  '
              xload(1,j)=0.d0
            endif
          enddo
        elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
          read(textpart(i)(11:90),'(a80)') amplitude
          do j=nam,1,-1
            if(amname(j).eq.amplitude) then
              iamptemp=j
              exit
            endif
          enddo
          if(j.eq.0) then
            write(*,*)
     &           '*ERROR reading *RADIATE: nonexistent amplitude'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
          iamptemp=j
        elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
          if(idelay1.ne.0) then
            write(*,*) 
     &           '*ERROR reading *RADIATE: the parameter TIME DELAY'
            write(*,*) '       is used twice in the same keyword'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          else
            idelay1=1
          endif
          nam=nam+1
          if(nam.gt.nam_) then
            write(*,*) '*ERROR reading *RADIATE: increase nam_'
            ier=1
            return
          endif
          amname(nam)='
     &'
          if(iamptemp.eq.0) then
            write(*,*) '*ERROR reading *RADIATE: time delay must be'
            write(*,*) '       preceded by the amplitude parameter'
            ier=1
            return
          endif
          namta(3,nam)=sign(iamptemp,namta(3,iamptemp))
          iamptemp=nam
c     if(nam.eq.1) then
c     namtot=0
c     else
c     namtot=namta(2,nam-1)
c     endif
          namtot=namtot+1
          if(namtot.gt.namtot_) then
            write(*,*) '*ERROR radiates: increase namtot_'
            ier=1
            return
          endif
          namta(1,nam)=namtot
          namta(2,nam)=namtot
c     call reorderampl(amname,namta,nam)
          read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &         amta(1,namtot)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
        elseif(textpart(i)(1:19).eq.'RADIATIONAMPLITUDE=') then
          read(textpart(i)(20:99),'(a80)') amplitude
          do j=nam,1,-1
            if(amname(j).eq.amplitude) then
              iampradi=j
              exit
            endif
          enddo
          if(j.eq.0) then
            write(*,*)
     &           '*ERROR reading *RADIATE: nonexistent amplitude'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
          iampradi=j
        elseif(textpart(i)(1:19).eq.'RADIATIONTIMEDELAY=') THEN
          if(idelay2.ne.0) then
            write(*,*) 
     &           '*ERROR reading *RADIATE: the parameter RADIATION'
            write(*,*) '       TIME DELAY is used twice in the'
            write(*,*) '       same keyword; '
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          else
            idelay2=1
          endif
          nam=nam+1
          if(nam.gt.nam_) then
            write(*,*) '*ERROR reading *RADIATE: increase nam_'
            ier=1
            return
          endif
          amname(nam)='
     &'
          if(iampradi.eq.0) then
            write(*,*) 
     &           '*ERROR reading *RADIATE: radiation time delay'
            write(*,*) '       must be preceded by the radiation'
            write(*,*) '       amplitude parameter'
            ier=1
            return
          endif
          namta(3,nam)=sign(iampradi,namta(3,iampradi))
          iampradi=nam
c     if(nam.eq.1) then
c     namtot=0
c     else
c     namtot=namta(2,nam-1)
c     endif
          namtot=namtot+1
          if(namtot.gt.namtot_) then
            write(*,*) '*ERROR radiates: increase namtot_'
            ier=1
            return
          endif
          namta(1,nam)=namtot
          namta(2,nam)=namtot
c     call reorderampl(amname,namta,nam)
          read(textpart(i)(20:39),'(f20.0)',iostat=istat) 
     &         amta(1,namtot)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
        elseif(textpart(i)(1:7).eq.'ENVNODE') THEN
          environmentnode=.true.
        elseif(textpart(i)(1:7).eq.'CAVITY=') THEN
          read(textpart(i)(8:10),'(a3)',iostat=istat) cavlabel
        else
          write(*,*) 
     &         '*WARNING reading *RADIATE: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*RADIATE%")
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
        label(18:20)=cavlabel
!     
!     compatibility with ABAQUS for shells
!     
        if(label(2:4).eq.'NEG') label(2:4)='1  '
        if(label(2:4).eq.'POS') label(2:4)='2  '
!     
!     for plane stress elements: 'N' and 'P' are converted
!     into '5' and '6' and farther down in '1' and '2'
!     
        if(label(2:2).eq.'N') label(2:2)='5'
        if(label(2:2).eq.'P') label(2:2)='6'
!     
!     reference temperature and radiation coefficient
!     (for non uniform loading: use user routine radiation.f)
!     
        if((label(3:4).ne.'NU').and.(label(5:5).ne.'N')) then
          if(environmentnode) then
            read(textpart(3)(1:10),'(i10)',iostat=istat) node
          else
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmagtemp
            node=0
          endif
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
          read(textpart(4)(1:20),'(f20.0)',iostat=istat) xmagradi
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
        else
          if(environmentnode) then
            read(textpart(3)(1:10),'(i10)',iostat=istat) node
          else
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmagtemp
            node=0
          endif
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RADIATE%",ier)
            return
          endif
        endif
        if(((label(1:2).ne.'R1').and.(label(1:2).ne.'R2').and.
     &       (label(1:2).ne.'R0').and.
     &       (label(1:2).ne.'R3').and.(label(1:2).ne.'R4').and.
     &       (label(1:2).ne.'R5').and.(label(1:2).ne.'R6')).or.
     &       ((label(3:5).ne.'   ').and.(label(3:5).ne.'NU ').and.
     &       (label(3:5).ne.'CR ').and.(label(3:5).ne.'CRN'))) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*RADIATE%",ier)
          return
        endif
!     
        read(textpart(1)(1:10),'(i10)',iostat=istat) l
        if(istat.eq.0) then
          if(l.gt.ne) then
            write(*,*) '*ERROR reading *RADIATE: element ',l
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
            if(label(1:2).eq.'R1') then
              label(1:2)='R3'
            elseif(label(1:2).eq.'R2') then
              label(1:2)='R4'
            elseif(label(1:2).eq.'R3') then
              label(1:2)='R5'
            elseif(label(1:2).eq.'R4') then
              label(1:2)='R6'
            elseif(label(1:2).eq.'R5') then
              label(1:2)='R1'
            elseif(label(1:2).eq.'R6') then
              label(1:2)='R2'
            endif
          endif
          call loadaddt(l,label,xmagradi,xmagtemp,nelemload,sideload,
     &         xload,nload,nload_,iamload,iamptemp,iampradi,nam,node,
     &         iload)
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
              write(*,*) '*ERROR reading *RADIATE: element set '
              write(*,*) '       or facial surface ',elset
              write(*,*) '       has not yet been defined. '
              call inputerror(inpc,ipoinpc,iline,
     &             "*RADIATE%",ier)
              return
            endif
          endif
!     
          l=ialset(istartset(i))
          if(surface) then
            write(label(2:2),'(i1)') l-10*(l/10)
            l=l/10
          endif
          if((lakon(l)(1:2).eq.'CP').or.
     &         (lakon(l)(2:2).eq.'A').or.
     &         (lakon(l)(7:7).eq.'E').or.
     &         (lakon(l)(7:7).eq.'S').or.
     &         (lakon(l)(7:7).eq.'A')) then
            if(label(1:2).eq.'R1') then
              label(1:2)='R3'
            elseif(label(1:2).eq.'R2') then
              label(1:2)='R4'
            elseif(label(1:2).eq.'R3') then
              label(1:2)='R5'
            elseif(label(1:2).eq.'R4') then
              label(1:2)='R6'
            elseif(label(1:2).eq.'R5') then
              label(1:2)='R1'
            elseif(label(1:2).eq.'R6') then
              label(1:2)='R2'
            endif
          endif
!     
          do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
              l=ialset(j)
              if(surface) then
                write(label(2:2),'(i1)') l-10*(l/10)
                l=l/10
              endif
              call loadaddt(l,label,xmagradi,xmagtemp,nelemload,
     &             sideload,xload,nload,nload_,iamload,
     &             iamptemp,iampradi,nam,node,iload)
            else
              l=ialset(j-2)
              do
                l=l-ialset(j)
                if(l.ge.ialset(j-1)) exit
                call loadaddt(l,label,xmagradi,xmagtemp,nelemload,
     &               sideload,xload,nload,nload_,iamload,
     &               iamptemp,iampradi,nam,node,iload)
              enddo
            endif
          enddo
        endif
      enddo
!     
      return
      end
