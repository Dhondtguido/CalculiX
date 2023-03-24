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
      subroutine dloads(inpc,textpart,set,istartset,iendset,
     &     ialset,nset,nelemload,sideload,xload,nload,nload_,
     &     ielmat,iamload,amname,nam,lakon,ne,dload_flag,istep,
     &     istat,n,iline,ipol,inl,ipoinp,inp,cbody,ibody,xbody,nbody,
     &     nbody_,xbodyold,iperturb,physcon,nam_,namtot_,namta,amta,
     &     nmethod,ipoinpc,maxsectors,mi,idefload,idefbody,ipkon,
     &     thicke,iamplitudedefault,namtot,ier)
!     
!     reading the input deck: *DLOAD
!     
      implicit none
!     
      logical dload_flag,submodel,edgeload,surface
!     
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 sideload(*),label
      character*80 amname(*),amplitude
      character*81 set(*),elset,cbody(*)
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),nelemload(2,*),mi(*),
     &     ielmat(mi(3),*),nset,nload,nload_,istep,istat,n,i,j,l,key,
     &     iamload(2,*),nam,iamplitude,ipos,ne,iline,ipol,iperturb(*),
     &     inl,ipoinp(2,*),inp(3,*),ibody(3,*),nbody,nbody_,nam_,namtot,
     &     namtot_,namta(3,*),idelay,nmethod,lc,isector,node,id,
     &     ipoinpc(0:*),maxsectors,jsector,iglobstep,idefload(*),
     %     idefbody(*),ipkon(*),k,indexe,iamplitudedefault,ier
!     
      real*8 xload(2,*),xbody(7,*),xmagnitude,dd,p1(3),p2(3),bodyf(3),
     &     xbodyold(7,*),physcon(*),amta(2,*),xxmagnitude,
     &     thicke(mi(3),*),thickness
!     
      iamplitude=iamplitudedefault
      idelay=0
      lc=1
      isector=0
      submodel=.false.
      iglobstep=0
      edgeload=.false.
      surface=.false.
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *DLOAD: *DLOAD should only be used'
        write(*,*) '  within a STEP'
        ier=1
        return
      endif
!     
      do i=2,n
        if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.dload_flag)) then
          do j=1,nload
            if(sideload(j)(1:1).eq.'P') then
              sideload(j)(3:4)='  '
              xload(1,j)=0.d0
            endif
          enddo
          do j=1,nbody
            xbody(1,j)=0.d0
          enddo
        elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
          read(textpart(i)(11:90),'(a80)') amplitude
          do j=1,nam
            if(amname(j).eq.amplitude) then
              iamplitude=j
              exit
            endif
          enddo
          if(j.gt.nam) then
            write(*,*)'*ERROR reading *DLOAD: nonexistent amplitude'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          endif
          iamplitude=j
        elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
          if(idelay.ne.0) then
            write(*,*) 
     &           '*ERROR reading *DLOAD: the parameter TIME DELAY'
            write(*,*) '       is used twice in the same keyword'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          else
            idelay=1
          endif
          nam=nam+1
          if(nam.gt.nam_) then
            write(*,*) '*ERROR reading *DLOAD: increase nam_'
            ier=1
            return
          endif
          amname(nam)='
     &'
          if(iamplitude.eq.0) then
            write(*,*) '*ERROR reading *DLOAD: time delay must be'
            write(*,*) '       preceded by the amplitude parameter'
            ier=1
            return
          endif
          namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
          iamplitude=nam
c     if(nam.eq.1) then
c     namtot=0
c     else
c     namtot=namta(2,nam-1)
c     endif
          namtot=namtot+1
          if(namtot.gt.namtot_) then
            write(*,*) '*ERROR dloads: increase namtot_'
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
     &           "*DLOAD%",ier)
            return
          endif
        elseif(textpart(i)(1:9).eq.'LOADCASE=') then
          read(textpart(i)(10:19),'(i10)',iostat=istat) lc
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          endif
          if(nmethod.ne.5) then
            write(*,*) 
     &           '*ERROR reading *DLOAD: the parameter LOAD CASE'
            write(*,*) '       is only allowed in STEADY STATE'
            write(*,*) '       DYNAMICS calculations'
            ier=1
            return
          endif
        elseif(textpart(i)(1:7).eq.'SECTOR=') then
          read(textpart(i)(8:17),'(i10)',iostat=istat) isector
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          endif
          if((nmethod.le.3).or.(iperturb(1).gt.1)) then
            write(*,*) '*ERROR reading *DLOAD: the parameter SECTOR'
            write(*,*) '       is only allowed in MODAL DYNAMICS or'
            write(*,*) '       STEADY STATE DYNAMICS calculations'
            ier=1
            return
          endif
          if(isector.gt.maxsectors) then
            write(*,*) '*ERROR reading *DLOAD: sector ',isector
            write(*,*) '       exceeds number of sectors'
            ier=1
            return
          endif
          isector=isector-1
        elseif(textpart(i)(1:8).eq.'SUBMODEL') then
          submodel=.true.
        elseif(textpart(i)(1:5).eq.'STEP=') then
          read(textpart(i)(6:15),'(i10)',iostat=istat) iglobstep
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          endif
        elseif(textpart(i)(1:8).eq.'DATASET=') then
          read(textpart(i)(9:18),'(i10)',iostat=istat) iglobstep
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*DLOAD%",ier)
            return
          endif
!     
!     the mode number for submodels
!     is stored as a negative global step
!     
          iglobstep=-iglobstep
        else
          write(*,*) 
     &         '*WARNING reading *DLOAD: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DLOAD%")
        endif
      enddo
!     
!     check for the presence of an amplitude in submodel cases
!     
      if(submodel) then
        if(iamplitude.ne.0) then
          write(*,*) '*WARNING reading *DSLOAD:'
          write(*,*) '         no amplitude definition is allowed'
          write(*,*) '         in combination with a submodel'
        endif
      endif
!     
!     check whether global step was specified for submodel
!     
      if((submodel).and.(iglobstep.eq.0)) then
        write(*,*) '*ERROR reading *DLOAD: no global step'
        write(*,*) '       step specified for the submodel'
        call inputerror(inpc,ipoinpc,iline,
     &       "*DLOAD%",ier)
        return
      endif
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
!     
        read(textpart(2)(1:20),'(a20)',iostat=istat) label
!     
!     for submodels the load label is modified and the global
!     step is stored in iamload(1,*)
!     
        if(submodel) then
          label(3:4)='SM'
          iamplitude=iglobstep
        endif
!     
        if(label(3:4).ne.'NP') then
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmagnitude
        else
          read(textpart(3)(1:10),'(i10)',iostat=istat) node
        endif
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DLOAD%",ier)
          return
        endif
        if(label(1:7).eq.'CENTRIF') then
          do i=1,3
            read(textpart(i+3)(1:20),'(f20.0)',iostat=istat) p1(i)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DLOAD%",ier)
              return
            endif
          enddo
          do i=1,3
            read(textpart(i+6)(1:20),'(f20.0)',iostat=istat) p2(i)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DLOAD%",ier)
              return
            endif
          enddo
          dd=dsqrt(p2(1)**2+p2(2)**2+p2(3)**2)
          do i=1,3
            p2(i)=p2(i)/dd
          enddo
        elseif(label(1:4).eq.'GRAV') then
          do i=1,3
            read(textpart(i+3)(1:20),'(f20.0)',iostat=istat) bodyf(i)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DLOAD%",ier)
              return
            endif
          enddo
        elseif(label(1:6).eq.'NEWTON') then
          if(iperturb(1).le.1) then
            write(*,*) '*ERROR reading *DLOAD: NEWTON gravity force'
            write(*,*) '       can only be used in a nonlinear'
            write(*,*) '       procedure'
            ier=1
            return
          endif
          if(physcon(3).le.0.d0) then
            write(*,*) '*ERROR reading *DLOAD: NEWTON gravity force'
            write(*,*) '       requires the definition of a'
            write(*,*) '       positive gravity constant with'
            write(*,*) '       a *PHYSICAL CONSTANTS card'
            ier=1
            return
          endif
        elseif(((label(1:2).ne.'P1').and.(label(1:2).ne.'P2').and.
     &         (label(1:2).ne.'P3').and.(label(1:2).ne.'P4').and.
     &         (label(1:2).ne.'P5').and.(label(1:2).ne.'P6').and.
     &         (label(1:2).ne.'P ').and.(label(1:2).ne.'BX').and.
     &         (label(1:2).ne.'BY').and.(label(1:2).ne.'BZ').and.
c     BernhardiStart
     &         (label(1:2).ne.'ED')).or.
     &         ((label(3:6).ne.'NOR1').and.(label(3:6).ne.'NOR2').and.
     &         (label(3:6).ne.'NOR3').and.(label(3:6).ne.'NOR4')).and.
c     BernhardiEnd
     &         ((label(3:4).ne.'  ').and.(label(3:4).ne.'NU').and.
     &         (label(3:4).ne.'NP').and.(label(3:4).ne.'SM'))) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*DLOAD%",ier)
          return
        endif
!     
        read(textpart(1)(1:10),'(i10)',iostat=istat) l
        if(istat.eq.0) then
          if(l.gt.ne) then
            write(*,*) '*ERROR reading *DLOAD: element ',l
            write(*,*) '       is not defined'
            ier=1
            return
          endif
          if((label(1:7).eq.'CENTRIF').or.(label(1:4).eq.'GRAV').or.
     &         (label(1:6).eq.'NEWTON')) then
            elset(1:80)=textpart(1)(1:80)
            elset(81:81)=' '
            call bodyadd(cbody,ibody,xbody,nbody,nbody_,elset,label,
     &           iamplitude,xmagnitude,p1,p2,bodyf,xbodyold,lc,idefbody)
          else
            xxmagnitude=xmagnitude
            if((lakon(l)(1:2).eq.'CP').or.
     &           (lakon(l)(2:2).eq.'A').or.
     &           (lakon(l)(7:7).eq.'E').or.
     &           (lakon(l)(7:7).eq.'S').or.
     &           (lakon(l)(7:7).eq.'A')) then
              if(label(1:2).eq.'P1') then
                label(1:2)='P3'
              elseif(label(1:2).eq.'P2') then
                label(1:2)='P4'
              elseif(label(1:2).eq.'P3') then
                label(1:2)='P5'
              elseif(label(1:2).eq.'P4') then
                label(1:2)='P6'
              endif
            elseif((lakon(l)(1:1).eq.'B').or.
     &             (lakon(l)(7:7).eq.'B')) then
              if(label(1:2).eq.'P2') label(1:2)='P5'
            elseif((lakon(l)(1:1).eq.'S').or.
     &             (lakon(l)(7:7).eq.'L')) then
c     BernhardiStart
              if(label(1:6).eq.'EDNOR1') then
                label(1:2)='P3'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR2') then
                label(1:2)='P4'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR3') then
                label(1:2)='P5'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR4') then
                label(1:2)='P6'
                edgeload=.true.
              else
                label(1:2)='P1'
              endif
!     
!     EDNOR is an edge load
!     
              if(edgeload) then
                indexe=ipkon(l)
                thickness=0.d0
                do k=1,mi(3)
                  if(ielmat(k,l).ne.0) then
                    thickness=thickness+thicke(k,indexe+1)
                  else
                    exit
                  endif
                enddo
                xxmagnitude=xmagnitude/thickness
              endif
c     BernhardiEnd
            endif
            if(lc.ne.1) then
              jsector=isector+maxsectors
            else
              jsector=isector
            endif
            if(label(3:4).ne.'NP') then
              call loadadd(l,label,xxmagnitude,nelemload,sideload,
     &             xload,nload,nload_,iamload,iamplitude,
     &             nam,jsector,idefload)
            else
              call loadaddp(l,label,nelemload,sideload,
     &             xload,nload,nload_,iamload,iamplitude,
     &             nam,node)
            endif
          endif
        else
          read(textpart(1)(1:80),'(a80)',iostat=istat) elset
          elset(81:81)=' '
          ipos=index(elset,' ')
!     
!     check for element set
!     
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
              write(*,*) '*ERROR reading *DLOAD: element set '
              write(*,*) '       or facial surface ',elset
              write(*,*) '       has not yet been defined. '
              call inputerror(inpc,ipoinpc,iline,
     &             "*DLOAD%",ier)
              return
            endif
          endif
!     
          if((label(1:7).eq.'CENTRIF').or.(label(1:4).eq.'GRAV').or.
     &         (label(1:6).eq.'NEWTON')) then
            call bodyadd(cbody,ibody,xbody,nbody,nbody_,elset,label,
     &           iamplitude,xmagnitude,p1,p2,bodyf,xbodyold,lc,idefbody)
          else
            l=ialset(istartset(i))
            if(surface) then
              write(label(2:2),'(i1)') l-10*(l/10)
              l=l/10
            endif
            if((lakon(l)(1:2).eq.'CP').or.
     &           (lakon(l)(2:2).eq.'A').or.
     &           (lakon(l)(7:7).eq.'E').or.
     &           (lakon(l)(7:7).eq.'S').or.
     &           (lakon(l)(7:7).eq.'A')) then
              if(label(1:2).eq.'P1') then
                label(1:2)='P3'
              elseif(label(1:2).eq.'P2') then
                label(1:2)='P4'
              elseif(label(1:2).eq.'P3') then
                label(1:2)='P5'
              elseif(label(1:2).eq.'P4') then
                label(1:2)='P6'
              endif
            elseif((lakon(l)(1:1).eq.'B').or.
     &             (lakon(l)(7:7).eq.'B')) then
              if(label(1:2).eq.'P2') label(1:2)='P5'
            elseif((lakon(l)(1:1).eq.'S').or.
     &             (lakon(l)(7:7).eq.'L')) then
c     BernhardiStart
              if(label(1:6).eq.'EDNOR1') then
                label(1:2)='P3'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR2') then
                label(1:2)='P4'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR3') then
                label(1:2)='P5'
                edgeload=.true.
              elseif(label(1:6).eq.'EDNOR4') then
                label(1:2)='P6'
                edgeload=.true.
              else
                label(1:2)='P1'
              endif
c     BernhardiEnd
            endif
!     
            do j=istartset(i),iendset(i)
              if(ialset(j).gt.0) then
                l=ialset(j)
                if(surface) then
                  write(label(2:2),'(i1)') l-10*(l/10)
                  l=l/10
                endif
                xxmagnitude=xmagnitude
!     
!     EDNOR is an edge load
!     
                if(edgeload) then
                  indexe=ipkon(l)
                  thickness=0.d0
                  do k=1,mi(3)
                    if(ielmat(k,l).ne.0) then
                      thickness=thickness+thicke(k,indexe+1)
                    else
                      exit
                    endif
                  enddo
                  xxmagnitude=xmagnitude/thickness
                endif
!     
                if(lc.ne.1) then
                  jsector=isector+maxsectors
                else
                  jsector=isector
                endif
                if(label(3:4).ne.'NP') then
                  call loadadd(l,label,xxmagnitude,nelemload,
     &                 sideload,xload,nload,nload_,iamload,
     &                 iamplitude,nam,jsector,idefload)
                else
                  call loadaddp(l,label,nelemload,
     &                 sideload,xload,nload,nload_,iamload,
     &                 iamplitude,nam,node)
                endif
              else
                l=ialset(j-2)
                do
                  l=l-ialset(j)
                  if(l.ge.ialset(j-1)) exit
                  xxmagnitude=xmagnitude
!     
!     EDNOR is an edge load
!     
                  if(edgeload) then
                    indexe=ipkon(l)
                    thickness=0.d0
                    do k=1,mi(3)
                      if(ielmat(k,l).ne.0) then
                        thickness=thickness+thicke(k,indexe+1)
                      else
                        exit
                      endif
                    enddo
                    xxmagnitude=xmagnitude/thickness
                  endif
!     
                  if(lc.ne.1) then
                    jsector=isector+maxsectors
                  else
                    jsector=isector
                  endif
                  if(label(3:4).ne.'NP') then
                    call loadadd(l,label,xxmagnitude,nelemload,
     &                   sideload,xload,nload,nload_,
     &                   iamload,iamplitude,nam,jsector,idefload)
                  else
                    call loadaddp(l,label,nelemload,
     &                   sideload,xload,nload,nload_,
     &                   iamload,iamplitude,nam,node)
                  endif
                enddo
              endif
            enddo
          endif
        endif
      enddo
!     
      return
      end
