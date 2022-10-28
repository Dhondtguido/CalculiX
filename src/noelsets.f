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
      subroutine noelsets(inpc,textpart,set,istartset,iendset,ialset,
     &     nset,nset_,nalset,nalset_,nk,ne,irstrt,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
!     reading the input deck: *NSET and *ELSET
!     
      implicit none
!     
      logical igen
!     
      character*1 inpc(*)
      character*81 set(*),noelset
      character*132 textpart(16)
!     
      integer nset,nset_,nalset,nalset_,istep,istat,n,key,i,nk,ne,
     &     kode,ipos,j,k,m,iset,nn,irstrt(*),istartset(*),iendset(*),
     &     ialset(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),
     &     ier,id
!     
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
        write(*,*) 
     &       '*ERROR reading *NSET/ELSET: *NSET/*ELSET should be placed'
        write(*,*) '  before all step definitions'
        ier=1
        return
      endif
!     
      igen=.false.
!     
!     reading the name of the set
!     
      if(textpart(1)(1:5).eq.'*NSET') then
        do i=2,n
          if(textpart(i)(1:5).eq.'NSET=') then
            noelset(1:80)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
              write(*,*) 
     &             '*ERROR reading *NSET/ELSET: set name too long'
              write(*,*) '       (more than 80 characters)'
              write(*,*) '       set name:',textpart(2)(1:132)
              ier=1
              return
            endif
            noelset(81:81)=' '
            ipos=index(noelset,' ')
            noelset(ipos:ipos)='N'
            kode=0
          elseif(textpart(i)(1:8).eq.'GENERATE') then
            igen=.true.
          else
            write(*,*) 
     &        '*WARNING reading *NSET/ELSET: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &           "*NSET or *ELSET%")
          endif
        enddo
      else
        do i=2,n
          if(textpart(i)(1:6).eq.'ELSET=') then
            noelset(1:80)=textpart(i)(7:86)
            if(textpart(i)(87:87).ne.' ') then
              write(*,*) 
     &             '*ERROR reading *NSET/ELSET: set name too long'
              write(*,*) '       (more than 80 characters)'
              write(*,*) '       set name',textpart(2)(1:132)
              ier=1
              return
            endif
            noelset(81:81)=' '
            ipos=index(noelset,' ')
            noelset(ipos:ipos)='E'
            kode=1
          elseif(textpart(i)(1:8).eq.'GENERATE') then
            igen=.true.
          else
            write(*,*) 
     &         '*WARNING reading *NSET/ELSET: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &           "*NSET or *ELSET%")
          endif
        enddo
      endif
!     
!     check whether new set or old set 
!     
ccc   to remove start      
c     do iset=1,nset
c     if(set(iset).eq.noelset) then
ccc   to remove end
      iset=0
      call cident81(set,noelset,nset,id)
      if(id.gt.0) then
        if(set(id).eq.noelset) then
          iset=id
!     
!     existent set
!     
          if(iendset(iset).ne.nalset) then
!     
!     rearranging set information towards the end
!     
            nn=iendset(iset)-istartset(iset)+1
            if(nalset+nn.gt.nalset_) then
              write(*,*)
     &             '*ERROR reading *NSET/ELSET: increase nalset_'
              ier=1
              return
            endif
            do k=1,nn
              ialset(nalset+k)=ialset(istartset(iset)+k-1)
            enddo
            if(nn.gt.0) then
              do k=istartset(iset),nalset
                ialset(k)=ialset(k+nn)
              enddo
              do k=1,nset
                if(istartset(k).gt.iendset(iset)) then
                  istartset(k)=istartset(k)-nn
                  iendset(k)=iendset(k)-nn
                endif
              enddo
            endif
            istartset(iset)=nalset-nn+1
            iendset(iset)=nalset
          endif
        endif
      endif
ccc   to remove start
c     enddo
c     if(iset.gt.nset) then
ccc   to remove end
      if(iset.eq.0) then
        nset=nset+1
        if(nset.gt.nset_) then
          write(*,*) '*ERROR reading *NSET/ELSET: increase nset_'
          ier=1
          return
        endif
ccc   to remove start          
c     set(nset)=noelset
c     istartset(nset)=nalset+1
c     iendset(nset)=0
c     iset=nset
ccc   to remove end
        do j=nset,id+2,-1
          istartset(j)=istartset(j-1)
          iendset(j)=iendset(j-1)
          set(j)=set(j-1)
        enddo
        set(id+1)=noelset
        istartset(id+1)=nalset+1
        iendset(id+1)=0
        iset=id+1
      endif
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) then
ccc   to remove start           
c     if(iendset(nset).eq.0) then
ccc   to remove end
          if(iendset(iset).eq.0) then
            do j=iset+1,nset
              istartset(j-1)=istartset(j)
              iendset(j-1)=iendset(j)
              set(j-1)=set(j)
            enddo
            nset=nset-1
          endif
          return
        endif
        if(igen) n=3
        if(nalset+n.gt.nalset_) then
          write(*,*) '*ERROR reading *NSET/ELSET: increase nalset_'
          ier=1
          return
        endif
!     
        if(igen) then
          if(textpart(3)(1:1).eq.' ') then
            textpart(3)='1
     &
     &'
          endif
          do i=1,3
            read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &           ialset(nalset+i)
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*NSET or *ELSET%",ier)
              return
            endif
          enddo
          if(kode.eq.0) then
            if(ialset(nalset+1).gt.nk) then
              write(*,*) 
     &             '*ERROR reading *NSET/ELSET: starting value in'
              write(*,*) '       set ',
     &             set(iset)(1:index(set(iset),' ')-2),' > nk'
              ier=1
              return
            elseif(ialset(nalset+2).gt.nk) then
              write(*,*) 
     &             '*WARNING reading *NSET/ELSET: end value in'
              write(*,*) '         set ',
     &             set(iset)(1:index(set(iset),' ')-2),' > nk;'
              write(*,*) '         replaced by nk'
              ialset(nalset+2)=nk
            elseif(ialset(nalset+3).le.0) then
              write(*,*) '*ERROR reading *NSET/ELSET: increment in'
              write(*,*) '       set ',
     &             set(iset)(1:index(set(iset),' ')-2),' <=0'
              ier=1
              return
            endif
          else
            if(ialset(nalset+1).gt.ne) then
              write(*,*) 
     &             '*ERROR reading *NSET/ELSET: starting value in'
              write(*,*) '       set ',
     &             set(iset)(1:index(set(iset),' ')-2),' > ne'
              ier=1
              return
            elseif(ialset(nalset+2).gt.ne) then
              write(*,*) 
     &             '*WARNING reading *NSET/ELSET: end value in'
              write(*,*) '         set ',
     &             set(iset)(1:index(set(iset),' ')-2),' > ne;'
              write(*,*) '         replaced by ne'
              ialset(nalset+2)=nk
            elseif(ialset(nalset+3).le.0) then
              write(*,*) '*ERROR reading *NSET/ELSET: increment in'
              write(*,*) '       set ',
     &             set(iset)(1:index(set(iset),' ')-2),' <=0'
              ier=1
              return
            endif
          endif
          if(ialset(nalset+1).eq.ialset(nalset+2)) then
            ialset(nalset+2)=0
            ialset(nalset+3)=0
            nalset=nalset+1
          else
            ialset(nalset+3)=-ialset(nalset+3)
            nalset=nalset+3
          endif
          iendset(iset)=nalset
        else
          do i=1,n
            read(textpart(i)(1:10),'(i10)',iostat=istat) 
     &           ialset(nalset+1)
            if(istat.gt.0) then
!     
!     set name
!     
              noelset=textpart(i)(1:80)
              noelset(81:81)=' '
              ipos=index(noelset,' ')
              if(kode.eq.0) then
                noelset(ipos:ipos)='N'
              else
                noelset(ipos:ipos)='E'
              endif
ccc   to remove start                   
c     do j=1,nset
c     if(j.eq.iset)cycle
c     if(noelset.eq.set(j)) then
c     m=iendset(j)-istartset(j)+1
c     do k=1,m
c     ialset(nalset+k)=ialset(istartset(j)+k-1)
c     enddo
c     nalset=nalset+m
c     exit
c     endif
c     enddo
ccc   to remove end
              j=0
              call cident81(set,noelset,nset,id)
              if(id.gt.0) then
                if(set(id).eq.noelset) then
                  if(id.ne.iset) then
                    m=iendset(id)-istartset(id)+1
                    do k=1,m
                      ialset(nalset+k)=ialset(istartset(id)+k-1)
                    enddo
                    nalset=nalset+m
                  endif
                  j=id
                endif
              endif
              if(j.eq.0) then
                noelset(ipos:ipos)=' '
                if(kode.eq.0) then
                  write(*,*) 
     &                 '*ERROR reading *NSET/ELSET: node set ',
     &                 noelset
                else
                  write(*,*) 
     &                 '*ERROR reading *NSET/ELSET: element set ',
     &                 noelset
                endif
                write(*,*) '       has not been defined yet'
                ier=1
                return
              endif
            else
!     
!     node or element number
!     
              if(kode.eq.0) then
                if(ialset(nalset+1).gt.nk) then
                  write(*,*) 
     &                 '*WARNING reading *NSET/ELSET: value ',
     &                 ialset(nalset+1)
                  write(*,*) '         in set ',
     &                 set(iset)(1:index(set(iset),' ')-2),' > nk'
                else
                  nalset=nalset+1
                endif
              else
                if(ialset(nalset+1).gt.ne) then
                  write(*,*) 
     &                 '*WARNING reading *NSET/ELSET: value ',
     &                 ialset(nalset+1)
                  write(*,*) '         in set ',
     &                 set(iset)(1:index(set(iset),' ')-2),' > ne;'
                  write(*,*) '         This is only allowed for'
                  write(*,*) 
     &                 '         global elsets in combination'
                  write(*,*) '         with submodels'
c     else
c     nalset=nalset+1
                endif
                nalset=nalset+1
              endif
            endif
          enddo
          iendset(iset)=nalset
        endif
      enddo
!     
      return
      end
