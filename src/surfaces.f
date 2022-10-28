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
      subroutine surfaces(inpc,textpart,set,istartset,iendset,ialset,
     &     nset,nset_,nalset,nalset_,nk,ne,istep,istat,n,iline,ipol,
     &     inl,ipoinp,inp,lakon,ipoinpc,ier)
!     
!     reading the input deck: *SURFACE
!     
      implicit none
!     
      character*1 type,inpc(*)
      character*8 lakon(*)
      character*20 label,newlabel
      character*81 set(*),noset,elset,noelset
      character*132 textpart(16)
!     
      integer nset,nset_,nalset,nalset_,istep,istat,n,key,i,nk,ne,
     &     j,istartset(*),iendset(*),ialset(*),ipos,iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),iside,l,k,kstart,kend,ipoinpc(0:*),
     &     iset,nn,kincrement,ier,id
!     
      if(istep.gt.0) then
        write(*,*) '*ERROR reading *SURFACE: *SURFACE should be placed'
        write(*,*) '       before all step definitions'
        ier=1
        return
      endif
!     
      kstart=0
      kend=0
!     
      type='T'
!     
      do i=2,n
        if(textpart(i)(1:5).eq.'NAME=')
     &       then
          noelset(1:80)=textpart(i)(6:85)
          noelset(81:81)=' '
          if(textpart(i)(86:86).ne.' ') then
            write(*,*) 
     &           '*ERROR reading *SURFACE: surface name too long'
            write(*,*) '       (more than 80 characters)'
            write(*,*) '       surface name:',textpart(i)(1:132)
            ier=1
            return
          endif
        elseif(textpart(i)(1:5).eq.'TYPE=') then
          if(textpart(i)(6:12).eq.'ELEMENT') then
            type='T'
          elseif(textpart(i)(6:9).eq.'NODE') then
            type='S'
          else
            write(*,*) 
     &           '*ERROR reading *SURFACE: unknown surface type'
            ier=1
            return
          endif
        else
          write(*,*) 
     &         '*WARNING reading *SURFACE: parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*SURFACE%")
        endif
      enddo
!     
      ipos=index(noelset,' ')
      if(ipos.eq.1) then
        write(*,*) '*ERROR reading *SURFACE: no name specified'
        ier=1
        return
      endif
      noelset(ipos:ipos)=type
!     
!     check whether new set or old set (a *SURFACE can be used to
!     extend an already existing surface)
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
c            exit
c          else
!     
!     rearranging set information towards the end
!     
            nn=iendset(iset)-istartset(iset)+1
            if(nalset+nn.gt.nalset_) then
              write(*,*)'*ERROR reading *SURFACE: increase nalset_'
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
          write(*,*) '*ERROR reading *SURFACE: increase nset_'
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
      if(type.eq.'S') then
!     
!     node surface
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
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
          if(n.gt.1) then
            write(*,*) '*ERROR reading *SURFACE: only one entry per'
            write(*,*) '       line allowed'
            call inputerror(inpc,ipoinpc,iline,
     &           "*SURFACE%",ier)
            return
          endif
!     
          if(nalset+1.gt.nalset_) then
            write(*,*) '*ERROR reading *SURFACE: increase nalset_'
            ier=1
            return
          endif
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat)ialset(nalset+1)
          if(istat.gt.0) then
            noset=textpart(1)(1:80)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
ccc   to remove start
c     do i=1,nset
c     if(set(i).eq.noset) then
ccc   to remove end
            i=0
            call cident81(set,noset,nset,id)
            if(id.gt.0) then
              if(set(id).eq.noset) then
                i=id
              endif
            endif
            if(i.gt.0) then
              do j=istartset(i),iendset(i)
                if(ialset(j).gt.0) then
                  nalset=nalset+1
                  if(nalset.gt.nalset_) then
                    write(*,*) 
     &                   '*ERROR reading *SURFACE: increase nalset_'
                    ier=1
                    return
                  endif
                  ialset(nalset)=ialset(j)
                else
                  kstart=ialset(nalset-1)
                  kend=ialset(nalset)
                  nalset=nalset-1
                  kincrement=-ialset(j)
                  do k=kstart+kincrement,kend,kincrement
                    nalset=nalset+1
                    if(nalset.gt.nalset_) then
                      write(*,*) 
     &                     '*ERROR reading *SURFACE: increase nalset_'
                      ier=1
                      return
                    endif
                    ialset(nalset)=k
                  enddo
                endif
              enddo
              iendset(iset)=nalset
c              exit
ccc   to remove start
c     endif
c     enddo
c     if(i.gt.nset) then
ccc   to remove end
            else
              noset(ipos:ipos)=' '
              write(*,*) '*ERROR reading *SURFACE: node set ',noset
              write(*,*) '       does not exist'
              ier=1
              return
            endif
          else
            if(ialset(nalset+1).gt.nk) then
              write(*,*) '*WARNING reading *SURFACE: value ',
     &             ialset(nalset+1)
              write(*,*) '         in set ',set(iset),' > nk'
            else
              nalset=nalset+1
              iendset(iset)=nalset
            endif
          endif
        enddo
!     
      else
!     
!     element surface
!     
        do
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
          if((istat.lt.0).or.(key.eq.1)) then
            if(iendset(nset).eq.0) then
              nset=nset-1
            endif
            return
          endif
          if(nalset+1.gt.nalset_) then
            write(*,*) '*ERROR reading *SURFACE: increase nalset_'
            ier=1
            return
          endif
!     
          read(textpart(2)(1:20),'(a20)',iostat=istat) label
!     
          if(label(2:4).eq.'NEG') then
            label(2:4)='1  '
          elseif(label(2:4).eq.'POS') then
            label(2:4)='2  '
          endif
!     
!     for plane stress elements: 'N' and 'P' are converted
!     into '5' and '6' and farther down in '1' and '2'
!     
          if(label(2:2).eq.'N') then
            label(2:2)='5'
          elseif(label(2:2).eq.'P') then
            label(2:2)='6'
          endif
!     
          if((label(1:2).ne.'S1').and.(label(1:2).ne.'S2').and.
     &         (label(1:2).ne.'S3').and.(label(1:2).ne.'S4').and.
     &         (label(1:2).ne.'S5').and.(label(1:2).ne.'S6')) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SURFACE%",ier)
            return
          endif
!     
          read(textpart(1)(1:10),'(i10)',iostat=istat)l
          if(istat.gt.0) then
            elset=textpart(1)(1:80)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
ccc   to remove start
c     do i=1,nset
c     if(set(i).eq.elset) then
ccc   to remove end
            i=0
            call cident81(set,elset,nset,id)
            if(id.gt.0) then
              if(set(id).eq.elset) then
                i=id
              endif
            endif
            if(i.gt.0) then
              do j=istartset(i),iendset(i)
                l=ialset(j)
                if(l.gt.0) then
                  kstart=kend
                  kend=l
                  nalset=nalset+1
                  if(nalset.gt.nalset_) then
                    write(*,*) 
     &                   '*ERROR reading *SURFACE: increase nalset_'
                    ier=1
                    return
                  endif
                  newlabel=label
                  if((lakon(l)(1:2).eq.'CP').or.
     &                 (lakon(l)(2:2).eq.'A')) then
                    if(label(1:2).eq.'S1') then
                      newlabel(1:2)='S3'
                    elseif(label(1:2).eq.'S2') then
                      newlabel(1:2)='S4'
                    elseif(label(1:2).eq.'S3') then
                      newlabel(1:2)='S5'
                    elseif(label(1:2).eq.'S4') then
                      newlabel(1:2)='S6'
                    elseif(label(1:2).eq.'S5') then
                      newlabel(1:2)='S1'
                    elseif(label(1:2).eq.'S6') then
                      newlabel(1:2)='S2'
                    endif
                  endif
                  read(newlabel(2:2),'(i1)',iostat=istat) iside
                  ialset(nalset)=iside+10*l
                else
                  kstart=kstart
                  nalset=nalset-1
                  kincrement=-ialset(j)
                  do l=kstart+kincrement,kend,kincrement
                    nalset=nalset+1
                    if(nalset.gt.nalset_) then
                      write(*,*) 
     &                     '*ERROR reading *SURFACE: increase nalset_'
                      ier=1
                      return
                    endif
                    newlabel=label
                    if((lakon(l)(1:2).eq.'CP').or.
     &                   (lakon(l)(2:2).eq.'A')) then
                      if(label(1:2).eq.'S1') then
                        newlabel(1:2)='S3'
                      elseif(label(1:2).eq.'S2') then
                        newlabel(1:2)='S4'
                      elseif(label(1:2).eq.'S3') then
                        newlabel(1:2)='S5'
                      elseif(label(1:2).eq.'S4') then
                        newlabel(1:2)='S6'
                      elseif(label(1:2).eq.'S5') then
                        newlabel(1:2)='S1'
                      elseif(label(1:2).eq.'S6') then
                        newlabel(1:2)='S2'
                      endif
                    endif
                    read(newlabel(2:2),'(i1)',iostat=istat) 
     &                   iside
                    ialset(nalset)=iside+10*l
                  enddo
                endif
              enddo
              iendset(iset)=nalset
c              exit
ccc   to remove start
c     endif
c     enddo
c     if(i.gt.nset) then
ccc   to remove end
            else
              elset(ipos:ipos)=' '
              write(*,*) '*ERROR reading *SURFACE: element set ',
     &             elset
              write(*,*) '       does not exist'
              ier=1
              return
            endif
          else
            if(l.gt.ne) then
              write(*,*) '*WARNING reading *SURFACE: element ',
     &             l
              write(*,*) '         in set ',set(iset),' > ne'
            else
              newlabel=label
              if((lakon(l)(1:2).eq.'CP').or.
     &             (lakon(l)(2:2).eq.'A')) then
                if(label(1:2).eq.'S1') then
                  newlabel(1:2)='S3'
                elseif(label(1:2).eq.'S2') then
                  newlabel(1:2)='S4'
                elseif(label(1:2).eq.'S3') then
                  newlabel(1:2)='S5'
                elseif(label(1:2).eq.'S4') then
                  newlabel(1:2)='S6'
                elseif(label(1:2).eq.'S5') then
                  newlabel(1:2)='S1'
                elseif(label(1:2).eq.'S6') then
                  newlabel(1:2)='S2'
                endif
              endif
              read(newlabel(2:2),'(i1)',iostat=istat) iside
              nalset=nalset+1
              ialset(nalset)=iside+10*l
              iendset(iset)=nalset
            endif
          endif
        enddo
      endif
!     
      return
      end

