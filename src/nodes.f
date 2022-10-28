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
      subroutine nodes(inpc,textpart,co,nk,nk_,set,istartset,
     &     iendset,ialset,nset,nset_,nalset,nalset_,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
!     reading the input deck: *NODE
!     
      implicit none
!     
      character*1 inpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!     
      integer nk,nk_,nset,nset_,nalset,nalset_,istep,istat,n,key,id,
     &     i,j,js,k,nn,inoset,ipos,istartset(*),iendset(*),ialset(*),
     &     iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!     
      real*8 co(3,*)
!     
      inoset=0
!     
!     checking for set definition
!     
      loop: do i=2,n
      if(textpart(i)(1:5).eq.'NSET=') then
        noset=textpart(i)(6:85)
        if(textpart(i)(86:86).ne.' ') then
          write(*,*) '*ERROR reading *NODE: set name too long'
          write(*,*) '       (more than 80 characters)'
          write(*,*) '       set name:',textpart(i)(1:132)
          ier=1
          return
        endif
        noset(81:81)=' '
        ipos=index(noset,' ')
        noset(ipos:ipos)='N'
        inoset=1
c        do js=1,nset
c          if(set(js).eq.noset) then
        call cident81(set,noset,nset,js)
        if(js.gt.0) then
          if(noset.eq.set(js)) then
!     
!     existent set
!     
            if(iendset(js).eq.nalset) then
              cycle loop
            else
              nn=iendset(js)-istartset(js)+1
              if(nalset+nn.gt.nalset_) then
                write(*,*) 
     &               '*ERROR reading *NODE: increase nalset_'
                ier=1
                return
              endif
              do k=1,nn
                ialset(nalset+k)=ialset(istartset(js)+k-1)
              enddo
              if(nn.gt.0) then
                do k=istartset(js),nalset
                  ialset(k)=ialset(k+nn)
                enddo
                do k=1,nset
                  if(istartset(k).gt.iendset(js)) then
                    istartset(k)=istartset(k)-nn
                    iendset(k)=iendset(k)-nn
                  endif
                enddo
              endif
              istartset(js)=nalset-nn+1
              iendset(js)=nalset
              cycle loop
            endif
          endif
        endif
c        enddo
!     
!     new set
!     
cc        call cident81(set,noset,nset,id)
        nset=nset+1
        if(nset.gt.nset_) then
          write(*,*) '*ERROR reading *NODE: increase nset_'
          ier=1
          return
        endif
        do j=nset,js+2,-1
          istartset(j)=istartset(j-1)
          iendset(j)=iendset(j-1)
          set(j)=set(j-1)
        enddo
        js=js+1
ccc   to remove start        
c     js=nset
ccc   to remove end        
        istartset(js)=nalset+1
        iendset(js)=nalset
        set(js)=noset
        exit
      else
        write(*,*) 
     &       '*WARNING reading *NODE: parameter not recognized:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputwarning(inpc,ipoinpc,iline,
     &       "*NODE%")
      endif
      enddo loop
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
        read(textpart(1)(1:10),'(i10)',iostat=istat) i
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*NODE%",ier)
          return
        endif
        if(n.eq.1) then
          co(1,i)=0.d0
        else
          read(textpart(2)(1:20),'(f20.0)',iostat=istat) co(1,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        if(n.le.2) then
          co(2,i)=0.d0
        else
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) co(2,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        if(n.le.3) then
          co(3,i)=0.d0
        else
          read(textpart(4)(1:20),'(f20.0)',iostat=istat) co(3,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        nk=max(nk,i)
        if(nk.gt.nk_) then
          write(*,*) '*ERROR reading *NODE: increase nk_'
          ier=1
          return
        endif
!     
!     assigning node to set
!     
        if(inoset.eq.1) then
          if(nalset+1.gt.nalset_) then
            write(*,*) '*ERROR reading *NODE: increase nalset_'
            ier=1
            return
          endif
          nalset=nalset+1
          ialset(nalset)=i
          iendset(js)=nalset
        endif
!     
      enddo
!     
      return
      end
