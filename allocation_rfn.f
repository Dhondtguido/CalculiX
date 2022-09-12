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
      subroutine allocation_rfn(nk_,ne_,nkon_,ipoinp,ipoinpc,inpc,inp)
!     
!     calculates a conservative estimate of the size of to be allocated
!     
      implicit none
!
      character*1 inpc(*)
      character*8 label
      character*132 textpart(16)
!     
      integer nk_,ne_,nkon_,ipoinp(2,*),ipoinpc(0:*),inp(3,*),ier,i,
     &     nteller,nopeexp,nope,nentries,n,key,istat,ipol,inl,iline
!     
      parameter(nentries=18)
!     
      ier=0
!     
!     initialisation of ipoinp
!     
      do i=1,nentries
        if(ipoinp(1,i).ne.0) then
          ipol=i
          inl=ipoinp(1,i)
          iline=inp(1,inl)-1
          exit
        endif
      enddo
!     
      istat=0
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      loop: do
        if(istat.lt.0) then
          exit
        endif
!     
        if(textpart(1)(1:8).eq.'*ELEMENT') then
!     
          loop1: do i=2,n
            if(textpart(i)(1:5).eq.'TYPE=') then
              read(textpart(i)(6:13),'(a8)') label
              if(label.eq.'        ') then
                write(*,*) 
     &               '*ERROR in allocation: element type is lacking'
                write(*,*) '       '
                call inputerror(inpc,ipoinpc,iline,
     &               "*ELEMENT or *ELEMENT OUTPUT%",ier)
                exit
              endif
!     
              nopeexp=0
!     
              if(label.eq.'C3D10   ') then
                nope=10
                nopeexp=10
              elseif(label.eq.'C3D4    ') then
                nope=4
                nopeexp=4
              endif
            endif
          enddo loop1
!     
          loop2:do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(1)(1:10),'(i10)',iostat=istat) i
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*ELEMENT%",ier)
              exit
            endif
            nteller=n-1
            if(nteller.lt.nope) then
              do
                call getnewline(inpc,textpart,istat,n,key,iline,
     &               ipol,inl,ipoinp,inp,ipoinpc)
                if((istat.lt.0).or.(key.eq.1)) exit loop2
                if(nteller+n.gt.nope) n=nope-nteller
                nteller=nteller+n
                if(nteller.eq.nope) exit
              enddo
            endif
            ne_=max(ne_,i)
            nkon_=nkon_+nopeexp
          enddo loop2
        elseif(textpart(1)(1:5).eq.'*NODE') then
!     
          do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(1)(1:10),'(i10)',iostat=istat) i
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*NODE%",ier)
              exit
            endif
            nk_=max(nk_,i)
          enddo
        else
!     
          call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
        endif
      enddo loop
!     
      if(ier.eq.1) then
        write(*,*) '*ERROR in allocation: at least one fatal'
        write(*,*) '       error message while reading the'
        write(*,*) '       input deck: CalculiX stops.'
        write(*,*)
        call exit(201)
      endif
!     
      return
      end
