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
      subroutine clearances(inpc,textpart,tieset,istat,n,iline,
     &           ipol,inl,ipoinp,inp,ntie,ipoinpc,istep,tietol,irstrt,
     &           ier)
!
!     reading the input deck: *CLEARANCE
!
      implicit none
!
      logical contactpair
!
      character*1 inpc(*)
      character*81 tieset(3,*)
      character*132 textpart(16),master,slave
!
      integer istat,n,i,j,key,ipos,iline,ipol,inl,ipoinp(2,*),irstrt(*),
     &  inp(3,*),ntie,ipoinpc(0:*),iposslave,iposmaster,itie,istep,
     &  ier
!
      real*8 tietol(4,*),value
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *CLEARANCE: *CLEARANCE should be'
         write(*,*) '       placed before all step definitions'
         ier=1
         return
      endif
!
      contactpair=.false.
!
      do i=2,n
         if(textpart(i)(1:7).eq.'MASTER=') then
            iposmaster=index(textpart(i),' ')
            master(1:iposmaster-8)=textpart(i)(8:iposmaster-1)
            do j=iposmaster-7,132
               master(j:j)=' '
            enddo
         elseif(textpart(i)(1:6).eq.'SLAVE=') then
            iposslave=index(textpart(i),' ')
            slave(1:iposslave-7)=textpart(i)(7:iposslave-1)
            do j=iposslave-6,132
               slave(j:j)=' '
            enddo
         elseif(textpart(i)(1:6).eq.'VALUE=') then
            read(textpart(i)(7:26),'(f20.0)',iostat=istat) value
         else
            write(*,*) 
     &       '*WARNING reading *CLEARANCE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CLEARANCE%")
         endif
      enddo
!
!     selecting the appropriate action
!
      iposslave=index(slave(1:80),' ')
      iposmaster=index(master(1:80),' ')
      do i=1,ntie
         if((tieset(1,i)(81:81).ne.'C').and.
     &      (tieset(1,i)(81:81).ne.'-')) cycle
         ipos=index(tieset(2,i),' ')-1
         if(ipos.ne.iposslave) cycle
         if(tieset(2,i)(1:ipos-1).ne.slave(1:ipos-1)) cycle
         ipos=index(tieset(3,i),' ')-1
         if(ipos.ne.iposmaster) cycle
         if(tieset(3,i)(1:ipos-1).ne.master(1:ipos-1)) cycle
         itie=i
         exit
      enddo
!
      if(i.gt.ntie) then
         write(*,*) '*ERROR reading *CLEARANCE: no such contact pair'
         call inputerror(inpc,ipoinpc,iline,
     &        "*CLEARANCE%",ier)
         return
      endif
      tietol(3,i)=value
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



