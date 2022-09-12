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
      subroutine depvars(inpc,textpart,nelcon,nmat,
     &        nstate_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &        ncocon,ipoinpc,ier)
!
!     reading the input deck: *DEPVAR
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,istep,nstate_,ncocon(2,*),ipoinpc(0:*),
     &  n,key,istat,nstate,irstrt(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),i,ier
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *DEPVAR: *DEPVAR should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *DEPVAR: *DEPVAR should be preceded'
         write(*,*) '  by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *DEPVAR: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*DEPVAR%")
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *DEPVAR: incomplete definition'
         ier=1
         return
      endif
      read(textpart(1)(1:10),'(i10)',iostat=istat) nstate
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*DEPVAR%",ier)
         return
      endif
      nstate_=max(nstate_,nstate)
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

