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
      subroutine specificgasconstants(inpc,textpart,shcon,nshcon,
     &  nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *SPECIFIC GAS CONSTANT
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nshcon(*),nmat,ntmat_,istep,istat,n,ipoinpc(0:*),
     &  key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),i,ier
!
      real*8 shcon(0:3,ntmat_,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &   '*ERROR reading *SPECIFIC GAS CONSTANT: *SPECIFIC GAS CONSTANT'
         write(*,*) '  should be placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &   '*ERROR reading *SPECIFIC GAS CONSTANT: *SPECIFIC GAS CONSTANT'
         write(*,*) '  should be preceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &     '*WARNING in specificgasconstants: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*SPECIFIC GAS CONSTANT%")
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &        shcon(3,1,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SPECIFIC GAS CONSTANT%",ier)
            return
         endif
      enddo
!
      return
      end

