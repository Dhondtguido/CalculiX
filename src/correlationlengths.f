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
      subroutine correlationlengths(inpc,textpart,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,physcon,ier)        
!     
!     reading the input deck: *CORRELATION LENGTH
!     
      implicit none
!     
      character*1 inpc(*)
      character*132 textpart(16)
!     
      integer istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),ipoinpc(0:*),ier
!     
      real*8 physcon(*),corrlen
!
!
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *CORRELATION LENGTH:'
        write(*,*) '       *CORRELATION LENGTH can only'
        write(*,*) '       only be used within a'
        write(*,*) '       *ROBUST DESIGN step'
        call inputerror(inpc,ipoinpc,iline,
     &       "*CORRELATION LENGTH%",ier)
        return
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
!     reading in the correlation length
!     
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) corrlen
      if(istat.gt.0) then
        call inputerror(inpc,ipoinpc,iline,
     &       "*CORRELATION LENGTH%",ier)
        return
      endif
      if(corrlen.lt.0.d0) then
        write(*,*) '*ERROR reading *CORRELATION LENGTH'
        write(*,*) '      Correlation length for computation'
        write(*,*) '      of the random fields cannot'
        write(*,*) '      be less than 0' 
        write(*,*)  
        call inputerror(inpc,ipoinpc,iline,
     &       "*CORRELATION LENGTH%",ier)
        return
      endif
      physcon(12)=corrlen
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end
      
      
