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
      subroutine feasibledirections(inpc,textpart,istat,n,key,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,nmethod,objectset,nobject,istep,
     &     ier)         
!     
!     reading the input deck: *FEASIBLE DIRECTION
!     
      implicit none
!     
      character*1 inpc(*)
      character*81 objectset(5,*)
      character*132 textpart(16)
!     
      integer nmethod,nobject,i,n,key,istat,istep,iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *FEASIBLE DIRECTION:'
        write(*,*) '       *FEASIBLE DIRECTION can only be used'
        write(*,*) '       within a STEP'    
        ier=1
        return
      endif
!     
      nmethod=16 
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)     
!     
      return
      end
      
      
