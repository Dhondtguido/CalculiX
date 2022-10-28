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
      subroutine noanalysiss(inpc,textpart,nmethod,iperturb,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,tper,ier)
!
!     reading the input deck: *NO ANALYSIS
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmethod,iperturb(*),istep,istat,n,key,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!
      real*8 tper
!
      if(istep.lt.1) then
         write(*,*)
     &      '*ERROR reading *NO ANALYSIS: *NO ANALYSIS can only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
      endif
!
      write(*,*) '*WARNING: no analysis option was chosen'
!
      nmethod=0
      iperturb(1)=0
      tper=1.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

