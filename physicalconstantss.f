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
      subroutine physicalconstantss(inpc,textpart,physcon,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *PHYSICAL CONSTANTS
!
!     physcon(1): absolute zero temperature
!            (2): Stefan-Boltzmann constant
!            (3): Newton gravity constant
!            (4): Static temperature at infinity (CFD)
!            (5): Norm of the velocity vector at infinity (CFD)
!            (6): Static pressure at infinity (CFD)
!            (7): Density at infinity (CFD)
!            (8): Length of the computational domain (CFD)
!            (9): perturbation flag (CFD)
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer i,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*),ier
!
      real*8 physcon(*)
!
      if(istep.gt.0) then
         write(*,*) 
     &      '*ERROR reading *PHYSICAL CONSTANTS: *PHYSICAL CONSTANTS'
         write(*,*) '        should only be used before the first STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:13).eq.'ABSOLUTEZERO=') then
            read(textpart(i)(14:33),'(f20.0)',iostat=istat) physcon(1)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*PHYSICAL CONSTANTS%",ier)
               return
            endif
         elseif(textpart(i)(1:16).eq.'STEFANBOLTZMANN=') then
            read(textpart(i)(17:36),'(f20.0)',iostat=istat) physcon(2)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*PHYSICAL CONSTANTS%",ier)
               return
            endif
         elseif(textpart(i)(1:14).eq.'NEWTONGRAVITY=') then
            read(textpart(i)(15:24),'(f20.0)',iostat=istat) physcon(3)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*PHYSICAL CONSTANTS%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING in physicalconstants: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*PHYSICAL CONSTANTS%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end







