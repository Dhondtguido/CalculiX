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
      subroutine randomfields(inpc,textpart,istep,istat,n,iline,
     &        ipol,inl,ipoinp,inp,ipoinpc,nener,physcon,ier)        
!
!     reading the input deck: *RANDOM FIELD
!
!     characterized by standarddeviation and correlation length
!       
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*),nener,k,ipos,neigenvectors,ier
!
      real*8 physcon(*),dummy
!      
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *RANDOM FIELD: *RANDOM FIELD can'
         write(*,*) '       only be used within a *SENSITIVITY step'
         ier=1
         return
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!  
!     reading the standard deviation and the correlation length
!    
!     Number of eigenvectors used for the creation of the random field
!
      read(textpart(1)(1:20),'(i10)',iostat=istat) neigenvectors
      physcon(11)=1.d0*neigenvectors
!
!     Standard deviation
!
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) physcon(12)
!
!     Correlation length
!
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) physcon(13)
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end
      
      
