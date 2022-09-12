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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine writeview(ntr,adview,auview,fenv,nzsrad,
     &  jobnamef)
!     
!     writing the viewfactors to file
!
      implicit none
!
      character*80 version
      character*132 jobnamef(*),fnvw
!     
      integer ntr,nzsrad,i,k
!
      real*8 adview(*),auview(*),fenv(*)
!     
      write(*,*) 'Writing the viewfactors to file'
      write(*,*)
!     
      if(jobnamef(3)(1:1).eq.' ') then
         do i=1,132
            if(jobnamef(1)(i:i).eq.' ') exit
         enddo
         i=i-1
         fnvw=jobnamef(1)(1:i)//'.vwf'
      else
         fnvw=jobnamef(3)
      endif
      open(10,file=fnvw,status='unknown',form='unformatted',
     &     access='sequential',err=10)
!
      call getversion(version)
!     
      write(10) version
      write(10) (adview(k),k=1,ntr)
      write(10) (auview(k),k=1,2*nzsrad)
      write(10)(fenv(k),k=1,ntr)
      close(10)
!     
      return
!
 10   write(*,*) '*ERROR in writeview: could not open file ',fnvw
      call exit(201)
      end
      
