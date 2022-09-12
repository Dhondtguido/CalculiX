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
      subroutine readview(ntr,adview,auview,fenv,nzsrad,ithermal,
     &  jobnamef)
!     
!     reading the viewfactors from file
!
      implicit none
!     
      logical exi
!
      character*80 versionvwf,version
      character*132 jobnamef(*),fnvw
!     
      integer ntr,nzsrad,ithermal(*),i,k,length,lengthvwf
!
      real*8 adview(*),auview(*),fenv(*)
!     
      if(ithermal(1).eq.3) then
         write(*,*) '*WARNING in readview: viewfactors are being'
         write(*,*) '         read from file for a thermomechani-'
         write(*,*) '         cal calculation: they will not be '
         write(*,*) '         recalculated in every iteration.'
      endif
!     
      write(*,*) 'Reading the viewfactors from file'
      write(*,*)
!     
      if(jobnamef(2)(1:1).eq.' ') then
         do i=1,132
            if(jobnamef(1)(i:i).eq.' ') exit
         enddo
         i=i-1
         fnvw=jobnamef(1)(1:i)//'.vwf'
      else
         fnvw=jobnamef(2)
      endif
      inquire(file=fnvw,exist=exi)
      if(exi) then
         open(10,file=fnvw,status='old',form='unformatted',
     &        access='sequential',err=10)
      else
         write(*,*) '*ERROR in readview: viewfactor file ',fnvw
         write(*,*) 'does not exist'
         call exit(201)
      endif
!     
      read(10) versionvwf
      read(10) (adview(k),k=1,ntr)
      read(10) (auview(k),k=1,2*nzsrad)
      read(10)(fenv(k),k=1,ntr)
!     
      close(10)
!     
!     check whether the CalculiX version in the viewfactor file corresponds
!     with the actual version
!
!     string includes the last nonblank character
!
      do i=80,1,-1
         if(versionvwf(i:i).ne.' ') then
            lengthvwf=i
            exit
         endif
      enddo
!
!     string stops before the first "p" or "_" character (from
!     patch or _i8)    
!
      do i=1,lengthvwf
         if((versionvwf(i:i).eq.'p').or.
     &      (versionvwf(i:i).eq.'_')) then
            lengthvwf=i-1
            exit
         endif
      enddo
!      
!     string includes the last nonblank character
!
      call getversion(version)
      do i=80,1,-1
         if(version(i:i).ne.' ') then
            length=i
            exit
         endif
      enddo
!
!     string stops before the first "p" or "_" character (from
!     patch or _i8)    
!
      do i=1,length
         if((version(i:i).eq.'p').or.
     &      (version(i:i).eq.'_')) then
            length=i-1
            exit
         endif
      enddo
!
      if(versionvwf(1:lengthvwf).ne.version(1:length)) then
         write(*,*) '*ERROR in readview: CalculiX ',
     &      versionvwf(1:lengthvwf)
         write(*,*) '       in viewfactor file ',
     %           fnvw(1:index(fnvw,' ')-1),' does not'
         write(*,*) '       correspond to the actual CalculiX ',
     &         version(1:length)
         write(*,*)
         call exit(201)
      endif
!
      return
!
 10   write(*,*) '*ERROR in readview: could not open file ',fnvw
      call exit(201)
      end
      
