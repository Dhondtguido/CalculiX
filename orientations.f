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
      subroutine orientations(inpc,textpart,orname,orab,norien,
     &  norien_,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *ORIENTATION
!
      implicit none
!
      logical distribution
!
      character*1 inpc(*)
      character*80 orname(*),distname,oriename
      character*132 textpart(16)
!
      integer norien,norien_,istep,istat,n,key,i,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),iaxis,j,ier,idis,iorien
!
      real*8 orab(7,*),a(3,3),c(3,3),angle,p(3),dc,ds,pi,system,ab(6)
!
      if(istep.gt.0) then
         write(*,*) 
     &       '*ERROR reading *ORIENTATION: *ORIENTATION should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      distribution=.false.
!
!     rectangular coordinate system: orab(7,norien)=1
!     cylindrical coordinate system: orab(7,norien)=-1
!     default is rectangular
!
      system=1.d0
      iaxis=0
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            oriename=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR reading *ORIENTATION: name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       orientation name:',textpart(i)(1:132)
               ier=1
               return
            endif
         elseif(textpart(i)(1:7).eq.'SYSTEM=') then
            if(textpart(i)(8:8).eq.'C') then
               system=-1.d0
            endif
         else
            write(*,*) 
     &        '*WARNING reading *ORIENTATION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*ORIENTATION%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)
     &      '*ERROR reading *ORIENTATION: definition of the following'
         write(*,*) '  orientation is not complete: ',orname(norien)
         ier=1
         return
      endif
!
      do i=1,6
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) ab(i)
         if(istat.gt.0) then
            distribution=.true.
            read(textpart(1)(1:80),'(a80)',iostat=istat) distname
            exit
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      if((istat.ge.0).and.(key.ne.1)) then
         read(textpart(1)(1:10),'(i10)',iostat=istat) iaxis
         if(istat.gt.0) then
            write(*,*) '*ERROR reading *ORIENTATION: expected',
     &        'axis of rotation'            
            call inputerror(inpc,ipoinpc,iline,
     &           "*ORIENTATION%",ier)
            return
         endif
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) angle
         if(istat.gt.0) then
            write(*,*) '*ERROR reading *ORIENTATION: expected',
     &        'angle of rotation'            
            call inputerror(inpc,ipoinpc,iline,
     &           "*ORIENTATION%",ier)
            return
         endif
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!
      if(.not.distribution) then
!     
         norien=norien+1
         if(norien.gt.norien_) then
            write(*,*) '*ERROR reading *ORIENTATION: increase norien_'
            ier=1
            return
         endif
         do i=1,6
            orab(i,norien)=ab(i)
         enddo
         orab(7,norien)=system
         orname(norien)=oriename
         if(iaxis.eq.0) return
!
!        additional rotation about an angle only for rectangular
!        coordinate systems
!
         if(orab(7,norien).lt.0.d0) then
            write(*,*) '*ERROR reading *ORIENTATION'
            write(*,*) '       additional rotation about an angle'
            write(*,*) '       is only allowed for rectangular systems'
            ier=1
            return
         endif
!     
         call transformatrix(orab(1,norien),p,a)
!     
!        vector on the rotation axis
!     
         do i=1,3
            p(i)=a(i,iaxis)
         enddo
!     
!        rotation matrix
!     
         pi=4.d0*datan(1.d0)
         angle=angle*pi/180.d0
         dc=dcos(angle)
         ds=dsin(angle)
         c(1,1)=dc+(1.d0-dc)*p(1)*p(1)
         c(1,2)=-ds*p(3)+(1.d0-dc)*p(1)*p(2)
         c(1,3)=ds*p(2)+(1.d0-dc)*p(1)*p(3)
         c(2,1)=ds*p(3)+(1.d0-dc)*p(2)*p(1)
         c(2,2)=dc+(1.d0-dc)*p(2)*p(2)
         c(2,3)=-ds*p(1)+(1.d0-dc)*p(2)*p(3)
         c(3,1)=-ds*p(2)+(1.d0-dc)*p(3)*p(1)
         c(3,2)=ds*p(1)+(1.d0-dc)*p(3)*p(2)
         c(3,3)=dc+(1.d0-dc)*p(3)*p(3)
!     
!        rotate vector along the local x-axis and store 
!        as first point in orab
!     
         do i=1,3
            orab(i,norien)=0.d0
            do j=1,3
               orab(i,norien)=orab(i,norien)+c(i,j)*a(j,1)
            enddo
         enddo
!     
!        rotate vector along the local y-axis and store as 
!        second point in orab
!     
         do i=1,3
            orab(i+3,norien)=0.d0
            do j=1,3
               orab(i+3,norien)=orab(i+3,norien)+c(i,j)*a(j,2)
            enddo
         enddo
      else
!
!        distribution
!
         idis=0
         do iorien=1,norien
            if(orname(iorien).eq.distname) then
               idis=idis+1
               orab(7,iorien)=system
               orname(iorien)=oriename
               if(iaxis.ne.0) then
!
!     additional rotation about an angle only for rectangular
!     coordinate systems
!     
                  if(orab(7,iorien).lt.0.d0) then
                     write(*,*) '*ERROR reading *ORIENTATION'
                     write(*,*)
     &                '       additional rotation about an angle'
                     write(*,*) 
     &                '       is only allowed for rectangular systems'
                     ier=1
                     return
                  endif
!     
                  call transformatrix(orab(1,iorien),p,a)
!     
!     vector on the rotation axis
!     
                  do i=1,3
                     p(i)=a(i,iaxis)
                  enddo
!     
!     rotation matrix
!     
                  pi=4.d0*datan(1.d0)
                  angle=angle*pi/180.d0
                  dc=dcos(angle)
                  ds=dsin(angle)
                  c(1,1)=dc+(1.d0-dc)*p(1)*p(1)
                  c(1,2)=-ds*p(3)+(1.d0-dc)*p(1)*p(2)
                  c(1,3)=ds*p(2)+(1.d0-dc)*p(1)*p(3)
                  c(2,1)=ds*p(3)+(1.d0-dc)*p(2)*p(1)
                  c(2,2)=dc+(1.d0-dc)*p(2)*p(2)
                  c(2,3)=-ds*p(1)+(1.d0-dc)*p(2)*p(3)
                  c(3,1)=-ds*p(2)+(1.d0-dc)*p(3)*p(1)
                  c(3,2)=ds*p(1)+(1.d0-dc)*p(3)*p(2)
                  c(3,3)=dc+(1.d0-dc)*p(3)*p(3)
!     
!     rotate vector along the local x-axis and store 
!     as first point in orab
!     
                  do i=1,3
                     orab(i,iorien)=0.d0
                     do j=1,3
                        orab(i,iorien)=orab(i,iorien)+c(i,j)*a(j,1)
                     enddo
                  enddo
!     
!     rotate vector along the local y-axis and store as 
!     second point in orab
!     
                  do i=1,3
                     orab(i+3,iorien)=0.d0
                     do j=1,3
                        orab(i+3,iorien)=orab(i+3,iorien)+c(i,j)*a(j,2)
                     enddo
                  enddo
!     
               endif
            endif
         enddo
!
!        if no corresponding distribution is found: error
!
         if(idis.eq.0) then
            write(*,*) '*ERROR reading *ORIENTATION: distribution',
     &          distname,' has not been defined'
            call inputerror(inpc,ipoinpc,iline,
     &           "*ORIENTATION%",ier)
            return
         endif
       endif
c       write(*,*) (orab(i,norien),i=1,7)
!     
      return
      end
      
      
