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
      subroutine initialstrainincreases(inpc,textpart,prestr,iprestr,
     &  mi,istep,istat,n,iline,ipol,inl,ipoinp,inp,ne,ipoinpc,ier)
!
!     reading the input deck: *INITIAL STRAIN INCREASE
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer iprestr,istep,istat,n,i,j,k,l,key,mi(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ne,ipoinpc(0:*),ier
!
      real*8 beta(8),prestr(6,mi(1),*)
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *INITIAL STRAIN INCREASE:'
         write(*,*) 
     &     '       *INITIAL STRAIN INCREASE should only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
      endif
!
      if(iprestr.ne.2) then
         write(*,*) '*ERROR reading *INITIAL STRAIN INCREASE:'
         write(*,*) '       a strain increase is only allowed'
         write(*,*) '       in an input deck with either a'
         write(*,*) '       *INITIAL CONDITIONS,TYPE=PLASTIC STRAIN'
         write(*,*) '       card, a *MODEL CHANGE,ADD card or a '
         write(*,*) '       *MODEL CHANGE,ADD=STRAIN FREE card.'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) '*WARNING reading *INITIAL STRAIN INCREASE:'
         write(*,*) '         parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*INITIAL STRAIN INCREASE%")
      enddo
!
!     storing the increase of strain as negative initial strains
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         do j=1,6
            read(textpart(j+2)(1:20),'(f20.0)',iostat=istat) 
     &           beta(j)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*INITIAL STRAIN INCREASE%",ier)
               return
            endif
         enddo
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.ne.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*INITIAL STRAIN INCREASE%",ier)
            return
         endif
         if(l.gt.ne) then
            write(*,*) 
     &           '*WARNING reading *INITIAL STRAIN INCREASE: element ',l
            write(*,*)'          exceeds the largest defined ',
     &           'element number'
            cycle
         endif
         read(textpart(2)(1:10),'(i10)',iostat=istat) k
         if(istat.eq.0) then
            do j=1,6
               prestr(j,k,l)=prestr(j,k,l)+beta(j)
            enddo
         else
            call inputerror(inpc,ipoinpc,iline,
     &           "*INITIAL STRAIN INCREASE%",ier)
            return
         endif
      enddo
!     
      return
      end
      
