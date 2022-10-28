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
      subroutine frictions(inpc,textpart,elcon,nelcon,
     &  imat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,nstate_,ichangefriction,mortar,ier)
!
!     reading the input deck: *FRICTION
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),imat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ncmat_,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  nstate_,ichangefriction,mortar,ier
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0).and.
     &   (ichangefriction.eq.0)) then
         write(*,*) '*ERROR reading *FRICTION:'
         write(*,*) '       *FRICTION should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(imat.eq.0) then
         write(*,*) '*ERROR reading *FRICTION:'
         write(*,*) '       *FRICTION should be preceded'
         write(*,*) '       by a *SURFACE INTERACTION card'
         ier=1
         return
      endif
!
      nstate_=max(nstate_,9)
!
c      if(nelcon(1,imat).gt.0) nelcon(1,imat)=max(nelcon(1,imat),7)
!
!     "8" is for Mortar contact
!
      if(nelcon(1,imat).ne.-51) nelcon(1,imat)=max(nelcon(1,imat),8)
      nelcon(2,imat)=1
!
!     no temperature dependence allowed; last line is decisive
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
c         do i=1,3
         do i=1,2
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(5+i,1,imat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*FRICTION%",ier)
               return
            endif
         enddo
         if(elcon(6,1,imat).le.0.d0) then
            write(*,*) '*ERROR reading *FRICTION: friction coefficient'
            write(*,*) '       must be strictly positive'
            call inputerror(inpc,ipoinpc,iline,
     &           "*FRICTION%",ier)
            return
         endif
         if(elcon(7,1,imat).le.0.d0) then
            write(*,*) '*WARNING reading *FRICTION: stick slope'
            write(*,*) '         must be strictly positive'
            write(*,*) 
     &       '         the following default will be used:',
     &                  elcon(1,1,1)/2.d0
            write(*,*) 
     &       '         the user is advised to analyze the results' 
            write(*,*) 
     &       '         carefully and, if possible, to come up with' 
            write(*,*) 
     &       '         a experimentally based stick slope' 
            call inputwarning(inpc,ipoinpc,iline,
     &"*FRICTION%")
            elcon(7,1,imat)=elcon(1,1,1)/2.d0
         endif
         elcon(0,1,imat)=0.d0
      enddo
!     
      return
      end

