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
      subroutine changesurfacebehaviors(inpc,textpart,matname,nmat,
     &  nmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  imat,ier)
!
!     reading the input deck: *CHANGE SURFACE BEHAVIOR
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*),interactionname
      character*132 textpart(16)
!
      integer nmat,nmat_,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),imat,ier
!
      if(istep.eq.0) then
         write(*,*) '*ERROR reading *CHANGE SURFACE BEHAVIOR:'
         write(*,*) '       *CHANGE SURFACE BEHAVIOR'
         write(*,*) '       cannot be used before the first step'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:12).eq.'INTERACTION=') then
            interactionname=textpart(i)(13:92)
         else
            write(*,*) '*WARNING reading *CHANGE SURFACE BEHAVIOR:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CHANGE SURFACE BEHAVIOR%")
         endif
      enddo
!
!     check whether the interaction exists
!
      imat=0
      do i=1,nmat
         if(matname(i).eq.interactionname) then
            imat=i
            exit
         endif
      enddo
!
      if(imat.eq.0) then
         write(*,*) '*ERROR reading *CHANGE SURFACE BEHAVIOR:',
     &           interactionname
         write(*,*) '       is a nonexistent interaction'
         ier=1
         return
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

