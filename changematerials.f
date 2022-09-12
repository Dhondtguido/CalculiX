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
      subroutine changematerials(inpc,textpart,matname,nmat,nmat_,
     &  irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  imat,ier)
!
!     reading the input deck: *CHANGE MATERIAL
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*),materialname
      character*132 textpart(16)
!
      integer nmat,nmat_,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),imat,ier
!
      if(istep.eq.0) then
         write(*,*) '*ERROR reading *CHANGE MATERIAL: *CHANGE MATERIAL'
         write(*,*) '       cannot be used before the first step'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            materialname=textpart(i)(6:85)
         else
            write(*,*) 
     &    '*WARNING reading *CHANGE MATERIAL: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CHANGE MATERIAL%")
         endif
      enddo
!
!     check whether the material exists
!
      imat=0
      do i=1,nmat
         if(matname(i).eq.materialname) then
            imat=i
            exit
         endif
      enddo
!
      if(imat.eq.0) then
         write(*,*) '*ERROR reading *CHANGE MATERIAL:',materialname
         write(*,*) '       is a nonexistent material'
         ier=1
         return
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

