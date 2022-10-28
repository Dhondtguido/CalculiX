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
      subroutine spcmatch(xboun,nodeboun,ndirboun,nboun,xbounold,
     &   nodebounold,ndirbounold,nbounold,ikboun,ilboun,vold,reorder,
     &   nreorder,mi,typeboun)
!
!     matches SPC boundary conditions of one step with those of
!     the previous step
!
      implicit none
!
      character*1 typeboun(*)
!
      integer nodeboun(*),ndirboun(*),nboun,nodebounold(*),ilboun(*),
     &  ndirbounold(*),nbounold,i,kflag,idof,id,nreorder(*),ikboun(*),
     &  mi(*)
!
      real*8 xboun(*),xbounold(*),vold(0:mi(2),*),reorder(*)
!
      kflag=2
!
      do i=1,nboun
         nreorder(i)=0
      enddo
!
!     determining for which new SPC there exists an old value
!     for new SPC i the value of reorder(i) is set to the
!     old value and nreorder(i) is set to 1 (default for
!     both fields is 0)
!
      do i=1,nbounold
         idof=8*(nodebounold(i)-1)+ndirbounold(i)
         if(nboun.gt.0) then
            call nident(ikboun,idof,nboun,id)
         else
            id=0
         endif
         if((id.gt.0).and.(ikboun(id).eq.idof)) then
            reorder(ilboun(id))=xbounold(i)
            nreorder(ilboun(id))=1
         endif
      enddo
!
!     for new SPC's which do not have an equivalent SPC
!     the last value of the corresponding field is stored
!     in reorder (not for fluid SPC's since the facial
!     values of the field are not stored in the restart
!     file)
!
      do i=1,nboun
         if(typeboun(i).eq.'F') cycle
         if(nreorder(i).eq.0) then
            if(ndirboun(i).gt.4) then
               reorder(i)=0.d0
            else
               reorder(i)=vold(ndirboun(i),nodeboun(i))
            endif
         endif
      enddo
!
      do i=1,nboun
         if(typeboun(i).eq.'F') cycle
         xbounold(i)=reorder(i)
      enddo
!
      return
      end

