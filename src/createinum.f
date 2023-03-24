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
      subroutine createinum(ipkon,inum,kon,lakon,nk,ne,cflag,nelemload,
     &  nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
     &  ielprop,prop)
!
!     determines inum in case no extrapolation is requested in the
!     input deck (e.g. only nodal variables are requested)
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,nope,nfield,mi(*),
     &  nk,i,j,nelemload(2,*),nload,node,nboun,nlayer,nopeexp,
     &  nodeboun(*),ndirboun(*),ithermal(*),ielmat(mi(3),*),
     &  ielprop(*)
!
      real*8 yn,co(3,*),vold(0:mi(2),*),prop(*)
!
      force=.false.
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') then
            nlayer=0
            do j=1,mi(3)
               if(ielmat(j,i).gt.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               nopeexp=28
            elseif(lakonl(4:5).eq.'15') then
               nopeexp=21
            endif
         endif
!
         if(lakonl(1:1).eq.'F') then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:1).eq.'E').and.
     &          ((lakon(i)(7:7).eq.'A').or.
     &           (lakon(i)(7:7).eq.'2'))) then
            inum(kon(indexe+1))=inum(kon(indexe+1))+1
            inum(kon(indexe+2))=inum(kon(indexe+2))+1
            cycle
         elseif(lakonl(1:7).eq.'ESPRNGF') then
            read(lakonl(8:8),'(i1)') nope
            nope=nope+1
            inum(kon(indexe+nope))=-1
            cycle
          elseif(lakonl(1:1).eq.'U') then
            nope=ichar(lakonl(8:8))
         else
            cycle
         endif
!
!        counting the number of elements a node belongs to
!
         if(lakonl(7:8).ne.'LC') then
            do j=1,nope
               inum(kon(indexe+j))=inum(kon(indexe+j))+1
            enddo
         else
            do j=1,nope*nlayer
               inum(kon(indexe+nopeexp+j))=inum(kon(indexe+nopeexp+j))+1
            enddo
         endif
c     Bernhardi start
c        incompatible modes elements
c         if(lakonl(1:5).eq.'C3D8I') then
c            do j=1,3
c               inum(kon(indexe+nope+j))=inum(kon(indexe+nope+j))+1
c            enddo
c         endif
c     Bernhardi end
!
      enddo
!
!     for 1d and 2d elements only:
!     finding the solution in the original nodes
!
      if((cflag.ne.' ').and.(cflag.ne.'E')) then
         nfield=0
         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,ne,cflag,co,
     &         vold,force,mi,ielprop,prop)
      endif
!
      return
      end
