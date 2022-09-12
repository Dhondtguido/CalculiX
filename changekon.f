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
      subroutine changekon(ne,ipkon,lakon,mi,nkon,thicke,ielmat,kon)
!
!     for composites the connectivity has to be changed to 
!     allow for the multiple expansions of the composite elements.
!
!     for a 8-node composite shell the connectivity is set up as follows:
!     1..20: simple expansion of the shell
!     21..28: original connectivity of the element
!     29..48: expansion of the first layer of the shell
!     49..68: expansion of the second layer of the shell
!     ..
!     29+(n-1)*20..28+n*20: expansion of the nth layer of the shell
!
!     similar for 6-node composite shell elements; at the start of the present
!     routine the connectivity does not provide for the expansion of the
!     layers yet. At the end of the routine appropriate space has been
!     provided and the regular connectivity of the elements appropriately
!     moved
!
      implicit none
!
      character*8 lakon(*)
!
      integer ne,ipkon(*),mi(*),nkon,ielmat(mi(3),*),nkondiff,i,j,k,
     &  kon(*),nexp,nopeexp,nlayer,ipointer
!
      real*8 thicke(mi(3),*)
!
      integer,dimension(:),allocatable::koncp
      real*8,dimension(:,:),allocatable::thickecp
!
!     calculate the extra space needed in kon
!
      nkondiff=0
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'S') cycle
         if(lakon(i)(8:8).ne.'C') cycle
         do  j=1,mi(3)
            if(ielmat(j,i).ne.0) then
               if(lakon(i)(2:2).eq.'8') then
                  nkondiff=nkondiff+20
               elseif(lakon(i)(2:2).eq.'6') then
                  nkondiff=nkondiff+15
               endif
            endif
         enddo
      enddo
!
!     move the topology in order to create appropriate space
!
      nkon=nkon+nkondiff
      ipointer=nkon
!
      allocate(koncp(nkon))
      allocate(thickecp(mi(3),nkon))
!
      do i=ne,1,-1
         if(ipkon(i).lt.0) cycle
!
!        calculating the size of the expanded connectivity
!
         if(lakon(i)(1:5).eq.'C3D8I') then
            nopeexp=11
         elseif(lakon(i)(4:5).eq.'20') then
            nopeexp=20
         elseif(
     &      (lakon(i)(1:4).eq.'CPE8').or.(lakon(i)(1:4).eq.'CPS8').or.
     &      (lakon(i)(1:4).eq.'CAX8').or.(lakon(i)(1:2).eq.'S8')) then
            nopeexp=28
         elseif(
     &      (lakon(i)(1:4).eq.'CPE6').or.(lakon(i)(1:4).eq.'CPS6').or.
     &      (lakon(i)(1:4).eq.'CAX6').or.(lakon(i)(1:2).eq.'S6')) then
            nopeexp=21
         elseif(lakon(i)(1:1).eq.'B') then
            nopeexp=23
         elseif(lakon(i)(4:4).eq.'8') then
            nopeexp=8
         elseif(lakon(i)(4:5).eq.'10') then
            nopeexp=10
         elseif(lakon(i)(4:4).eq.'4') then
            nopeexp=4
         elseif(lakon(i)(4:5).eq.'15') then
            nopeexp=15
         elseif(lakon(i)(4:4).eq.'6') then
            nopeexp=6
         elseif(lakon(i)(1:8).eq.'DASHPOTA') then
            nopeexp=2
         elseif(lakon(i)(1:1).eq.'D') then
            nopeexp=3
         elseif(lakon(i)(1:1).eq.'G') then
            nopeexp=2
         elseif(lakon(i)(1:7).eq.'SPRINGA') then
            nopeexp=2
         else
            write(*,*) '*ERROR in changekon: element type unknown:',
     &                 ' element: ',i,' type: ',lakon(i)
            call exit(201)
         endif
!
         nlayer=0
!
         if(lakon(i)(8:8).eq.'C') then
            do j=1,mi(3)
               if(ielmat(j,i).ne.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
         endif
!        
         nexp=0
!
         if(lakon(i)(2:2).eq.'8') then
            nexp=20
         elseif(lakon(i)(2:2).eq.'6') then
            nexp=15
         endif
!
         ipointer=ipointer-nopeexp-nlayer*nexp
!
         do j=nopeexp,1,-1
            koncp(ipointer+j)=kon(ipkon(i)+j)
            do k=1,mi(3)
               thickecp(k,ipointer+j)=thicke(k,ipkon(i)+j)
            enddo
         enddo
         ipkon(i)=ipointer
c         write(*,*) 'changekon '
c         write(*,*) i,ipkon(i)
c         write(*,*) (kon(ipointer+j),j=1,nopeexp+nlayer*nexp)
      enddo
!
      do i=1,nkon
         kon(i)=koncp(i)
         do j=1,mi(3)
            thicke(j,i)=thickecp(j,i)
         enddo
      enddo
!
      deallocate(koncp)
      deallocate(thickecp)
!     
      return
      end


