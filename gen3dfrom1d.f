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
      subroutine gen3dfrom1d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &  thicke,ntrans,inotr,trab,nk,nk_,co,offset,mi)
!
!     expands 1d element i into a 3d element
!
      implicit none
!
      character*8 lakon(*)
!
      integer mi(*),i,kon(*),ipkon(*),ne,iponor(2,*),knor(*),ntrans,
     &  inotr(2,*),nk,nk_,indexe,j,nodel(8),indexx,indexk,k,nodeb(8,3),
     &  nope,ishift
!
      real*8 xnor(*),thicke(mi(3),*),trab(7,*),co(3,*),offset(2,*),
     &  thickb(2,3),xnorb(6,3),sc
!
      indexe=ipkon(i)
!
!     check whether linear or quadratic
!
      if((lakon(i)(3:3).eq.'1').or.(lakon(i)(4:4).eq.'2')) then
         nope=2
         if(lakon(i)(4:4).eq.'R') then
            ishift=8
         else
            ishift=11
         endif
      elseif((lakon(i)(3:3).eq.'2').or.(lakon(i)(4:5).eq.'3')) then
         nope=3
         ishift=20
      endif
!
!     localizing the nodes, thicknesses and normals for the
!     beam element
!            
c      do j=1,3
      do j=1,nope
         nodel(j)=kon(indexe+j)
         kon(indexe+ishift+j)=nodel(j)
         indexx=iponor(1,indexe+j)
         indexk=iponor(2,indexe+j)
         thickb(1,j)=thicke(1,indexe+j)
         thickb(2,j)=thicke(2,indexe+j)
         do k=1,6
            xnorb(k,j)=xnor(indexx+k)
         enddo
         do k=1,8
            nodeb(k,j)=knor(indexk+k)
         enddo
         if(ntrans.gt.0) then
            do k=1,8
               inotr(1,nodeb(k,j))=inotr(1,nodel(j))
            enddo
         endif
      enddo
!
!     generating the 3-D element topology for beam elements
!
!     rectangular cross section (or parent section)
!
      if(lakon(i)(8:8).ne.'C') then
         kon(indexe+1)=nodeb(1,1)
         do j=1,3
            co(j,nodeb(1,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &           +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
         enddo
         kon(indexe+2)=nodeb(1,nope)
         do j=1,3
            co(j,nodeb(1,nope))=co(j,nodel(nope))
     &           -thickb(1,nope)*xnorb(j,nope)*(.5d0+offset(1,i))
     &           +thickb(2,nope)*xnorb(j+3,nope)*(.5d0-offset(2,i))
         enddo
         kon(indexe+3)=nodeb(2,nope)
         do j=1,3
            co(j,nodeb(2,nope))=co(j,nodel(nope))
     &           -thickb(1,nope)*xnorb(j,nope)*(.5d0+offset(1,i))
     &           -thickb(2,nope)*xnorb(j+3,nope)*(.5d0+offset(2,i))
         enddo
         kon(indexe+4)=nodeb(2,1)
         do j=1,3
            co(j,nodeb(2,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
         enddo
         kon(indexe+5)=nodeb(4,1)
         do j=1,3
            co(j,nodeb(4,1))=co(j,nodel(1))
     &           +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &           +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
         enddo
         kon(indexe+6)=nodeb(4,nope)
         do j=1,3
            co(j,nodeb(4,nope))=co(j,nodel(nope))
     &           +thickb(1,nope)*xnorb(j,nope)*(.5d0-offset(1,i))
     &           +thickb(2,nope)*xnorb(j+3,nope)*(.5d0-offset(2,i))
         enddo
         kon(indexe+7)=nodeb(3,nope)
         do j=1,3
            co(j,nodeb(3,nope))=co(j,nodel(nope))
     &           +thickb(1,nope)*xnorb(j,nope)*(.5d0-offset(1,i))
     &           -thickb(2,nope)*xnorb(j+3,nope)*(.5d0+offset(2,i))
         enddo
         kon(indexe+8)=nodeb(3,1)
         do j=1,3
            co(j,nodeb(3,1))=co(j,nodel(1))
     &           +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
         enddo
!
!        middle nodes for quadratic elements
!
         if(nope.eq.3) then
            kon(indexe+9)=nodeb(1,2)
            do j=1,3
               co(j,nodeb(1,2))=co(j,nodel(2))
     &              -thickb(1,2)*xnorb(j,2)*(.5d0+offset(1,i))
     &              +thickb(2,2)*xnorb(j+3,2)*(.5d0-offset(2,i))
            enddo
            kon(indexe+10)=nodeb(5,3)
            do j=1,3
               co(j,nodeb(5,3))=co(j,nodel(3))
     &              -thickb(1,3)*xnorb(j,3)*(.5d0+offset(1,i))
     &              -thickb(2,3)*xnorb(j+3,3)*offset(2,i)
            enddo
            kon(indexe+11)=nodeb(2,2)
            do j=1,3
               co(j,nodeb(2,2))=co(j,nodel(2))
     &              -thickb(1,2)*xnorb(j,2)*(.5d0+offset(1,i))
     &              -thickb(2,2)*xnorb(j+3,2)*(.5d0+offset(2,i))
            enddo
            kon(indexe+12)=nodeb(5,1)
            do j=1,3
               co(j,nodeb(5,1))=co(j,nodel(1))
     &              -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &              -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
            enddo
            kon(indexe+13)=nodeb(4,2)
            do j=1,3
               co(j,nodeb(4,2))=co(j,nodel(2))
     &              +thickb(1,2)*xnorb(j,2)*(.5d0-offset(1,i))
     &              +thickb(2,2)*xnorb(j+3,2)*(.5d0-offset(2,i))
            enddo
            kon(indexe+14)=nodeb(7,3)
            do j=1,3
               co(j,nodeb(7,3))=co(j,nodel(3))
     &              +thickb(1,3)*xnorb(j,3)*(.5d0-offset(1,i))
     &              -thickb(2,3)*xnorb(j+3,3)*offset(2,i)
            enddo
            kon(indexe+15)=nodeb(3,2)
            do j=1,3
               co(j,nodeb(3,2))=co(j,nodel(2))
     &              +thickb(1,2)*xnorb(j,2)*(.5d0-offset(1,i))
     &              -thickb(2,2)*xnorb(j+3,2)*(.5d0+offset(2,i))
            enddo
            kon(indexe+16)=nodeb(7,1)
            do j=1,3
               co(j,nodeb(7,1))=co(j,nodel(1))
     &              +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &              -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
            enddo
            kon(indexe+17)=nodeb(8,1)
            do j=1,3
               co(j,nodeb(8,1))=co(j,nodel(1))
     &              -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &              +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
            enddo
            kon(indexe+18)=nodeb(8,3)
            do j=1,3
               co(j,nodeb(8,3))=co(j,nodel(3))
     &              -thickb(1,3)*xnorb(j,3)*offset(1,i)
     &              +thickb(2,3)*xnorb(j+3,3)*(.5d0-offset(2,i))
            enddo
            kon(indexe+19)=nodeb(6,3)
            do j=1,3
               co(j,nodeb(6,3))=co(j,nodel(3))
     &              -thickb(1,3)*xnorb(j,3)*offset(1,i)
     &              -thickb(2,3)*xnorb(j+3,3)*(.5d0+offset(2,i))
            enddo
            kon(indexe+20)=nodeb(6,1)
            do j=1,3
               co(j,nodeb(6,1))=co(j,nodel(1))
     &              -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &              -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
            enddo
!
!           generating coordinates for the expanded nodes which
!           are not used by the C3D20(R) element (needed for the
!           determination of the knot dimension)
!
            do j=1,3
               co(j,nodeb(5,2))=co(j,nodel(2))
     &              -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &              -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
            enddo
            do j=1,3
               co(j,nodeb(7,2))=co(j,nodel(2))
     &              +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &              -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
            enddo
            do j=1,3
               co(j,nodeb(8,2))=co(j,nodel(2))
     &              -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &              +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
            enddo
            do j=1,3
               co(j,nodeb(6,2))=co(j,nodel(2))
     &              -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &              -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
            enddo
         endif
      else
!
!                 circular cross section
!
         if(nope.eq.2) then
            write(*,*) '*ERROR in gen3dfrom1d: element ',i,
     &        'is a linear beam element with circular cross section'
            write(*,*) '       Please use quadratic elements for beams
     &with circular cross section.'
            call exit(201)
         endif
!
         sc=.5d0/dsqrt(2.d0)
         kon(indexe+1)=nodeb(1,1)
         do j=1,3
            co(j,nodeb(1,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*(sc+offset(1,i))
     &           +thickb(2,1)*xnorb(j+3,1)*(sc-offset(2,i))
         enddo
         kon(indexe+2)=nodeb(1,3)
         do j=1,3
            co(j,nodeb(1,3))=co(j,nodel(3))
     &           -thickb(1,3)*xnorb(j,3)*(sc+offset(1,i))
     &           +thickb(2,3)*xnorb(j+3,3)*(sc-offset(2,i))
         enddo
         kon(indexe+3)=nodeb(2,3)
         do j=1,3
            co(j,nodeb(2,3))=co(j,nodel(3))
     &           -thickb(1,3)*xnorb(j,3)*(sc+offset(1,i))
     &           -thickb(2,3)*xnorb(j+3,3)*(sc+offset(2,i))
         enddo
         kon(indexe+4)=nodeb(2,1)
         do j=1,3
            co(j,nodeb(2,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*(sc+offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*(sc+offset(2,i))
         enddo
         kon(indexe+5)=nodeb(4,1)
         do j=1,3
            co(j,nodeb(4,1))=co(j,nodel(1))
     &           +thickb(1,1)*xnorb(j,1)*(sc-offset(1,i))
     &           +thickb(2,1)*xnorb(j+3,1)*(sc-offset(2,i))
         enddo
         kon(indexe+6)=nodeb(4,3)
         do j=1,3
            co(j,nodeb(4,3))=co(j,nodel(3))
     &           +thickb(1,3)*xnorb(j,3)*(sc-offset(1,i))
     &           +thickb(2,3)*xnorb(j+3,3)*(sc-offset(2,i))
         enddo
         kon(indexe+7)=nodeb(3,3)
         do j=1,3
            co(j,nodeb(3,3))=co(j,nodel(3))
     &           +thickb(1,3)*xnorb(j,3)*(sc-offset(1,i))
     &           -thickb(2,3)*xnorb(j+3,3)*(sc+offset(2,i))
         enddo
         kon(indexe+8)=nodeb(3,1)
         do j=1,3
            co(j,nodeb(3,1))=co(j,nodel(1))
     &           +thickb(1,1)*xnorb(j,1)*(sc-offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*(sc+offset(2,i))
         enddo
         kon(indexe+9)=nodeb(1,2)
         do j=1,3
            co(j,nodeb(1,2))=co(j,nodel(2))
     &           -thickb(1,2)*xnorb(j,2)*(sc+offset(1,i))
     &           +thickb(2,2)*xnorb(j+3,2)*(sc-offset(2,i))
         enddo
         kon(indexe+10)=nodeb(5,3)
         do j=1,3
            co(j,nodeb(5,3))=co(j,nodel(3))
     &           -thickb(1,3)*xnorb(j,3)*(.5d0+offset(1,i))
     &           -thickb(2,3)*xnorb(j+3,3)*offset(2,i)
         enddo
         kon(indexe+11)=nodeb(2,2)
         do j=1,3
            co(j,nodeb(2,2))=co(j,nodel(2))
     &           -thickb(1,2)*xnorb(j,2)*(sc+offset(1,i))
     &           -thickb(2,2)*xnorb(j+3,2)*(sc+offset(2,i))
         enddo
         kon(indexe+12)=nodeb(5,1)
         do j=1,3
            co(j,nodeb(5,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
         enddo
         kon(indexe+13)=nodeb(4,2)
         do j=1,3
            co(j,nodeb(4,2))=co(j,nodel(2))
     &           +thickb(1,2)*xnorb(j,2)*(sc-offset(1,i))
     &           +thickb(2,2)*xnorb(j+3,2)*(sc-offset(2,i))
         enddo
         kon(indexe+14)=nodeb(7,3)
         do j=1,3
            co(j,nodeb(7,3))=co(j,nodel(3))
     &           +thickb(1,3)*xnorb(j,3)*(.5d0-offset(1,i))
     &           -thickb(2,3)*xnorb(j+3,3)*offset(2,i)
         enddo
         kon(indexe+15)=nodeb(3,2)
         do j=1,3
            co(j,nodeb(3,2))=co(j,nodel(2))
     &           +thickb(1,2)*xnorb(j,2)*(sc-offset(1,i))
     &           -thickb(2,2)*xnorb(j+3,2)*(sc+offset(2,i))
         enddo
         kon(indexe+16)=nodeb(7,1)
         do j=1,3
            co(j,nodeb(7,1))=co(j,nodel(1))
     &           +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
         enddo
         kon(indexe+17)=nodeb(8,1)
         do j=1,3
            co(j,nodeb(8,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &           +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
         enddo
         kon(indexe+18)=nodeb(8,3)
         do j=1,3
            co(j,nodeb(8,3))=co(j,nodel(3))
     &           -thickb(1,3)*xnorb(j,3)*offset(1,i)
     &           +thickb(2,3)*xnorb(j+3,3)*(.5d0-offset(2,i))
         enddo
         kon(indexe+19)=nodeb(6,3)
         do j=1,3
            co(j,nodeb(6,3))=co(j,nodel(3))
     &           -thickb(1,3)*xnorb(j,3)*offset(1,i)
     &           -thickb(2,3)*xnorb(j+3,3)*(.5d0+offset(2,i))
         enddo
         kon(indexe+20)=nodeb(6,1)
         do j=1,3
            co(j,nodeb(6,1))=co(j,nodel(1))
     &           -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &           -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
         enddo
!
!           generating coordinates for the expanded nodes which
!           are not used by the C3D20(R) element (needed for the
!           determination of the knot dimension)
!
         do j=1,3
            co(j,nodeb(5,2))=co(j,nodel(2))
     &           -thickb(1,1)*xnorb(j,1)*(.5d0+offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
         enddo
         do j=1,3
            co(j,nodeb(7,2))=co(j,nodel(2))
     &           +thickb(1,1)*xnorb(j,1)*(.5d0-offset(1,i))
     &           -thickb(2,1)*xnorb(j+3,1)*offset(2,i)
         enddo
         do j=1,3
            co(j,nodeb(8,2))=co(j,nodel(2))
     &           -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &           +thickb(2,1)*xnorb(j+3,1)*(.5d0-offset(2,i))
         enddo
         do j=1,3
            co(j,nodeb(6,2))=co(j,nodel(2))
     &           -thickb(1,1)*xnorb(j,1)*offset(1,i)
     &           -thickb(2,1)*xnorb(j+3,1)*(.5d0+offset(2,i))
         enddo
      endif
!
      return
      end


