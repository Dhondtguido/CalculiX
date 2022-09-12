!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine extrapolatefem(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc)
!
!     extrapolates field values at the integration points to the 
!     nodes
!
!     the number of internal state variables is limited to 999
!     (cfr. array field)
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,indexe,nope,
     &  nfield,nk,i,j,k,l,ndim,
     &  iorienloc,iorien,ielorien(mi(3),*),konl,
     &  mint3d,m,iflag
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),field(999,20*mi(3)),a8(8,8),
     &  a4(4,4),a2(6,2),orab(7,*),co(3,*),
     &  coords(3,27),xi,et,ze,xl(3,20),xsj,shp(4,20),weight,
     &  yiloc(6,27),a(3,3),b(3,3),c(3,3)
!
      include "gauss.f"
!
      data a2 /  1.1455,-0.1455,1.1455,-0.1455,1.1455,-0.1455,
     &           -0.1455,1.1455,-0.1455,1.1455,-0.1455,1.1455/
      data a4 /  1.92705, -0.30902, -0.30902, -0.30902,
     &          -0.30902,  1.92705, -0.30902, -0.30902,
     &          -0.30902, -0.30902,  1.92705, -0.30902,
     &          -0.30902, -0.30902, -0.30902,  1.92705/
!
!     extrapolation from a 2x2x2=8 integration point scheme in a hex to
!     the vertex nodes
!    
      data a8 /2.549,-.683,.183,-.683,-.683,.183,
     &        -.04904,.183,-.683,2.549,-.683,.183,
     &        .183,-.683,.183,-.04904,-.683,.183,
     &        -.683,2.549,.183,-.04904,.183,-.683,
     &        .183,-.683,2.549,-.683,-.04904,.183,
     &        -.683,.183,-.683,.183,-.04904,.183,
     &        2.549,-.683,.183,-.683,.183,-.683,
     &        .183,-.04904,-.683,2.549,-.683,.183,
     &        .183,-.04904,.183,-.683,-.683,.183,
     &        -.683,2.549,-.04904,.183,-.683,.183,
     &        .183,-.683,2.549,-.683/  
!
      data iflag /1/
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,nk
         do j=1,nfield
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).le.-1) cycle
         indexe=ipkon(i)
!
         lakonl=lakon(i)
!
         if(lakonl(1:1).ne.'F') then
            cycle
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         else
            cycle
         endif
!
!     storage in local coordinates
!
!     calculation of the integration point coordinates for
!     output in the local system
!
         if((iorienloc.ne.0).and.(ielorien(1,i).ne.0)) then
!
           iorien=ielorien(1,i)
!     
            if(lakon(i)(4:5).eq.'8R') then
               mint3d=1
            elseif(lakon(i)(4:4).eq.'8') then
                  mint3d=8
            elseif(lakon(i)(4:4).eq.'4') then
               mint3d=1
            elseif(lakon(i)(4:4).eq.'6') then
               mint3d=2
            endif
!
            do j=1,nope
               konl=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl)
               enddo
            enddo
!
            do j=1,mint3d
               if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,j)
                  et=gauss3d1(2,j)
                  ze=gauss3d1(3,j)
                  weight=weight3d1(j)
               elseif(lakon(i)(4:4).eq.'8')then
                     xi=gauss3d2(1,j)
                     et=gauss3d2(2,j)
                     ze=gauss3d2(3,j)
                     weight=weight3d2(j)
               elseif(lakon(i)(4:4).eq.'4') then
                  xi=gauss3d4(1,j)
                  et=gauss3d4(2,j)
                  ze=gauss3d4(3,j)
                  weight=weight3d4(j)
               elseif(lakon(i)(4:4).eq.'6') then
                  xi=gauss3d7(1,j)
                  et=gauss3d7(2,j)
                  ze=gauss3d7(3,j)
                  weight=weight3d7(j)
               endif
!
               if(nope.eq.8) then
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.4) then
                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
               endif
!
               if(iorien.eq.0) then
                  if(nfield.eq.3) then
                     do k=1,3
                        yiloc(k,j)=yi(k,j,i)
                     enddo
                  elseif(nfield.eq.6) then
                     do k=1,6
                        yiloc(k,j)=yi(k,j,i)
                     enddo
                  endif
                  cycle
               endif
!
               do k=1,3
                  coords(k,j)=0.d0
                  do l=1,nope
                     coords(k,j)=coords(k,j)+xl(k,l)*shp(4,l)
                  enddo
               enddo
!
               if(nfield.eq.3) then
                  call transformatrix(orab(1,iorien),coords(1,j),a)
                  yiloc(1,j)=yi(1,j,i)*a(1,1)+yi(2,j,i)*a(2,1)+
     &                     yi(3,j,i)*a(3,1)
                  yiloc(2,j)=yi(1,j,i)*a(1,2)+yi(2,j,i)*a(2,2)+
     &                     yi(3,j,i)*a(3,2)
                  yiloc(3,j)=yi(1,j,i)*a(1,3)+yi(2,j,i)*a(2,3)+
     &                     yi(3,j,i)*a(3,3)
               elseif(nfield.eq.6) then
                  call transformatrix(orab(1,iorien),coords(1,j),a)
                  b(1,1)=yi(1,j,i)
                  b(2,2)=yi(2,j,i)
                  b(3,3)=yi(3,j,i)
                  b(1,2)=yi(4,j,i)
                  b(1,3)=yi(5,j,i)
                  b(2,3)=yi(6,j,i)
                  b(2,1)=b(1,2)
                  b(3,1)=b(1,3)
                  b(3,2)=b(2,3)
                  do k=1,3
                     do l=1,3
                        c(k,l)=0.d0
                        do m=1,3
                           c(k,l)=c(k,l)+b(k,m)*a(m,l)
                        enddo
                     enddo
                  enddo
                  do k=1,3
                     do l=k,3
                        b(k,l)=0.d0
                        do m=1,3
                           b(k,l)=b(k,l)+a(m,k)*c(m,l)
                        enddo
                     enddo
                  enddo
                  yiloc(1,j)=b(1,1)
                  yiloc(2,j)=b(2,2)
                  yiloc(3,j)=b(3,3)
                  yiloc(4,j)=b(1,2)
                  yiloc(5,j)=b(1,3)
                  yiloc(6,j)=b(2,3)
               endif
            enddo
!
            if(lakonl(4:5).eq.'8 ') then
                  do j=1,8
                     do k=1,nfield
                        field(k,j)=0.d0
                        do l=1,8
                           field(k,j)=field(k,j)+a8(j,l)*yiloc(k,l)
                        enddo
                     enddo
                  enddo
            elseif(lakonl(4:4).eq.'8') then
               do j=1,8
                  do k=1,nfield
                     field(k,j)=yiloc(k,1)
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'4') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=yiloc(k,1)
                  enddo
               enddo
            else
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,2
                        field(k,j)=field(k,j)+a2(j,l)*yiloc(k,l)
                     enddo
                  enddo
               enddo
            endif
         else
!
!        storage in global coordinates
!
!        determining the field values in the vertex nodes
!        for C3D20R and C3D8: trilinear extrapolation (= use of the
!                             C3D8 linear interpolation functions)
!        for C3D8R: constant field value in each element
!        for C3D10: use of the C3D4 linear interpolation functions
!        for C3D4: constant field value in each element
!        for C3D15: use of the C3D6 linear interpolation functions
!        for C3D6: use of a linear interpolation function
!
            if(lakonl(4:5).eq.'8 ') then
                  do j=1,8
                     do k=1,nfield
                        field(k,j)=0.d0
                        do l=1,8
                           field(k,j)=field(k,j)+a8(j,l)*yi(k,l,i)
                        enddo
                     enddo
                  enddo
            elseif(lakonl(4:4).eq.'8') then
               do j=1,8
                  do k=1,nfield
                     field(k,j)=yi(k,1,i)
                  enddo
               enddo
            elseif(lakonl(4:4).eq.'4') then
               do j=1,4
                  do k=1,nfield
                     field(k,j)=yi(k,1,i)
                  enddo
               enddo
            else
               do j=1,6
                  do k=1,nfield
                     field(k,j)=0.d0
                     do l=1,2
                        field(k,j)=field(k,j)+a2(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            endif
         endif
!
!        transferring the field values into yn
!
            do j=1,nope
               do k=1,nfield
                  yn(k,kon(indexe+j))=yn(k,kon(indexe+j))+
     &                 field(k,j)
               enddo
               inum(kon(indexe+j))=inum(kon(indexe+j))+1
            enddo
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!
      return
      end
