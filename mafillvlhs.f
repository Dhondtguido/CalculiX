!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine mafillvlhs(nk,kon,ipkon,lakon,ne,icolv,jqv,irowv,
     &     nzsv,adbv,aubv,ipvar,var)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer kon(*),icolv(*),jqv(*),nzsv,konl(8),irowv(*),ipkon(*),
     &     ipvar(*),nk,ne,nzlv,i,j,l,jdof1,jdof2,indexe,nope
!     
      real*8 sm(8,8),adbv(*),aubv(*),var(*)
!     
!     determining nzlv
!     
      nzlv=0
      do i=nk,1,-1
        if(icolv(i).gt.0) then
          nzlv=i
          exit
        endif
      enddo
!     
      do i=1,nk
        adbv(i)=0.d0
      enddo
      do i=1,nzsv
        aubv(i)=0.d0
      enddo
!     
!     loop over all fluid elements
!     
      do i=1,ne
!     
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'F') cycle
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'4') then
          nope=4
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
!     
        call e_c3d_vlhs(lakon(i),sm,i,ipvar,var)
!     
        do j=1,nope
!     
          jdof1=kon(indexe+j)
!     
          do l=j,nope
!     
            jdof2=kon(indexe+l)
!     
            call add_sm_fl(aubv,adbv,jqv,irowv,jdof1,jdof2,
     &           sm(j,l),j,l)
          enddo
        enddo
      enddo
!     
      return
      end
