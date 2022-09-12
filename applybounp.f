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
      subroutine applybounp(nodeboun,ndirboun,nboun,xbounact,nk,vold,
     &     v,ipompc,nodempc,coefmpc,nmpc,inomat,mi)
!     
!     applies velocity boundary conditions
!     
      implicit none
!     
      integer mi(*),nodeboun(*),ndirboun(*),nk,i,nboun,node,imat,
     &     index,ipompc(*),nodempc(3,*),nmpc,ist,ndir,inomat(*),ndiri,
     &     nodei
!     
      real*8 vold(0:mi(2),*),xbounact(*),v(nk,0:mi(2)),coefmpc(*),residu
!     
!     inserting the pressure boundary conditions
!     
      do i=1,nboun
        if(ndirboun(i).ne.4) cycle
!     
        node=nodeboun(i)
        v(node,4)=xbounact(i)-vold(4,node)
      enddo
!     
!     inserting the pressure MPC conditions
!     
      do i=1,nmpc
        ist=ipompc(i)
        ndir=nodempc(2,ist)
        if(ndir.ne.4) cycle
        node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
        imat=inomat(node)
!     
        index=nodempc(3,ist)
        residu=0.d0
        if(index.ne.0) then
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
!     
            residu=residu+coefmpc(index)*
     &           (vold(ndiri,nodei)+v(nodei,ndiri))
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        endif
!     
        v(node,ndir)=-residu/coefmpc(ist)-vold(ndir,node)
      enddo
!     
      return
      end
