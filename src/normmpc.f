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
      subroutine normmpc(nmpc,ipompc,nodempc,coefmpc,inomat,coefmodmpc,
     &     ikboun,nboun)
!     
!     normalizing the coefficients of MPC's for fluid applications
!     (CFD)
!     
      implicit none
!     
      integer i,node,index,nmpc,inomat(*),nodempc(3,*),ipompc(*),
     &     ikboun(*),nboun,nodei,ndiri,idof,id
!     
      real*8 coefmpc(*),size,coefmodmpc(*)
!     
!     normalizing the MPC-coefficients
!     
      do i=1,nmpc
        index=ipompc(i)
!     
!     check whether fluid node
!     
        node=nodempc(1,index)
        if(inomat(node).eq.0) cycle
!     
!     calculating sum of the square of the MPC coefficients
!     
        size=coefmpc(index)**2
        coefmodmpc(index)=coefmpc(index)
        do
          index=nodempc(3,index)
          if(index.eq.0) exit
!     
!     check whether term is subject to SPC
!     
          nodei=nodempc(1,index)
          ndiri=nodempc(2,index)
          idof=8*(nodei-1)+ndiri
          call nident(ikboun,idof,nboun,id)
          if(id.gt.0) then
            if(ikboun(id).eq.idof) then
              coefmodmpc(index)=0.d0
              cycle
            endif
          endif
          coefmodmpc(index)=coefmpc(index)
!
          size=size+coefmpc(index)**2
        enddo
!
        size=dsqrt(size)
!     
!     normalizing all terms of the MPC
!     
        index=ipompc(i)
        do
          coefmodmpc(index)=coefmodmpc(index)/size
          coefmpc(index)=coefmpc(index)/size
          index=nodempc(3,index)
          if(index.eq.0) exit
        enddo
      enddo
!     
      return
      end
      
