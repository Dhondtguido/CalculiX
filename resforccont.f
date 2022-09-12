!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine resforccont(vold,nk,mi,aubi,irowbi,jqbi,neqtot,ktot,
     &     fext,gapdisp,auib,irowib,jqib,nactdof,volddof,neq,qik_kbi)
!     
      implicit none
!     
      integer irowbi(*),jqbi(*),nk,neqtot,irowib(*),jqib(*),itranspose,
     &     mi(*),nactdof(0:mi(2),*),i,j,neq(*),ktot(*)
!     
      real*8 vold(0:mi(2),*),fext(*),qik_kbi(*),
     &     auib(*),aubi(*),gapdisp(*),volddof(*)
!     
!     create field volddof, from vold (displacements)
!     from sorting nodes to DOF
!     
      do i=1,nk
        do j=1,3
          if(nactdof(j,i).gt.0) then
            volddof(nactdof(j,i))=vold(j,i)
          endif
        enddo
      enddo
!     
!     We compute g as volddof=(Kbi*volddof)+(Kib*volddof) in qik_kbi
!     to account for the missing terms due to the low triangle structure
!     of the matrices
!     
!     calculate Kbi*volddof
!     
      itranspose=0
      call mulmatvec_asym(aubi,jqbi,irowbi,neq(1),volddof,qik_kbi,
     &     itranspose)
!     
!     calculate Kib^T*volddof and add to g.
!     transposed multiplication
!     
      itranspose=1
      call mulmatvec_asym(auib,jqib,irowib,neqtot,volddof,qik_kbi,
     &     itranspose)
!     
!     add external force of BOUNDARY DOF
!     
      do i=1,neqtot
        gapdisp(i)=fext(ktot(i))-qik_kbi(i)
      enddo
!     
      return
      end
