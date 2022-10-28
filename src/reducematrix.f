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
!     MASSLESS DYNAMIC CONTACT
!     
C     >    *MASSLESS DYNAMIC CONTACT*: Reducing the global stiffness matrix to the ii matrix by setting contact diagonal entries to 1 and contact subdiagonal entries to 0

C     > Submatrices:
C     > + the bb matrix (neqtot x neqtot)
C     > + the bi matrix (neqtot x neq(1))
C     > + the ib matrix (neq(1) x neqtot)
C     > + the ii matrix (neq(1) x neq(1))
C     >    - this matrix has the same
C     >     size of the global matrix, in which the bi-entries,
C     >     ib-entries, bb-entries (off-diagonal) are set to zero
C     >     the diagonal bb-entries are set to one

C     >     @param au        LOWER triangle of STIFFNESS matrix of size: neq[0]
C     >     @param ad        DIAGONAL       of STIFFNESS matrix of size: neq[0]
C     >     @param jq        Location in field **irow** of the first subdiagonal nonzero in column i (only for symmetric matrices)
C     >     @param irow      Row of element i in field au (i.e. au(i)) 
C     >     @param  neq        NTOT of equations: 
C     >                            + neq[0]: mechanical TOTDOF, 
C     >                            + neq[1]: TOTDOF_mech + TOTDOF_thermal
C     >                            + neq[2]: neq[1] + # of single point constraints (only for modal calculations)
C     >     @param neqtot    the param
C     >     @param ktot       slave and master dofs list
C     
C     >    @see massless.c
C     >    @author Guido Dhondt
!     
      subroutine reducematrix(au,ad,jq,irow,neq,neqtot,ktot)
!     
      implicit none
!     
      integer jq(*),irow(*),neq(*),neqtot,kbb,i,j,ktot(*),id
!     
      real*8 au(*),ad(*)
!     
!     treating column by column in the global stiffness matrix
!     
      kbb=1
      do i=1,neq(1)
        if(i.ne.ktot(kbb)) then
!
!     i column: set contact dof entries to zero
!
          do j=jq(i),jq(i+1)-1
            call nident(ktot,irow(j),neqtot,id)
            if(id.gt.0) then
              if(ktot(id).eq.irow(j)) then
                au(j)=0.d0
              endif
            endif
          enddo
          cycle
        else
!
!     b column: set all off-diagonal entries to zero and
!               diagonal entry to one.
!
          do j=jq(i),jq(i+1)-1
            au(j)=0.d0
          enddo
          ad(i)=1.d0
          if(kbb.lt.neqtot) kbb=kbb+1
        endif
      enddo
!     
      return
      end
