!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
!     Authors: Carlo Monjaraz and Guido Dhondt  
!
!     extract matrices bb, bi and ib from the stiffness matrix
!     setting those entries to zero (off-diagonal terms) and one
!     (diagonal terms) to get ii
!
!     all matrices only contain the subdiagonal terms
!
!     to obtain the complete bi-matrix it has to be complemented by
!     the transpose of the ib-matrix
!
      subroutine extract_matrices(au,ad,jq,irow,neq,aubb,adbb,jqbb,
     &     irowbb,neqtot,nzsbb,aubi,
     &     jqbi,irowbi,nzsbi,auib,jqib,irowib,nzsib,ktot,icolbb)
!      
      implicit none
!     
      integer jq(*),irow(*),neq(*),jqbb(*),irowbb(*),neqtot,
     &     nzsbb,jqbi(*),irowbi(*),nzsbi,jqib(*),irowib(*),nzsib,
     &     kbb,i,j,ktot(*),id,lbb, icolbb(*)
!     
      real*8 au(*),ad(*),aubb(*),adbb(*),aubi(*),auib(*)
!     
!     number of potential nonzero subdiagonal entries
!     
      nzsbb=0
      nzsbi=0
      nzsib=0
!     
!     treating column by column in the global stiffness matrix
!
!     column counter for the contact dofs
!
      kbb=1
!      
      do i=1,neq(1)
!     
!     filling submatrix bi
!     
        if(i.ne.ktot(kbb)) then
!
!         i-column
!
          jqbi(i)=nzsbi+1
!
          do j=jq(i),jq(i+1)-1
            call nident(ktot,irow(j),neqtot,id)
            if(id.gt.0) then
              if(ktot(id).eq.irow(j)) then
!
!                b-row: copy entry to bi-matrix; remove from ii-matrix.
!
                nzsbi=nzsbi+1
                aubi(nzsbi)=au(j)
                au(j)=0.d0
                irowbi(nzsbi)=id
              endif
            endif
          enddo
          cycle
        endif
!        
!       b-column
!
!       row counter for the contact dofs
!
        lbb=1
!
        jqbi(i)=nzsbi+1
        jqbb(kbb)=nzsbb+1
        jqib(kbb)=nzsib+1
!     
        loop: do j=jq(i),jq(i+1)-1
          if(lbb.gt.neqtot) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
            nzsib=nzsib+1
            auib(nzsib)=au(j)
            au(j)=0.d0
            irowib(nzsib)=irow(j)
            cycle
          elseif(irow(j).lt.ktot(lbb)) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!           
            nzsib=nzsib+1
            auib(nzsib)=au(j)
            au(j)=0.d0
            irowib(nzsib)=irow(j)
            cycle
          elseif(irow(j).eq.ktot(lbb)) then
!
!           b-row: copy to bb-matrix and remove from ii-matrix          
!           
            nzsbb=nzsbb+1
            aubb(nzsbb)=au(j)
            au(j)=0.d0
            irowbb(nzsbb)=lbb
            lbb=lbb+1
            cycle
          else
!     
!     actual row entry irow(j) exceeds actual contact dof ktot(lbb).
!     Look for the next contact dof exceeding the actual row entry
!     (the procedure with kbb could be replaced by a call to nident;
!     however, it is assumed that the kbb-procedure is faster)
!     
            do
              lbb=lbb+1
              if(lbb.gt.neqtot) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
                nzsib=nzsib+1
                auib(nzsib)=au(j)
                au(j)=0.d0
                irowib(nzsib)=irow(j)
                cycle loop
              elseif(irow(j).lt.ktot(lbb)) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
                nzsib=nzsib+1
                auib(nzsib)=au(j)
                au(j)=0.d0
                irowib(nzsib)=irow(j)
                cycle loop
              elseif(irow(j).eq.ktot(lbb)) then
!
!           b-row: copy to bb-matrix and remove from ii-matrix          
!
                nzsbb=nzsbb+1
                aubb(nzsbb)=au(j)
                au(j)=0.d0
                irowbb(nzsbb)=lbb
                lbb=lbb+1
                cycle loop
              endif
            enddo
          endif
        enddo loop
!
!       b-diagonal row: copy to bb-matrix and set entry in
!       ii-matrix to 1
!
        adbb(kbb)=ad(i)
        ad(i)=1.d0
!        
        if(kbb.lt.neqtot) kbb=kbb+1
      enddo
!
!     finalizing the jq-entries
!
      jqbb(neqtot+1)=nzsbb+1
      jqbi(neq(1)+1)=nzsbi+1
      jqib(neqtot+1)=nzsib+1
!     
      do i=1,neqtot
        icolbb(i)=jqbb(i+1)-jqbb(i)
      enddo
!
      return
      end
