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
!
!     Expand contact directions matrix W_b
!     from nodal size into DOF size
!      
      subroutine expand_auw(auw,jqw,iroww,nslavs,auwnew,jqwnew,
     &     irowwnew,nactdof,mi,ktot,neqtot,islavnode,imastnode)
!     
      implicit none
!     
      integer nzswnew,jqwnew(*),jqw(*),iroww(*),nslavs,
     &     irowwnew(*),i,idof,nodes,j,
     &     length,id,mi(*),nactdof(0:mi(2),*),ktot(*),neqtot,
     &     nodem,islavnode(*),imastnode(*),kflag,
     &     jdofm,idirs,idirm,nodemrel,filedebug
!     
      real*8 auw(*),auwnew(*)
!     
      filedebug=0
!     
      nzswnew=0
!     
!     loop over all columns of auwnew 
!     
      do i=1,3*nslavs
        jqwnew(i)=nzswnew+1
        nodes=islavnode((i-1)/3+1)
!
!     check for zero column in W-matrix (no opposite master face
!     to the slave node)
!
        if(jqw(i+1).eq.jqw(i)) cycle
!     
!     first node is the slave node; identify its row
!     position
!     
        do idirs=1,3
          idof=nactdof(idirs,nodes)
          call nident(ktot,idof,neqtot,id)
          if(id.gt.0) then
            if(ktot(id).eq.idof) then
              nzswnew=nzswnew+1
              auwnew(nzswnew)=auw(jqw(i)-1+idirs)
              irowwnew(nzswnew)=id
            endif
          endif
        enddo
!     
!     independent (master) dofs
!     
        do j=jqw(i)+3,jqw(i+1)-1
          jdofm=iroww(j)-3*nslavs
          nodemrel=(jdofm-1)/3+1 ! position of master node in field imastnode
          nodem=imastnode(nodemrel) ! master node
          idirm=jdofm-3*(nodemrel-1) ! coordinate direction
          idof=nactdof(idirm,nodem)
          call nident(ktot,idof,neqtot,id)
          if(id.gt.0) then
            if(ktot(id).eq.idof) then
              nzswnew=nzswnew+1
              auwnew(nzswnew)=auw(j)
              irowwnew(nzswnew)=id
              cycle
            endif
          endif
        enddo
      enddo
!     
      jqwnew(3*nslavs+1)=nzswnew+1
!     
!     sorting the column in auwnew
!     
      kflag=2
      do i=1,3*nslavs
        length=jqwnew(i+1)-jqwnew(i)
        call isortid(irowwnew(jqwnew(i)),auwnew(jqwnew(i)),length,
     &       kflag)
      enddo
!
      return
      end

