!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine projectgrad(vector,ndesi,nodedesi,dgdxglob,nactive,
     &     nobject,nnlconst,ipoacti,nk,rhs,iconst,objectset,xlambd,xtf,
     &     objnorm,gradproj,g0)         
!     
!     calculates the projected gradient, the lagrange multipliers and 
!     the correction term due to nonlinear constraints
!     
      implicit none
!     
      character*81 objectset(5,*)
!     
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &     ipoacti(*),nk,ipos,node,iconst,i
!     
      real*8 dgdxglob(2,nk,nobject),vector(ndesi),rhs(*),scalar,dd,
     &     len,xlambd(*),xtf(*),objnorm(*),gradproj(3,*),g0(*),val
!     
!     calculation of the vector
!     
      do irow=1,ndesi
        do icol=1,nactive
          if(icol.le.nnlconst) then
            ipos=ipoacti(icol)
            node=nodedesi(irow)
            vector(irow)=vector(irow)
     &           +rhs(icol)*dgdxglob(2,node,ipos)
          else
            if(ipoacti(icol).eq.irow) then
              vector(irow)=vector(irow)+rhs(icol)
            endif
          endif
        enddo
      enddo
!     
!     calculation of the scalar value
!     
      scalar=0.d0
      if(iconst.le.nnlconst) then
        do irow=1,ndesi
          ipos=ipoacti(iconst)      
          node=nodedesi(irow)
          scalar=scalar+dgdxglob(2,node,1)*dgdxglob(2,node,ipos)
        enddo
      else
        node=nodedesi(ipoacti(iconst))
        scalar=dgdxglob(2,node,1)
      endif            
!     
!     multiplication of scalar and vector
!     
      do irow=1,ndesi
        node=nodedesi(irow)
        gradproj(2,node)=gradproj(2,node)+vector(irow)*scalar
      enddo
!     
!     calculation of lagrange multiplier
!     
      do icol=1,nactive
        xlambd(icol)=xlambd(icol)+(-1)*rhs(icol)*xtf(iconst)
      enddo
!     
!     calculation of correction term
!     
      if(objnorm(iconst).gt.0.d0) then
        if(iconst.le.nnlconst) then
          do irow=1,ndesi
            node=nodedesi(irow)     
            gradproj(1,node)=gradproj(1,node)
     &           +vector(irow)*objnorm(iconst)
          enddo
        else
          irow=ipoacti(iconst)
          node=nodedesi(irow)   
          val=vector(irow)*objnorm(iconst)
          gradproj(1,node)=gradproj(1,node)
     &         +vector(irow)*objnorm(iconst)
        endif
      endif
!     
      return        
      end
