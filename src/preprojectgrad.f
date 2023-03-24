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
      subroutine preprojectgrad(vector,ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,rhs,objectset,xtf)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node
!
      real*8 dgdxglob(2,nk,nobject),vector(ndesi),rhs(*),scalar,xtf(*)
!
!     calculate the second part of xlambd
!     
      do icol=1,nactive
         if(icol.le.nnlconst) then   
            do irow=1,ndesi      
               ipos=ipoacti(icol)   
               node=nodedesi(irow)
               xtf(icol)=xtf(icol)+dgdxglob(2,node,1)
     &                   *dgdxglob(2,node,ipos)
            enddo
         else
            ipos=ipoacti(icol)
            node=nodedesi(ipos)
            xtf(icol)=dgdxglob(2,node,1)
         endif
      enddo
!
      return        
      end
