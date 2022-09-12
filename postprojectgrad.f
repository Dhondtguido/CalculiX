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
      subroutine postprojectgrad(ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,objectset,inameacti)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node,iconst,i,m,inameacti(*) 
!
      real*8 dgdxglob(2,nk,nobject),scalar,dd,len
!
!     calculation of final projected gradient 
!     in case of an active constraint
!
      if(nactive.gt.0) then
         do irow=1,ndesi
            node=nodedesi(irow)
            dgdxglob(2,node,nobject)=dgdxglob(2,node,1)
     &          -dgdxglob(2,node,nobject)
         enddo
         objectset(1,nobject)(1:11)='PROJECTGRAD'        
!
         write(*,*)
         write(*,*) '*INFO: at least 1 constraint active.'
         write(*,*) '       projected gradient calculated'   
!
!     prepare output of objective sensitivity as feasible direction
!
      else
         objectset(1,1)(1:11)='PROJECTGRAD' 
         do i=12,20
            objectset(1,1)(i:i)=' '
         enddo
         do i=1,ndesi
            dgdxglob(1,nodedesi(i),1)=0.d0
         enddo
         write(*,*)
         write(*,*) '*INFO: no constraint active'    
         write(*,*) '       no projected gradient calculated'
         write(*,*) '       senstivity of the objective function' 
         write(*,*) '       taken as feasible direction'
         write(*,*)
      endif    
!
      return        
      end

