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
      subroutine nmatrix(ad,au,jqs,irows,ndesi,nodedesi,dgdxglob,
     &   nactive,nobject,nnlconst,ipoacti,nk)         
!
!     calculates the values of the expression: N^(T)N
!
      implicit none
!
      integer jqs(*),irows(*),ndesi,nodedesi(*),idof,i,j,jdof,
     &   nactive,nobject,nnlconst,ipos,jpos,ipoacti(*),nk,node 
!
      real*8 ad(*),au(*),dgdxglob(2,nk,nobject)
!
      do idof=1,nactive
         if(idof.le.nnlconst) then
            ipos=ipoacti(idof)
            do i=1,ndesi
               node=nodedesi(i)
               ad(idof)=ad(idof)+dgdxglob(2,node,ipos)**2
            enddo
            do i=jqs(idof),jqs(idof+1)-1
               jdof=irows(i)
               if(jdof.le.nnlconst) then
                  jpos=ipoacti(i)
                  do j=1,ndesi
                     node=nodedesi(j)
                     au(i)=au(i)+dgdxglob(2,node,ipos)
     &                          *dgdxglob(2,node,jpos)
                  enddo
               else
                  node=nodedesi(ipoacti(i))
                  au(i)=dgdxglob(2,node,ipos)
               endif
            enddo
         else
            ad(idof)=1 
         endif
      enddo
!
      return        
      end




