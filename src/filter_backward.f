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
      subroutine filter_backward(au,jq,irow,icol,ndesi,nodedesi,
     &     dgdxglob,nobject,nk,nobjectstart,objectset)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      character*81 objectset(5,*)
!
      integer jq(*),irow(*),icol(*),ndesi,nodedesi(*),nobject,ipos,
     &   inode1,inode2,kk,jj,ii,nk,nobjectstart,i
!     
      real*8 au(*),dgdxglob(2,nk,*)
!  
!     loop over all colums
      do kk=1,ndesi
         inode1=nodedesi(kk)
!
!        loop over all rows
         do jj=jq(kk),jq(kk+1)-1
            ipos=irow(jj)
            inode2=nodedesi(ipos)
!
!           filtered sensitivities
            do ii=1+nobjectstart,nobject
               if(objectset(1,ii)(4:13).eq.'MEMBERSIZE') cycle 
               if(objectset(1,ii)(4:9).eq.'PACKAGING') cycle            
               if(objectset(1,ii)(1:9).eq.'MAXGROWTH') cycle             
               if(objectset(1,ii)(1:12).eq.'MAXSHRINKAGE') cycle
               dgdxglob(2,inode2,ii)=dgdxglob(2,inode2,ii)+au(jj)
     &            *dgdxglob(1,inode1,ii)
            enddo
         enddo    
      enddo  
!     
      return        
      end
