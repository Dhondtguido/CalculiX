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
      subroutine clonesensitivities(nobject,nk,objectset,g0,dgdxglob)
!
      implicit none
!      
      character*81 objectset(5,*),rpname,entrybuffer(5)
!              
      integer nobject,nk,i,j,k,l,objective_id
!      
      real*8 g0(*),dgdxglob(2,nk,*),buffer        
!
!     as design responses are stored first inside objectset,
!     we start at the end and look for an entry with the same name.
!     If there is a previous entry with the same name, we copy the entries
!     within dgdxglob and g0. We only copy from the very first occurrence of 
!     the name.
!
      do i=nobject,1,-1
          rpname=objectset(5,i)(1:80)
          if(objectset(5,i)(81:81).eq.'O') then
             objective_id=i
          endif
          do j=1,i-1
             if(objectset(5,i)(81:81).ne.'G') then
                if(objectset(5,j)(1:80).eq.rpname(1:80)) then
                   g0(i)=g0(j)         
                   do k=1, nk
                      dgdxglob(1,k,i)=dgdxglob(1,k,j)
                      dgdxglob(2,k,i)=dgdxglob(2,k,j)
                   enddo            
                   exit
                endif
             endif
          enddo
      enddo    
!
!     its important for the objective to be placed first inside the
!     objectset. if objective_id is not equal to 1, we swap the entries
!     to the first place
!     
      if(objective_id.eq.1) then
         return
      endif
!      
      do i=1,5
!
!        we need to be careful when copying the objectset entries
!        as we don't want to move the filters away from i=2 but move
!        the kreisselmeier-steinhauser values.
!
         if(i.ne.2) then 
            entrybuffer(i)(1:81)=objectset(i,1)(1:81)
            objectset(i,1)(1:81)=objectset(i,objective_id)(1:81)
            objectset(i,objective_id)(1:81)=entrybuffer(i)(1:81)
         else
            entrybuffer(i)(40:81)=objectset(i,1)(40:81)
            objectset(i,1)(40:81)=objectset(i,objective_id)(40:81)
            objectset(i,objective_id)(40:81)=entrybuffer(i)(40:81)
         endif
      enddo
!      
      buffer=g0(1)
      g0(1)=g0(objective_id)
      g0(objective_id)=buffer
!         
      do k=1,nk
         buffer=dgdxglob(1,k,1)
         dgdxglob(1,k,1)=dgdxglob(1,k,objective_id)
         dgdxglob(1,k,objective_id)=buffer
!         
         buffer=dgdxglob(2,k,1)
         dgdxglob(2,k,1)=dgdxglob(2,k,objective_id)
         dgdxglob(2,k,objective_id)=buffer
      enddo
!     
      return
      end
         
