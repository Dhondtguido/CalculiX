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
      subroutine calcfeasibledirection_gd(ndesi,nodedesi,dgdxglob,
     &   nactive,nobject,nk,gradproj,objectset)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(5,*)
!
      integer ndesi,nodedesi(*),nactive,nobject,nk,node,i 
!
      real*8 dgdxglob(2,nk,*),gradproj(3,*),dd,xi,cosphi,lambda1,
     &   lambda2,primaleig,dualeig,cosalpha1,cosalpha2
!
      xi=0.98
!
!     copy sensitivities of objective function to field gradproj(2,*)
!
      do i=1,ndesi
         node=nodedesi(i)
         gradproj(2,node)=dgdxglob(2,node,1)
      enddo  
!
!     Assembly of feasible direction
!
      do i=1,ndesi
         node=nodedesi(i)
         if(nobject.gt.1) then
            gradproj(3,node)=gradproj(2,node)-xi*gradproj(1,node)
         else
            gradproj(3,node)=gradproj(2,node)
         endif  
      enddo
!
!     Normalization of feasible direction
!
      if(nobject.gt.1) then
         dd=0.d0
         do i=1,ndesi
            node=nodedesi(i)
            dd=dd+gradproj(3,node)**2
         enddo
         if(dd.le.0.d0) then
            dd=1.d0
         endif
         dd=dsqrt(dd) 
         do i=1,ndesi
            node=nodedesi(i)
            gradproj(3,node)=gradproj(3,node)/dd
        enddo
      endif
!
!     calculation of coeffitions alpha1 and alpha2
!
      cosphi=0
      do i=1,ndesi    
         node=nodedesi(i)
         cosphi=cosphi+gradproj(1,node)*gradproj(2,node)
      enddo
      lambda1=1-cosphi
      lambda2=1+cosphi
      if(lambda1.lt.1.0e-10) then
         lambda1=0.d0
      endif
      if(lambda2.lt.1.0e-10) then
         lambda2=0.d0
      endif
      primaleig=dsqrt(lambda1)
      dualeig=dsqrt(lambda2)
      cosalpha1=primaleig/sqrt(2.0)
      cosalpha2=dualeig/sqrt(2.0)
!      
      write(5,*) ''
      write(5,*) ''
      write(5,*) '  #######################################
     &#########################'
      write(5,*) '  S I N G U L A R   V A L U E   
     &D E C O M P O S I T I O N'
      write(5,*) ''
      write(5,'(3x,a18,e14.7)') 'PRIMAL EIGENVALUE: ', primaleig
      write(5,'(3x,a18,e14.7)') 'DUAL EIGENVALUE:   ', dualeig
      write(5,'(3x,a18,e14.7)') 'COS ALPHA1:        ', cosalpha1
      write(5,'(3x,a18,e14.7)') 'COS ALPHA2:        ', cosalpha2
      write(5,*) ''
      write(5,*) '  #######################################
     &#########################'
!
      return        
      end
