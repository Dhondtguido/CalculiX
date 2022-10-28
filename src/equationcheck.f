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
!     Checking for zero columns in the gas matrix
!     
      subroutine equationcheck(ac,nteq,nactdog,itg,ntg,nacteq,
     &  network)
!     
      implicit none
!     
      integer nactdog(0:3,*),nteq,itg(*),ntg,i,j,k,l,node,
     &  nacteq(0:3,*),network
!     
      real*8 ac(nteq,*)
!   
!     checking the columns
!
      loop1: do j=1,nteq
         do i=1,nteq
            if(dabs(ac(i,j)).gt.0.d0) cycle loop1
         enddo
         do k=1,ntg
            node=itg(k)
            do l=0,2
               if(nactdog(l,node).eq.j) then
                  if(l.eq.0) then
                     write(*,*)      '*INFO in equationcheck: temperatur
     &e in node ',node,     ' cannot be determined: probably no incoming
     & mass flow'
                  elseif(l.eq.2) then
                     write(*,*)      '*INFO in equationcheck: pressure i 
     &n node ',node,     ' cannot be determined: all incoming elements a
     &re probably critical'
                  endif
                  cycle loop1
               endif
            enddo
         enddo
      enddo loop1
!  
!     checking the rows
! 
      loop2: do i=1,nteq
         do j=1,nteq
            if(dabs(ac(i,j)).gt.0.d0) cycle loop2
         enddo
         do k=1,ntg
            node=itg(k)
            do l=0,2
               if(nacteq(l,node).eq.i) then
                  if(l.eq.0) then
                     write(*,*)      '*INFO in equationcheck: energy equ
     &ation in node ',node,     ' is identically zero: probably no incom
     &ing mass flow'
                  elseif(l.eq.2) then
                     write(*,*)      '*INFO in equationcheck: element eq 
     &uation in node',node, ' is identically zero: the element is probab
     &ly critical'
                  endif
                  cycle loop2
               endif
            enddo
         enddo
      enddo loop2
!     
      if(network.le.2) then
         write(*,*) '*ERROR in equationcheck: singular system in'
         write(*,*) '       thermal network'
         call exit(201)
      endif
!
      return
      end
      
      
