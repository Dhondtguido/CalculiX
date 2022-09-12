!
!      CalculiX - A 3-dimensional finite element program
!               Copyright (C) 1998-2015 Guido Dhondt
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU General Public License as
!      published by the Free Software Foundation(version 2);
!     
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of 
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!      GNU General Public License for more details.
!
!      You should have received a copy of the GNU General Public License
!      along with this program; if not, write to the Free Software
!      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
       subroutine writemaccs(mac,nev,nm)
!
!      writes the results of MAC-caculation in *_mac.dat
!
!      nm is the nodal diameter in case of complex frequency
!      nev is the number of eigenvectors
!      mac contains the MAC-Values
!
       implicit none
!
       integer i,j,nev,nm(*)
       real*8 mac(nev,*)
!      
       write(5,*)
       write(5,*)'    Modal Assurance Criterium'
       write(5,*)'       Nodal Diameter',nm(1)
!
       do i=1,nev
          write(5,100) (mac(i,j),j=1,nev)
       enddo
!
 100   format(15(1x,e11.4))
       return
       end
