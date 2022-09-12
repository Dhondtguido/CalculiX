!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine checkintegrity(nktet,ipofa,ifac,itetfa,ifatet)
!
      implicit none
!
      integer i,nktet,ipofa(*),ifac(4,*),indexef,itetfa(2,*),
     &  ifatet(4,*),iel1,iel2,j1,j2,j
!
      do i=1,nktet
        indexef=ipofa(i)
        do
          if(indexef.eq.0) exit
          if(itetfa(2,indexef).ne.0) then
            iel1=itetfa(1,indexef)
            iel2=itetfa(2,indexef)
            do j=1,4
              if(abs(ifatet(j,iel1)).eq.indexef) then
                j1=j
                exit
              endif
            enddo
            do j=1,4
              if(abs(ifatet(j,iel2)).eq.indexef) then
                j2=j
                exit
              endif
            enddo
            if(ifatet(j1,iel1).ne.-ifatet(j2,iel2)) then
              write(*,*) '*ERROR in checkintegrity:'
              write(*,*) '       element 1 ',iel1,' face 1 ',
     &             ifatet(j1,iel1)
              write(*,*) '       element 2 ',iel2,' face 2 ',
     &             ifatet(j2,iel2)
              call exit(201)
            endif
          endif
          indexef=ifac(4,indexef)
        enddo
      enddo
!
      return
      end
