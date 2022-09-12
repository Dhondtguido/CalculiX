!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine resultsk(nk,nactdok,v,solk,solt,ipompc,nodempc,
     &  coefmpc,nmpc,mi)
!
!     calculates the turbulence correction (STEP 5) in the nodes
!
      implicit none
!
      integer ipompc(*),nodempc(3,*),nmpc,nk,nactdok(*),i,ist,
     &  node,ndir,index,mi(*)
!
      real*8 coefmpc(*),solk(*),v(0:mi(2),*),fixed_dispk,fixed_dispt,
     &  solt(*)
!
!     extracting the turbulence correction from the solution
!
      do i=1,nk
         if(nactdok(i).gt.0) then
            v(5,i)=solk(nactdok(i))
            v(6,i)=solt(nactdok(i))
         else
            v(5,i)=0.d0
            v(6,i)=0.d0
         endif
      enddo
!     
!     inserting the mpc information: it is assumed that the
!     temperature MPC's also apply to the turbulence
!     
c      do i=1,nmpc
c         ist=ipompc(i)
c         node=nodempc(1,ist)
c         ndir=nodempc(2,ist)
c         if(ndir.ne.0) cycle
c         index=nodempc(3,ist)
c         fixed_dispk=0.d0
c         fixed_dispt=0.d0
c         if(index.ne.0) then
c            do
c               fixed_dispk=fixed_dispk-coefmpc(index)*
c     &              vtu(1,nodempc(1,index))
c               fixed_dispt=fixed_dispt-coefmpc(index)*
c     &              vtu(2,nodempc(1,index))
c               index=nodempc(3,index)
c               if(index.eq.0) exit
c            enddo
c         endif
c         fixed_dispk=fixed_dispk/coefmpc(ist)
c         vtu(1,node)=fixed_dispk
c         fixed_dispt=fixed_dispt/coefmpc(ist)
c         vtu(2,node)=fixed_dispt
c      enddo
!
      return
      end
