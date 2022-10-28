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
!     calculating h0 on the interface between A or A-V domains and 
!     the phi-domain. At the start of the routine h0 is available in
!     the complete phi-domain. Using the MPC's developed in tiedcontact
!     to tie the phi-values at the border of the A or A-V domains to
!     those in the phi-domain, the h0 values are calculated in a
!     similar way
!
      subroutine calch0interface(nmpc,ipompc,nodempc,coefmpc,h0)
!
      implicit none
!
      integer i,j,nmpc,ist,ipompc(*),ndir,nodempc(3,*),node,index
!
      real*8 h0(3,*),coefmpc(*),h0l(3)
!
      do i=1,nmpc
         ist=ipompc(i)
         if(ist.gt.0) then
            ndir=nodempc(2,ist)
!
!           looking for MPC's tying phi between the A or A-V
!           domains and the phi-domain
!
            if(ndir.eq.5) then
               node=nodempc(1,ist)
               index=nodempc(3,ist)
               do j=1,3
                  h0l(j)=0.d0
               enddo
               if(index.ne.0) then
                  do
                     do j=1,3
                        h0l(j)=h0l(j)-coefmpc(index)*
     &                       h0(j,nodempc(1,index))
                     enddo
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               do j=1,3
                  h0(j,node)=h0l(j)/coefmpc(ist)
               enddo
            endif
         endif
      enddo
!
      return
      end
