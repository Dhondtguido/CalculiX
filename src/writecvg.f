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
      subroutine writecvg(istep,iinc,icutb,iit,ne,ne0,ram,qam,cam,uam,
     &  ithermal)
!
      implicit none
!
!     writes convergence information in the .cvg-file
!
      integer istep,iinc,iit,ne,ne0,ithermal(*),icutb
!
      real*8 ram(*),qam(*),cam(*),uam(*),residforce,corrdisp,
     &  residflux,corrtemp
!
      if(ithermal(1).eq.2) then
         residforce=0.d0
         corrdisp=0.d0
      else
         if(dabs(qam(1)).lt.1.d-30) then
            if(dabs(ram(1)).lt.1.d-30) then
               residforce=1.d-30
            else
               residforce=1.d30
            endif
         else
            residforce=ram(1)/qam(1)*100.d0
         endif
!
         if(dabs(uam(1)).lt.1.d-30) then
            if(dabs(cam(1)).lt.1.d-30) then
               corrdisp=1.d-30
            else
               corrdisp=1.d30
            endif
         else
            corrdisp=cam(1)/uam(1)*100.d0
         endif
      endif
!     
      if(ithermal(1).le.1) then
         residflux=0.d0
         corrtemp=0.d0
      else
         if(dabs(qam(2)).lt.1.d-30) then
            if(dabs(ram(2)).lt.1.d-30) then
               residflux=1.d-30
            else
               residflux=1.d30
            endif
         else
            residflux=ram(2)/qam(2)*100.d0
         endif
         if(dabs(uam(2)).lt.1.d-30) then
            if(dabs(cam(2)).lt.1.d-30) then
               corrtemp=1.d-30
            else
               corrtemp=1.d30
            endif
         else
            corrtemp=cam(2)/uam(2)*100.d0
         endif
      endif
!     
      write(11,'(2x,i4,2x,i4,2x,i4,2x,i4,2x,i7,4(1x,e11.4))') istep,
     &  iinc,icutb+1,iit,ne-ne0,residforce,corrdisp,residflux,corrtemp
!
      flush(11)
!
      return
      end
