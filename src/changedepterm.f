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
      subroutine changedepterm(ikmpc,ilmpc,nmpc,mpc,idofrem,idofins)
!
!     changes the dependent term in ikmpc and ilmpc for MPC mpc.
!
      implicit none
!
      integer ikmpc(*),ilmpc(*),nmpc,idofrem,idofins,id,k,mpc
!
!     remove MPC from ikmpc
!
      call nident(ikmpc,idofrem,nmpc,id)
      if(id.gt.0) then
         if(ikmpc(id).eq.idofrem) then
            do k=id+1,nmpc
               ikmpc(k-1)=ikmpc(k)
               ilmpc(k-1)=ilmpc(k)
            enddo
         else
            write(*,*) '*ERROR in changedepterm'
            write(*,*) '       ikmpc database corrupted'
            call exit(201)
         endif
      else
         write(*,*) '*ERROR in changedepterm'
         write(*,*) '       ikmpc database corrupted'
         call exit(201)
      endif
!
!     insert new MPC
!
      call nident(ikmpc,idofins,nmpc-1,id)
      if((id.gt.0).and.(ikmpc(id).eq.idofins)) then
         write(*,*) '*ERROR in changedepterm: dependent DOF'
         write(*,*) '       of nonlinear MPC cannot be changed'
         write(*,*) '       since new dependent DOF is already'
         write(*,*) '       used in another MPC'
         call exit(201)
      else
         do k=nmpc,id+2,-1
            ikmpc(k)=ikmpc(k-1)
            ilmpc(k)=ilmpc(k-1)
         enddo
         ikmpc(id+1)=idofins
         ilmpc(id+1)=mpc
      endif
!
      return
      end











