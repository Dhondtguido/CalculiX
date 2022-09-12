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
      subroutine materialdata_em(elcon,nelcon,alcon,nalcon,
     &  imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!
      implicit none
!
!     determines the electric conductance and the magnetic
!     permeability for temperature t1l
!
      integer nelcon(2,*),nalcon(2,*),imat,k,ntmat_,nelconst,
     &  ncmat_,id,seven
!
      real*8 elcon(0:ncmat_,ntmat_,*),alcon(0:6,ntmat_,*),alpha(6),t1l,
     &  elconloc(21)
!
      seven=7
!     
      nelconst=nelcon(1,imat)
!     
!     calculating the electric conductance
!     
      call ident2(alcon(0,1,imat),t1l,nalcon(2,imat),seven,id)
      if(nalcon(2,imat).eq.0) then
         do k=1,6
            alpha(k)=0.d0
         enddo
         continue
      elseif(nalcon(2,imat).eq.1) then
         do k=1,nalcon(1,imat)
            alpha(k)=alcon(k,1,imat)
         enddo
      elseif(id.eq.0) then
         do k=1,nalcon(1,imat)
            alpha(k)=alcon(k,1,imat)
         enddo
      elseif(id.eq.nalcon(2,imat)) then
         do k=1,nalcon(1,imat)
            alpha(k)=alcon(k,id,imat)
         enddo
      else
         do k=1,nalcon(1,imat)
            alpha(k)=(alcon(k,id,imat)+
     &           (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &           (t1l-alcon(0,id,imat))/
     &           (alcon(0,id+1,imat)-alcon(0,id,imat)))
         enddo
      endif
!     
!     calculating the permeability
!     
      call ident2(elcon(0,1,imat),t1l,nelcon(2,imat),ncmat_+1,id)
      if(nelcon(2,imat).eq.0) then
         continue
      elseif(nelcon(2,imat).eq.1) then
         do k=1,nelconst
            elconloc(k)=elcon(k,1,imat)
         enddo
      elseif(id.eq.0) then
         do k=1,nelconst
            elconloc(k)=elcon(k,1,imat)
         enddo
      elseif(id.eq.nelcon(2,imat)) then
         do k=1,nelconst
            elconloc(k)=elcon(k,id,imat)
         enddo
      else
         do k=1,nelconst
            elconloc(k)=elcon(k,id,imat)+
     &           (elcon(k,id+1,imat)-elcon(k,id,imat))*
     &           (t1l-elcon(0,id,imat))/
     &           (elcon(0,id+1,imat)-elcon(0,id,imat))
         enddo
      endif
!     
      return
      end
