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
      subroutine materialramping(nelcon,elcon,ncmat_,ntmat_,nmat,iramp,
     &     xramp,idel,xdel,nk,mi,nactdof,b,iponoel,inoel,idivergence,
     &     ipkon)
!     
      implicit none
!     
      integer nelcon(2,*),ncmat_,ntmat_,nmat,i,j,k,iramp,idel,
     &     imat,nk,mi(*),nactdof(0:mi(2),*),iponoel(*),inoel(2,*),
     &     idivergence,index,nelem,ipkon(*)
!     
      real*8 elcon(0:ncmat_,ntmat_,*),b(*),xdel,xramp
!
      if(idivergence.eq.1) then
!
!       decreasing the elastic material properties
!
        if(idel.eq.0) then
          do imat=1,nmat
            if(nelcon(1,imat).eq.2) then
              do j=1,nelcon(2,imat)
                elcon(1,j,imat)=elcon(1,j,imat)/xramp
              enddo
            elseif(nelcon(1,imat).gt.0) then
              do j=1,nelcon(2,imat)
                do k=1,nelcon(1,imat)
                  elcon(k,j,imat)=elcon(k,j,imat)/xramp
                enddo
              enddo
            endif
          enddo
        else
!
!       resetting the elastic material properties;
!
          do imat=1,nmat
            if(nelcon(1,imat).eq.2) then
              do j=1,nelcon(2,imat)
                elcon(1,j,imat)=elcon(1,j,imat)*xramp**iramp
              enddo
            elseif(nelcon(1,imat).gt.0) then
              do j=1,nelcon(2,imat)
                do k=1,nelcon(1,imat)
                  elcon(k,j,imat)=elcon(k,j,imat)*xramp**iramp
                enddo
              enddo
            endif
          enddo
!
!       deleting elements belonging to nodes with a too large residual force
!
          do i=1,nk
            do j=0,mi(2)
              if(nactdof(j,i).gt.0) then
                if(dabs(b(nactdof(j,i))).gt.xdel) then
                  index=iponoel(i)
                  do
                    nelem=inoel(1,index)
                    if(ipkon(nelem).gt.-1) then
                      ipkon(nelem)=-2-ipkon(nelem)
                    endif
                    index=inoel(2,index)
                    if(index.eq.0) exit
                  enddo
                endif
              endif
            enddo
          enddo
        endif
      else
!
!       resetting the elastic material properties;
!
        do imat=1,nmat
          if(nelcon(1,imat).eq.2) then
            do j=1,nelcon(2,imat)
              elcon(1,j,imat)=elcon(1,j,imat)*xramp
            enddo
          elseif(nelcon(1,imat).gt.0) then
            do j=1,nelcon(2,imat)
              do k=1,nelcon(1,imat)
                elcon(k,j,imat)=elcon(k,j,imat)*xramp
              enddo
            enddo
          endif
        enddo
        iramp=iramp-1
      endif
!     
      return
      end
