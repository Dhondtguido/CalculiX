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
      subroutine correctem(stx,stn,
     &     prlab,nprint,ne,ipkon,lakon,elcon,ncmat_,ntmat_,nk,om,
     &     filab,mi,ielmat)
!
!     correct the real and imaginary part of the electric field
!     for a harmonic current      
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*),lakonl
      character*87 filab(*)
!     
      integer nprint,ne,ipkon(*),ncmat_,ntmat_,nk,i,ii,mint3d,k,
     &     kk,mi(*),imat,ielmat(mi(3),*)
!     
      real*8 stx(6,mi(1),*),stn(6,*),elcon(0:ncmat_,ntmat_,*),om,temp
!     
!     correction of the electric field (real and imaginary part)
!     at the integration points (only needed if the Joule heating
!     power is prequested by a *el print card)
!     
      do ii=1,nprint
        if(prlab(ii)(1:4).eq.'EBHE') then
          do i=1,ne
!     
            if(ipkon(i).lt.0) cycle
!
            imat=ielmat(1,i)
            if(int(elcon(2,1,imat)).ne.2) cycle
!     
            lakonl(1:8)=lakon(i)(1:8)
!     
            if(lakonl(4:5).eq.'8R') then
              mint3d=1
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R')) then
              mint3d=8
            elseif(lakonl(4:4).eq.'2') then
              mint3d=27
            elseif(lakonl(4:5).eq.'10') then
              mint3d=4
            elseif(lakonl(4:4).eq.'4') then
              mint3d=1
            elseif(lakonl(4:5).eq.'15') then
              mint3d=9
            elseif(lakonl(4:4).eq.'6') then
              mint3d=2
            endif
!     
            do kk=1,mint3d
!     
!     correct the electric field E in A-V domain
!     
              do k=1,3
                temp=stx(k,kk,i)
                stx(k,kk,i)=-om*stx(k,kk,i+ne)
                stx(k,kk,i+ne)=om*temp
              enddo
!     
            enddo
          enddo
          exit
        endif
      enddo
!     
!     correction of the electric field (real and imaginary part)
!     at the nodes (only needed of the electric or induction
!     field is requested by *el file or *element output)
!     
      if((filab(44)(1:4).eq.'EMFE').or.(filab(45)(1:4).eq.'EMFB'))
     &     then
        do i=1,nk
          do k=1,3
            temp=stn(k,i)
            stn(k,i)=-om*stn(k,i+nk)
            stn(k,i+nk)=om*temp
          enddo
        enddo
      endif
!     
      return
      end
