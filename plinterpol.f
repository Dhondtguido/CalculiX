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
      subroutine plinterpol(plcon,nplcon,itemp,f,df,npmat_,ntmat_,
     &  imat,nelem,epl)
!
      implicit none
!
!     interpolation of isotropic or kinematic hardening data
!     input: hardening data plcon and nplcon, temperature itemp,
!       size parameters npmat_ and ntmat_, material number imat
!       and equivalent plastic strain at which the coefficients
!       are to be determined
!     output: hardening coefficient and its local derivative f and df
!
      integer npmat_,ntmat_,nplcon(0:ntmat_,*),itemp,ndata,imat,j,
     &  nelem
!
      real*8 plcon(0:2*npmat_,ntmat_,*),f,df,epl
!
!
!
      ndata=nplcon(itemp,imat)
!
      do j=1,ndata
         if(epl.lt.plcon(2*j,itemp,imat)) exit
      enddo
!
      if((j.eq.1).or.(j.gt.ndata)) then
         if(j.eq.1) then
            f=plcon(1,itemp,imat)
            df=0.d0
         else
            f=plcon(2*ndata-1,itemp,imat)
            df=0.d0
         endif
         write(*,*) '*WARNING in plinterpol: plastic strain ',epl
         write(*,*) '         outside material plastic strain range'
         write(*,*) '         in element ',nelem,' and material ',imat
         write(*,*) '         for temperature ',plcon(0,itemp,imat)
      else
         df=(plcon(2*j-1,itemp,imat)-plcon(2*j-3,itemp,imat))/
     &      (plcon(2*j,itemp,imat)-plcon(2*j-2,itemp,imat))
         f=plcon(2*j-3,itemp,imat)+
     &      df*(epl-plcon(2*j-2,itemp,imat))
      endif
!
      return
      end
