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
      subroutine skip(nset,nalset,nload,nbody,
     &     nforc,nboun,nk,ne,nkon,
     &     mi,nmpc,mpcend,nmat,ntmat_,npmat_,ncmat_,norien,ntrans,nam,
     &     nprint,nlabel,ncs_,ne1d,ne2d,infree,nmethod,
     &     iperturb,nener,ithermal,nstate_,iprestr,mcs,ntie,
     &     nslavs,nprop,mortar,ifacecount,nintpoint,nef,nheading_,
     &     nfc,ndc)
!     
      implicit none
!     
      integer nset,nalset,nload,nforc,nboun,nk,ne,nkon,mi(*),
     &     nmpc,mpcend,nmat,ntmat_,npmat_,ncmat_,norien,ntrans,nam,
     &     nprint,nlabel,ncs_,ne1d,ne2d,infree(4),i,mt,nprop,mortar,
     &     nmethod,iperturb(*),nener,ithermal(*),nstate_,iprestr,i4,
     &     maxamta,mcs,ntie,nbody,nslavs,nintpoint,ifacecount,nef,
     &     nheading_,nfc,ndc
!     
      character*1 c1
      character*3 c3
      character*6 c6
      character*8 c8
      character*20 c20
      character*66 c66
      character*80 c80
      character*81 c81
      character*87 c87
!     
      real*8 r8
!     
      mt=mi(2)+1
!     
!     skipping the next entries
!     
!     
!     sets
!     
      read(15)(c81,i=1,nset)
      read(15)(i4,i=1,nset)
      read(15)(i4,i=1,nset)
      do i=1,nalset
        read(15)i4
      enddo
!
!     header lines
!
      read(15)(c66,i=1,nheading_)
!     
!     mesh
!     
      read(15)(r8,i=1,3*nk)
      read(15)(i4,i=1,nkon)
      read(15)(i4,i=1,ne)
      read(15)(c8,i=1,ne)
!     
!     single point constraints
!     
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(c1,i=1,nboun)
      read(15)(r8,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      if(nam.gt.0) read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(r8,i=1,nboun)
!     
!     multiple point constraints
!     
      read(15)(i4,i=1,nmpc)
      read(15)(c20,i=1,nmpc)
      read(15)(i4,i=1,nmpc)
      read(15)(i4,i=1,nmpc)
      read(15)(r8,i=1,nmpc)
      read(15)(i4,i=1,3*mpcend)
      read(15)(r8,i=1,mpcend)
!     
!     force constraints
!
      if(nfc.gt.0) then
        read(15)(r8,i=1,7*nfc)
        read(15)(i4,i=1,ndc)
        read(15)(r8,i=1,12*ndc)
      endif
!     
!     point forces
!     
      read(15)(i4,i=1,2*nforc)
      read(15)(i4,i=1,nforc)
      read(15)(r8,i=1,nforc)
      read(15)(i4,i=1,nforc)
      read(15)(i4,i=1,nforc)
      if(nam.gt.0) read(15)(i4,i=1,nforc)
      read(15)(r8,i=1,nforc)
!     
!     distributed loads
!     
      read(15)(i4,i=1,2*nload)
      read(15)(c20,i=1,nload)
      read(15)(r8,i=1,2*nload)
      if(nam.gt.0) read(15)(i4,i=1,2*nload)
      read(15)(r8,i=1,2*nload)
      read(15)(c81,i=1,nbody)
      read(15)(i4,i=1,3*nbody)
      read(15)(r8,i=1,7*nbody)
      read(15)(r8,i=1,7*nbody)
!     
!     prestress
!     
      if(iprestr.gt.0) read(15) (r8,i=1,6*mi(1)*ne)
!     
!     labels
!     
      read(15)(c6,i=1,nprint)
      read(15)(c81,i=1,nprint)
      read(15)(c87,i=1,nlabel)
!     
!     elastic constants
!     
      read(15)(r8,i=1,(ncmat_+1)*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
!     
!     density
!     
      read(15)(r8,i=1,2*ntmat_*nmat)
      read(15)(i4,i=1,nmat)
!     
!     specific heat
!     
      read(15)(r8,i=1,4*ntmat_*nmat)
      read(15)(i4,i=1,nmat)
!     
!     conductivity
!     
      read(15)(r8,i=1,7*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
!     
!     expansion coefficients
!     
      read(15)(r8,i=1,7*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
      read(15)(r8,i=1,nmat)
!     
!     physical constants
!     
      read(15)(r8,i=1,14)
!     
!     plastic data
!     
      if(npmat_.ne.0)then
        read(15)(r8,i=1,(2*npmat_+1)*ntmat_*nmat)
        read(15)(i4,i=1,(ntmat_+1)*nmat)
        read(15)(r8,i=1,(2*npmat_+1)*ntmat_*nmat)
        read(15)(i4,i=1,(ntmat_+1)*nmat)
      endif
!     
!     material orientation
!     
      if(norien.ne.0)then
        read(15)(c80,i=1,norien)
        read(15)(r8,i=1,7*norien)
        read(15)(i4,i=1,mi(3)*ne)
      endif
!     
!     fluid section properties
!     
      if(nprop.ne.0) then
        read(15)(i4,i=1,ne)
        read(15)(r8,i=1,nprop)
      endif
!     
!     transformations
!     
      if(ntrans.ne.0)then
        read(15)(r8,i=1,7*ntrans)
        read(15)(i4,i=1,2*nk)
      endif
!     
!     amplitudes
!     
      if(nam.gt.0)then
        read(15)(c80,i=1,nam)
        read(15)(i4,i=1,3*nam-1)
        maxamta=2*i4
        read(15)i4
        read(15)(r8,i=1,maxamta)
      endif
!     
!     temperatures
!     
      if(ithermal(1).gt.0)then
        read(15)(r8,i=1,nk)
        read(15)(r8,i=1,nk)
        if((ne1d.gt.0).or.(ne2d.gt.0))then
          read(15)(r8,i=1,2*nk)
          read(15)(r8,i=1,2*nk)
        endif
        if(nam.gt.0) read(15)(i4,i=1,nk)
        read(15)(r8,i=1,nk)
      endif
!     
!     materials
!     
      read(15)(c80,i=1,nmat)
      read(15)(i4,i=1,mi(3)*ne)
!     
!     temperature, displacement, static pressure, velocity and acceleration
!     
      read(15)(r8,i=1,mt*nk)
      if((nmethod.eq.4).or.((nmethod.eq.1).and.(iperturb(1).ge.2))) 
     &     then
        read(15)(r8,i=1,mt*nk)
      endif
!     
!     CFD results at the element centers
!     
      if(nef.gt.0) then
        read(15)(r8,i=1,8*nef)
        read(15)(r8,i=1,8*nef)
        read(15)(r8,i=1,8*nef)
      endif
!     
!     1d and 2d elements
!     
      if((ne1d.gt.0).or.(ne2d.gt.0))then
        read(15)(i4,i=1,2*nkon)
        read(15)(r8,i=1,infree(1))
        read(15)(i4,i=1,infree(2))
        read(15)(r8,i=1,mi(3)*nkon)
        read(15)(r8,i=1,2*ne)
        read(15)(i4,i=1,infree(4))
        read(15)(i4,i=1,3*(infree(3)-1))
        read(15)(i4,i=1,infree(4))
        read(15)(i4,i=1,2*infree(4))
      endif
!     
!     tie constraints
!     
      if(ntie.gt.0) then
        read(15)(c81,i=1,3*ntie)
        read(15)(r8,i=1,4*ntie)
      endif
!     
!     cyclic symmetry
!     
      if(ncs_.gt.0)then
        read(15)(i4,i=1,ncs_)
      endif
      if(mcs.gt.0) then
        read(15)(r8,i=1,17*mcs)
      endif
!     
!     integration point variables
!     
      read(15)(r8,i=1,6*mi(1)*ne)
      read(15)(r8,i=1,6*mi(1)*ne)
      if(nener.eq.1) then
        if(mortar.ne.1) then
          read(15)(r8,i=1,2*mi(1)*(ne+nslavs))
        else
          read(15)(r8,i=1,2*mi(1)*ne)
        endif
      endif
      if(nstate_.gt.0) then
        if(mortar.ne.1) then
          read(15)(r8,i=1,nstate_*mi(1)*(ne+nslavs))
        else
          read(15)(r8,i=1,nstate_*mi(1)*(ne+nintpoint))
        endif
      endif
!     
!     face-to-face penalty contact variables
!     
      if(mortar.eq.1) then
        read(15) (i4,i=1,2*ifacecount+2)
        read(15) (r8,i=1,3*nintpoint)
        read(15) (r8,i=1,3*9*ifacecount)
      endif
!     
!     control parameters
!     
      read(15) (r8,i=1,52)
      read(15) (r8,i=1,2)
      read(15) c3
      read(15) r8
!     
!     restart parameters
!     
      read(15) (i4,i=1,2)
!     
      return
      end
