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
      subroutine restartshort(nset,nload,nbody,nforc,nboun,nk,ne,
     &     nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &     ntrans,ncs_,namtot,ncmat_,memmpc_,ne1d,ne2d,nflow,
     &     set,meminset,rmeminset,jobnamec,irestartstep,icntrl,ithermal,
     &     nener,nstate_,ntie,nslavs,nkon,mcs,nprop,mortar,
     &     ifacecount,nintpoint,infree,nef,mpcend,nheading_,network,
     &     nfc,ndc,iprestr)
!     
!     istartset := meminset
!     iendset := rmeminset
!     
      implicit none
!     
      character*80 version
      character*81 set(*)
      character*132 fnrstrt,jobnamec(*)
!     
      integer istep,nset,nload,nforc,nboun,nk,ne,nmpc,nalset,
     &     nmat,ntmat_,npmat_,norien,nam,nprint,mi(*),ntrans,ncs_,
     &     nprop,namtot,ncmat_,memmpc_,ne1d,ne2d,nflow,infree(4),mortar,
     &     nmethod,iperturb(2),meminset(*),rmeminset(*),nintpoint,
     &     i,j,k,ipos,icntrl,nener,irestartstep,im0,im1,im2,mem,iact,
     &     istat,nkon,nlabel,iplas,ithermal(*),nstate_,iprestr,
     &     mcs,ntie,nbody,nslavs,ifacecount,iversion,nef,mpcend,
     &     maxlenmpc,nheading_,network,nfc,ndc
!     
      if(icntrl.eq.0) then
!     
!     this branch is called from readinput.c
!     its purpose is to read the value of nset
!     
!     determining the name of the restart file
!     
        ipos=index(jobnamec(1),char(0))
        fnrstrt(1:ipos-1)=jobnamec(1)(1:ipos-1)
        fnrstrt(ipos:ipos+3)=".rin"
        do i=ipos+4,132
          fnrstrt(i:i)=' '
        enddo
!     
!     opening the restart file
!     
        open(15,file=fnrstrt,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     &       err=15)
!     
        iversion=0
!     
        do
!     
          read(15,iostat=istat) version
          if(istat.lt.0) then
            if(irestartstep.eq.0) then
!     
!     reading the last step
!     
              irestartstep=istep
              close(15)
              open(15,file=fnrstrt,ACCESS='SEQUENTIAL',
     &             FORM='UNFORMATTED',err=15)
              read(15) version
            else
              write(*,*) '*ERROR in restartshort: requested step'
              write(*,*) '       is not in the restart file'
              call exit(201)
            endif
          endif
!     
          if(iversion.eq.0) then
            write(*,*)
            write(*,*) '*INFO: restart file ',fnrstrt
            write(*,*) '       has been opened for reading.'
            write(*,*) '       it was created with CalculiX ',version
            iversion=1
          endif
!     
          read(15)istep
!     
!     reading the number of sets
!     
          read(15)nset
!     
          if(istep.eq.irestartstep) exit
!     
          read(15)nalset
!     
!     load size
!     
          read(15)nload
          read(15)nbody
          read(15)nforc
          read(15)nboun
          read(15)nflow
!     
!     mesh size
!     
          read(15)nk
          read(15)ne
          read(15)nef
          read(15)nkon
          read(15)(mi(i),i=1,3)
!     
!     constraint size
!     
          read(15)nmpc
          read(15)mpcend
          read(15)maxlenmpc
          read(15)memmpc_
!     
!     coupling, distributing
!     
          read(15)nfc
          read(15)ndc
!     
!     material size
!     
          read(15)nmat
          read(15)ntmat_
          read(15)npmat_
          read(15)ncmat_
!     
!     property info
!     
          read(15)nprop
!     
!     transformation size
!     
          read(15)norien
          read(15)ntrans
!     
!     amplitude size
!     
          read(15)nam
          read(15)namtot
!     
!     print size
!     
          read(15)nprint
          read(15)nlabel
!     
!     tie size
!     
          read(15)ntie
!     
!     cyclic symmetry size
!     
          read(15)ncs_
          read(15)mcs
!     
!     1d and 2d element size
!     
          read(15)ne1d 
          read(15)ne2d 
          read(15)(infree(i),i=1,4)
!     
!     procedure info
!     
          read(15)nmethod
          read(15)(iperturb(i),i=1,2)
          read(15)nener
          read(15)iplas
          read(15)ithermal(1)
          read(15)nstate_
          read(15)nslavs
          read(15)iprestr
          read(15)mortar
          if(mortar.eq.1) then
            read(15)ifacecount
            read(15)nintpoint
          endif
          read(15)nheading_
          read(15)network
!     
!     skipping the next entries
!     
          call skip(nset,nalset,nload,nbody,
     &         nforc,nboun,nk,ne,nkon,
     &         mi,nmpc,mpcend,nmat,ntmat_,npmat_,ncmat_,norien,
     &         ntrans,nam,nprint,nlabel,ncs_,ne1d,ne2d,infree,
     &         nmethod,iperturb,nener,ithermal,nstate_,iprestr,
     &         mcs,ntie,nslavs,nprop,mortar,ifacecount,nintpoint,
     &         nef,nheading_,nfc,ndc)
!     
        enddo
!     
        close(15)
!     
        return
      endif
!     
!     determining the name of the restart file
!     
      ipos=index(jobnamec(1),char(0))
      fnrstrt(1:ipos-1)=jobnamec(1)(1:ipos-1)
      fnrstrt(ipos:ipos+3)=".rin"
      do i=ipos+4,132
        fnrstrt(i:i)=' '
      enddo
!     
!     opening the restart file
!     
      open(15,file=fnrstrt,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     &     err=15)
!     
      do
!     
        read(15,iostat=istat) version
        if(istat.lt.0) then
          if(irestartstep.eq.0) then
!     
!     reading the last step
!     
            irestartstep=istep
            close(15)
            open(15,file=fnrstrt,ACCESS='SEQUENTIAL',
     &           FORM='UNFORMATTED',err=15)
            read(15) version
          else
            write(*,*) '*ERROR in restartshort: requested step'
            write(*,*) '       is not in the restart file'
            call exit(201)
          endif
        endif
        read(15)istep
!     
!     set size
!     
        read(15)nset
        read(15)nalset
!     
!     load size
!     
        read(15)nload
        read(15)nbody
        read(15)nforc
        read(15)nboun
        read(15)nflow
!     
!     mesh size
!     
        read(15)nk
        read(15)ne
        read(15)nef
        read(15)nkon
        read(15)(mi(i),i=1,3)
!     
!     constraint size
!     
        read(15)nmpc
        read(15)mpcend
        read(15)maxlenmpc
        read(15)memmpc_
!     
!     coupling, distributing
!     
        read(15)nfc
        read(15)ndc
!     
!     material size
!     
        read(15)nmat
        read(15)ntmat_
        read(15)npmat_
        read(15)ncmat_
!     
!     property info
!     
        read(15)nprop
!     
!     transformation size
!     
        read(15)norien
        read(15)ntrans
!     
!     amplitude size
!     
        read(15)nam
        read(15)namtot
!     
!     print size
!     
        read(15)nprint
        read(15)nlabel
!     
!     tie size
!     
        read(15)ntie
!     
!     cyclic symmetry size
!     
        read(15)ncs_
        read(15)mcs
!     
!     1d and 2d element size
!     
        read(15)ne1d 
        read(15)ne2d 
        read(15)(infree(i),i=1,4)
!     
!     procedure info
!     
        read(15)nmethod
        read(15)(iperturb(i),i=1,2)
        read(15)nener
        read(15)iplas
        read(15)ithermal(1)
        read(15)nstate_
        read(15)nslavs
        read(15)iprestr
        read(15)mortar
        if(mortar.eq.1) then
          read(15)ifacecount
          read(15)nintpoint
        endif
        read(15)nheading_
        read(15)network
!     
        if(istep.eq.irestartstep) exit
!     
!     skipping the next entries
!     
        call skip(nset,nalset,nload,nbody,nforc,nboun,nk,ne,
     &       nkon,mi,nmpc,mpcend,nmat,ntmat_,npmat_,ncmat_,norien,
     &       ntrans,nam,nprint,nlabel,ncs_,ne1d,ne2d,infree,nmethod,
     &       iperturb,nener,ithermal,nstate_,iprestr,mcs,ntie,
     &       nslavs,nprop,mortar,ifacecount,nintpoint,nef,nheading_,
     &       nfc,ndc)
!     
      enddo
!     
!     sets
!     
      read(15)(set(i),i=1,nset)
!     
!     the contents of istartset is temporarily stored in meminset
!     
      read(15)(meminset(i),i=1,nset)
!     
!     the contents of iendset is temporarily stored in rmeminset
!     
      read(15)(rmeminset(i),i=1,nset)
!     
!     reordering the information of istartset, iendset and ialset
!     into meminset and rmeminset
!     
      iact=0
      do j=1,nalset
        if(iact.eq.0) then
          do k=1,nset
            if(meminset(k).eq.j) then
              meminset(k)=0
              mem=0
              iact=1
              exit
            endif
          enddo
          if(k.gt.nset) cycle
        endif
        mem=mem+1
        im2=im1
        im1=im0
        read(15) im0
        if(im0.gt.0) then
          meminset(k)=meminset(k)+1
        else
!     
!     im0<0 and two elements are already stored
!     
          meminset(k)=meminset(k)+(im2-im1)/im0-1
        endif
        if(rmeminset(k).eq.j) then
          iact=0
          rmeminset(k)=mem
!     
!     make set k ineligible in further iterations
!     
          meminset(k)=-meminset(k)
        endif
      enddo
!     
!     restore the sign of meminset
!     
      do k=1,nset
        meminset(k)=-meminset(k)
      enddo
!     
      close(15)
!     
      return
!     
 15   write(*,*) '*ERROR in restartshort: could not open file ',fnrstrt
      call exit(201)
      end
