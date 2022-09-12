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
      subroutine calinput_rfn(co,filab,set,istartset,iendset,ialset,
     &     nset,nset_,nalset,nalset_,mi,kon,ipkon,lakon,nkon,ne,ne_,
     &     iponor,xnor,istep,ipoinp,inp,iaxial,ipoinpc,network,
     &     nlabel,iuel,nuel_,ielmat,inpc,iperturb,iprestr,nk,nk_,ntie,
     &     tieset,iparentel,tietol)
!     
      implicit none
!
      logical solid,out3d
!     
      character*1 inpc(*)
      character*8 lakon(*)
      character*81 set(*),tieset(3,*)
      character*87 filab(*)
      character*132 textpart(16)
!     
      integer kon(*),i,istartset(*),iendset(*),ialset(*),nset,nset_,
     &     nalset,nalset_,mi(*),ipkon(*),nkon,ne,ne_,ixfree,ielmat(*),
     &     iponor(2,*),istep,istatrfn,n,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),iaxial,ipoinpc(0:*),network,nlabel,iuel,nuel_,ier,
     &     iperturb(*),iprestr,key,nk,nk_,ntie,iparentel(*)
!     
      real*8 co(3,*),xnor(*),tietol(4,*)
!     
      integer nentries
      parameter(nentries=18)
!
      ier=0
!     
      do i=1,nentries
        if(ipoinp(1,i).ne.0) then
          ipol=i
          inl=ipoinp(1,i)
          iline=inp(1,inl)-1
          exit
        endif
      enddo
!     
      call getnewline(inpc,textpart,istatrfn,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      loop: do
!     
        if(istatrfn.lt.0) exit
!     
        if(textpart(1)(1:8).eq.'*ELEMENT') then
          call elements(inpc,textpart,kon,ipkon,lakon,nkon,
     &         ne,ne_,set,istartset,iendset,ialset,nset,nset_,nalset,
     &         nalset_,mi(1),ixfree,iponor,xnor,istep,istatrfn,n,iline,
     &         ipol,inl,ipoinp,inp,iaxial,ipoinpc,solid,
     &         network,filab,nlabel,out3d,iuel,nuel_,ier,iparentel)
!     
        elseif(textpart(1)(1:12).eq.'*MODELCHANGE') then
          call modelchanges(inpc,textpart,tieset,istatrfn,n,iline,
     &         ipol,inl,ipoinp,inp,ntie,ipoinpc,istep,ipkon,nset,
     &         istartset,iendset,set,ialset,ne,mi,ielmat,iprestr,
     &         iperturb,ier,tietol)
!     
        elseif(textpart(1)(1:5).eq.'*NODE') then
          call nodes(inpc,textpart,co,nk,nk_,set,istartset,iendset,
     &         ialset,nset,nset_,nalset,nalset_,istep,istatrfn,n,iline,
     &         ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
!     check for zero-character-lines?
!     
        elseif(ipoinpc(iline-1).eq.ipoinpc(iline)) then
          call getnewline(inpc,textpart,istatrfn,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
!     
        else
          write(*,*) '*WARNING in calinput_rfn. Card image cannot be',
     &         ' interpreted:'
          call inputwarning(inpc,ipoinpc,iline,
     &         "the input file%")
          call getnewline(inpc,textpart,istatrfn,n,key,iline,ipol,inl,
     &         ipoinp,inp,ipoinpc)
        endif
!     
        if(ier.eq.1) then
          do
            call getnewline(inpc,textpart,istatrfn,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if(key.eq.1) exit
          enddo
          ier=2
        endif
!     
      enddo loop
!     
      if(ier.ge.1) then
        write(*,*) '*ERROR in calinput: at least one fatal'
        write(*,*) '       error message while reading the'
        write(*,*) '       input deck: CalculiX stops.'
        write(*,*)
        call exit(201)
      endif
!     
      return
      end
