!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine changecontacttypes(inpc,textpart,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,iperturb,ipoinpc,mortar,ier,iexpl,
     &     nmethod)
!
!     reading the input deck: *CHANGE CONTACT TYPE
!
      implicit none
!     
      character*1 inpc(*)
      character*132 textpart(16)
!     
      integer istep,istat,n,i,key,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),iperturb(*),ipoinpc(0:*),iexpl,nmethod,
     &     mortar,ier
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *CHANGE CONTACT TYPE:'
        write(*,*) '       *CHANGE CONTACT TYPE can only be used'
        write(*,*) '       within a STEP'
        ier=1
        return
      endif
!
      if((nmethod.ne.4).or.(iperturb(1).le.1)) then
        write(*,*) '*ERROR reading *CHANGE CONTACT TYPE:'
        write(*,*) '       *CHANGE CONTACT TYPE can only be used'
        write(*,*) '       in a nonmodal dynamic step'
        ier=1
        return
      endif
!     
      do i=2,n
        if(textpart(i)(1:15).eq.'TONODETOSURFACE') then
          mortar=0
        elseif(textpart(i)(1:12).eq.'TOMASSLESS') then
          if(iexpl.le.1) then
            write(*,*) '*ERROR reading *CHANGE CONTACT TYPE:'
            write(*,*) '       *CHANGE CONTACT TYPE,TO MASSLESS'
            write(*,*) '       can only be used in explicit dynamics'
            ier=1
            return
          endif
          mortar=-1
        else
          write(*,*) '*WARNING reading *CHANGE CONTACT TYPE:'
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*CONTACT PAIR%")
        endif
      enddo
!
      write(*,*) '*INFO reading *CHANGE CONTACT TYPE:'
      write(*,*) '      actual contact type:'
      if(mortar.eq.-1) then
        write(*,*) '      MASSLESS'
      elseif(mortar.eq.0) then
        write(*,*) '      NODE TO SURFACE'
      elseif(mortar.eq.1) then
        write(*,*) '      SURFACE TO SURFACE'
      elseif(mortar.eq.2) then
        write(*,*) '      MORTAR'
      elseif(mortar.eq.3) then
        write(*,*) '      LINMORTAR'
      elseif(mortar.eq.4) then
        write(*,*) '      PGLINMORTAR'
      elseif(mortar.eq.5) then
        write(*,*) '      PGMORTAR'
      endif
      write(*,*)
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end
