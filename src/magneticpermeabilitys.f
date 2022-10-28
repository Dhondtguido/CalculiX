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
      subroutine magneticpermeabilitys(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *MAGNETIC PERMEABILITY
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ityp,ncmat_,irstrt(*),iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),idomain,ier
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *MAGNETIC PERMEABILITY:'
         write(*,*) '       *MAGNETIC PERMEABILITY should be placed'
         write(*,*) '        before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *MAGNETIC PERMEABILITY:'
         write(*,*) '       *MAGNETIC PERMEABILITY should be preceded'
         write(*,*) '       by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=2
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
               ityp=2
            endif
            exit
         else
            write(*,*) '*WARNING reading *MAGNETIC PERMEABILITY:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MAGNETIC PERMEABILITY%")
         endif
      enddo
!
      nelcon(1,nmat)=ityp
!
      if(ityp.eq.2) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nelcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *MAGNETIC PERMEABILITY:'
               write(*,*) '       increase ntmat_'
               ier=1
               return
            endif
!
            read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &           elcon(1,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MAGNETIC PERMEABILITY%",ier)
               return
            endif
!
            read(textpart(2)(1:10),'(i10)',iostat=istat) idomain
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MAGNETIC PERMEABILITY%",ier)
               return
            endif
            elcon(2,ntmat,nmat)=idomain+0.5d0
!
            if(textpart(3)(1:1).ne.' ') then
               read(textpart(3)(1:20),'(f20.0)',iostat=istat)
     &                   elcon(0,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*MAGNETIC PERMEABILITY%",ier)
                  return
               endif
            else
               elcon(0,ntmat,nmat)=0.d0
            endif
         enddo
      else
         write(*,*) '*ERROR reading *MAGNETIC PERMEABILITY:'
         write(*,*) '       no anisotropy allowed'
         ier=1
         return
      endif
!
      if(ntmat.eq.0) nelcon(1,nmat)=0
!
      return
      end

