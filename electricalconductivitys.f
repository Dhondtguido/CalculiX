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
      subroutine electricalconductivitys(inpc,textpart,alcon,nalcon,
     &  alzero,nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *ELECTRICAL CONDUCTIVITY
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nalcon(2,*),nmat,ntmat,ntmat_,istep,istat,n,
     &  ipoinpc(0:*),ier,
     &  i,ityp,key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 alcon(0:6,ntmat_,*),alzero(*)
!
      ntmat=0
      alzero(nmat)=0.d0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *ELECTRICAL CONDUCTIVITY:'
         write(*,*) '       *ELECTRICAL CONDUCTIVITY should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *ELECTRICAL CONDUCTIVITY:'
         write(*,*) '       *ELECTRICAL CONDUCTIVITY should be preceded'
         write(*,*) '       by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=1
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
               ityp=1
            endif
         else
            write(*,*) '*WARNING reading *ELECTRICAL CONDUCTIVITY:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*ELECTRICAL CONDUCTIVITY%")
         endif
      enddo
!
      nalcon(1,nmat)=ityp
!
      if(ityp.eq.1) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            nalcon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) '*ERROR reading *ELECTRICAL CONDUCTIVITY:'
               write(*,*) '       increase ntmat_'
               ier=1
               return
            endif
            do i=1,1
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 alcon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*ELECTRICAL CONDUCTIVITY%",ier)
                  return
               endif
            enddo
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &                 alcon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*ELECTRICAL CONDUCTIVITY%",ier)
               return
            endif
         enddo
      else
         write(*,*) '*ERROR reading *ELECTRICAL CONDUCTIVITY:'
         write(*,*) '       no anisotropy allowed'
         ier=1
         return
      endif
!
      return
      end

