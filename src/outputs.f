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
      subroutine outputs(inpc,textpart,jout,itpamp,istep,istat,n,iline,
     &  ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the *OUTPUT card in the input deck
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,ii,jout(2),joutl,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),itpamp,ier
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OUTPUT'
         write(*,*) '       *OUTPUT'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         ier=1
         return
      endif
!
      do ii=2,n
        if(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*OUTPUT %",ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(1)=joutl
              itpamp=0
           endif
        elseif(textpart(ii)(1:11).eq.'FREQUENCYF=') then
           read(textpart(ii)(12:21),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*OUTPUT %",ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(2)=joutl
              itpamp=0
           endif
        else
            write(*,*) 
     &             '*WARNING reading *NODE FILE or *EL FILE:' 
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*OUTPUT %")
        endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end
