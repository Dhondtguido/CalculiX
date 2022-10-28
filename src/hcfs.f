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
      subroutine hcfs(inpc,textpart,istep,istat,n,iline,ipol,inl,ipoinp,
     &     inp,ipoinpc,ier,jobnamec,mei,tincf)
!     
!     reading the input deck: *HCF
!     
      implicit none
!     
      logical input
!     
      character*1 inpc(*)
      character*132 textpart(16),jobnamec(*)
!     
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &     ipoinpc(0:*),ier,j,k,l,mei(4)
!     
      real*8 tincf
!     
      if(istep.ne.1) then
        write(*,*) '*ERROR reading *HCF:'
        write(*,*) '       *HCF can only be used'
        write(*,*) '       within the first STEP'
        ier=1
        return
      endif
!     
!     defaults 
!     
      mei(1)=0
      mei(2)=0
      mei(3)=0
      tincf=0.d0
      input=.false.
!     
      do i=2,n
        if(textpart(i)(1:5).eq.'MODE=') then
          read(textpart(i)(6:15),'(i10)',iostat=istat) mei(1)
        elseif(textpart(i)(1:12).eq.'MISSIONSTEP=') then
          read(textpart(i)(13:22),'(i10)',iostat=istat) mei(2)
        elseif(textpart(i)(1:9).eq.'MAXCYCLE=') then
          read(textpart(i)(10:19),'(i10)',iostat=istat) mei(3)
        elseif(textpart(i)(1:8).eq.'SCALING=') then
          read(textpart(i)(9:28),'(f20.0)',iostat=istat) tincf
        elseif(textpart(i)(1:6).eq.'INPUT=') then
          input=.true.
          jobnamec(5)(1:126)=textpart(i)(7:132)
          jobnamec(5)(127:132)='      '
          loop1: do j=1,126
          if(jobnamec(5)(j:j).eq.'"') then
            do k=j+1,126
              if(jobnamec(5)(k:k).eq.'"') then
                do l=k-1,126
                  jobnamec(5)(l:l)=' '
                  exit loop1
                enddo
              endif
              jobnamec(5)(k-1:k-1)=jobnamec(5)(k:k)
            enddo
            jobnamec(5)(126:126)=' '
          endif
        enddo loop1
      else
        write(*,*) 
     &       '*WARNING reading *HCF: parameter not recognized:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputwarning(inpc,ipoinpc,iline,
     &       "*HCF%")
      endif
      enddo
!     
!     check for the INPUT parameter
!
      if(.not.input) then
        write(*,*) 
     &       '*ERROR reading *HCF: no input file specified:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputerror(inpc,ipoinpc,iline,
     &       "*HCF%")
      endif
!     
      if(mei(1).eq.0) then
        write(*,*) 
     &       '*ERROR reading *HCF: no mode specified:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputerror(inpc,ipoinpc,iline,
     &       "*HCF%")
      endif
!     
      if(mei(2).eq.0) then
        write(*,*) 
     &       '*ERROR reading *HCF: no mission step specified:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputerror(inpc,ipoinpc,iline,
     &       "*HCF%")
      endif
!     
      if(tincf.eq.0.d0) then
        write(*,*) '*WARNING reading *HCF: no scaling specified;'
        write(*,*) '         a default of 1. is taken:'
        write(*,*) '         ',
     &       textpart(i)(1:index(textpart(i),' ')-1)
        call inputwarning(inpc,ipoinpc,iline,
     &       "*HCF%")
        tincf=1.d0
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

