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
      subroutine timepointss(inpc,textpart,amname,amta,namta,nam,
     &  nam_,namtot_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,namtot,ier)
!
!     reading the input deck: *AMPLITUDE
!
      implicit none
!
      character*1 inpc(*)
      character*80 amname(*)
      character*132 textpart(16)
!
      integer namta(3,*),nam,nam_,istep,istat,n,key,i,namtot,
     &  namtot_,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipos,
     &  ipoinpc(0:*),nttp,ier
!
      logical igen
!
!
      real*8 amta(2,*),x,tpmin,tpmax,tpinc
!
      igen=.false.

      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *TIME POINTS: *AMPLITUDE should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      nam=nam+1
      if(nam.gt.nam_) then
         write(*,*) '*ERROR reading *TIME POINTS: increase nam_'
         ier=1
         return
      endif
      namta(3,nam)=nam
      amname(nam)='
     &                           '
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            amname(nam)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*)
     &           '*ERROR reading *TIME POINTS: amplitude name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       amplitude name:',textpart(i)(1:132)
               ier=1
               return
            endif
         elseif(textpart(i)(1:14).eq.'TIME=TOTALTIME') then
            namta(3,nam)=-nam
         elseif(textpart(i)(1:8).eq.'GENERATE') then
            igen=.true.
         else
            write(*,*) 
     &        '*WARNING reading *TIME POINTS: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*TIME POINTS%")
         endif
      enddo
!
      if(amname(nam).eq.'                                               
     &                                 ') then
         write(*,*) '*ERROR reading *TIME POINTS: Amplitude has no name'
         call inputerror(inpc,ipoinpc,iline,
     &        "*TIME POINTS%",ier)
         return
      endif
!
c      if(nam.eq.1) then
c         namtot=0
c      else
c         namtot=namta(2,nam-1)
c      endif
      namta(1,nam)=namtot+1
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         if(.not.igen)then
            do i=1,8
               if(textpart(i)(1:1).ne.' ') then  
                  namtot=namtot+1
                  if(namtot.gt.namtot_) then
                     write(*,*) 
     &                '*ERROR reading *TIME POINTS: increase namtot_'
                     ier=1
                     return
                  endif
                  read(textpart(i),'(f20.0)',iostat=istat) x
                  if(istat.gt.0) then
                     call inputerror(inpc,ipoinpc,iline,
     &                    "*TIME POINTS%",ier)
                     return
                  endif
                  amta(1,namtot)=x
                  namta(2,nam)=namtot
               else
                  exit
               endif
            enddo
         else
            read(textpart(1)(1:20),'(f20.0)',iostat=istat) tpmin
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TIME POINTS%",ier)
               return
            endif
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) tpmax
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TIME POINTS%",ier)
               return
            endif
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) tpinc
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*TIME POINTS%",ier)
               return
            endif
!            
            nttp=INT((tpmax-tpmin)/tpinc)
!            
            if(namtot+2+nttp.gt.namtot_) then
               write(*,*) '*ERROR in timepoints: increase namtot_'
               ier=1
               return
            endif
            amta(1,namtot+1)=tpmin
            do i=1,nttp
               amta(1,namtot+1+i)=tpmin+(i*tpinc)
            enddo
            namtot=namtot+2+nttp
            amta(1,namtot)=tpmax
            namta(2,nam)=namtot
         endif
         if(textpart(9)(1:1).ne.' ') then 
            write(*,*) '*WARNING reading *TIME POINTS:'
            write(*,*) '         only 8 entries per line allowed'
            write(*,*) '         9th entry and above will be discarded'
            call inputwarning(inpc,ipoinpc,iline,
     &"*TIME POINTS%")
         endif
      enddo
!
      if(namta(1,nam).gt.namta(2,nam)) then
         ipos=index(amname(nam),' ')
         write(*,*) 
     &    '*WARNING reading *TIME POINTS: *TIME POINTS definition ',
     &        amname(nam)(1:ipos-1)
         write(*,*) '         has no data points'
         nam=nam-1
c      else
c         call reorderampl(amname,namta,nam)
      endif
!
      return
      end

