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
      subroutine viewfactors(textpart,iviewfile,istep,inpc,
     &     istat,n,key,iline,ipol,inl,ipoinp,inp,jobnamec,ipoinpc,ier,
     &     irestartread)
!
!     reading the input deck: *VIEWFACTOR 
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16),jobnamec(*)
!
      integer i,iviewfile,istep,n,istat,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),key,j,k,l,ipoinpc(0:*),ier,irestartread
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *VIEWFACTOR: *VIEWFACTOR can '
         write(*,*) '       only be used within a STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:4).eq.'READ') then
            if(iviewfile.eq.0) then
               iviewfile=-1
            elseif(iviewfile.gt.0) then
               write(*,*) '*ERROR reading *VIEWFACTOR: READ and WRITE/'
               write(*,*) '       WRITE ONLY are mutually exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%",ier)
               return
            endif
         elseif(textpart(i)(1:8).eq.'NOCHANGE') then
            if(istep.eq.1) then
               write(*,*) '*ERROR reading *VIEWFACTOR: NO CHANGE cannot'
               write(*,*) '       be used in the first step'
               call inputwarning(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%")
            elseif(irestartread.eq.1) then
               write(*,*) '*ERROR reading *VIEWFACTOR: NO CHANGE cannot'
               write(*,*) '       be used in the first step of a'
               write(*,*) '       restart calculation'
               call inputwarning(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%")
            elseif(iviewfile.le.0) then
               iviewfile=-2
            elseif(iviewfile.gt.0) then
               write(*,*) '*ERROR reading *VIEWFACTOR: NO CHANGE and'
               write(*,*) '       WRITE/WRITE ONLY are mutually'
               write(*,*) '       exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%",ier)
               return
            endif
         elseif(textpart(i)(1:9).eq.'WRITEONLY') then
            if(iviewfile.eq.0) then
               iviewfile=3
            elseif(iviewfile.lt.0) then
               write(*,*) '*ERROR reading *VIEWFACTOR: '
               write(*,*) '       WRITE ONLY and READ/NO CHANGE'
               write(*,*) '       are mutually exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%",ier)
               return
            endif
         elseif(textpart(i)(1:5).eq.'WRITE') then
            if(iviewfile.eq.0) then
               iviewfile=2
            elseif(iviewfile.lt.0) then
               write(*,*) '*ERROR reading *VIEWFACTOR: WRITE'
               write(*,*) '       and READ/NO CHANGE'
               write(*,*) '       are mutually exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*VIEWFACTOR%",ier)
               return
            endif
         elseif(textpart(i)(1:6).eq.'INPUT=') then
            jobnamec(2)(1:126)=textpart(i)(7:132)
            jobnamec(2)(127:132)='      '
            loop1: do j=1,126
               if(jobnamec(2)(j:j).eq.'"') then
                  do k=j+1,126
                     if(jobnamec(2)(k:k).eq.'"') then
                        do l=k-1,126
                           jobnamec(2)(l:l)=' '
                           exit loop1
                        enddo
                     endif
                     jobnamec(2)(k-1:k-1)=jobnamec(2)(k:k)
                  enddo
                  jobnamec(2)(126:126)=' '
               endif
            enddo loop1
         elseif(textpart(i)(1:7).eq.'OUTPUT=') then
            jobnamec(3)(1:125)=textpart(i)(8:132)
            jobnamec(3)(126:132)='      '
            loop2: do j=1,125
               if(jobnamec(3)(j:j).eq.'"') then
                  do k=j+1,125
                     if(jobnamec(3)(k:k).eq.'"') then
                        do l=k-1,125
                           jobnamec(3)(l:l)=' '
                           exit loop2
                        enddo
                     endif
                     jobnamec(3)(k-1:k-1)=jobnamec(3)(k:k)
                  enddo
                  jobnamec(3)(125:125)=' '
               endif
            enddo loop2
         else
            write(*,*) 
     &        '*WARNING reading *VIEWFACTOR: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*VIEWFACTOR%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

