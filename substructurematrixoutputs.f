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
      subroutine substructurematrixoutputs(textpart,istep,
     &  inpc,istat,n,key,iline,ipol,inl,ipoinp,inp,jobnamec,ipoinpc,
     &  ier)
!
!     reading the input deck: *SUBSTRUCTURE MATRIX OUTPUT
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16),jobnamec(*)
!
      integer i,istep,n,istat,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),key,j,k,l,ipoinpc(0:*),ier
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *SUBSTRUCTURE MATRIX OUTPUT:'
         write(*,*) '       *SUBSTRUCTURE MATRIX OUTPUT can '
         write(*,*) '       only be used within a STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:13).eq.'STIFFNESS=YES') then
         elseif(textpart(i)(1:22).eq.'OUTPUTFILE=USERDEFINED') then
         elseif(textpart(i)(1:9).eq.'FILENAME=') then
            jobnamec(5)(1:123)=textpart(i)(10:132)
            jobnamec(5)(124:132)='      '
            loop2: do j=1,123
               if(jobnamec(5)(j:j).eq.'"') then
                  do k=j+1,123
                     if(jobnamec(5)(k:k).eq.'"') then
                        do l=k-1,123
                           jobnamec(5)(l:l)=' '
                           exit loop2
                        enddo
                     endif
                     jobnamec(5)(k-1:k-1)=jobnamec(5)(k:k)
                  enddo
                  jobnamec(5)(123:123)=' '
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

