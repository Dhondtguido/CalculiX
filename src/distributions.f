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
      subroutine distributions(inpc,textpart,orname,orab,norien,
     &     norien_,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier,
     &     set,istartset,iendset,ialset,nset,ne,ielorien,mi)
!
!     reading the input deck: *DISTRIBUTION
!
      implicit none
!
      character*1 inpc(*)
      character*80 orname(*),distname
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer norien,norien_,istep,istat,n,key,i,iline,ipol,inl,id,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),j,ier,istartset(*),
     &     iendset(*),ialset(*),nset,l,ne,mi(*),ielorien(mi(3),*),ipos
!
      real*8 orab(7,*)
!
!
!
      if(istep.gt.0) then
         write(*,*) 
     &       '*ERROR reading *DISTRIBUTION: *DISTRIBUTION should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            distname=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR reading *DISTRIBUTION: name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       distribution name:',textpart(i)(1:132)
               ier=1
               return
            endif
         elseif(textpart(i)(1:9).eq.'LOCATION=') then
         elseif(textpart(i)(1:6).eq.'TABLE=') then
         else
            write(*,*) 
     &       '*WARNING reading *DISTRIBUTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*DISTRIBUTION%")
         endif
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         norien=norien+1
         if(norien.gt.norien_) then
            write(*,*) '*ERROR reading *DISTRIBUTION: increase norien_'
            ier=1
            return
         endif
         orname(norien)=distname
!
!        reading the element number or element set
!        for the default (first line underneath *DISTRIBUTION) the
!        element number is zero
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if(l.gt.ne) then
               write(*,*) '*ERROR reading *DISTRIBUTION: element ',l
               write(*,*) '       is not defined'
               ier=1
               return
            endif
            if(l.gt.0) ielorien(1,l)=-norien
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
!
!           check for element set
!
            elset(ipos:ipos)='E'
c            do i=1,nset
c               if(set(i).eq.elset) exit
c            enddo
            call cident81(set,elset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(elset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
               elset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *DISTRIBUTION: element set '
               write(*,*) '       or facial surface ',elset
               write(*,*) '       has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTION%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  l=ialset(j)
                  ielorien(1,l)=-norien
               else
                  l=ialset(j-2)
                  do
                     l=l-ialset(j)
                     if(l.ge.ialset(j-1)) exit
                     ielorien(1,l)=-norien
                  enddo
               endif
            enddo
         endif
!
         do i=1,6
            read(textpart(i+1)(1:20),'(f20.0)',iostat=istat)
     &             orab(i,norien)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTION%",ier)
               return
            endif
         enddo
      enddo
!
      return
      end

