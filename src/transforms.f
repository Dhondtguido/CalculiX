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
      subroutine transforms(inpc,textpart,trab,ntrans,ntrans_,
     &     inotr,set,istartset,iendset,ialset,nset,istep,istat,
     &     n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *TRANSFORM
!
      implicit none
!
      real*8 trab(7,*)
!
      character*1 inpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer ntrans,ntrans_,istep,istat,n,key,i,j,k,inotr(2,*),
     &  istartset(*),iendset(*),ialset(*),nset,ipos,iline,ipol,
     &  inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier,id
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *TRANSFORM: *TRANSFORM should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      ntrans=ntrans+1
      if(ntrans.gt.ntrans_) then
         write(*,*) '*ERROR reading *TRANSFORM: increase ntrans_'
         ier=1
         return
      endif
!
      ipos=1
      noset(1:1)=' '
!
!     rectangular coordinate system: trab(7,norien)=1
!     cylindrical coordinate system: trab(7,norien)=-1
!     default is rectangular
!
      trab(7,ntrans)=1.d0
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NSET=') then
            noset=textpart(i)(6:85)
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
         elseif(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:6).eq.'C') then
               trab(7,ntrans)=-1.d0
            endif
         else
            write(*,*) 
     &        '*WARNING reading *TRANSFORM: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*TRANSFORM%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR reading *TRANSFORM: definition of a'
         write(*,*) '  transformation is not complete'
         call inputerror(inpc,ipoinpc,iline,
     &        "*TRANSFORM%",ier)
         return
      endif
!
      do i=1,6
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) trab(i,ntrans)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*TRANSFORM%",ier)
            return
         endif
      enddo
!
      if(noset(1:1).eq.' ') then
         write(*,*) '*ERROR reading *TRANSFORM: no node set defined'
         ier=1
         return
      endif
!         
      call cident81(set,noset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(noset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         noset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *TRANSFORM: node set ',noset
         write(*,*) '       has not yet been defined.'
         ier=1
         return
      endif
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            inotr(1,ialset(j))=ntrans
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               inotr(1,k)=ntrans
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

