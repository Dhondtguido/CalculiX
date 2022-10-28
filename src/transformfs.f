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
      subroutine transformfs(inpc,textpart,trab,ntrans,ntrans_,
     &     set,istartset,iendset,ialset,nset,istep,istat,
     &     n,iline,ipol,inl,ipoinp,inp,ipoinpc,xload,sideload,
     &     nelemload,idefload,nload,nload_,ne,nam,iamload,ier)
!
!     reading the input deck: *TRANSFORMF
!
      implicit none
!
      real*8 trab(7,*)
!
      character*1 inpc(*)
      character*20 sideload(*),label
      character*81 set(*),surfaceset
      character*132 textpart(16)
!
      integer ntrans,ntrans_,istep,istat,n,key,i,j,id,
     &  istartset(*),iendset(*),ialset(*),nset,ipos,iline,ipol,
     &  inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),nelemload(2,*),ne,
     &  nam,nload,l,iamload(*),iamplitude,idefload(*),
     &  nload_,iset,nelem,ifacel,ier
!
      real*8 xload(2,*),xmagnitude
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *TRANSFORMF: *TRANSFORMF should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      ntrans=ntrans+1
      if(ntrans.gt.ntrans_) then
         write(*,*) '*ERROR reading *TRANSFORMF: increase ntrans_'
         ier=1
         return
      endif
!
      ipos=1
      xmagnitude=0.d0
      iamplitude=0
!
!     rectangular coordinate system: trab(7,norien)=1
!     cylindrical coordinate system: trab(7,norien)=-1
!     default is rectangular
!
      trab(7,ntrans)=1.d0
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:6).eq.'C') then
               trab(7,ntrans)=-1.d0
            endif
         elseif(textpart(i)(1:8).eq.'SURFACE=') then
            surfaceset(1:80)=textpart(i)(9:88)
            surfaceset(81:81)=' '
            ipos=index(surfaceset,' ')
            surfaceset(ipos:ipos)='T'
c            do iset=1,nset
c               if(set(iset).eq.surfaceset) exit
c            enddo
            call cident81(set,surfaceset,nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(surfaceset.eq.set(id)) then
                iset=id
              endif
            endif
            if(iset.gt.nset) then
               write(*,*) 
     &             '*WARNING reading *TRANSFORMF: element surface ',
     &              surfaceset(1:ipos-1),' does not exist'
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *TRANSFORMF: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*TRANSFORMF%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR reading *TRANSFORMF: definition of a'
         write(*,*) '  transformation is not complete'
         call inputerror(inpc,ipoinpc,iline,
     &        "*TRANSFORMF%",ier)
         return
      endif
!
      do i=1,6
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) trab(i,ntrans)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*TRANSFORMF%",ier)
            return
         endif
      enddo
!
      label(1:20)='T                   '
      do j=istartset(iset),iendset(iset)
         l=ialset(j)
         nelem=int(l/10.d0)
         ifacel=l-10*nelem
         write(label(2:2),'(i1)') ifacel
         call loadadd(nelem,label,xmagnitude,nelemload,
     &        sideload,xload,nload,nload_,iamload,
     &        iamplitude,nam,ntrans,idefload)
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

