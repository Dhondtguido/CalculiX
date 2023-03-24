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
      subroutine massflows(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nelemload,sideload,xload,nload,nload_,iamload,
     &  lakon,ne,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ipoinpc,idefload,nam,ier)
!
!     reading the input deck: *MASS FLOW
!
      implicit none
!
      logical surface
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 sideload(*),label
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nelemload(2,*),id,
     &  nset,nload,nload_,istep,istat,n,i,j,l,key,idefload(*),
     &  iamload(2,*),nam,iamplitude,ipos,ne,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),idelay,isector,ipoinpc(0:*),ier
!
      real*8 xload(2,*),xmagnitude
!
      iamplitude=0
      idelay=0
      isector=0
      surface=.false.
!
      if(istep.ne.1) then
         write(*,*) 
     &     '*ERROR reading *MASS FLOW: *MASS FLOW should only be used'
         write(*,*) '  in the first STEP'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *MASS FLOW: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &        "*MASS FLOW%")
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:20),'(a20)',iostat=istat) label
!
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmagnitude
         if(xmagnitude.ne.0.d0) then
            write(*,*) '*WARNING reading *MASS FLOW:'
            write(*,*) '         magnitude for label: ',label
            write(*,*) '         is not zero but'
            write(*,*) '         takes the value: ',xmagnitude
            write(*,*) '         it is set to zero'
            xmagnitude=0.d0
         endif
!
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*MASS FLOW%",ier)
            return
         endif
         if((label(1:2).ne.'M1').and.(label(1:2).ne.'M2').and.
     &           (label(1:2).ne.'M ').and.
     &           (label(1:2).ne.'M3').and.(label(1:2).ne.'M4').and.
     &           (label(1:2).ne.'M5').and.(label(1:2).ne.'M6')) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*MASS FLOW%",ier)
            return
         endif
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if(l.gt.ne) then
               write(*,*) '*ERROR reading *MASS FLOW: element ',l
               write(*,*) '       is not defined'
               ier=1
               return
            endif
!
            if(lakon(l)(1:1).ne.'F') then
               write(*,*) '*ERROR reading *MASS FLOW: element ',l
               write(*,*) '       is not a fluid element*'
               ier=1
               return
            endif
            call loadadd(l,label,xmagnitude,nelemload,sideload,
     &           xload,nload,nload_,iamload,iamplitude,
     &           nam,isector,idefload)
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) elset
            elset(81:81)=' '
            ipos=index(elset,' ')
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
!
!              check for facial surface
!
               surface=.true.
               elset(ipos:ipos)='T'
c               do i=1,nset
c                  if(set(i).eq.elset) exit
c               enddo
               call cident81(set,elset,nset,id)
               i=nset+1
               if(id.gt.0) then
                 if(elset.eq.set(id)) then
                   i=id
                 endif
               endif
               if(i.gt.nset) then
                  elset(ipos:ipos)=' '
                  write(*,*) '*ERROR reading *MASS FLOW: element set '
                  write(*,*) '       or facial surface ',elset
                  write(*,*) '       has not yet been defined. '
                  call inputerror(inpc,ipoinpc,iline,
     &                                  "*MASS FLOW%",ier)
                  return
               endif
            endif
!
            l=ialset(istartset(i))
            if(.not.surface) then
               if(lakon(l)(1:1).ne.'F') then
                  write(*,*) '*ERROR reading *MASS FLOW: element ',l
                  write(*,*) '       is not a fluid element*'
                  ier=1
                  return
               endif
            else
               if(lakon(l/10)(1:1).ne.'F') then
                  write(*,*) '*ERROR reading *MASS FLOW: element ',l/10
                  write(*,*) '       is not a fluid element*'
                  ier=1
                  return
               endif
            endif
!
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  l=ialset(j)
                  if(surface) then
                     write(label(2:2),'(i1)') l-10*(l/10)
                     l=l/10
                  endif
                  call loadadd(l,label,xmagnitude,nelemload,sideload,
     &                 xload,nload,nload_,iamload,iamplitude,
     &                 nam,isector,idefload)
               else
                  l=ialset(j-2)
                  do
                     l=l-ialset(j)
                     if(l.ge.ialset(j-1)) exit
                     call loadadd(l,label,xmagnitude,nelemload,
     &                    sideload,xload,nload,nload_,
     &                    iamload,iamplitude,nam,isector,idefload)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

