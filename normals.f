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
      subroutine normals(inpc,textpart,iponor,xnor,ixfree,
     &  ipkon,kon,nk,nk_,ne,lakon,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *NORMAL
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*132 textpart(16)
!
      integer iponor(2,*),ixfree,ipkon(*),kon(*),nk,ipoinpc(0:*),
     &  nk_,ne,istep,istat,n,ielement,node,j,indexe,i,ier,
     &  key,iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 xnor(*),x,y,z,dd
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *NORMAL: *NORMAL should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *NORMAL: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*NORMAL%")
      enddo
!
      loop:do
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) ielement
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NORMAL%",ier)
            return
         endif
         read(textpart(2)(1:10),'(i10)',iostat=istat) node
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NORMAL%",ier)
            return
         endif
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) x
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NORMAL%",ier)
            return
         endif
         if(n.le.3) then
            y=0.d0
         else
            read(textpart(4)(1:20),'(f20.0)',iostat=istat) y
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*NORMAL%",ier)
               return
            endif
         endif
         if(n.le.4) then
            z=0.d0
         else
            read(textpart(5)(1:20),'(f20.0)',iostat=istat) z
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*NORMAL%",ier)
               return
            endif
         endif
!
!        normalizing the normal
!
         dd=dsqrt(x*x+y*y+z*z)
         x=x/dd
         y=y/dd
         z=z/dd
!
         if(ielement.gt.ne) then
            write(*,*) '*ERROR reading *NORMAL: element number',ielement
            write(*,*) '       exceeds ne'
            ier=1
            return
         endif
!
         indexe=ipkon(ielement)
         do j=1,8
            if(kon(indexe+j).eq.node) then
               iponor(1,indexe+j)=ixfree
               if(lakon(ielement)(1:1).eq.'B') then
                  xnor(ixfree+4)=x
                  xnor(ixfree+5)=y
                  xnor(ixfree+6)=z
                  ixfree=ixfree+6
               elseif(lakon(ielement)(1:2).ne.'C3') then
                  xnor(ixfree+1)=x
                  xnor(ixfree+2)=y
                  xnor(ixfree+3)=z
                  ixfree=ixfree+3
               else
                  write(*,*) 
     &               '*WARNING reading *NORMAL: specifying a normal'
                  write(*,*) '         3-D element does not make sense'
               endif
               cycle loop
            endif
         enddo
         write(*,*) '*WARNING: node ',node,' does not belong to'
         write(*,*) '          element ',ielement
         write(*,*) '          normal definition discarded'
!
      enddo loop
!
      return
      end










