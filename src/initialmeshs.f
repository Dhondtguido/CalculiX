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
      subroutine initialmeshs(inpc,textpart,coini,nk,nk_,set,istat,n,
     &     iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!     
!     reading the input deck: *INITIAL MESH
!     
      implicit none
!     
      character*1 inpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!     
      integer nk,nk_,nkini,istat,n,key,id,i,j,js,k,nn,inoset,ipos,
     &     iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!     
      real*8 coini(3,*)
!     
!     read in nodal coordinates
!
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
        read(textpart(1)(1:10),'(i10)',iostat=istat) i
        if(istat.gt.0) then
          call inputerror(inpc,ipoinpc,iline,
     &         "*NODE%",ier)
          return
        endif
        if(n.eq.1) then
          coini(1,i)=0.d0
        else
          read(textpart(2)(1:20),'(f20.0)',iostat=istat) coini(1,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        if(n.le.2) then
          coini(2,i)=0.d0
        else
          read(textpart(3)(1:20),'(f20.0)',iostat=istat) coini(2,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        if(n.le.3) then
          coini(3,i)=0.d0
        else
          read(textpart(4)(1:20),'(f20.0)',iostat=istat) coini(3,i)
          if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODE%",ier)
            return
          endif
        endif
        nkini=max(nk,i)
        if(nkini.gt.nk_) then
          write(*,*) '*ERROR reading *INITIAL MESH:  Current '
          write(*,*) '       mesh and initial mesh do not match'
          write(*,*) '       w.r.t. highest node number'
          ier=1
          return
        endif   
      enddo
!     
      return
      end
