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
      subroutine geometrictolerances(inpc,textpart,set,istartset,
     &     iendset,ialset,nset,nk,istep,istat,n,iline,ipol,inl,ipoinp,
     &     inp,ipoinpc,ier,irobustdesign,irandomtype,randomval)
!     
!     reading the input deck: *GEOMETRIC TOLERANCE
!     
      implicit none
!     
      character*1 inpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!     
      integer istartset(*),iendset(*),ialset(*),nset,istep,istat,
     &     n,i,j,k,key,ipos,nk,iline,ipol,inl,ipoinp(2,*),id,
     &     inp(3,*),ipoinpc(0:*),ier,itype,inode,irandomtype(*),
     &     irobustdesign(*),ndesi  
!     
      real*8 meanval,stddev,randomval(2,*)
!     
      if(istep.lt.1) then
        write(*,*) '*ERROR reading *GEOMETRIC TOLERANCE:' 
        write(*,*) '       *GEOMETRIC TOLERANCE should only be used'
        write(*,*) '       within a *ROBUST DESIGN STEP'
        ier=1
        return
      endif
!     
!     irobustdesign(2)=1 --> homogeneous randomfield
!         --> randomfield can be characterized by
!             a single mean value and a standard deviation
!     irobustdesign(2)=0 --> different mean and standard deviation at
!             every node
!      
!     irobustdesign(3)=1 --> constrained random field: boundaries of the 
!             design space are set to zero (default: irobustdesign(3)=0)
!     
      do i=1,n
        if(textpart(i)(6:11).eq.'NORMAL') then
          itype=1
        elseif(textpart(i)(1:11).eq.'CONSTRAINED') then
          irobustdesign(3)=1
        endif
      enddo
!
!     it is always assumed that the randomfield is nonhomogeneous
!     due to the fact that a homogeneous is just a special case
!     of a nonhomogeneous randomfield
!
      irobustdesign(2)=0
!     
!     normal distribution is defined as default
!      
      if(itype.eq.0) then       
        itype=1 
      endif
!     
      ndesi=0
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((istat.lt.0).or.(key.eq.1)) return
!     
        read(textpart(2)(1:20),'(f20.0)',iostat=istat) meanval
        read(textpart(3)(1:20),'(f20.0)',iostat=istat) stddev
!     
        read(textpart(1)(1:10),'(i10)',iostat=istat) inode
        if(istat.eq.0) then
          if(inode.gt.nk) then
            write(*,*) '*ERROR reading *GEOMETRIC TOLERANCE:'
            write(*,*) '       node ',inode,' is not defined'
            ier=1
            call inputerror(inpc,ipoinpc,iline,
     &           "*GEOMETRIC TOLERANCE%",ier)
            return
          endif
          ndesi=ndesi+1
          irandomtype(inode)=itype
          randomval(1,inode)=meanval
          randomval(2,inode)=stddev
        else
          read(textpart(1)(1:80),'(a80)',iostat=istat) noset
          noset(81:81)=' '
          ipos=index(noset,' ')
          noset(ipos:ipos)='N'
c     do i=1,nset
c     if(set(i).eq.noset) exit
c     enddo
          call cident81(set,noset,nset,id)
          i=nset+1
          if(id.gt.0) then
            if(noset.eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
            noset(ipos:ipos)=' '
            write(*,*) '*ERROR reading *GEOMETRIC TOLERANCE:'
            write(*,*) '       node set ',noset,' has not'
            write(*,*) '       yet been defined.'
            call inputerror(inpc,ipoinpc,iline,
     &           "*GEOMETRIC TOLERANCE%",ier)
            return
          endif
          do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
              k=ialset(j)
              ndesi=ndesi+1
              irandomtype(k)=itype
              randomval(1,k)=meanval
              randomval(2,k)=stddev
            else
              k=ialset(j-2)
              do
                k=k-ialset(j)
                if(k.ge.ialset(j-1)) exit
                ndesi=ndesi+1
                irandomtype(k)=itype
                randomval(1,k)=meanval
                randomval(2,k)=stddev
              enddo
            endif
          enddo
        endif
      enddo
!
!     check if some variables are defined
!     
      if(ndesi.eq.0) then
        write(*,*) '*ERROR reading *GEOMETRIC TOLERANCE:'
        write(*,*) '     no nodes for the compuation of '
        write(*,*) '     the random field have been defined.'
        call inputerror(inpc,ipoinpc,iline,
     &       "*GEOMETRIC TOLERANCE%",ier)
        return
      endif
!     
      return
      end
