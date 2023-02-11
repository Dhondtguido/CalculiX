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
      subroutine geometricconstraints(inpc,textpart,istat,n,
     &     iline,ipol,inl,ipoinp,inp,ipoinpc,nobject,objectset,
     &     set,nset,ier,nmethod)     
!     
!     reading the input deck: *GEOMETRIC CONSTRAINT
!     
      implicit none
!     
      character*1 inpc(*), settype
      character*132 textpart(16)
      character*81 objectset(5,*),set(*),drname
!     
      integer istat,n,key,i,iline,ipol,inl,ipoinp(2,*),nset,id,m,
     &     inp(3,*),ipoinpc(0:*),nobject,k,ipos,ier,nmethod,nsets
!     
      real*8 absval
!     
!     geometric constraints can only be defined 
!     within feasible direction steps.
!     
      if(nmethod.ne.16) then
        write(*,*) '*ERROR reading *GEOMETRIC CONSTRAINT'
        write(*,*) '       *GEOMETRIC CONSTRAINT can only be specified'
        write(*,*) '       within a *FEASIBILE DIRECTION step.'   
        call inputerror(inpc,ipoinpc,iline,
     &       "*GEOMETRIC CONSTRAINT%",ier)
        return
      endif
!     
      do
!     
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
!        
        if((istat.lt.0).or.(key.eq.1)) exit
        drname=textpart(1)(1:80)
!     
!       store the constraint (a geometric constraint is not 
!       connected to an existing design response
!     
        nobject=nobject+1
        do i=1,5
          objectset(i,nobject)(1:81)=' '
        enddo
!        
        objectset(5,nobject)(81:81)='G'
        if(textpart(1)(1:12).eq.'MAXSHRINKAGE') then
          objectset(1,nobject)(1:12)='MAXSHRINKAGE'
c          do k=13,20
c            objectset(1,nobject)(k:k)=' '
c          enddo
          settype='N'
          objectset(1,nobject)(19:20)='LE'
          nsets=1
        elseif(textpart(1)(1:9).eq.'MAXGROWTH') then
          objectset(1,nobject)(1:9)='MAXGROWTH'
c          do k=10,20
c            objectset(1,nobject)(k:k)=' '
c          enddo
          settype='N'
          objectset(1,nobject)(19:20)='LE'
          nsets=1
        elseif (textpart(1)(1:13).eq.'MAXMEMBERSIZE') then
          objectset(1,nobject)(1:13)='MAXMEMBERSIZE'
c          do k=14,20
c            objectset(1,nobject)(k:k)=' '
c          enddo
          settype='N'
          objectset(1,nobject)(19:20)='LE'
          nsets=2
        elseif (textpart(1)(1:13).eq.'MINMEMBERSIZE') then
          objectset(1,nobject)(1:13)='MINMEMBERSIZE'
c          do k=14,20
c            objectset(1,nobject)(k:k)=' '
c          enddo
          settype='N'
          objectset(1,nobject)(19:20)='GE'
          nsets=2
        elseif (textpart(1)(1:9).eq.'PACKAGING') then
          objectset(1,nobject)(1:9)='PACKAGING'
c          do k=9,20
c            objectset(1,nobject)(k:k)=' '
c          enddo
          settype='N'
          objectset(1,nobject)(19:20)='GE'
          nsets=2
        else
          write(*,*) '*ERROR reading *GEOMETRIC CONSTRAINT'
          write(*,*) '       given constraint type is not a'
          write(*,*) '       valid option.'
          call inputerror(inpc,ipoinpc,iline,
     &         "*GEOMETRIC CONSTRAINT%",ier)
          return 
        endif
!     
!       reading the sets needed for the geometric constraint
!     
        do m=1,nsets
          objectset(2+m,nobject)(1:80)=textpart(1+m)(1:80)
          ipos=index(objectset(2+m,nobject),' ')
          if(n.lt.m+1) then
            write(*,*)'*ERROR reading *GEOMETRIC CONSTRAINT'
            write(*,*)'       set ',m,' is lacking'
            call inputerror(inpc,ipoinpc,iline,
     &           "*GEOMETRIC CONSTRAINT%",ier)
            return
          endif
          objectset(2+m,nobject)(ipos:ipos)=settype
!          
c          do i=1,nset
c            if(set(i).eq.objectset(2+l,nobject)) exit
c          enddo
          call cident81(set,objectset(2+m,nobject),nset,id)
          i=nset+1
          if(id.gt.0) then
            if(objectset(2+m,nobject).eq.set(id)) then
              i=id
            endif
          endif
          if(i.gt.nset) then
            write(*,*) '*ERROR reading *GEOMETRIC CONSTRAINT'
            write(*,*) '       unknown set name: '
            write(*,*) objectset(2+m,nobject)
            call inputerror(inpc,ipoinpc,iline,
     &           "*GEOMETRIC CONSTRAINT%",ier)
            return
          endif
        enddo
!     
!     except of the PACKAGING constraint (does not need any value) 
!     all the other geometric constraints always take ONLY a single 
!     absolute value and no relative values! 
!     
        if(objectset(1,nobject)(1:9).ne.'PACKAGING') then
          if(n.ge.(2+nsets)) then
            read(textpart(2+nsets)(1:20),'(f20.0)',iostat=istat) absval
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*GEOMETRIC CONSTRAINT%",ier)
              return
            endif
            if(istat.le.0) then
              objectset(1,nobject)(61:80)=textpart(2+nsets)(1:20)
            endif
          else
            write(*,*) '*ERROR reading *GEOMETRIC CONSTRAINT'
            write(*,*) '       no absolute value for geometric'
            write(*,*) '       constraint specified.'
            call inputerror(inpc,ipoinpc,iline,
     &             "*GEOMETRIC CONSTRAINT%",ier)
            return
          endif
        endif
      enddo
!      
      return
      end
