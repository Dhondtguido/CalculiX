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
      subroutine refinemeshs(inpc,textpart,filab,istep,istat,n,iline,
     &     ipol,inl,ipoinp,inp,ipoinpc,ier,set,nset)
!     
!     reading *REFINE MESH in the input deck
!     
      implicit none
!     
      character*1 inpc(*)
      character*81 set(*),elset
      character*87 filab(*)
      character*132 textpart(16)
!     
      integer istep,istat,n,key,i,ii,ier,iline,ipol,inl,ipoinp(2,*),
     &     inp(3,*),ipoinpc(0:*),nset,ipos,id
!     
      if(istep.lt.1) then
        write(*,*) 
     &       '*ERROR reading *REFINE MESH: *REFINE MESH'
        write(*,*) '       should only be used within a *STEP' 
        write(*,*) '       definition'
        ier=1
        return
      endif
!     
!     reset the refinement controlling field
!     
      filab(48)(1:6)='RM    '
!
      ipos=0
      do ii=2,n
        if(textpart(ii)(1:6).eq.'LIMIT=') then
          filab(48)(7:26)=textpart(ii)(7:26)
c        elseif(textpart(ii)(1:4).eq.'USER') then
c          filab(48)(1:6)='RMUSER'
        elseif(textpart(ii)(1:6).eq.'ELSET=') THEN
          elset(1:60)=textpart(ii)(7:66)
          if(textpart(ii)(67:67).ne.' ') then
            write(*,*) 'set name too long (more than 60 characters)'
            call inputerror(inpc,ipoinpc,iline,
     &           "*REFINE MESH%",ier)
            return
          endif
          do i=61,81
            elset(i:i)=' '
          enddo
          ipos=index(elset,' ')
          elset(ipos:ipos)='E'
        else
          write(*,*) 
     &         '*WARNING reading *REFINE MESH:' 
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(ii)(1:index(textpart(ii),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*REFINE MESH %")
        endif
      enddo
!
!     check whether a set was defined, and if any, whether it exists
!      
      if(ipos.gt.0) then
c        do i=1,nset
c          if(set(i).eq.elset) exit
c        enddo
        call cident81(set,elset,nset,id)
        i=nset+1
        if(id.gt.0) then
          if(elset.eq.set(id)) then
            i=id
          endif
        endif
        if(i.gt.nset) then
          elset(ipos:ipos)=' '
          write(*,*) '*ERROR reading *REFINE MESH: element set',elset
          write(*,*) '       has not been defined'
          call inputerror(inpc,ipoinpc,iline,
     &         "*REFINE MESH%",ier)
          return
        endif
        write(filab(48)(27:87),'(a61)') elset(1:61)
      else
        do i=27,87
          filab(48)(i:i)=' '
        enddo
      endif
!     
      do
        call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &       ipoinp,inp,ipoinpc)
        if((key.eq.1).or.(istat.lt.0)) return
!     
        if(n.gt.1) then
          write(*,*) '*WARNING reading *REFINE MESH'
          write(*,*) '         only one refinement controlling'
          write(*,*) '         field allowed'
          call inputwarning(inpc,ipoinpc,iline,
     &         "*REFINE MESH %")
          return
        endif
!     
        do ii=1,n
          if((textpart(ii)(1:4).eq.'ERR ').or.
     &         (textpart(ii)(1:4).eq.'U   ').or.
     &         (textpart(ii)(1:4).eq.'NT  ').or.
     &         (textpart(ii)(1:4).eq.'S   ').or.
     &         (textpart(ii)(1:4).eq.'E   ').or.
     &         (textpart(ii)(1:4).eq.'ME  ').or.
     &         (textpart(ii)(1:4).eq.'PEEQ').or.
     &         (textpart(ii)(1:4).eq.'ENER').or.
     &         (textpart(ii)(1:4).eq.'HFL')) then
            filab(48)(3:6)=textpart(ii)(1:4)
          else
            write(*,*) 
     &           '*WARNING reading *REFINE MESH: label not applicable'
            write(*,*) '         or unknown; '
            call inputwarning(inpc,ipoinpc,iline,
     &           "*REFINE MESH %")
          endif
        enddo
      enddo
!     
      return
      end
