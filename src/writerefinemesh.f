!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine writerefinemesh(kontet,netet_,cotet,nktet,jobnamec,
     &     iquad,iedtet,iedgmid,number,jfix,iparentel,nk,iwrite,
     &     maxnnewnodes,kontetor)
!
      implicit none
!
      character*10 elestr
      character*132 fnrfn,jobnamec(*),el_header
      character*256 fn
!
      integer kontet(4,*),netet_,i,j,k,nktet,node,iquad,iedtet(6,*),
     &     iedgmid(*),number(*),nk,jfix(*),iparentel(*),iwrite,ilen,
     &     nodeold,nodenew,kontetor(6,*),maxnnewnodes
!
      real*8 cotet(3,*)
!
!     give nodes of the unrefined mesh which were not fixed
!     a new node number in order to avoid collisions with the
!     refined mesh
!      
c      do i=1,netet_
c        if(kontet(1,i).ne.0) then
c          do j=1,4
c            node=kontet(j,i)
c            if((jfix(node).ne.1).and.(node.le.nk)) then
c              if(number(node).ne.0) then
c                kontet(j,i)=number(node)
c              else
c                nktet=nktet+1
c                number(node)=nktet
c                kontet(j,i)=nktet
c                do k=1,3
c                  cotet(k,nktet)=cotet(k,node)
c                enddo
c              endif
c            endif
c          enddo
c        endif
c      enddo
!
!     stores the refined mesh in input format
!
      do i=1,132
         if(ichar(jobnamec(1)(i:i)).eq.0) exit
      enddo
      if(i.gt.125) then
         write(*,*) '*ERROR in writerefinemesh'
         write(*,*) '       jobname has more than 124 characters'
         call exit(201)
      endif
      fnrfn(1:i+7)=jobnamec(1)(1:i-1)//'.rfn.inp'
!
!     storing the mesh in input format
!
      if(maxnnewnodes.gt.0) then
!
!     new nodes were created: do not keep the original numbering
!     of the midnodes
!
        open(2,file=fnrfn(1:i+7),status='unknown',position='append')
!     
!     storing the nodes
!     
        write(2,102)
 102    format('*NODE')
        do i=1,nktet
!     
!     setting too small numbers to zero (else the exponent in the
!     output contains 3 digits and the letter "D" is omitted)
!     
          do j=1,3
            if(dabs(cotet(j,i)).lt.1.d-99) cotet(j,i)=0.d0
          enddo
          write(2,100) i,(cotet(j,i),j=1,3)
        enddo
!     
!     storing the tetrahedral elements
!     
        if(iquad.eq.0) then
          do i=1,netet_
            if(kontet(1,i).ne.0) then
!     
!     keyword card
!     
              write(elestr,'(i10)') iparentel(i)
              do k=1,10
                if(elestr(k:k).ne.' ') exit
              enddo
              el_header='*ELEMENT,PARENT='//elestr(k:10)//
     &             ',TYPE=C3D4'               
              write(2,*) el_header(1:36-k+1)
!     
!     topology
!     
              write(2,101) i,(kontet(j,i),j=1,4)
            endif
          enddo
        else
          do i=1,netet_
            if(kontet(1,i).ne.0) then
!     
!     keyword card
!     
              write(elestr,'(i10)') iparentel(i)
              do k=1,11
                if(elestr(k:k).ne.' ') exit
              enddo
              el_header='*ELEMENT,PARENT='//elestr(k:10)//
     &             ',TYPE=C3D10'               
              write(2,*) el_header(1:37-k+1)
!     
!     topology
!     
              write(2,101) i,(kontet(j,i),j=1,4),
     &             (iedgmid(iedtet(j,i)),j=1,6)
            endif
          enddo
        endif
!     
        close(2)
!
      else
!
!       no new nodes were created: keep the original numbering
!       of the midnodes (only smoothing)
!
        open(2,file=fnrfn(1:i+7),status='unknown',err=51)
        close(2,status='delete',err=52)
        open(2,file=fnrfn(1:i+7),status='unknown',err=51)
!
        write(2,102)
        if(iquad.eq.0) then
          do i=1,nktet
!     
!           setting too small numbers to zero (else the exponent in the
!           output contains 3 digits and the letter "D" is omitted)
!     
            do j=1,3
              if(dabs(cotet(j,i)).lt.1.d-99) cotet(j,i)=0.d0
            enddo
            write(2,100) i,(cotet(j,i),j=1,3)
          enddo
        else
!     
!         the modified midnodes (modified w.r.t. their location) got
!         new node numbers. If only smoothing was performed, these
!         modified coordinates should be attached to the orginal midnode
!         numbers
!     
          do i=1,netet_
            if(kontet(1,i).ne.0) then
c              write(*,*) i,(kontet(j,i),j=1,10)
              do j=1,6
                nodeold=kontetor(j,i)
                nodenew=iedgmid(iedtet(j,i))
c                write(*,*) 'writerefinemesh',nodeold,nodenew
                do k=1,3
                  cotet(k,nodeold)=cotet(k,nodenew)
                enddo
              enddo
            endif
          enddo
!     
!         writing the coordinates to file     
!     
          do i=1,nk
            write(2,100) i,(cotet(j,i),j=1,3)
          enddo
        endif
!     
        close(2)
!     
      endif
!     
 100  format(i10,',',e20.13,',',e20.13,',',e20.13)
 101  format(11(i10,','))
!
      if(iwrite.eq.1) then
        ilen=index(jobnamec(1),char(0))-1
        fn=jobnamec(1)(1:ilen)//'_WarnNodeNotProjected.nam'
        write(*,*) '*INFO in writerefinemesh:'
        write(*,*) '      not (completely) projected nodes'
        write(*,*) '      are stored in file'
        write(*,*) '      ',fn(1:ilen+25)
        write(*,*) '      This file can be loaded into'
        write(*,*) '      an active cgx-session by typing'
        write(*,*) 
     &       '      read ',fn(1:ilen+25),' inp'
        write(*,*)
        close(40)
      else
        close(40,status='delete')
      endif
!     
      return
!
 51   write(*,*) '*ERROR in openfile: could not open file ',fnrfn(1:i+7)
      call exit(201)
 52   write(*,*) '*ERROR in openfile: could not delete file ',
     &  fnrfn(1:i+7)
      call exit(201)
!
      end
