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
     &     iquad,iedtet,iedgmid,number,jfix,iparentel,nk,iwrite)
!
      implicit none
!
      character*10 elestr
      character*132 fnrfn,jobnamec(*),el_header
      character*256 fn
!
      integer kontet(4,*),netet_,i,j,k,nktet,node,iquad,iedtet(6,*),
     &     iedgmid(*),number(*),nk,jfix(*),iparentel(*),iwrite,ilen
!
      real*8 cotet(3,*)
!
!     give nodes of the unrefined mesh which were not fixed
!     a new node number in order to avoid collisions with the
!     refined mesh
!      
      do i=1,netet_
        if(kontet(1,i).ne.0) then
          do j=1,4
            node=kontet(j,i)
            if((jfix(node).ne.1).and.(node.le.nk)) then
              if(number(node).ne.0) then
                kontet(j,i)=number(node)
              else
                nktet=nktet+1
                number(node)=nktet
                kontet(j,i)=nktet
                do k=1,3
                  cotet(k,nktet)=cotet(k,node)
                enddo
              endif
            endif
          enddo
        endif
      enddo
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
      open(2,file=fnrfn(1:i+7),status='unknown',position='append')
!
!     storing the nodes
!
      write(2,102)
 102  format('*NODE')
      do i=1,nktet
!
!        setting too small numbers to zero (else the exponent in the
!        output contains 3 digits and the letter "D" is omitted)
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
     &           ',TYPE=C3D4'               
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
     &           ',TYPE=C3D10'               
            write(2,*) el_header(1:37-k+1)
!     
!     topology
!     
            write(2,101) i,(kontet(j,i),j=1,4),
     &           (iedgmid(iedtet(j,i)),j=1,6)
          endif
        enddo
      endif
!     
      close(2)
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
      end
