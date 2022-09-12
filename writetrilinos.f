!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine writetrilinos(jq,irow,ad,au,b,neq,nzs,symmetryflag,
     &     inputformat,co,nk,nactdof,jobnamef,mi)
!     
      implicit none
!
      character*132 jobnamef,fnlhs,fnrhs,fnrig
!     
      integer jq(*),irow(*),neq,nzs,symmetryflag,inputformat,nk,
     &     mi(*),nactdof(0:mi(2),*),iperm(3*nk),i,j,k
!     
      real*8 ad(*),au(*),b(*),co(3,*),null,one
!
      fnlhs=jobnamef(1:index(jobnamef,' ')-1)//'.lhs'
      fnrhs=jobnamef(1:index(jobnamef,' ')-1)//'.rhs'
      fnrig=jobnamef(1:index(jobnamef,' ')-1)//'.rig'
!
      null=0.d0
      one=1.d0
!
      open(50,file=fnlhs,status='unknown')
!
!     left hand side: fill-in equations
!
      k=0
      do i=1,nk
        if((nactdof(1,i).le.0).and.(nactdof(2,i).le.0).and.
     &       (nactdof(3,i).le.0)) cycle
        do j=1,3
          k=k+1
          if(nactdof(j,i).gt.0) then
            iperm(nactdof(j,i))=k
          else
            write(50,*) k,k,one
          endif
        enddo
      enddo
!
!     left hand side: significant contributions
!
      if(symmetryflag.eq.0) then
!
!       symmetric
!
        do i=1,neq
          write(50,*) iperm(i),iperm(i),ad(i)
        enddo
        do i=1,neq
          do j=jq(i),jq(i+1)-1
            write(50,*) iperm(irow(j)),iperm(i),au(j)
            write(50,*) iperm(i),iperm(irow(j)),au(j)
          enddo
        enddo
      elseif(inputformat.eq.1) then
!
!        not symmetric
!
        do i=1,neq
          write(50,*) iperm(i),iperm(i),ad(i)
        enddo
        do i=1,neq
          do j=jq(i),jq(i+1)-1
            write(50,*) iperm(irow(j)),iperm(i),au(j)
            write(50,*) iperm(i),iperm(irow(j)),au(j+nzs)
          enddo
        enddo
      else
        write(*,*) '*ERROR in writetrilinos: input format'
        write(*,*) '       not known'
        call exit(201)
      endif
      close(50)
!
!     right hand side
!      
      open(50,file=fnrhs,status='unknown')
      do i=1,neq
        write(50,*) iperm(i),b(i)
      enddo
      close(50)
!
!     rigid body modes
!      
      open(50,file=fnrig,status='unknown')
      do i=1,nk
        if((nactdof(1,i).le.0).and.(nactdof(2,i).le.0).and.
     &       (nactdof(3,i).le.0)) cycle
        write(50,*) one,null,null,null,co(3,i),-co(2,i)
        write(50,*) null,one,null,-co(3,i),null,co(1,i)
        write(50,*) null,null,one,co(2,i),-co(1,i),null
      enddo
      close(50)
!     
      stop
      end

