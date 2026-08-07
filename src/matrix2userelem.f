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
      subroutine matrix2userelem(textpart,n,iuel,nuel,inpc,ipoinpc,
     &     iline,ier,ipoinp,inp,inl,ipol,nk_,mi,icfd,ne_,nkon_,nmat_)
!     
!     reading the input deck: *MATRIX ASSEMBLE: 
!     creating a user element
!     1) determining the name
!     2) opening the stiffness file and determining the number of nodes
!     
      implicit none
!
      logical input,name
!     
      character*1 inpc(*)
      character*80 filename
      character*132 textpart(16)
!     
      integer n,iuel(4,*),nuel,i,j,k,l,istat,number,ipoinpc(0:*),iline,
     &     four,nope,intpoints,maxdof,id,ier,key,ipoinp(2,*),inp(3,*),
     &     inl,ipol,node,nk_,mi(*),icfd,ne_,nkon_,nmat_,idof,ndof
!     
      integer,dimension(:),allocatable::inode
!     
      four=4
      name=.false.
      input=.false.
      ndof=0
!     
      do i=2,n
        if(textpart(i)(1:5).eq.'NAME=') then
          name=.true.
          number=ichar(textpart(i)(6:6))*256**3+
     &         ichar(textpart(i)(7:7))*256**2+
     &         ichar(textpart(i)(8:8))*256+
     &         ichar(textpart(i)(9:9))
        elseif(textpart(i)(1:14).eq.'STIFFNESSFILE=') then
          input=.true.
          filename(1:80)=textpart(i)(15:94)
          loop1: do j=1,80
            if(filename(j:j).eq.'"') then
              do k=j+1,80
                if(filename(k:k).eq.'"') then
                  do l=k-1,80
                    filename(l:l)=' '
                    exit loop1
                  enddo
                endif
                filename(k-1:k-1)=filename(k:k)
              enddo
              filename(80:80)=' '
            endif
          enddo loop1
        endif
      enddo
!     
!     check for the NAME=U parameter
!
      if(.not.name) then
        write(*,*) 
     &   '*ERROR reading *MATRIX ASSEMBLE: no name specified:'
        call inputerror(inpc,ipoinpc,iline,
     &       "*MATRIX ASSEMBLE%",ier)
      endif
!     
!     check for the INPUT parameter
!
      if(.not.input) then
        write(*,*) 
     &   '*ERROR reading *MATRIX ASSEMBLE: no stiffness file specified:'
        call inputerror(inpc,ipoinpc,iline,
     &       "*MATRIX ASSEMBLE%",ier)
      endif
!     
!     determine the number of nodes in the stiffness matrix
!     
      allocate(inode(nk_))
      open(20,file=filename,status='old',err=3)
      nope=0
      do
        read(20,*,end=1,err=2) node,idof
        ndof=max(ndof,idof)
        call nident(inode,node,nope,id)
        if(id.gt.0) then
          if(inode(id).eq.node) cycle
        endif
        nope=nope+1
        do j=nope,id+2,-1
          inode(j)=inode(j-1)
        enddo
        inode(id+1)=node
      enddo
 1    close(20)
      deallocate(inode)
!     
!     check range
!     
      if(nope.gt.255) then
        write(*,*) '*ERROR reading *MATRIX ASSEMBLE'
        write(*,*) '       number of nodes ',nope,' exceeds 255'
        ier=1
        return
      endif
!     
!     storing the element information in iuel
!     
      call nidentk(iuel,number,nuel,id,four)
!     
      if(id.gt.0) then
        if(iuel(1,id).eq.number) then
          write(*,*) '*ERROR reading *MATRIX ASSEMBLE'
          write(*,*) '       name was already used for'
          write(*,*) '       another *MATRIX ASSEMBLE or a'
          write(*,*) '       user element'
          ier=1
          return
        endif
      endif
!
      if((ndof.ne.3).and.(ndof.ne.6)) then
        write(*,*) '*ERROR in matrix2userelem'
        write(*,*) '       a substructure (=user element) should'
        write(*,*) '       either have 3 or 6 dofs;'
        write(*,*) '       actual number of dofs: ',ndof
        ier=1
        return
      endif
!     
      nuel=nuel+1
      do i=nuel,id+2,-1
        do j=1,4
          iuel(j,i)=iuel(j,i-1)
        enddo
      enddo
      iuel(1,id+1)=number
      iuel(2,id+1)=0
c      iuel(3,id+1)=3
      iuel(3,id+1)=ndof
      iuel(4,id+1)=nope
!
!     3 dofs per node assumed
!
c      mi(2)=max(mi(2),3)
      mi(2)=max(mi(2),ndof)
!
!     2 materials (matname is used to store the stiffness file and
!     the mass file)
!
      mi(3)=max(mi(3),2)
      nmat_=nmat_+2
!
      if(icfd.eq.-1) then
        icfd=0
      elseif(icfd.eq.1) then
        icfd=2
      endif
!
      nkon_=nkon_+nope
      ne_=ne_+1
!     
      return
 2    write(*,*)
      write(*,*) '*ERROR in stiffness matrix file ',filename
      write(*,*) '       incorrect format'
      call exit(201)
 3    write(*,*)
      write(*,*) '*ERROR reading stiffness matrix file: ',filename
      write(*,*) '       does not exist'
      call exit(201)
      end







