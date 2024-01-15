!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine matrixassembles(textpart,n,iuel,nuel_,inpc,ipoinpc,
     &     iline,ier,ipoinp,inp,inl,ipol,lakon,ipkon,kon,nkon,ne,ne_,
     &     ielmat,mi,matname,nmat,nmat_,irstrt,istep)
!     
!     reading the input deck: *MATRIX ASSEMBLE: 
!     creating a user element
!     1) creating the name
!     2) storing the topology in kon(*)
!     
      implicit none
!
      logical stiffness,mass
!     
      character*1 inpc(*)
      character*8 lakon(*),label
      character*80 matname(*),filestiff,filemass
      character*132 textpart(16)
!     
      integer n,iuel(4,*),nuel_,i,j,k,l,istat,number,ipoinpc(0:*),iline,
     &     four,nodes,intpoints,id,ier,key,ipoinp(2,*),inp(3,*),
     &     inl,ipol,node,ipkon(*),kon(*),ne,ne_,nkon,indexe,mi(*),
     &     ielmat(mi(3),*),ndof,nmat,nmat_,nope,irstrt(*),istep
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &       '*ERROR reading *MATRIX ASSEMBLE: *MATRIX ASSEMBLE should'
         write(*,*) '       be placed before all step definitions'
         ier=1
         return
      endif
!     
      four=4
      stiffness=.false.
      mass=.false.
!     
      do i=2,n
        if(textpart(i)(1:5).eq.'NAME=') then
          number=ichar(textpart(i)(6:6))*256**3+
     &         ichar(textpart(i)(7:7))*256**2+
     &         ichar(textpart(i)(8:8))*256+
     &         ichar(textpart(i)(9:9))
          label(1:1)='U'
          label(2:5)=textpart(i)(6:9)
        elseif(textpart(i)(1:14).eq.'STIFFNESSFILE=') then
          filestiff(1:80)=textpart(i)(15:94)
          loop1: do j=1,80
            if(filestiff(j:j).eq.'"') then
              do k=j+1,80
                if(filestiff(k:k).eq.'"') then
                  do l=k-1,80
                    filestiff(l:l)=' '
                    exit loop1
                  enddo
                endif
                filestiff(k-1:k-1)=filestiff(k:k)
              enddo
              filestiff(80:80)=' '
            endif
          enddo loop1
          stiffness=.true.
        elseif(textpart(i)(1:9).eq.'MASSFILE=') then
          filemass(1:80)=textpart(i)(10:89)
          loop2: do j=1,80
            if(filemass(j:j).eq.'"') then
              do k=j+1,80
                if(filemass(k:k).eq.'"') then
                  do l=k-1,80
                    filemass(l:l)=' '
                    exit loop2
                  enddo
                endif
                filemass(k-1:k-1)=filemass(k:k)
              enddo
              filemass(80:80)=' '
            endif
          enddo loop2
          mass=.true.
        endif
      enddo
!
      if(.not.stiffness) then
        write(*,*) '*ERROR reading *MATRIX ASSEMBLE: no stiffness'
        write(*,*) '       matrix given'
        call inputerror(inpc,ipoinpc,iline,
     &       "*MATRIX ASSEMBLE%",ier)
        return
      endif
!     
!     determine the number of nodes in the stiffness matrix
!     
      call nidentk(iuel,number,nuel_,id,four)
      intpoints=iuel(2,id)
      ndof=iuel(3,id)
      nope=iuel(4,id)
      write(label(6:6),'(a1)') char(intpoints)
      write(label(7:7),'(a1)') char(ndof)
      write(label(8:8),'(a1)') char(nope)
!
!     new element
!
      ne=ne+1
      if(ne.gt.ne_) then
        write(*,*) '*ERROR reading *ELEMENT: increase ne_'
        call exit(201)
      endif
      ipkon(ne)=nkon
      lakon(ne)=label
      indexe=nkon
!
      nkon=nkon+nope
!
      open(20,file=filestiff,status='old')
      nodes=0
      do
        read(20,*,end=1) node
        call nident(kon(indexe+1),node,nodes,id)
        if(id.gt.0) then
          if(kon(indexe+id).eq.node) cycle
        endif
        nodes=nodes+1
        do j=nodes,id+2,-1
          kon(indexe+j)=kon(indexe+j-1)
        enddo
        kon(indexe+id+1)=node
      enddo
 1    close(20)
!
!     storing the stiffness file name and mass file name as
!     material names
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
        write(*,*) '*ERROR reading *MATRIX ASSEMBLE: increase nmat_'
        ier=1
        return
      endif
      matname(nmat)=filestiff
      ielmat(1,ne)=nmat
!
      if(mass) then
        nmat=nmat+1
        if(nmat.gt.nmat_) then
          write(*,*) '*ERROR reading *MATRIX ASSEMBLE: increase nmat_'
          ier=1
          return
        endif
        matname(nmat)=filemass
        ielmat(2,ne)=nmat
      else
        ielmat(2,ne)=0
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end







