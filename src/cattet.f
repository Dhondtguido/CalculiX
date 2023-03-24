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
      subroutine cattet(kontet,netet_,ifac,ne,ipkon,kon,ifatet,ifreetet,
     &     bc,itetfa,ifreefa,planfa,ipofa,cotet,cg,ipoeln,ieln,ifreeln,
     &     lakon,kontetor,iquad,istartset,iendset,ialset,set,nset,filab,
     &     jfix,iparentel,jobnamec)
!     
!     catalogueing the tetrahedral elements of the mesh
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 set(*),elset
      character*87 filab(*)
      character*132 fnrfn,jobnamec(*)
      character*256 fn
!     
      integer kontet(4,*),netet_,i,ifac(4,*),ne,ipkon(*),kon(*),index1,
     &     j,nodes(4),ifatet(4,*),ifreetet,itetfa(2,*),ifreefa,ipofa(*),
     &     ipoeln(*),ieln(2,*),node,ifreeln,kontetor(6,*),iquad,nset,
     &     istartset(*),iendset(*),ialset(*),indexe,k,nope,jfix(*),
     &     iparentel(*),id,ilen
!     
      real*8 bc(4,*),planfa(4,*),cotet(3,*),cg(3,*)
!
!     open a file for nodes which are not properly projected on the
!     free surface      
!
      ilen=index(jobnamec(1),char(0))-1
      fn=jobnamec(1)(1:ilen)//'_WarnNodeNotProjected.nam'
      open(40,file=fn,status='unknown')
      write(40,*) '*NSET,NSET=WarnNodeNotProjected'
!
!     read the element set name to refine, if any
!     default: all tetrahedral elements are refined      
!
      read(filab(48)(27:87),'(a61)') elset(1:61)
!
      do i=62,81
        elset(i:i)=' '
      enddo
!
!     determining the number of the set
!
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
!
!     if element set detected:
!
      if(i.le.nset) then
!
!       replace 'C' by 'A' in all tet elements
!
        do j=1,ne
          if((lakon(j)(1:4).eq.'C3D4').or.
     &         (lakon(j)(1:5).eq.'C3D10')) then
            lakon(j)(1:1)='A'
          endif
        enddo
!
!     replace 'A' by 'C' in all tet elements belonging to the set
!     to refine        
!
        do j=istartset(i),iendset(i)
          if(ialset(j).gt.0) then
            k=ialset(j)
            lakon(k)(1:1)='C'
          else
            k=ialset(j-2)
            do
              k=k-ialset(j)
              if(k.ge.ialset(j-1)) exit
              lakon(k)(1:1)='C'
            enddo
          endif
        enddo
!
!     setting field jfix in all VERTEX nodes belonging to
!     elements not to be refined to 1 (including the common
!     boundaries between the set to be refined and its complement)        
!
        do i=1,ne
          if(lakon(i)(1:5).eq.'C3D8I') then
            nope=8
          elseif(lakon(i)(4:5).eq.'20') then
            nope=8
          elseif(lakon(i)(4:4).eq.'8') then
            nope=8
          elseif(lakon(i)(1:5).eq.'C3D10') then
            cycle
          elseif(lakon(i)(1:4).eq.'C3D4') then
            cycle
          elseif(lakon(i)(4:5).eq.'15') then
            nope=6
          elseif(lakon(i)(4:4).eq.'6') then
            nope=6
          elseif(lakon(i)(1:2).eq.'ES') then
            if(lakon(i)(7:7).eq.'C') then
              cycle
            else
              nope=ichar(lakon(i)(8:8))-47
            endif
          elseif(lakon(i)(1:4).eq.'MASS') then
            nope=1
          elseif(lakon(i)(1:1).eq.'A') then
            nope=4
          else
            cycle
          endif
!
          indexe=ipkon(i)
          do j=1,nope
            jfix(kon(indexe+j))=1
          enddo
        enddo
      endif
!     
!     determine the first unused element
!     
      do i=1,ne
        if(lakon(i)(1:1).eq.char(0)) exit
      enddo
!      
      ifreetet=i
!     
      do 
        do j=i+1,netet_
!     
!     element number supersedes largest one
!     
          if(j.gt.ne) then
            kontet(4,i)=j
            i=j
            exit
          endif
!     
!     unused element
!     
          if(lakon(j)(1:1).eq.char(0)) then
            kontet(4,i)=j
            i=j
            exit
          endif
        enddo
        if(j.eq.netet_) exit
      enddo
      kontet(4,netet_)=0
!     
!     initialization of ipofa and ifac
!     
      do i=1,4*netet_
        ifac(4,i)=i+1
      enddo
      ifac(4,4*netet_)=0
!     
!     initialization of ieln
!     
      do i=1,4*netet_
        ieln(2,i)=i+1
      enddo
      ieln(2,4*netet_)=0
!     
!     adding the tetrahedral elements one by one in kontet
!     tagging the original elements in kon to be removed
!     by using a model change in input deck jobname.rfn.inp
!     
!     creating the file name
!
      do i=1,132
        if(ichar(jobnamec(1)(i:i)).eq.0) exit
      enddo
      if(i.gt.125) then
        write(*,*) '*ERROR in cattet'
        write(*,*) '       jobname has more than 124 characters:'
        write(*,*) jobnamec(1)(1:132)
        call exit(201)
      endif
      fnrfn(1:i+7)=jobnamec(1)(1:i-1)//'.rfn.inp'
!
!     opening the file
!
      open(2,file=fnrfn(1:i+7),status='unknown')
!
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if((lakon(i)(1:4).ne.'C3D4').and.
     &       (lakon(i)(1:5).ne.'C3D10')) cycle
        index1=ipkon(i)
        do j=1,4
          nodes(j)=kon(index1+j)
        enddo
        write(2,101)
        write(2,102) i
!     
!     if C3D10: store the middle nodes
!     If there is at least one C3D10 element iquad is set to 1
!     which means that the refined mesh will be fully quadratic
!     
        if(lakon(i)(4:4).eq.'1') then
          iquad=1
          do j=1,6
            kontetor(j,ifreetet)=kon(index1+4+j)
          enddo
        endif
!
!       defining the initial parent element (= element itself)
!        
        iparentel(ifreetet)=i
        call generatetet_refine(kontet,ifatet,ifreetet,bc,ifac,itetfa,
     &       ifreefa,planfa,ipofa,nodes,cotet,cg)
      enddo
!
 101  format('*MODEL CHANGE,TYPE=ELEMENT,REMOVE')
 102  format(i10)
      close(2)
!     
!     generating the element per node relationship
!     
      do j=1,netet_
        if(kontet(1,j).eq.0) cycle
        do i=1,4
          node=kontet(i,j)
          index1=ifreeln
          ieln(1,index1)=j
          ifreeln=ieln(2,index1)
          if(ifreeln.eq.0) then
            write(*,*) '*ERROR in cattet: increase the'
            write(*,*) '       dimension of ieln'
            call exit(201)
          endif
          ieln(2,index1)=ipoeln(node)
          ipoeln(node)=index1
        enddo
      enddo
!     
      return
      end

