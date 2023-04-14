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
      subroutine normalsforequ_se(nk,co,iponoelfa,inoelfa,konfa,
     &     ipkonfa,lakonfa,ne,iponor,xnor,nodedesiinv,jobnamef,
     &     iponexp,nmpc,labmpc,ipompc,nodempc,ipretinfo,kon,ipkon,lakon,
     &     iponoel,inoel,iponor2d,knor2d,nod2nd3rd,ipoface,nodface)
!     
!     calculates normals on surface for mesh modification
!     purposes in an optimization loop
!     
!     during optimization the coordinates of the design variables
!     are changed leading to a changed geometry. In order to keep
!     a good quality mesh the other nodes may have to be moved as
!     well. The external shape in these nodes has to be kept, which
!     can be guaranteed by MPC's. These MPC's are based on the local
!     normal(s) in a node. At sharp corners more than one normal 
!     may be necessary.
!     
!     the equations are stored in file jobname.equ
!     
!     the user can use this file for the appropriate mesh
!     modifications. 
!     
      implicit none
!     
      character*132 jobnamef,fnequ
      character*8 lakonfa(*)
      character*20 labmpc(*)
      character*8 lakon(*)
!     
      integer nk,iponoelfa(*),inoelfa(3,*),konfa(*),ipkonfa(*),ne,
     &     i,ndepnodes,index,nexp,nel,ielem,indexe,j,iel(100),
     &     jl(100),ial(100),k,l,nemin,jact,ixfree,
     &     node,iponor(*),nodedesiinv(*),len,ndet(3),nsort(3),two,
     &     three,iponexp(2,*),nmpc,ipompc(*),nodempc(3,*),indexf,
     &     node1,node2,node3,ipretinfo(*),ieq,pretflag,inoel(2,*),nope,
     &     nodepret,ixfreei,ixfreej,kon(*),ipkon(*),iponoel(*),iface,
     &     inode,ifaceq(8,6),ifacew(8,5),iposn,iponor2d(2,*),flag2d,
     &     knor2d(*),node2d,nod2nd3rd(2,*),nopesurf(8),ipoface(*),
     &     nodface(5,*),konl(20),nopem,ifaceqmid(6),ifacewmid(5),node3d
!     
      real*8 co(3,*),xnor(*),xno(3,100),xi,et,coloc6(2,6),coloc8(2,8),
     &     xl(3,8),dd,xnoref(3),dot,xnorloc(3,3),det(3),sort(3),xdir,
     &     ydir,zdir
!     
!     In this routine the faces at the free surface play an
!     important role. They are considered to be like a layer of
!     shell elements. Therefore, the term "shell elements" in this
!     routine is basically equivalent to "external faces"
!     
      data coloc6 /0.d0,0.d0,1.d0,0.d0,0.d0,1.d0,0.5d0,0.d0,
     &     0.5d0,0.5d0,0.d0,0.5d0/
      data coloc8 /-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &     0.d0,-1.d0,1.d0,0.d0,0.d0,1.d0,-1.d0,0.d0/
!     
      data ifaceqmid /0,
     &     0,
     &     5,
     &     6,
     &     7,
     &     8/
!     
      data ifacewmid /0,
     &     0,
     &     4,
     &     5,
     &     6/
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
      two=2
      three=3
!     
      do len=1,132
        if(jobnamef(len:len).eq.' ') exit
      enddo
      len=len-1
      
      fnequ=jobnamef(1:len)//'.equ'
      open(20,file=fnequ(1:len+4),status='unknown',err=100)
      close(20,status='delete',err=101)
      open(20,file=fnequ(1:len+4),status='unknown',err=100)
      write(20,102)
!     write(20,103)
 102  format('**SUMMARY OF EQUATIONS FOR MESH-UPDATE')
 103  format('*EQUATION')
!     
      ixfree=0
!     
      do i=1,nk
        ndepnodes=0
        index=iponoelfa(i)
        if(index.eq.0) cycle
!     
!     nexp indicates how many different normals there are in the node
!     
        nexp=0
!     
!     locating all external faces to which node i belongs
!     
        nel=0
        do
          if(index.eq.0) exit
          nel=nel+1
          if(nel.gt.100) then
            write(*,*) '*ERROR in normalsforequ_se: more '
            write(*,*) '  than 100 shell elements '
            write(*,*) '  share the same node'
            call exit(201)
          endif
          jl(nel)=inoelfa(2,index)
          iel(nel)=inoelfa(1,index)
          index=inoelfa(3,index)
        enddo
!     
        if(nel.gt.0) then
          do j=1,nel
            ial(j)=0
          enddo
!     
!     estimate the normal
!     
          do j=1,nel
            indexf=ipkonfa(iel(j))
!     
!     local normal on the element (Jacobian)
!     
            if(lakonfa(iel(j))(2:2).eq.'3') then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,3
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell3(xi,et,xl,xno(1,j))
            elseif(lakonfa(iel(j))(2:2).eq.'4') then
              xi=coloc8(1,jl(j))
              et=coloc8(2,jl(j))
              do k=1,4
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell4(xi,et,xl,xno(1,j))
            elseif(lakonfa(iel(j))(2:2).eq.'6') then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,6
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell6(xi,et,xl,xno(1,j))
            elseif(lakonfa(iel(j))(2:2).eq.'8') then
              xi=coloc8(1,jl(j))
              et=coloc8(2,jl(j))
              do k=1,8
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell8(xi,et,xl,xno(1,j))
            endif
!     
            dd=dsqrt(xno(1,j)**2+xno(2,j)**2+xno(3,j)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in normalsforequ_se: size '
              write(*,*) '       of estimatedshell normal in 
     &node ',i,' element ',iel(j)
              write(*,*) '       is smaller than 1.e-10'
              call exit(201)
            endif
            do k=1,3
              xno(k,j)=xno(k,j)/dd
            enddo
          enddo
!     
          do
!     
!     determining a fixed normal which was not treated yet,
!     or, if none is left, the minimum element number of all
!     elements containing node i and for which no normal was
!     determined yet
!     
!     if ial(j)=0: the normal on this element has not been
!     treated yet
!     if ial(j)=2: normal has been treated
!     
            nemin=ne+1
            do j=1,nel
              if(ial(j).eq.0) then
                if(iel(j).lt.nemin) then
                  nemin=iel(j)
                  jact=j
                endif
              endif
            enddo
            if(nemin.eq.ne+1) exit
!     
            do j=1,3
              xnoref(j)=xno(j,jact)
            enddo
!     
!     determining all elements whose normal in node i makes an
!     angle smaller than 20 degrees with the reference normal
!     
!     if ial(j)=1: normal on element is being treated now
!     
            do j=1,nel
              if(ial(j).eq.2) cycle
              if(j.eq.jact) then
                ial(jact)=1
              else
                dot=xno(1,j)*xnoref(1)+xno(2,j)*xnoref(2)+
     &               xno(3,j)*xnoref(3)
                if(dot.gt.0.939693d0) ial(j)=1
              endif
            enddo
!     
!     determining the mean normal for the selected elements
!     
            do j=1,3
              xnoref(j)=0.d0
            enddo
            do j=1,nel
              if(ial(j).eq.1) then
                do k=1,3
                  xnoref(k)=xnoref(k)+xno(k,j)
                enddo
              endif
            enddo
            dd=dsqrt(xnoref(1)**2+xnoref(2)**2+xnoref(3)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in normalsforequ_se: size of'
              write(*,*) '        estimated face normal is'
              write(*,*) '        smaller than 1.e-10'
              call exit(201)
            endif
            do j=1,3
              xnoref(j)=xnoref(j)/dd
            enddo
!     
!     updating the pointers iponor
!     
            nexp=nexp+1
            do j=1,nel
              if(ial(j).eq.1) then
                ial(j)=2
                iponor(ipkonfa(iel(j))+jl(j))=ixfree
              endif
            enddo
!     
!     storing the normal in xnor
!     
            do j=1,3
              xnor(ixfree+j)=xnoref(j)
            enddo
            ixfree=ixfree+3
!     
          enddo
        endif
!     
!     save nexp (number of normals in node i) and ixfree (pointer to
!     normals for node i+1
!     
        iponexp(1,i)=nexp
        iponexp(2,i)=ixfree
!     
      enddo     
!     
!     find nodes created by "*PRETENSION SECTION"
!     
!     find pretension node if existing
!     
      pretflag=0
      do i=1,nmpc
        if(labmpc(i)(1:10).eq.'PRETENSION') then
          pretflag=1
          index=ipompc(i)
          index=nodempc(3,index)
          index=nodempc(3,index)
          nodepret=nodempc(1,index)
          exit
        endif
      enddo
!     
      if(pretflag.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:11).eq.'THERMALPRET') cycle
!     
          ieq=0
          index=ipompc(i)
          if(index.eq.0) cycle      
          node1=nodempc(1,index)          
          index=nodempc(3,index)
          node2=nodempc(1,index)               
          index=nodempc(3,index)
          node3=nodempc(1,index)
          if(node3.eq.nodepret) then
            ipretinfo(node2)=node1 
            ipretinfo(node1)=-1      
          endif        
        enddo
      endif
!     
!     correct nodes on free pretension surface
!     
      do i=1,nk
        if(ipretinfo(i).le.0) cycle 
!     
        nexp=iponexp(1,i)
        ixfreei=iponexp(2,i)
        ixfreej=iponexp(2,ipretinfo(i))
!     
        do j=1,nexp
          k=j*3-3
          zdir=xnor(ixfreei+1-1-k)+xnor(ixfreej+1-1-k)
          ydir=xnor(ixfreei+1-2-k)+xnor(ixfreej+1-2-k)
          xdir=xnor(ixfreei+1-3-k)+xnor(ixfreej+1-3-k)      
          dd=(xdir)**2+(ydir)**2+(zdir)**2
!     
          if(dd.gt.1.0e-12) then      
            ipretinfo(i)=0
          endif
!     
        enddo   
!     
      enddo
!     
!---------------------------------------------------------------------------
!     
!     write equations in file "jobname.equ"
!     in case of a 2D model just write the node numbers in the file
!     
      do i=1,nk
        flag2d=0
!     
!     check for additional pretension nodes
!     
        if(ipretinfo(i).ne.0) cycle
!     
!     check if node is a designvariable     
!     
        if(nodedesiinv(i).eq.0) then   
!     
!     consideration of plain stress/strain 2d-elements
!     and rotational symmetry elements        
!     
          if(iponoel(i).eq.0) cycle
          ielem=inoel(1,iponoel(i))
          if((lakon(ielem)(7:7).eq.'A').or.
     &         (lakon(ielem)(7:7).eq.'S').or.
     &         (lakon(ielem)(7:7).eq.'E')) then
!     
            if(lakon(ielem)(4:5).eq.'20') then
              nope=20
            elseif (lakon(ielem)(4:4).eq.'8') then
              nope=8
            elseif (lakon(ielem)(4:5).eq.'15') then
              nope=15
            else
              cycle
            endif
!     
            indexe=ipkon(ielem)
            do inode=1,nope
              if(i.eq.kon(indexe+inode)) then
                exit
              endif
            enddo
            if(lakon(ielem)(4:5).eq.'20') then
              if((inode.ne.17).and.
     &             (inode.ne.18).and.
     &             (inode.ne.19).and.
     &             (inode.ne.20)) cycle
!     
!     replace 3D node number by 2D node number     
!     
              node=kon(indexe+inode+4)
              flag2d=1         
            elseif(lakon(ielem)(4:5).eq.'15') then
              if((inode.ne.13).and.
     &             (inode.ne.14).and.
     &             (inode.ne.15)) cycle
!     
!     replace 3D node number by 2D node number
!     
              node=kon(indexe+inode+3)
              flag2d=1
            else
              cycle
            endif
          elseif(lakon(ielem)(7:7).eq.'L') then
!     
!     no output for shell elements necessary
!     
            cycle
          else
!     
!     in case of a 3D model no change of node number     
            node=i
          endif
!     
!     write equations in case nexp is greater or equal 3
!     
          nexp=iponexp(1,i)
          ixfree=iponexp(2,i)
!     
          if((nexp.ge.3).and.(flag2d.eq.0)) then
            do j=1,3
              write(20,106) 1
              write(20,105) node,j,1
            enddo
!     
!     write equations in case nexp is 1
!     
          elseif((nexp.eq.1).and.(flag2d.eq.0)) then
            j=1
            do l=1,3
              xnorloc(4-l,j)=xnor(ixfree+1-l)
              sort(4-l)=dabs(xnor(ixfree+1-l))
              nsort(4-l)=4-l            
            enddo
            call dsort(sort,nsort,three,two)
            write(20,106) 3  
            write(20,104) node,nsort(3),xnorloc(nsort(3),1),
     &           node,nsort(2),xnorloc(nsort(2),1),
     &           node,nsort(1),xnorloc(nsort(1),1)
!     
!     write equations in case nexp is 2
!     
          elseif((nexp.eq.2).and.(flag2d.eq.0)) then
            do j=1,nexp
              k=j*3-3
              do l=1,3
                xnorloc(4-l,j)=xnor(ixfree+1-l-k)
              enddo
            enddo
            ndet(1)=1
            ndet(2)=2
            ndet(3)=3
            det(1)=dabs(xnorloc(1,1)*xnorloc(2,2)-
     &           xnorloc(1,2)*xnorloc(2,1))
            det(2)=dabs(xnorloc(1,1)*xnorloc(3,2)-
     &           xnorloc(1,2)*xnorloc(3,1))
            det(3)=dabs(xnorloc(2,1)*xnorloc(3,2)-
     &           xnorloc(2,2)*xnorloc(3,1))
            call dsort(det,ndet,three,two)
            
            if(ndet(3).eq.1) then
              if((dabs(xnorloc(1,1)).gt.1.d-5).and.
     &             (dabs(xnorloc(2,2)).gt.1.d-5)) then
                write(20,106) 3  
                write(20,104) node,1,xnorloc(1,1),
     &               node,2,xnorloc(2,1),node,3,xnorloc(3,1)
                write(20,106) 3  
                write(20,104) node,2,xnorloc(2,2),
     &               node,1,xnorloc(1,2),node,3,xnorloc(3,2)
              else
                write(20,106) 3  
                write(20,104) node,2,xnorloc(2,1),
     &               node,1,xnorloc(1,1),node,3,xnorloc(3,1)
                write(20,106) 3  
                write(20,104) node,1,xnorloc(1,2),
     &               node,2,xnorloc(2,2),node,3,xnorloc(3,2)
              endif
            elseif(ndet(3).eq.2) then
              if((dabs(xnorloc(1,1)).gt.1.d-5).and.
     &             (dabs(xnorloc(3,2)).gt.1.d-5)) then
                write(20,106) 3  
                write(20,104) node,1,xnorloc(1,1),
     &               node,3,xnorloc(3,1),node,2,xnorloc(2,1)
                write(20,106) 3  
                write(20,104) node,3,xnorloc(3,2),
     &               node,1,xnorloc(1,2),node,2,xnorloc(2,2)
              else
                write(20,106) 3  
                write(20,104) node,3,xnorloc(3,1),
     &               node,1,xnorloc(1,1),node,2,xnorloc(2,1)
                write(20,106) 3  
                write(20,104) node,1,xnorloc(1,2),
     &               node,3,xnorloc(3,2),node,2,xnorloc(2,2)
              endif
            elseif(ndet(3).eq.3) then
              if((dabs(xnorloc(2,1)).gt.1.d-5).and.
     &             (dabs(xnorloc(3,2)).gt.1.d-5)) then
                write(20,106) 3  
                write(20,104) node,2,xnorloc(2,1),
     &               node,3,xnorloc(3,1),node,1,xnorloc(1,1)
                write(20,106) 3  
                write(20,104) node,3,xnorloc(3,2),
     &               node,2,xnorloc(2,2),node,1,xnorloc(1,2)
              else
                write(20,106) 3  
                write(20,104) node,3,xnorloc(3,1),
     &               node,2,xnorloc(2,1),node,1,xnorloc(1,1)
                write(20,106) 3  
                write(20,104) node,2,xnorloc(2,2),
     &               node,3,xnorloc(3,2),node,1,xnorloc(1,2)
              endif     
            endif
!     
!     WORKAROUND: MPC's in combination with expanded 2D models does not work
!     in case of expanded 2D models create a set with all surface nodes
!     which are not in the designvariables set. These nodes are fully
!     constrained
!     
          elseif(flag2d.eq.1) then
            write(20,'(i10,a1)') node,','   
          endif          
        elseif(nodedesiinv(i).eq.1) then
          if(iponoel(i).eq.0) cycle
          ielem=inoel(1,iponoel(i))
          if((lakon(ielem)(7:7).eq.'A').or.
     &         (lakon(ielem)(7:7).eq.'S').or.
     &         (lakon(ielem)(7:7).eq.'L').or.
     &         (lakon(ielem)(7:7).eq.'E')) then
!     
            nodedesiinv(nod2nd3rd(1,i))=-1  
            nodedesiinv(nod2nd3rd(2,i))=-1 
          endif
        endif
!     
      enddo
!     
!-------------------------------------------------------------------------------
!     
!     in case of plain strain/stress/axi 2D models write midnodes belonging
!     to these 2D elements to file. This naturally only applies to
!     quadratic elements
!     
      do i=1,nk
        
        if(ipoface(i).eq.0) cycle
        indexf=ipoface(i)
!     
        do
          ielem=nodface(3,indexf)
          iface=nodface(4,indexf)
!     
          if((lakon(ielem)(7:7).eq.'A').or.
     &         (lakon(ielem)(7:7).eq.'S').or.
     &         (lakon(ielem)(7:7).eq.'E')) then
!     
!     faces in z-direction (expansion direction) do not play
!     a role in the optimization of plane stress/strain/axi
!     elements (corresponds to iface=1 and iface=2)
!     
            if(iface.gt.2) then
!     
              if(lakon(ielem)(4:5).eq.'20') then
                nope=20
                nopem=8
              elseif(lakon(ielem)(4:5).eq.'15') then
                nope=15
                nopem=8
              else
                indexf=nodface(5,indexf)
                if(indexf.eq.0) then
                  exit
                else
                  cycle
                endif
              endif
!     
!     node number and equation of the 2D node
!     
              if(nope.eq.20) then
                node2d=kon(ipkon(ielem)+nope+ifaceqmid(iface))
                iposn=iponor2d(2,ipkon(ielem)+ifaceqmid(iface))
              elseif(nope.eq.15) then
                node2d=kon(ipkon(ielem)+nope+ifacewmid(iface))
                iposn=iponor2d(2,ipkon(ielem)+ifacewmid(iface))
              endif
!     
!     3D-equivalent of the 2D-design variable
!     The user defines 2D nodes of plane stress/strain/axi
!     elements as design variables. Internally, these are
!     replaced by the first node in the 3D expansion
!     
              node3d=knor2d(iposn+1)
!     
!     write the 2D node to file if it is not a design variable
!     
              if(nodedesiinv(node3d).eq.0) then
                write(20,'(i10,a1)') node2d,','   
              endif         
            endif
          endif      
          indexf=nodface(5,indexf)
          if(indexf.eq.0) exit
!     
        enddo  
      enddo       
!     
      do i=1,nk
        if(nodedesiinv(i).eq.-1) then
          nodedesiinv(i)=0
        endif
      enddo
!     
      close(20)
      return
!     
 104  format(3(i10,",",i1,",",e20.13,","))
 105  format(1(i10,",",i1,",",i1,","))
 106  format(i1)
 107  format(2(i10,",",i1,",",e20.13,",")) 
!     
 100  write(*,*) '*ERROR in openfile: could not open file ',
     &     fnequ(1:len+4)
      call exit(201)
 101  write(*,*) '*ERROR in openfile: could not delete file ',
     &     fnequ(1:len+4) 
      call exit(201)
!     
      end
