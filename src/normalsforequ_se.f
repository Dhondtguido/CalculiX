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
      subroutine normalsforequ_se(nk,co,iponoelfa,inoelfa,konfa,
     &  ipkonfa,lakonfa,ne,iponor,xnor,nodedesiinv,jobnamef,
     &  iponexp,nmpc,labmpc,ipompc,nodempc,ipretinfo,kon,ipkon,lakon,
     &  iponoel,inoel,iponor2d,knor2d,iponod2dto3d,ipoface,nodface)
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
      character*8 lakonfa(*),lakonfaloc
      character*20 labmpc(*)
      character*8 lakon(*)
!
      integer nk,iponoelfa(*),inoelfa(3,*),konfa(*),ipkonfa(*),ne,
     &  i,ndepnodes,index,nexp,nel,ielem,indexe,j,iel(100),
     &  jl(100),ial(100),ifi(100),indexx,k,l,ifix,nemin,jact,ixfree,
     &  node,iponor(*),nodedesiinv(*),len,ndet(3),nsort(3),two,
     &  three,iponexp(2,*),nmpc,ipompc(*),nodempc(3,*),
     &  node1,node2,node3,ipretinfo(*),ieq,pretflag,inoel(2,*),nope,
     &  nodepret,ixfreei,ixfreej,kon(*),ipkon(*),iponoel(*),iface,
     &  inode,ifaceq(8,6),ifacew(8,5),iposn,iponor2d(2,*),flag2d,
     &  knor2d(*),node2d,iponod2dto3d(2,*),nopesurf(8),ipoface(*),
     &  nodface(5,*),konl(20),nopem,ifaceqmid(6),ifacewmid(5),node3d
!
      real*8 co(3,*),xnor(*),xno(3,100),xi,et,coloc6(2,6),coloc8(2,8),
     &  xl(3,8),dd,xnoref(3),dot,xnorloc(3,3),det(3),sort(3),xdir,
     &  ydir,zdir
!
!
!
!     In this routine the faces at the free surface play an
!     important role. They are considered to be like a layer of
!     shell elements. Therefore, the term "shell elements" in this
!     routine is basically equivalent to "external faces"
!
      
      data coloc6 /0.d0,0.d0,1.d0,0.d0,0.d0,1.d0,0.5d0,0.d0,
     &             0.5d0,0.5d0,0.d0,0.5d0/
      data coloc8 /-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0,1.d0,
     &            0.d0,-1.d0,1.d0,0.d0,0.d0,1.d0,-1.d0,0.d0/
!
      data ifaceqmid /0,
     &                0,
     &                5,
     &                6,
     &                7,
     &                8/
!
!
!
      data ifacewmid /0,
     &                0,
     &                4,
     &                5,
     &                6/
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
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
!      write(20,103)
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
!     nexp indicates how many times the node was expanded
!     
         nexp=0
!     
!     locating the shell elements to which node i belongs
!     
         nel=0
         do
            if(index.eq.0) exit
            ielem=inoelfa(1,index)
            indexe=ipkonfa(ielem)
            nel=nel+1
            if(nel.gt.100) then
               write(*,*) '*ERROR in normalsforequ_se: more '
               write(*,*) '  than 100 shell elements '
               write(*,*) '  share the same node'
               call exit(201)
            endif
            j=inoelfa(2,index)
            jl(nel)=j
            iel(nel)=ielem
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
               indexe=ipkonfa(iel(j))
               indexx=iponor(indexe+jl(j))
               if(indexx.ge.0) then
                  do k=1,3
                     xno(k,j)=xnor(indexx+k)
                  enddo
                  ifi(j)=1
                  cycle
               else
                  ifi(j)=0
               endif
!     
!     local normal on the element (Jacobian)
!     
               lakonfaloc=lakonfa(iel(j))
               if(lakonfa(iel(j))(2:2).eq.'3') then
                  xi=coloc6(1,jl(j))
                  et=coloc6(2,jl(j))
                  do k=1,3
                     node=konfa(indexe+k)
                     do l=1,3
                        xl(l,k)=co(l,node)
                     enddo
                  enddo
                  call norshell3(xi,et,xl,xno(1,j))
               elseif(lakonfa(iel(j))(2:2).eq.'4') then
                  xi=coloc8(1,jl(j))
                  et=coloc8(2,jl(j))
                  do k=1,4
                     node=konfa(indexe+k)
                     do l=1,3
                        xl(l,k)=co(l,node)
                     enddo
                  enddo
                  call norshell4(xi,et,xl,xno(1,j))
               elseif(lakonfa(iel(j))(2:2).eq.'6') then
                  xi=coloc6(1,jl(j))
                  et=coloc6(2,jl(j))
                  do k=1,6
                     node=konfa(indexe+k)
                     do l=1,3
                        xl(l,k)=co(l,node)
                     enddo
                  enddo
                  call norshell6(xi,et,xl,xno(1,j))
               elseif(lakonfa(iel(j))(2:2).eq.'8') then
                  xi=coloc8(1,jl(j))
                  et=coloc8(2,jl(j))
                  do k=1,8
                     node=konfa(indexe+k)
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
               ifix=0
               nemin=ne+1
               do j=1,nel
                  if(ial(j).ne.0) cycle
                  if(ifi(j).eq.1) then
                     jact=j
                     ifix=1
                     exit
                  endif
               enddo
               if(ifix.eq.0) then
                  do j=1,nel
                     if(ial(j).eq.0) then
                        if(iel(j).lt.nemin) then
                           nemin=iel(j)
                           jact=j
                        endif
                     endif
                  enddo
                  if(nemin.eq.ne+1) exit
               endif
!     
               do j=1,3
                  xnoref(j)=xno(j,jact)
               enddo
!     
!     determining all elements whose normal in node i makes an
!     angle smaller than 0.5 or 20 degrees with the reference normal,
!     depending whether the reference normal was given by the
!     user or is being calculated; the thickness and offset must
!     also fit.
!     
!     if ial(j)=1: normal on element is being treated now
!     
               do j=1,nel
                  if(ial(j).eq.2) cycle
                  if(j.eq.jact) then
                     ial(jact)=1
                  else
                     dot=xno(1,j)*xnoref(1)+xno(2,j)*xnoref(2)+
     &                    xno(3,j)*xnoref(3)
                     if(ifix.eq.0) then
                        if(dot.gt.0.939693d0)then
                           if((lakonfa(iel(j))(1:3).eq.
     &                          lakonfa(iel(jact))(1:3))
     &                          .or.
     &                          ((lakonfa(iel(j))(1:1).eq.'S').and.
     &                          (lakonfa(iel(jact))(1:1).eq.'S')))
     &                          ial(j)=1
                        else
                           if((lakonfa(iel(j))(1:1).eq.'S').and.
     &                          (lakonfa(iel(jact))(1:1).eq.'S')) then
!     
!     if the normals have the opposite
!     direction, the expanded nodes are on a
!     straight line
!     
                              if(dot.le.-0.999962d0) then
                                 write(*,*) '*INFO in gen3dnor: in some 
     &nodes opposite normals are defined'
                              endif
                           endif
                        endif
                     else
                        if(dot.gt.0.999962d0) then
                           if((lakonfa(iel(j))(1:3).eq.
     &                          lakonfa(iel(jact))(1:3))
     &                          .or.
     &                          ((lakonfa(iel(j))(1:1).eq.'S').and.
     &                          (lakonfa(iel(jact))(1:1).eq.'S')))
     &                          ial(j)=1
                        else
                           if((lakonfa(iel(j))(1:1).eq.'S').and.
     &                          (lakonfa(iel(jact))(1:1).eq.'S')) then
!     
!     if the normals have the opposite
!     direction, the expanded nodes are on a
!     straight line
!     
                              if(dot.le.-0.999962d0) then
                                 write(*,*) '*INFO in gen3dnor: in some
     &nodes opposite normals are defined'
                              endif
                           endif
                        endif
                     endif
                  endif
               enddo
!     
!     determining the mean normal for the selected elements
!     
               if(ifix.eq.0) then
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
                     write(*,*) '*ERROR in gen3dnor: size of'
                     write(*,*) '        estimated shell normal is'
                     write(*,*) '        smaller than 1.e-10'
                     call exit(201)
                  endif
                  do j=1,3
                     xnoref(j)=xnoref(j)/dd
                  enddo
               endif
!     
!     updating the pointers iponor
!     
               nexp=nexp+1
               do j=1,nel
                  if(ial(j).eq.1) then
                     ial(j)=2
                     if(ifix.eq.0) then
                        iponor(ipkonfa(iel(j))+jl(j))=ixfree
                     elseif(j.ne.jact) then
                        iponor(ipkonfa(iel(j))+jl(j))=
     &                       iponor(ipkonfa(iel(jact))+jl(jact)) 
                     endif
                  endif
               enddo
!     
!     storing the normal in xnor and generating 3 nodes
!     for knor
!     
               if(ifix.eq.0) then
                  do j=1,3
                     xnor(ixfree+j)=xnoref(j)
                  enddo
                  ixfree=ixfree+3
               endif
!     
            enddo
         endif
!
!     save nexp and ixfree
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
            pretflag=1
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
!     correct nodes on free pretension surface"
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
!        check for additional pretension nodes
!
         if(ipretinfo(i).ne.0) cycle
!
!        check if node is a designvariable     
!
         if(nodedesiinv(i).eq.0) then   
!
!           consideration of plain stress/strain 2d-elements
!           and rotational symmetry elements        
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
                  if((inode.ne.ifaceq(6,3)).and.
     &               (inode.ne.ifaceq(8,3)).and.
     &               (inode.ne.ifaceq(6,5)).and.
     &               (inode.ne.ifaceq(8,5))) cycle
!
!                 replace 3D node number by 2D node number     
!
c                  write(*,*) 'normalsforequ_se ',i,iponoelfa(i)
c                write(*,*) 'normalsforequ_se ',i,inoelfa(1,iponoelfa(i))
c          write(*,*) 'normalsforequ_se ',(konfa(ipkonfa(iface)+j),j=1,4)
c                  iface=inoelfa(1,iponoelfa(i))
c                  if(iface.le.2) cycle
                  node=kon(indexe+inode+4)
                  flag2d=1         
               elseif(lakon(ielem)(4:5).eq.'15') then
                  if((inode.ne.ifacew(6,3)).and.
     &               (inode.ne.ifacew(8,3)).and.
     &               (inode.ne.ifacew(6,4))) cycle
!
!                 replace 3D node number by 2D node number
!
c                  iface=inoelfa(1,iponoelfa(i))
c                  if(iface.le.2) cycle
                  node=kon(indexe+inode+3)
                  flag2d=1
               else
                  cycle
               endif
            elseif(lakon(ielem)(7:7).eq.'L') then
!
!              no output for shell elements necessary
!
               cycle
            else
!
!              in case of a 3D model no change of node number     
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
     &              node,nsort(2),xnorloc(nsort(2),1),
     &              node,nsort(1),xnorloc(nsort(1),1)
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
     &              xnorloc(1,2)*xnorloc(2,1))
               det(2)=dabs(xnorloc(1,1)*xnorloc(3,2)-
     &              xnorloc(1,2)*xnorloc(3,1))
               det(3)=dabs(xnorloc(2,1)*xnorloc(3,2)-
     &              xnorloc(2,2)*xnorloc(3,1))
               call dsort(det,ndet,three,two)
               
               if(ndet(3).eq.1) then
!     if((dabs(xnorloc(1,1)).ge.dabs(xnorloc(2,1))).and.
!     &               (dabs(xnorloc(2,2)).gt.1.d-5)) then
                  if((dabs(xnorloc(1,1)).gt.1.d-5).and.
     &                 (dabs(xnorloc(2,2)).gt.1.d-5)) then
                     write(20,106) 3  
                     write(20,104) node,1,xnorloc(1,1),
     &                    node,2,xnorloc(2,1),node,3,xnorloc(3,1)
                     write(20,106) 3  
                     write(20,104) node,2,xnorloc(2,2),
     &                    node,1,xnorloc(1,2),node,3,xnorloc(3,2)
                  else
                     write(20,106) 3  
                     write(20,104) node,2,xnorloc(2,1),
     &                    node,1,xnorloc(1,1),node,3,xnorloc(3,1)
                     write(20,106) 3  
                     write(20,104) node,1,xnorloc(1,2),
     &                    node,2,xnorloc(2,2),node,3,xnorloc(3,2)
                  endif
               elseif(ndet(3).eq.2) then
!     if((dabs(xnorloc(1,1)).ge.dabs(xnorloc(3,1))).and.
!     &                (dabs(xnorloc(3,2)).gt.1.d-5)) then
                  if((dabs(xnorloc(1,1)).gt.1.d-5).and.
     &                 (dabs(xnorloc(3,2)).gt.1.d-5)) then
                     write(20,106) 3  
                     write(20,104) node,1,xnorloc(1,1),
     &                    node,3,xnorloc(3,1),node,2,xnorloc(2,1)
                     write(20,106) 3  
                     write(20,104) node,3,xnorloc(3,2),
     &                    node,1,xnorloc(1,2),node,2,xnorloc(2,2)
                  else
                     write(20,106) 3  
                     write(20,104) node,3,xnorloc(3,1),
     &                    node,1,xnorloc(1,1),node,2,xnorloc(2,1)
                     write(20,106) 3  
                     write(20,104) node,1,xnorloc(1,2),
     &                    node,3,xnorloc(3,2),node,2,xnorloc(2,2)
                  endif
               elseif(ndet(3).eq.3) then
!     if((dabs(xnorloc(2,1)).ge.dabs(xnorloc(3,1))).and.
!     &                (dabs(xnorloc(3,2)).gt.1.d-5)) then
                  if((dabs(xnorloc(2,1)).gt.1.d-5).and.
     &                 (dabs(xnorloc(3,2)).gt.1.d-5)) then
                     write(20,106) 3  
                     write(20,104) node,2,xnorloc(2,1),
     &                    node,3,xnorloc(3,1),node,1,xnorloc(1,1)
                     write(20,106) 3  
                     write(20,104) node,3,xnorloc(3,2),
     &                    node,2,xnorloc(2,2),node,1,xnorloc(1,2)
                  else
                     write(20,106) 3  
                     write(20,104) node,3,xnorloc(3,1),
     &                    node,2,xnorloc(2,1),node,1,xnorloc(1,1)
                     write(20,106) 3  
                     write(20,104) node,2,xnorloc(2,2),
     &                    node,3,xnorloc(3,2),node,1,xnorloc(1,2)
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
               nodedesiinv(iponod2dto3d(1,i))=-1  
               nodedesiinv(iponod2dto3d(2,i))=-1 
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
         indexe=ipoface(i)
!
         do
            ielem=nodface(3,indexe)
            iface=nodface(4,indexe)
!           
            if((lakon(ielem)(7:7).eq.'A').or.
     &         (lakon(ielem)(7:7).eq.'S').or.
     &           (lakon(ielem)(7:7).eq.'E')) then
!
!              faces in z-direction (expansion direction) do not play
!              a role in the optimization of plane stress/strain/axi
!              elements (corresponds to iface=1 and iface=2)
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
                     indexe=nodface(5,indexe)
                     if(indexe.eq.0) then
                        exit
                     else
                        cycle
                     endif
                  endif
c!     
c!                 actual position of the nodes belonging to the
c!                 surface and surface normal               
c                  do j=1,nope
c                     konl(j)=kon(ipkon(ielem)+j)
c                  enddo
c!
c                  do j=1,nopem
c                     nopesurf(j)=konl(ifaceq(j,iface))
c                     do k=1,3
c                        xl(k,j)=co(k,nopesurf(j))
c                     enddo
c                  enddo
c                  xi=0.d0
c                  et=0.d0
c                  call norshell8(xi,et,xl,xno(1,1))
c                  dd=dsqrt(xno(1,1)**2+xno(2,1)**2+
c     &               xno(3,1)**2)
c                  if(dd.gt.0.d0) then
c                     do j=1,3
c                        xno(j,1)=xno(j,1)/dd
c                     enddo
c                  endif
!     
!                 node number and equation of the 2D node
!                  
                  if(nope.eq.20) then
                     node2d=kon(ipkon(ielem)+nope+ifaceqmid(iface))
                     iposn=iponor2d(2,ipkon(ielem)+ifaceqmid(iface))
                  elseif(nope.eq.15) then
                     node2d=kon(ipkon(ielem)+nope+ifacewmid(iface))
                     iposn=iponor2d(2,ipkon(ielem)+ifacewmid(iface))
                  endif
!
!                 3D-equivalent of the 2D-design variable
!                 The user defines 2D nodes of plane stress/strain/axi
!                 elements as design variables. Internally, these are
!                 replaced by the first node in the 3D expansion
!
                  node3d=knor2d(iposn+1)
!
!                 write the 2D node to file if it is not a design variable
!
                  if(nodedesiinv(node3d).eq.0) then
                     write(20,'(i10,a1)') node2d,','   
!                    write(20,106) 3  
!                    write(20,104) node2d,1,xno(1,1),
!     &                            node2d,2,xno(2,1),
!     &                            node2d,3,xno(3,1)
                  endif         
              endif
            endif      
            indexe=nodface(5,indexe)
            if(indexe.eq.0) exit
!               
         enddo  
      enddo       
!
!-------------------------------------------------------------------------------
! 
!     in case of 2D shell models fix all nodes which are no design variable
!    
!      do i=1,nk
!         if(ipoface(i).eq.0) cycle
!         if(nodedesiinv(i).eq.0) then
!            write(20,'(i10,a1)') i,','
!         endif
!      enddo
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
