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
      subroutine writeinputdeck(nk,co,iponoelfa,inoelfa,konfa,
     &     ipkonfa,lakonfa,nsurfs,iponor,xnor,nodedesiinv,jobnamef,
     &     iponexp,nmpc,labmpc,ipompc,nodempc,ipretinfo,kon,ipkon,lakon,
     &     iponoel,inoel,iponor2d,knor2d,ipoface,nodface,ne,x,y,z,
     &     xo,yo,zo,nx,ny,nz,nodes,dist,ne2d,nod1st,nod2nd3rd,
     &     extnor,nodedesi,ndesi)
!     
!     calculates the normal boundary conditions on the surface for
!     mesh modification purposes in an optimization loop
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
      logical out
!     
      character*8 lakonfa(*),lakon(*),label
      character*20 labmpc(*)
      character*132 jobnamef,fnequ
!     
      integer nk,iponoelfa(*),inoelfa(3,*),konfa(*),ipkonfa(*),nsurfs,
     &     i,index,nexp,nfa,ielem,indexe,j,ifa(100),nopeexp,ixfree1,
     &     jl(100),ial(100),k,l,nemin,jact,ixfree,ixfree3,six,iflag,
     &     node,iponor(*),nodedesiinv(*),len,nsort(6),two,ne,kneigh,
     &     three,iponexp(2,*),nmpc,ipompc(*),nodempc(3,*),indexf,
     &     ipretinfo(*),pretflag,inoel(2,*),nope,idummy,isix,nodeext,
     &     nodepret,kon(*),ipkon(*),iponoel(*),iface,n,nodes(*),
     &     ifaceq(8,6),ifacew(8,5),iposn,iponor2d(2,*),inor(3),
     &     knor2d(*),node2d,ipoface(*),node1,node2,node3,kflag,
     &     nodface(5,*),nopem,ifaceqmid(6),ifacewmid(5),node3d,
     &     nnor1,nnor2,inor1(3),inor2(3),nx(*),ny(*),nz(*),id,
     &     neigh(1),ne2d,nod1st(*),nod2nd3rd(2,*),one,nodedesi(*),
     &     ndesi,ifree,ndist
!
      integer,dimension(:),allocatable::ipo
      integer,dimension(:,:),allocatable::idiste
!     
      real*8 co(3,*),xnor(*),xno(3,100),coloc6(2,6),coloc8(2,8),
     &     xl(3,20),dd,xnoref(3),dot,xnorloc(6),sort(6),x(*),y(*),z(*),
     &     xo(*),yo(*),zo(*),xi,et,ze,shp(4,20),xsj,p(3),dist(*),
     &     distmax,e,emax,extnor(3,*),xnorloc1(3),xnorloc2(3),alpha1,
     &     alpha2
!     
      real*8,dimension(:),allocatable::xdist
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
      one=1
      two=2
      three=3
      six=6
!     
      do len=1,132
        if(jobnamef(len:len).eq.' ') exit
      enddo
      len=len-1
      
      fnequ=jobnamef(1:len)//'.equ'
      open(20,file=fnequ(1:len+4),status='unknown',err=100)
      close(20,status='delete',err=101)
      open(20,file=fnequ(1:len+4),status='unknown',err=100)
!     
      ixfree=0
!     
      do i=1,nk
        index=iponoelfa(i)
        if(index.eq.0) cycle
!     
!     nexp indicates how many different normals there are in the node
!     
        nexp=0
!     
!     locating all external faces to which node i belongs
!     
        nfa=0
        do
          if(index.eq.0) exit
          nfa=nfa+1
          if(nfa.gt.100) then
            write(*,*) '*ERROR in meshmodboun: more '
            write(*,*) '  than 100 shell elements '
            write(*,*) '  share the same node'
            call exit(201)
          endif
          jl(nfa)=inoelfa(2,index)
          ifa(nfa)=inoelfa(1,index)
          index=inoelfa(3,index)
        enddo
!     
        if(nfa.gt.0) then
          do j=1,nfa
            ial(j)=0
          enddo
!     
!     estimate the normal
!     
          do j=1,nfa
            indexf=ipkonfa(ifa(j))
!     
!     local normal on the element (Jacobian)
!     
            if(lakonfa(ifa(j))(2:2).eq.'3') then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,3
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell3(xi,et,xl,xno(1,j))
            elseif(lakonfa(ifa(j))(2:2).eq.'4') then
              xi=coloc8(1,jl(j))
              et=coloc8(2,jl(j))
              do k=1,4
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell4(xi,et,xl,xno(1,j))
            elseif(lakonfa(ifa(j))(2:2).eq.'6') then
              xi=coloc6(1,jl(j))
              et=coloc6(2,jl(j))
              do k=1,6
                node=konfa(indexf+k)
                do l=1,3
                  xl(l,k)=co(l,node)
                enddo
              enddo
              call norshell6(xi,et,xl,xno(1,j))
            elseif(lakonfa(ifa(j))(2:2).eq.'8') then
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
              write(*,*) '*ERROR in meshmodboun: size '
              write(*,*) '       of estimatedshell normal in 
     &node ',i,' element ',ifa(j)
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
            nemin=nsurfs+1
            do j=1,nfa
              if(ial(j).eq.0) then
                if(ifa(j).lt.nemin) then
                  nemin=ifa(j)
                  jact=j
                endif
              endif
            enddo
            if(nemin.eq.nsurfs+1) exit
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
            do j=1,nfa
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
            do j=1,nfa
              if(ial(j).eq.1) then
                do k=1,3
                  xnoref(k)=xnoref(k)+xno(k,j)
                enddo
              endif
            enddo
            dd=dsqrt(xnoref(1)**2+xnoref(2)**2+xnoref(3)**2)
            if(dd.lt.1.d-10) then
              write(*,*) '*ERROR in meshmodboun: size of'
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
            do j=1,nfa
              if(ial(j).eq.1) then
                ial(j)=2
                iponor(ipkonfa(ifa(j))+jl(j))=ixfree
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
        iponexp(2,i)=ixfree-3*nexp
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
!     initializing ipretinfo
!
      do i=1,nk
        ipretinfo(i)=i
      enddo
!
      if(pretflag.eq.1) then
        do i=1,nmpc
          if(labmpc(i)(1:11).eq.'THERMALPRET') cycle
!     
          index=ipompc(i)
          if(index.eq.0) cycle      
          node1=nodempc(1,index)          
          index=nodempc(3,index)
          node2=nodempc(1,index)               
          index=nodempc(3,index)
          node3=nodempc(1,index)
!     
!         the value of ipretinfo for newly generated nodes for 
!         pretension purposes points to the old node (which was duplicated),
!         for all other nodes the value is the node number itself   
!     
          if(node3.eq.nodepret) then
            ipretinfo(node1)=node2 
          endif        
        enddo
      endif
!     
!     write the coordinates in file "jobname.equ"
!
      write(20,102)
 102  format('*NODE,NSET=Nall')
      do i=1,nk
        if(ipretinfo(i).ne.i) cycle
        if(ne2d.ne.0) then
          if((nod2nd3rd(2,i).ne.0).or.(nod1st(i).ne.0)) cycle
        endif
        write(20,103) i,(co(j,i),j=1,3)
      enddo
 103  format(i10,3(',',e15.8))
!     
!     write the topology in file "jobname.equ"
!
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        indexe=ipkon(i)
        if((lakon(i)(7:7).eq.'A').or.(lakon(i)(7:7).eq.'S').or.
     &     (lakon(i)(7:7).eq.'E').or.(lakon(i)(7:7).eq.'L')) then
          if(lakon(i)(4:5).eq.'20') then
            nopeexp=20
            nope=8
            label='CPS8    '
          elseif(lakon(i)(4:5).eq.'15') then
            nopeexp=15
            nope=6
            label='CPS6    '
          elseif(lakon(i)(4:4).eq.'8') then
            nopeexp=8
            nope=4
            label='CPS4    '
          else
            nopeexp=6
            nope=3
            label='CPS3    '
          endif
        else
          nopeexp=0
          if(lakon(i)(4:5).eq.'20') then
            nope=20
            label='C3D20   '
          elseif(lakon(i)(4:5).eq.'15') then
            nope=15
            label='C3D15   '
          elseif(lakon(i)(4:5).eq.'10') then
            nope=10
            label='C3D10   '
          elseif(lakon(i)(4:4).eq.'8') then
            nope=8
            label='C3D8    '
          elseif(lakon(i)(4:4).eq.'6') then
            nope=6
            label='C3D6    '
          elseif(lakon(i)(4:4).eq.'4') then
            nope=4
            label='C3D4    '
          else
            cycle
          endif
        endif
        write(20,107) label
 107    format('*ELEMENT,TYPE=',a8)
        write(20,108) i,(ipretinfo(kon(indexe+nopeexp+j)),j=1,nope)
 108    format(16(i10,','))
      enddo
!
c!     catalogueing all external face nodes in field nodes(*)
c!
c      n=0
c      do i=1,nsurfs
c        do j=ipkonfa(i)+1,ipkonfa(i+1)
c          node=konfa(j)
c          call nident(nodes,node,n,id)
c          if(id.gt.0) then
c            if(nodes(id).eq.node) cycle
c          endif
c          n=n+1
c          do k=n,id+2,-1
c            nodes(k)=nodes(k-1)
c          enddo
c          nodes(id+1)=node
c        enddo
c      enddo
!
!     preparing fields for near3d: considering all design nodes
!
      kneigh=1
      do j=1,ndesi
        xo(j)=co(1,nodedesi(j))
        x(j)=xo(j)
        nx(j)=j
        yo(j)=co(2,nodedesi(j))
        y(j)=yo(j)
        ny(j)=j
        zo(j)=co(3,nodedesi(j))
        z(j)=zo(j)
        nz(j)=j
      enddo
      kflag=2
      call dsort(x,nx,ndesi,kflag)
      call dsort(y,ny,ndesi,kflag)
      call dsort(z,nz,ndesi,kflag)
!
!     determining the center of each element and determine the
!     distance from this center to the closest design node
!
      distmax=0.d0
      iflag=1
      kflag=-2
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        indexe=ipkon(i)
        if(lakon(i)(4:5).eq.'20') then
          nope=20
        elseif(lakon(i)(4:5).eq.'15') then
          nope=15
        elseif(lakon(i)(4:5).eq.'10') then
          nope=10
        elseif(lakon(i)(4:4).eq.'8') then
          nope=8
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
        else
          nope=4
        endif
!     
!     computation of the coordinates of the local nodes
!     
        do j=1,nope
          do k=1,3
            xl(k,j)=co(k,kon(indexe+j))
          enddo
        enddo
!     
!       determine the shape functions at the center of the element     
!     
        if(lakon(i)(4:5).eq.'20') then
          nope=20
          xi=0.d0
          et=0.d0
          ze=0.d0
          call shape20h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(i)(4:5).eq.'15') then
          nope=15
          xi=1.d0/3.d0
          et=1.d0/3.d0
          ze=0.d0
          call shape15w(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(i)(4:5).eq.'10') then
          nope=10
          xi=0.25d0
          et=0.25d0
          ze=0.25d0
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(i)(4:4).eq.'8') then
          nope=8
          xi=0.d0
          et=0.d0
          ze=0.d0
          call shape8h(xi,et,ze,xl,xsj,shp,iflag)
        elseif(lakon(i)(4:4).eq.'6') then
          nope=6
          xi=1.d0/3.d0
          et=1.d0/3.d0
          ze=0.d0
          call shape6w(xi,et,ze,xl,xsj,shp,iflag)
        else
          nope=4
          xi=0.25d0
          et=0.25d0
          ze=0.25d0
          call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
        endif
!
!       determine the global coordinates of the center
!
        p(1)=0.d0
        p(2)=0.d0
        p(3)=0.d0
        do j=1,nope
          do k=1,3
            p(k)=p(k)+shp(4,j)*xl(k,j)
          enddo
        enddo
!     
!     determining the neighboring design node
!     
        call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &       ndesi,neigh,kneigh)
!
        nodeext=nodedesi(neigh(1))
        dist(i)=dsqrt((p(1)-co(1,nodeext))**2+
     &                (p(2)-co(2,nodeext))**2+
     &                (p(3)-co(3,nodeext))**2)
        distmax=max(distmax,dist(i))
      enddo
!
!     determine the Young's modulus depending on the distance from
!     the closest design node; to reduce the number of materials
!     only ndist different E-moduli are taken. To this end the 
!     normalized distance from the nearest design node is sorted
!     into ndist intervals
!
      ndist=10
      allocate(xdist(ndist))
      do i=1,ndist
        xdist(i)=(i-1)/(1.d0*ndist)
      enddo
!
      allocate(ipo(ndist))
      do i=1,ndist
        ipo(i)=0
      enddo
      allocate(idiste(2,ne))
      ifree=0
      do i=1,ne
!
!       normalizing the distance (0 to 1)
!
        dist(i)=dist(i)/distmax
!
!       determining to which of the ndist intervals i belongs
!
        call ident(xdist,dist(i),ndist,id)
        ifree=ifree+1
        idiste(2,ifree)=ipo(id)
        idiste(1,ifree)=i
        ipo(id)=ifree
      enddo
!
!     fictitious Young's modulus
!
      emax=1.d0
      do i=1,ndist
        xi=xdist(i)
        e=(2.d0-xi)*emax/2.d0
        write(20,113) i
        index=ipo(i)
        do
          if(index.eq.0) exit
          write(20,108) idiste(1,index)
          index=idiste(2,index)
        enddo
        write(20,110) i
        write(20,111) e
        write(20,112) i,i
      enddo
 110  format('*MATERIAL,NAME=',i10,'M')
 111  format('*ELASTIC',/,e15.8,',0.3')
 112  format('*SOLID SECTION,ELSET=',i10,'E,MATERIAL=',i10,'M')
 113  format('*ELSET,ELSET=',i10,'E')
!
      deallocate(xdist)
      deallocate(ipo)
      deallocate(idiste)
!      
c!
c!     determine the Young's modulus depending on the distance from
c!     the free surface
c!
c      emax=1.d0
c      do i=1,ne
c        if(ipkon(i).lt.0) cycle
c        xi=dist(i)/distmax
cc        write(*,*) i,dist(i),distmax
c        e=(2.d0-xi)*emax/2.d0
c        write(20,110) i
c        write(20,111) e
c        write(20,112) i,i
c      enddo
c 110  format('*MATERIAL,NAME=',i10,'M')
c 111  format('*ELASTIC',/,e15.8,',0.3')
c 112  format('*SOLID SECTION,ELSET=',i10,'E,MATERIAL=',i10,'M')
!     
!     write equations in file "jobname.equ"
!     
      write(20,109)
 109  format('*EQUATION')
      do i=1,nk
        if((iponoel(i).eq.0).or.(ipretinfo(i).ne.i)) cycle
!     
!     check if node is a designvariable     
!     
        if(nodedesiinv(i).ge.0) then   
!     
!     consideration of plain stress/strain 2d-elements
!     and axisymmetric elements        
!     
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
            elseif (lakon(ielem)(4:4).eq.'6') then
              nope=6
            else
              cycle
            endif
!     
            indexe=ipkon(ielem)
            do j=1,nope
              if(i.eq.kon(indexe+j)) then
                exit
              endif
            enddo
!     
!     replace 3D node number by 2D node number     
!     the only 3D nodes to be replaced are the expanded node with   
!     the lowest node number (of the three expanded nodes; so nodes
!     with local numbers 1,2,3,4,9,10,11,12 for a C3D20(R) element). 
!     The other two are skipped by the statement nodedesiinv(.)=-1  
!     
            if(lakon(ielem)(4:5).eq.'20') then
              if(j.gt.4) j=j-4
              node=kon(indexe+nope+j)
            elseif(lakon(ielem)(4:5).eq.'15') then
              if(j.gt.3) j=j-3
              node=kon(indexe+nope+j)
            elseif(lakon(ielem)(4:4).eq.'8') then
              node=kon(indexe+nope+j)
            elseif(lakon(ielem)(4:5).eq.'6') then
              node=kon(indexe+nope+j)
            endif
!     
!     deactivate the other expansion nodes
!     
            nodedesiinv(i+1)=-1  
            nodedesiinv(i+2)=-1 
!     
!     taking the mean of the normals at expanded node 1
!     and 3
!     
            ixfree1=iponexp(2,i)
            ixfree3=iponexp(2,i+2)
!     
            do j=1,nexp
              do k=1,3
                xnor(ixfree1+3*(j-1)+k)=(xnor(ixfree1+3*(j-1)+k)+
     &               xnor(ixfree3+3*(j-1)+k))/2.d0
              enddo
              dd=dsqrt(xnor(ixfree1+3*(j-1)+1)**2+
     &             xnor(ixfree1+3*(j-1)+2)**2+
     &             xnor(ixfree1+3*(j-1)+3)**2)
              do k=1,3
                xnor(ixfree1+3*(j-1)+k)=xnor(ixfree1+3*(j-1)+k)/dd
              enddo
            enddo
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
          if(nexp.ge.3) then
c     write(*,*) nodedesiinv(i),j,i,extnor(j,i)
            do j=1,3
              if(nodedesiinv(i).eq.0) then
                out=.true.
              elseif(dabs(extnor(j,i)).le.1.d-10) then
                out=.true.
              else
                out=.false.
              endif
              if(out) then
c     do j=1,3
                write(20,106) one
                write(20,105) node,j,1
c     enddo
              endif
            enddo
!     
!     write equations in case nexp is 1
!     
          elseif(nexp.eq.1) then
            if(nodedesiinv(i).eq.0) then
              do j=1,3
                xnorloc(j)=xnor(ixfree+j)
                sort(j)=dabs(xnorloc(j))
                nsort(j)=j
              enddo
              call dsort(sort,nsort,three,two)
              write(20,106) three
              write(20,104) node,nsort(3),xnorloc(nsort(3)),
     &             node,nsort(2),xnorloc(nsort(2)),
     &             node,nsort(1),xnorloc(nsort(1))
            endif
!     
!     write equations in case nexp is 2
!     
          elseif(nexp.eq.2) then
            if(nodedesiinv(i).eq.0) then
!     
!     node is not a design variable: both normal directions     
!     are blocked
!     
              do j=1,6
                xnorloc(j)=xnor(ixfree+j)
                sort(j)=dabs(xnorloc(j))
                nsort(j)=j
              enddo
              call dsort(sort,nsort, six,two)
              nnor1=0
              nnor2=0
!     
!     sorting the two normals apart
!     
              do j=6,1,-1
                if(nsort(j).le.3) then
                  nnor1=nnor1+1
                  inor1(nnor1)=nsort(j)
                  if(j.eq.6) isix=1
                else
                  nnor2=nnor2+1
                  inor2(nnor2)=nsort(j)
                  if(j.eq.6) isix=2
                endif
              enddo
!     
!     check that the dependent dof in both normals is not identical
!     
              if(isix.eq.1) then
                if(inor2(1)-3.eq.inor1(1)) then
                  idummy=inor2(1)
                  inor2(1)=inor2(2)
                  inor2(2)=idummy
                endif
              else
                if(inor1(1)+3.eq.inor2(1)) then
                  idummy=inor1(1)
                  inor1(1)=inor1(2)
                  inor1(2)=idummy
                endif
              endif
!     
              write(20,106) three
              write(20,104) node,inor1(1),xnorloc(inor1(1)),
     &             node,inor1(2),xnorloc(inor1(2)),
     &             node,inor1(3),xnorloc(inor1(3))
!     
              write(20,106) three
              write(20,104) node,inor2(1)-3,xnorloc(inor2(1)),
     &             node,inor2(2)-3,xnorloc(inor2(2)),
     &             node,inor2(3)-3,xnorloc(inor2(3))
            else
!     
!             node is a design variable  
!             storing both normals     
!     
              do j=1,3
                xnorloc1(j)=xnor(ixfree+j)
                xnorloc2(j)=xnor(ixfree+j+3)
              enddo
!
!             scalar product of normals with sensitivity direction,
!             which is stored in extnor
!
              alpha1=extnor(1,i)*xnorloc1(1)+extnor(2,i)*xnorloc1(2)+
     &               extnor(3,i)*xnorloc1(3)
              alpha2=extnor(1,i)*xnorloc2(1)+extnor(2,i)*xnorloc2(2)+
     &             extnor(3,i)*xnorloc2(3)
!
!             sorting the indices for the sensitivity direction
!
              do j=1,3
                sort(j)=dabs(extnor(j,i))
                inor(j)=j
              enddo
              call dsort(sort,inor,three,kflag)
!
!             if the two largest values are equal, the entry with the
!             smallest index is taken as dependent entry
!
              if(dabs(sort(1)-sort(2)).lt.1.e-10) then
                if(inor(2).lt.inor(1)) then
                  idummy=inor(1)
                  inor(1)=inor(2)
                  inor(2)=idummy
                endif
              endif
!
!             determining the direction perpendicular to the
!             sensitivity direction
!
              do j=1,3
                xnorloc(j)=alpha2*xnorloc1(j)-alpha1*xnorloc2(j)
              enddo
              dd=dsqrt(xnorloc(1)**2+xnorloc(2)**2+xnorloc(3)**2)
              do j=1,3
                xnorloc(j)=xnorloc(j)/dd
              enddo
!
!             sorting the indices for the displacement restriction
!             perpendicular to the sensitivity direction
!
              do j=1,3
                sort(j)=dabs(xnorloc(j))
                inor1(j)=j
              enddo
              call dsort(sort,inor1,three,kflag)
!
!             looking for a dependent component which is different
!             from the dependent component (= inor(1)) for the condition in
!             sensitivity direction (which is writing in routine
!             writeinputdeck2.f from the feasible direction procedure)
!
              do j=1,3
                if(inor1(j).ne.inor(1)) exit
              enddo
!
              inor(1)=inor1(j)
              j=j+1
              if(j.gt.3) j=1
              inor(2)=inor1(j)
              j=j+1
              if(j.gt.3) j=1
              inor(3)=inor1(j)
!
              write(20,106) three
              write(20,104) node,inor(1),xnorloc(inor(1)),
     &             node,inor(2),xnorloc(inor(2)),
     &             node,inor(3),xnorloc(inor(3))
            endif
          endif
        endif
!     
      enddo
!     
      do i=1,nk
        if(nodedesiinv(i).eq.-1) then
          nodedesiinv(i)=0
        endif
      enddo
!     
c      close(20)
      return
!     
 100  write(*,*) '*ERROR in openfile: could not open file ',
     &     fnequ(1:len+4)
      call exit(201)
 101  write(*,*) '*ERROR in openfile: could not delete file ',
     &     fnequ(1:len+4) 
      call exit(201)
 104  format(3(i10,",",i1,",",e20.13,","))
 105  format(1(i10,",",i1,",",i1,","))
 106  format(i1)
!     
      end
