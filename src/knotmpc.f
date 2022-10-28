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
      subroutine knotmpc(ipompc,nodempc,coefmpc,irefnode,irotnode,
     &  iexpnode,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,nk,nk_,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,idepnodes,typeboun,co,xboun,istep,
     &  k,ndepnodes,idim,e1,e2,t1)
!
!     generates three knot MPC's for node "node" about reference
!     (translational) node irefnode and rotational node irotnode 
!
      implicit none
!
      character*1 typeboun(*)
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,nk,nk_,
     &  ikmpc(*),idepnodes(*),k,ndepnodes,idim,n,matz,ier,i,
     &  ilmpc(*),node,id,mpcfreeold,j,idof,l,nodeboun(*),
     &  ndirboun(*),ikboun(*),ilboun(*),nboun,nboun_,irefnode,
     &  irotnode,iexpnode,istep,ispcnode,inode
!
      real*8 coefmpc(*),co(3,*),xboun(*),e(3,3,3),dc(3,3,3),s(3,3),
     &  w(3),z(3,3),fv1(3),fv2(3),e1(3),e2(3),t1(3),sx,sy,sz,sxx,
     &  sxy,sxz,syy,syz,szz,u2(3,3),u3(3,3)
!
!     e_ijk symbol
!
      data e /0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,1.d0,0.d0,
     &        0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,-1.d0,0.d0,0.d0,
     &        0.d0,-1.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
!
!     dc_ijk=e_ikj
!
      data dc /0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,-1.d0,0.d0,
     &        0.d0,0.d0,-1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     &        0.d0,1.d0,0.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
!
c      idim=1
!
!     on entry: if idim=1: only 2-d elements
!               if idim=3: at least one 1-d element
!
      if((istep.gt.1).and.(k.eq.1)) then
!
!        determine the dimensionality of the knot
!        (only for step 2 or higher since in step 1 the geometry
!         is not known yet)
!
!        determine the area moments of inertia
!         
         sx=0.d0
         sy=0.d0
         sz=0.d0
         sxx=0.d0
         sxy=0.d0
         sxz=0.d0
         syy=0.d0
         syz=0.d0
         szz=0.d0
!
         do i=1,ndepnodes
            node=idepnodes(i)
            sx=sx+co(1,node)
            sy=sy+co(2,node)
            sz=sz+co(3,node)
            sxx=sxx+co(1,node)*co(1,node)
            sxy=sxy+co(1,node)*co(2,node)
            sxz=sxz+co(1,node)*co(3,node)
            syy=syy+co(2,node)*co(2,node)
            syz=syz+co(2,node)*co(3,node)
            szz=szz+co(3,node)*co(3,node)
         enddo
!
         sxx=sxx-sx*sx/ndepnodes
         sxy=sxy-sx*sy/ndepnodes
         sxz=sxz-sx*sz/ndepnodes
         syy=syy-sy*sy/ndepnodes
         syz=syz-sy*sz/ndepnodes
         szz=szz-sz*sz/ndepnodes
!
c         write(*,*) 'sxx...',sxx,sxy,sxz,syy,syz,szz
!
         s(1,1)=sxx
         s(1,2)=sxy
         s(1,3)=sxz
         s(2,1)=sxy
         s(2,2)=syy
         s(2,3)=syz
         s(3,1)=sxz
         s(3,2)=syz
         s(3,3)=szz
!
!        determining the eigenvalues
!
         n=3
         matz=1
         ier=0
         call rs(n,n,s,w,matz,z,fv1,fv2,ier)
         if(ier.ne.0) then
            write(*,*) '*ERROR in knotmpc while calculating the'
            write(*,*) '       eigenvalues/eigenvectors'
            call exit(201)
         endif
!
!        the eigenvalues are the moments of inertia w.r.t. the
!        plane orthogonal to the eigenvector
!
!        dimension=1 if the two lowest eigenvalues are zero
!        dimension=2 if only the lowest eigenvalue is zero
!        else dimension=3
!
c         write(*,*) 'eigenvalues ',w(1),w(2),w(3)
         if((w(1).lt.1.d-10).and.(w(2).lt.1.d-10)) then
            idim=min(idim,1)
c            idim=1
         elseif(w(1).lt.1.d-10) then
            idim=min(idim,2)
c            idim=2
         else
            idim=min(idim,1)
c            idim=3
         endif
c         write(*,*) 'knotmpc irefnode= ',irefnode,' idim= ',idim
!
!        defining a local coordinate system for idim=2
!
         if(idim.eq.2) then
            do i=1,3
               t1(i)=z(i,1)
               e2(i)=z(i,2)
               e1(i)=z(i,3)
            enddo
!
!           check whether e1-e2-t1 is a rhs system
!
            if(t1(1)*(e1(2)*e2(3)-e1(3)*e2(2))-
     &         t1(2)*(e1(1)*e2(3)-e1(3)*e2(1))+
     &         t1(3)*(e1(1)*e2(2)-e1(2)*e2(1)).lt.0.d0) then
               do i=1,3
                  t1(i)=-t1(i)
               enddo
            endif
c            write(*,*) 't1 ',t1(1),t1(2),t1(3)
c            write(*,*) 'e1 ',e1(1),e1(2),e1(3)
c            write(*,*) 'e2 ',e2(1),e2(2),e2(3)
!
!           storing t1 and e1 as coordinates of irotnode and
!           iexpnode, respectively
!
            do i=1,3
               co(i,irotnode)=t1(i)
               co(i,iexpnode)=e1(i)
            enddo
         endif
      endif
!
      inode=idepnodes(k)
!
      nk=nk+1
      if(nk.gt.nk_) then
         write(*,*) '*ERROR in knotmpc: increase nk_'
         call exit(201)
      endif
!
      ispcnode=nk
!
      if((idim.eq.1).or.(idim.eq.3)) then
!
!     knot nodes lie on a line
!
         do j=1,3
            idof=8*(inode-1)+j
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  cycle
               endif
            endif
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in knotmpc: increase nmpc_'
               call exit(201)
            endif
!     
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)='KNOT                '
            write(labmpc(nmpc)(5:5),'(i1)') idim
!     
            do l=nmpc,id+2,-1
               ikmpc(l)=ikmpc(l-1)
               ilmpc(l)=ilmpc(l-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            nodempc(1,mpcfree)=inode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
!     
!     translation term
!     
            nodempc(1,mpcfree)=irefnode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=-1.d0
            mpcfree=nodempc(3,mpcfree)
!     
!     expansion term
!     
            nodempc(1,mpcfree)=iexpnode
            nodempc(2,mpcfree)=1
            if(istep.gt.1) then
               coefmpc(mpcfree)=co(j,irefnode)-co(j,inode)
            endif
            mpcfree=nodempc(3,mpcfree)
!     
!     rotation terms
!     
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=1
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,1)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,1)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,1)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=2
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,2)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,2)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,2)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=3
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,3)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,3)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,3)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=ispcnode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=1.d0
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            nodempc(3,mpcfreeold)=0
            idof=8*(ispcnode-1)+j
            call nident(ikboun,idof,nboun,id)
            nboun=nboun+1
            if(nboun.gt.nboun_) then
               write(*,*) '*ERROR in knotmpc: increase nboun_'
               call exit(201)
            endif
            nodeboun(nboun)=ispcnode
            ndirboun(nboun)=j
            typeboun(nboun)='R'
            if(istep.gt.1) then
               xboun(nboun)=0.d0
            endif
            do l=nboun,id+2,-1
               ikboun(l)=ikboun(l-1)
               ilboun(l)=ilboun(l-1)
            enddo
            ikboun(id+1)=idof
            ilboun(id+1)=nboun
         enddo
      elseif(idim.eq.2) then
!
!        knot nodes lie in a plane
!
         do i=1,3
            do j=1,3
               u2(i,j)=2.d0*e1(i)*e1(j)
               u3(i,j)=2.d0*e2(i)*e2(j)
            enddo
         enddo
!
         do j=1,3
            idof=8*(inode-1)+j
            call nident(ikmpc,idof,nmpc,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  cycle
               endif
            endif
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) '*ERROR in knotmpc: increase nmpc_'
               call exit(201)
            endif
!     
            ipompc(nmpc)=mpcfree
            labmpc(nmpc)='KNOT2               '
!     
            do l=nmpc,id+2,-1
               ikmpc(l)=ikmpc(l-1)
               ilmpc(l)=ilmpc(l-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
!     
            nodempc(1,mpcfree)=inode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
!     
!           translation term
!     
            nodempc(1,mpcfree)=irefnode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=-1.d0
            mpcfree=nodempc(3,mpcfree)
!     
!           expansion terms
!
!           first term is "removed" (= amalgated with the second term)
!           since u1 is the null matrix
!     
            nodempc(1,mpcfree)=iexpnode
            nodempc(2,mpcfree)=2
            coefmpc(mpcfree)=0.d0
            mpcfree=nodempc(3,mpcfree)
!
            nodempc(1,mpcfree)=iexpnode
            nodempc(2,mpcfree)=2
            coefmpc(mpcfree)=u2(j,1)*(co(1,irefnode)-co(1,inode))+
     &           u2(j,2)*(co(2,irefnode)-co(2,inode))+
     &           u2(j,3)*(co(3,irefnode)-co(3,inode))
            mpcfree=nodempc(3,mpcfree)
!
            nodempc(1,mpcfree)=iexpnode
            nodempc(2,mpcfree)=3
            coefmpc(mpcfree)=u3(j,1)*(co(1,irefnode)-co(1,inode))+
     &           u3(j,2)*(co(2,irefnode)-co(2,inode))+
     &           u3(j,3)*(co(3,irefnode)-co(3,inode))
            mpcfree=nodempc(3,mpcfree)
!
!           rotation terms
!     
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=1
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,1)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,1)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,1)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=2
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,2)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,2)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,2)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=irotnode
            nodempc(2,mpcfree)=3
            if(istep.gt.1) then
               coefmpc(mpcfree)=dc(j,1,3)*(co(1,irefnode)-co(1,inode))+
     &              dc(j,2,3)*(co(2,irefnode)-co(2,inode))+
     &              dc(j,3,3)*(co(3,irefnode)-co(3,inode))
            endif
            mpcfree=nodempc(3,mpcfree)
            nodempc(1,mpcfree)=ispcnode
            nodempc(2,mpcfree)=j
            coefmpc(mpcfree)=1.d0
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            nodempc(3,mpcfreeold)=0
            idof=8*(ispcnode-1)+j
            call nident(ikboun,idof,nboun,id)
            nboun=nboun+1
            if(nboun.gt.nboun_) then
               write(*,*) '*ERROR in knotmpc: increase nboun_'
               call exit(201)
            endif
            nodeboun(nboun)=ispcnode
            ndirboun(nboun)=j
            typeboun(nboun)='R'
            if(istep.gt.1) then
               xboun(nboun)=0.d0
            endif
            do l=nboun,id+2,-1
               ikboun(l)=ikboun(l-1)
               ilboun(l)=ilboun(l-1)
            enddo
            ikboun(id+1)=idof
            ilboun(id+1)=nboun
         enddo
      else
!
!     to do
!
      endif
!     
      return
      end


