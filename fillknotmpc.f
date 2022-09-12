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
      subroutine fillknotmpc(co,ipompc,nodempc,coefmpc,labmpc,
     &  nmpc,nmpcold,mpcfree,idim,e1,e2,t1)
!
!     updates the coefficients in nonlinear MPC's
!
      implicit none
!
      character*20 labmpc(*)
!
      integer ipompc(*),nodempc(3,*),irefnode,irotnode,idir,idim,n,
     &  nmpc,index,ii,inode,nmpcold,iexpnode,irefnodeprev,i,ndepnodes,
     &  matz,ier,j,indexnext,node,mpcfree,nodeprev
!
      real*8 co(3,*),coefmpc(*),e(3,3,3),dc(3,3,3) ,sx,sy,sz,sxx,
     &  sxy,sxz,syy,syz,szz,s(3,3),w(3),z(3,3),fv1(3),fv2(3),e1(3),
     &  e2(3),t1(3),u2(3,3),u3(3,3)
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
      irefnodeprev=0
!
      do ii=nmpcold+1,nmpc
         if(labmpc(ii)(1:4).eq.'KNOT') then
!
!           from labmpc: if idim=1: only 2-d elements
!                        if idim=3: at least one 1-d element
!
            irefnode=nodempc(1,nodempc(3,ipompc(ii)))
!
            if(irefnode.ne.irefnodeprev) then
!
!              new knot
!
               irefnodeprev=irefnode
               read(labmpc(ii)(5:5),'(i1)') idim
!
!              determine the area moments of inertia
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
               ndepnodes=0
!
               nodeprev=0
               do i=ii,nmpc
                  if(labmpc(i)(1:4).eq.'KNOT') then
                     if(nodempc(1,nodempc(3,ipompc(i))).eq.irefnode)then
!     
!                       node belonging to the same knot
!
                        node=nodempc(1,ipompc(i))
!
                        if(node.ne.nodeprev) then
                           nodeprev=node
                           ndepnodes=ndepnodes+1
!     
                           sx=sx+co(1,node)
                           sy=sy+co(2,node)
                           sz=sz+co(3,node)
                           sxx=sxx+co(1,node)*co(1,node)
                           sxy=sxy+co(1,node)*co(2,node)
                           sxz=sxz+co(1,node)*co(3,node)
                           syy=syy+co(2,node)*co(2,node)
                           syz=syz+co(2,node)*co(3,node)
                           szz=szz+co(3,node)*co(3,node)
                        endif
                     else
                        exit
                     endif
                  else
                     exit
                  endif
               enddo
!
               sxx=sxx-sx*sx/ndepnodes
               sxy=sxy-sx*sy/ndepnodes
               sxz=sxz-sx*sz/ndepnodes
               syy=syy-sy*sy/ndepnodes
               syz=syz-sy*sz/ndepnodes
               szz=szz-sz*sz/ndepnodes
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
!              determining the eigenvalues
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
!              the eigenvalues are the moments of inertia w.r.t. the
!              plane orthogonal to the eigenvector
!
!              dimension=1 if the two lowest eigenvalues are zero
!              dimension=2 if only the lowest eigenvalue is zero
!              else dimension=3
!
!              the dimension cannot exceed the maximum dimension of
!              the elements connected to the knot (2d-element nodes: 
!              dimension 1, 1d-element nodes: dimension 3)
!
c               write(*,*) 'fillknotmpc eigenvalues ',w(1),w(2),w(3)
               if((w(1).lt.1.d-10).and.(w(2).lt.1.d-10)) then
                  idim=min(idim,1)
c                  idim=1
               elseif(w(1).lt.1.d-10) then
                  idim=min(idim,2)
c                  idim=2
               else
                  idim=min(idim,1)
c                  idim=3
               endif
c               write(*,*) 'fillknotmpc iref= ',irefnode,' idim= ',idim
!
!             defining a local coordinate system for idim=2
!
               if(idim.eq.2) then
                  do i=1,3
                     t1(i)=z(i,1)
                     e2(i)=z(i,2)
                     e1(i)=z(i,3)
                  enddo
!
!                 check whether e1-e2-t1 is a rhs system
!
                  if(t1(1)*(e1(2)*e2(3)-e1(3)*e2(2))-
     &                 t1(2)*(e1(1)*e2(3)-e1(3)*e2(1))+
     &                 t1(3)*(e1(1)*e2(2)-e1(2)*e2(1)).lt.0.d0) then
                     do i=1,3
                        t1(i)=-t1(i)
                     enddo
                  endif
c                  write(*,*) 't1 ',t1(1),t1(2),t1(3)
c                  write(*,*) 'e1 ',e1(1),e1(2),e1(3)
c                  write(*,*) 'e2 ',e2(1),e2(2),e2(3)
!
!                 storing t1 and e1 as coordinates of irotnode and
!                 iexpnode, respectively
!
                  iexpnode=nodempc(1,nodempc(3,nodempc(3,ipompc(ii))))
                  irotnode=
     &            nodempc(1,nodempc(3,nodempc(3,nodempc(3,ipompc(ii)))))
                  do i=1,3
                     co(i,irotnode)=t1(i)
                     co(i,iexpnode)=e1(i)
                  enddo
               endif
            endif
!
            if((idim.eq.1).or.(idim.eq.3)) then
!
!              knot on a line
!
               labmpc(ii)(5:5)='1'
!
!              dependent node
!
               index=ipompc(ii)
               inode=nodempc(1,index)
               idir=nodempc(2,index)
!
!              translation node
!
               index=nodempc(3,index)
               irefnode=nodempc(1,index)
!
!              expansion node
!
               index=nodempc(3,index)
               iexpnode=nodempc(1,index)
               coefmpc(index)=co(idir,irefnode)-co(idir,inode)
!     
!              rotation node
!
               index=nodempc(3,index)
               irotnode=nodempc(1,index)
!
!              determining the coefficients of the rotational degrees
!              of freedom
!
               coefmpc(index)=dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,1)*(co(3,irefnode)-co(3,inode))
!     
               index=nodempc(3,index)
               coefmpc(index)=dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,2)*(co(3,irefnode)-co(3,inode))
!     
               index=nodempc(3,index)
               coefmpc(index)=dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,3)*(co(3,irefnode)-co(3,inode))
!     
            elseif(idim.eq.2) then
!
!              nodes of knot lie in a plane
!
               labmpc(ii)(5:5)='2'
!
               do i=1,3
                  do j=1,3
                     u2(i,j)=2.d0*e1(i)*e1(j)
                     u3(i,j)=2.d0*e2(i)*e2(j)
                  enddo
               enddo
!
!              dependent node
!
               index=ipompc(ii)
               inode=nodempc(1,index)
               idir=nodempc(2,index)
!
!              translation node
!
               index=nodempc(3,index)
               irefnode=nodempc(1,index)
!
!              expansion node (first term is amalgated with the second
!              term since the coefficient matrix is zero)
!
               index=nodempc(3,index)
               iexpnode=nodempc(1,index)
               nodempc(2,index)=2
               coefmpc(index)=0.d0
               indexnext=nodempc(3,index)
               nodempc(3,index)=mpcfree
!
               nodempc(1,mpcfree)=iexpnode
               nodempc(2,mpcfree)=2
               coefmpc(mpcfree)=u2(idir,1)*(co(1,irefnode)-co(1,inode))+
     &              u2(idir,2)*(co(2,irefnode)-co(2,inode))+
     &              u2(idir,3)*(co(3,irefnode)-co(3,inode))
               mpcfree=nodempc(3,mpcfree)
!     
               nodempc(1,mpcfree)=iexpnode
               nodempc(2,mpcfree)=3
               coefmpc(mpcfree)=u3(idir,1)*(co(1,irefnode)-co(1,inode))+
     &              u3(idir,2)*(co(2,irefnode)-co(2,inode))+
     &              u3(idir,3)*(co(3,irefnode)-co(3,inode))
               index=mpcfree
               mpcfree=nodempc(3,mpcfree)
               nodempc(3,index)=indexnext
!
!              rotation node
!
               index=indexnext
               irotnode=nodempc(1,index)
!
!              determining the coefficients of the rotational degrees
!              of freedom
!
               coefmpc(index)=dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,1)*(co(3,irefnode)-co(3,inode))
!     
               index=nodempc(3,index)
               coefmpc(index)=dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,2)*(co(3,irefnode)-co(3,inode))
!     
               index=nodempc(3,index)
               coefmpc(index)=dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &              dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &              dc(idir,3,3)*(co(3,irefnode)-co(3,inode))
!
            endif
         endif
      enddo
!
      return
      end
