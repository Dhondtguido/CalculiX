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
      subroutine gen3dtruss(ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,ikmpc,ilmpc,labmpc,nk,ithermal,i,nodeboun,ndirboun,
     &  ikboun,ilboun,nboun,nboun_,typeboun,xboun,xta,jact,co,
     &  knor,ntrans,inotr,trab,vold,mi,nmethod,nk_,nam,iperturb,
     &  indexk,iamboun,iflagpl)
!
!     - connects the expanded nodes of a truss element to the 
!       original node
!     - sets the rotation about the truss axis to zero
!
      implicit none
!
      logical fixed
!
      character*1 type,typeboun(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,mi(*),
     &  ikmpc(*),ilmpc(*),i,j,idir,nk,newnode,idof,id,mpcfreenew,
     &  ithermal(*),jstart,jend,nodeboun(*),ndirboun(*),ikboun(*),
     &  ilboun(*),nboun,nboun_,jact,knor(*),ntrans,inotr(2,*),
     &  nnodes,nodeact,nmethod,nk_,k,iperturb(*),nam,indexk,
     &  iamplitude,idirref,iamboun(*),iflagpl
!
      real*8 coefmpc(*),xboun(*),xta(3,100),co(3,*),trab(7,*),
     &     vold(0:mi(2)),val
!
!
!
!     generating a hinge at a node of a truss element                
!
!     u(n_1)+u(n_2)+u(n_3)+u(n_4)=4*u(n)
!
      newnode=nk-7
!     
      if(ithermal(2).le.1) then
         jstart=1
         jend=3
      elseif(ithermal(2).eq.2) then
         jstart=0
         jend=0
      else
         jstart=0
         jend=3
      endif
!
      do idir=jstart,jend
         idof=8*(newnode-1)+idir
         call nident(ikmpc,idof,nmpc,id)
         if((id.le.0).or.(ikmpc(id).ne.idof)) then
            nmpc=nmpc+1
            if(nmpc.gt.nmpc_) then
               write(*,*) 
     &              '*ERROR in gen3dtruss: increase nmpc_'
               call exit(201)
            endif
            labmpc(nmpc)='                    '
            ipompc(nmpc)=mpcfree
            do j=nmpc,id+2,-1
               ikmpc(j)=ikmpc(j-1)
               ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id+1)=idof
            ilmpc(id+1)=nmpc
            nodempc(1,mpcfree)=newnode
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=1.d0
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) 
     &              '*ERROR in gen3dtruss: increase memmpc_'
               call exit(201)
            endif
            do k=2,4
               nodempc(1,mpcfree)=nk-8+k
               nodempc(2,mpcfree)=idir
               coefmpc(mpcfree)=1.d0
               mpcfree=nodempc(3,mpcfree)
               if(mpcfree.eq.0) then
                  write(*,*) 
     &                 '*ERROR in gen3dtruss: increase memmpc_'
                  call exit(201)
               endif
            enddo
            nodempc(1,mpcfree)=i
            nodempc(2,mpcfree)=idir
!
!           in the presence of 2D plane strain/stress/axi elements
!           (iflagpl=1) the displacements in z are to be fixed
!
            if((iflagpl.eq.1).and.(idir.eq.3)) then
               coefmpc(mpcfree)=0.d0
            else
               coefmpc(mpcfree)=-4.d0
            endif
            mpcfreenew=nodempc(3,mpcfree)
            if(mpcfreenew.eq.0) then
               write(*,*) 
     &              '*ERROR in gen3dtruss: increase memmpc_'
               call exit(201)
            endif
            nodempc(3,mpcfree)=0
            mpcfree=mpcfreenew
         endif
      enddo
!
!     mean rotation MPC to restrain rotation about the beam
!     axis
!
      label='MEANROTBS           '
!
!     axis of the beam is defined as x-axis in the local beam
!     system (only needed for printing in usermpc.f)
!
      idirref=1
      nnodes=0
      do j=4,1,-1
         nodeact=knor(indexk+j)
         do k=1,3
            nnodes=nnodes+1
            call usermpc(ipompc,nodempc,coefmpc,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &           nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &           nboun,nboun_,nnodes,nodeact,co,label,
     &           typeboun,iperturb,i,idirref,xboun)
         enddo
      enddo
!     
!     rotation value term
!     
      nodeact=nk+1
      do k=1,3
         co(k,nodeact)=xta(k,jact)
      enddo
      nnodes=nnodes+1
      call usermpc(ipompc,nodempc,coefmpc,
     &     labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &     nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &     nboun,nboun_,nnodes,nodeact,co,label,
     &     typeboun,iperturb,i,idirref,xboun)
!     
!     inhomogeneous term
!     
      nodeact=0
      call usermpc(ipompc,nodempc,coefmpc,
     &     labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &     nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &     nboun,nboun_,nnodes,nodeact,co,label,
     &     typeboun,iperturb,i,idirref,xboun)
!     
!     end meanrotationmpc
!     
!     SPC angle term
!     
      if(nodeact.ne.-1) then
         idir=1
         type='B'
         val=0.d0
         iamplitude=0
         fixed=.false.
         call bounadd(nk,idir,idir,val,nodeboun,
     &        ndirboun,xboun,nboun,nboun_,iamboun,
     &        iamplitude,nam,ipompc,nodempc,coefmpc,
     &        nmpc,nmpc_,mpcfree,inotr,trab,ntrans,
     &        ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &        type,typeboun,nmethod,iperturb,fixed,vold,
     &        nk,mi,label)
!     
!     storing the index of the SPC with the angle
!     value in ilboun(id)
!     
         ilboun(id)=nboun
      endif
!     
      return
      end
      
