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
      subroutine maxdesvardisp(nobject,nk,set,nset,istartset,iendset,
     &     ialset,iobject,nodedesiinv,dgdxglob,objectset,xdesi,coini,
     &     co,nodedesipos,ndesi,nodedesi,g0,extnorini)                       
!     
!     determination of the fixed nodes for the sensitivity analysis     
!     
      implicit none
!     
      character*81 objectset(5,*),set(*)
!     
      integer nk,nobject,nset,istartset(*),iendset(*),ialset(*),id,
     &     iobject,nodedesiinv(*),i,j,node,istat,nodedesipos(*),
     &     ndesi,nodedesi(*)
!     
      real*8 dgdxglob(2,nk,*),coini(3,*),xdesi(3,*),co(3,*),bound,
     &   scalprod,g0(*),dispvector(3),extnorini(3,*)
!
      real*4 actmove
!
!
!     determine position of designvaribles in nodedesi
!
      do i=1,ndesi
         node=nodedesi(i)
         nodedesipos(node)=i
      enddo
!     
!     determining the set of bounded design variables
!     
      call cident81(set,objectset(3,iobject),nset,id)
      i=nset+1
      if(id.gt.0) then
        if(objectset(3,iobject).eq.set(id)) then
          i=id
        endif
      endif
!     
      read(objectset(1,iobject)(61:80),'(f20.0)',iostat=istat) bound
!
      if(i.le.nset) then    
        do j=istartset(i),iendset(i)
          if(ialset(j).gt.0) then
            node=ialset(j)
            if(nodedesiinv(node).eq.1) then
!           
!              calculate actual movement of design variable
!
               dispvector(1)=co(1,node)-coini(1,node)
               dispvector(2)=co(2,node)-coini(2,node)
               dispvector(3)=co(3,node)-coini(3,node)
               actmove=dsqrt((dispvector(1)*extnorini(1,node))**2+
     &                       (dispvector(2)*extnorini(2,node))**2+
     &                       (dispvector(3)*extnorini(3,node))**2)
               if(actmove.lt.1.0e-8) then
                  actmove=0.d0
                  dispvector(1)=0.d0
                  dispvector(2)=0.d0
                  dispvector(3)=0.d0
               endif
               scalprod=dispvector(1)*xdesi(1,nodedesipos(node))+
     &                  dispvector(2)*xdesi(2,nodedesipos(node))+
     &                  dispvector(3)*xdesi(3,nodedesipos(node))
               if(scalprod.le.0.d0) then
                  dgdxglob(1,node,iobject)=-actmove
               else
                  dgdxglob(1,node,iobject)=actmove
               endif
               if((objectset(1,iobject)(1:9).eq.'MAXGROWTH').and.
     &            (scalprod.ge.0d0)) then
                  dgdxglob(2,node,iobject)=actmove-bound
               elseif((objectset(1,iobject)(1:9).eq.'MAXGROWTH').and.
     &            (scalprod.lt.0d0)) then
                  dgdxglob(2,node,iobject)=-actmove-bound
               elseif((objectset(1,iobject)(1:12).eq.'MAXSHRINKAGE')
     &            .and.(scalprod.le.0d0)) then
                  dgdxglob(2,node,iobject)=actmove-bound
               elseif((objectset(1,iobject)(1:12).eq.'MAXSHRINKAGE')
     &            .and.(scalprod.gt.0d0)) then
                  dgdxglob(2,node,iobject)=-actmove-bound
               endif
! 
!              count number of active nodes
!
               if(dgdxglob(2,node,iobject).ge.0.d0) then
                  g0(iobject)=g0(iobject)+1.d0
               endif
            endif
          else
            node=ialset(j-2)
            do
              node=node-ialset(j)
              if(node.ge.ialset(j-1)) exit
              if(nodedesiinv(node).eq.1) then
!           
!               calculate actual movement of design variable
!
                dispvector(1)=co(1,node)-coini(1,node)
                dispvector(2)=co(2,node)-coini(2,node)
                dispvector(3)=co(3,node)-coini(3,node)
                actmove=dsqrt((dispvector(1)*extnorini(1,node))**2+
     &                        (dispvector(2)*extnorini(2,node))**2+
     &                        (dispvector(3)*extnorini(3,node))**2)
                if(actmove.lt.1.0e-8) then
                   actmove=0.d0
                   dispvector(1)=0.d0
                   dispvector(2)=0.d0
                   dispvector(3)=0.d0
                endif
                scalprod=dispvector(1)*xdesi(1,nodedesipos(node))+
     &                   dispvector(2)*xdesi(2,nodedesipos(node))+
     &                   dispvector(3)*xdesi(3,nodedesipos(node))
                if(scalprod.le.0.d0) then
                   dgdxglob(1,node,iobject)=-actmove
                else
                   dgdxglob(1,node,iobject)=actmove
                endif
                if((objectset(1,iobject)(1:9).eq.'MAXGROWTH').and.
     &             (scalprod.ge.0d0)) then
                   dgdxglob(2,node,iobject)=actmove-bound
                elseif((objectset(1,iobject)(1:9).eq.'MAXGROWTH').and.
     &             (scalprod.lt.0d0)) then
                   dgdxglob(2,node,iobject)=-actmove-bound
                elseif((objectset(1,iobject)(1:12).eq.'MAXSHRINKAGE')
     &             .and.(scalprod.le.0d0)) then
                   dgdxglob(2,node,iobject)=actmove-bound
                elseif((objectset(1,iobject)(1:12).eq.'MAXSHRINKAGE')
     &             .and.(scalprod.gt.0d0)) then
                   dgdxglob(2,node,iobject)=-actmove-bound
                endif
! 
!               count number of active nodes
!
                if(dgdxglob(2,node,iobject).ge.0.d0) then
                   g0(iobject)=g0(iobject)+1.d0
                endif
              endif
            enddo
          endif
        enddo
      endif
!     
      return        
      end
