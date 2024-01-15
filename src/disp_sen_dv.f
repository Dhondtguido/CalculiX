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
      subroutine disp_sen_dv(nodeset,istartset,iendset,ialset,iobject,
     &     mi,nactdof,dgdu,vold,objectset,nactdofinv,neq,g0,nod1st,
     &     ne2d)
!     
!     calculates the sum of the square of the displacements of a node
!     set and its derivative w.r.t. the coordinates of the mesh
!     
      implicit none
!     
      character*81 objectset(5,*)
!     
      integer istartset(*),iendset(*),ialset(*),nodeset,idir,
     &     idof,iobject,mi(*),nactdof(0:mi(2),*),j,k,nactdofinv(*),
     &     inode,node,neq,mt,nod1st(*),ne2d
!     
      real*8 dgdu(*),vold(0:mi(2),*),g0(*)
!     
      mt=mi(2)+1
!     
!     check for the existence of a set, else take the complete mesh
!     
      if(nodeset.eq.0) then
        do idof=1,neq
          inode=nactdofinv(idof)               
          idir=inode-mt*(inode/mt);
          node=inode/mt+1;
          if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
            dgdu(idof)=vold(idir,node)/g0(iobject)
          elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
            if(idir.eq.1) then
              dgdu(idof)=vold(idir,node)/g0(iobject)
            endif   
          elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
            if(idir.eq.2) then
              dgdu(idof)=vold(idir,node)/g0(iobject)
            endif   
          elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
            if(idir.eq.3) then
              dgdu(idof)=vold(idir,node)/g0(iobject)
            endif   
          endif          
        enddo
      else
        do j=istartset(nodeset),iendset(nodeset)      
          if(ialset(j).gt.0) then
            node=ialset(j)
            if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
              do idir=1,3
                idof=nactdof(idir,node)
                if(idof.gt.0) then
                  dgdu(idof)=vold(idir,node)/g0(iobject)
                elseif(ne2d.ne.0) then
                  if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                    idof=nactdof(idir,node+1)
                    if(idof.gt.0) then
                      dgdu(idof)=vold(idir,node+1)/g0(iobject)
                    endif
                  endif
                endif
              enddo
            elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
              idof=nactdof(1,node)
              if(idof.gt.0) then
                dgdu(idof)=vold(1,node)/g0(iobject)
              elseif(ne2d.ne.0) then
                if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                  idof=nactdof(1,node+1)
                  if(idof.gt.0) then
                    dgdu(idof)=vold(1,node+1)/g0(iobject)
                  endif
                endif   
              endif
            elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
              idof=nactdof(2,node)
              if(idof.gt.0) then
                dgdu(idof)=vold(2,node)/g0(iobject)
              elseif(ne2d.ne.0) then
                if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                  idof=nactdof(2,node+1)
                  if(idof.gt.0) then
                    dgdu(idof)=vold(2,node+1)/g0(iobject)
                  endif
                endif   
              endif
            elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
              idof=nactdof(3,node)
              if(idof.gt.0) then
                dgdu(idof)=vold(3,node)/g0(iobject)
              elseif(ne2d.ne.0) then
                if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                  idof=nactdof(3,node+1)
                  if(idof.gt.0) then
                    dgdu(idof)=vold(3,node+1)/g0(iobject)
                  endif
                endif   
              endif
            endif    
          else
            k=ialset(j-2)
            do
              k=k-ialset(j)
              if(k.ge.ialset(j-1)) exit
              node=k
              if(objectset(1,iobject)(1:8).eq.'ALL-DISP') then
                do idir=1,3
                  idof=nactdof(idir,node)
                  if(idof.gt.0) then
                    dgdu(idof)=vold(idir,node)/g0(iobject)
                  elseif(ne2d.ne.0) then
                    if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                      idof=nactdof(idir,node+1)
                      if(idof.gt.0) then
                        dgdu(idof)=vold(idir,node+1)/g0(iobject)
                      endif
                    endif
                  endif
                enddo
              elseif(objectset(1,iobject)(1:6).eq.'X-DISP') then
                idof=nactdof(1,node)
                if(idof.gt.0) then
                  dgdu(idof)=vold(1,node)/g0(iobject)
                elseif(ne2d.ne.0) then
                  if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                    idof=nactdof(1,node+1)
                    if(idof.gt.0) then
                      dgdu(idof)=vold(1,node+1)/g0(iobject)
                    endif
                  endif
                endif     
              elseif(objectset(1,iobject)(1:6).eq.'Y-DISP') then
                idof=nactdof(2,node)
                if(idof.gt.0) then
                  dgdu(idof)=vold(2,node)/g0(iobject)
                elseif(ne2d.ne.0) then
                  if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                    idof=nactdof(2,node+1)
                    if(idof.gt.0) then
                      dgdu(idof)=vold(2,node+1)/g0(iobject)
                    endif
                  endif
                endif     
              elseif(objectset(1,iobject)(1:6).eq.'Z-DISP') then
                idof=nactdof(3,node)
                if(idof.gt.0) then
                  dgdu(idof)=vold(3,node)/g0(iobject)
                elseif(ne2d.ne.0) then
                  if(nod1st(node).eq.node) then
!     
!     2d-element: check node in the middle
!     
                    idof=nactdof(3,node+1)
                    if(idof.gt.0) then
                      dgdu(idof)=vold(3,node+1)/g0(iobject)
                    endif
                  endif     
                endif
              endif      
            enddo
          endif
        enddo  
      endif
!     
      return
      end
      
