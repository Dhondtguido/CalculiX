!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine mafillsmcsmatrix(ipompc,nodempc,coefmpc,nmpc,
     &     labmpc,ad,au,nactdof,jq,irow,mi,ner,
     &     k,m,node1,node2,jj,ll,val,i,mcs,ielcs,ics,cs)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     for a substructure (superelement) for cyclic symmetry calculations
!     
      implicit none
!     
      character*20 labmpc(*)
!     
      integer ipompc(*),nodempc(3,*),jq(*),mi(*),nactdof(0:mi(2),*),
     &     irow(*),nmpc,i,k,m,jj,icomplex,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,i0,icomplex1,icomplex2,
     &     ij,ilength,inode,inode1,inode2,lprev,mcs,ner,ielcs(*),ics(*)
!     
      real*8 coefmpc(*),val,ad(*),au(*),value,walue,tr,ti,cs(17,*)
!     
      if (mcs.gt.1)then
        if(ielcs(i).gt.0) then
          val=(cs(1,(ielcs(i)+1))/cs(1,1))*val
        endif
      endif
!     
      jdof1=nactdof(k,node1)
      jdof2=nactdof(m,node2)
!     
!     check whether one of the DOF belongs to a SPC or MPC
!     
      if((jdof1.gt.0).and.(jdof2.gt.0)) then
        call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
     &       val,jj,ll)
        call add_sm_st(au,ad,jq,irow,jdof1+ner,jdof2+ner,
     &       val,jj,ll)
      elseif((jdof1.gt.0).or.(jdof2.gt.0)) then
!     
!     idof1: genuine DOF
!     idof2: nominal DOF of the SPC/MPC
!     
        if(jdof1.le.0) then
          idof1=jdof2
          idof2=jdof1
        else
          idof1=jdof1
          idof2=jdof2
        endif
!     
        if(nmpc.gt.0) then
          if(idof2.ne.2*(idof2/2)) then
!     
!     regular DOF / MPC
!     
            id1=(-idof2+1)/2
            ist=ipompc(id1)
            index=nodempc(3,ist)
            if(index.eq.0) return
            do
              inode=nodempc(1,index)
              icomplex=0
              if(labmpc(id1)(1:6).eq.'CYCLIC') then
                read(labmpc(id1)(7:20),'(i14)') icomplex
              elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                do ij=1,mcs
                  ilength=int(cs(4,ij))
                  lprev=int(cs(14,ij))
                  call nident(ics(lprev+1),inode,ilength,id)
                  if(id.gt.0) then
                    if(ics(lprev+id).eq.inode) then
                      icomplex=ij
                      exit
                    endif
                  endif
                enddo
              endif
              idof2=nactdof(nodempc(2,index),inode)
              if(idof2.gt.0) then
                value=-coefmpc(index)*val/coefmpc(ist)
                if(idof1.eq.idof2) then
                  value=2.d0*value
                endif
                if(icomplex.eq.0) then
                  call add_sm_st(au,ad,jq,irow,
     &                 idof1,idof2,value,i0,i0)
                  call add_sm_st(au,ad,jq,irow,
     &                 idof1+ner,idof2+ner,value,i0,i0)
                else
                  walue=value*cs(15,icomplex)
                  call add_sm_st(au,ad,jq,irow,
     &                 idof1,idof2,walue,i0,i0)
                  call add_sm_st(au,ad,jq,irow,
     &                 idof1+ner,idof2+ner,walue,i0,i0)
                  if(idof1.ne.idof2) then
                    walue=value*cs(16,icomplex)
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1,idof2+ner,walue,i0,i0)
                    walue=-walue
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1+ner,idof2,walue,i0,i0)
                  endif
                endif
              endif
              index=nodempc(3,index)
              if(index.eq.0) exit
            enddo
            return
          endif
        endif
!     
      else
        idof1=jdof1
        idof2=jdof2
!     
        mpc1=0
        mpc2=0
        if(nmpc.gt.0) then
          if(idof1.ne.2*(idof1/2)) mpc1=1
          if(idof2.ne.2*(idof2/2)) mpc2=1
        endif
        if((mpc1.eq.1).and.(mpc2.eq.1)) then
          id1=(-idof1+1)/2
          id2=(-idof2+1)/2
          if(id1.eq.id2) then
!     
!     MPC id1 / MPC id1
!     
            ist=ipompc(id1)
            index1=nodempc(3,ist)
            if(index1.eq.0) return
            do
              inode1=nodempc(1,index1)
              icomplex1=0
              if(labmpc(id1)(1:6).eq.'CYCLIC') then
                read(labmpc(id1)(7:20),'(i14)') icomplex1
              elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                do ij=1,mcs
                  ilength=int(cs(4,ij))
                  lprev=int(cs(14,ij))
                  call nident(ics(lprev+1),inode1,
     &                 ilength,id)
                  if(id.gt.0) then
                    if(ics(lprev+id).eq.inode1) then
                      icomplex1=ij
                      exit
                    endif
                  endif
                enddo
              endif
              idof1=nactdof(nodempc(2,index1),inode1)
              index2=index1
              do
                inode2=nodempc(1,index2)
                icomplex2=0
                if(labmpc(id1)(1:6).eq.'CYCLIC') then
                  read(labmpc(id1)(7:20),'(i14)') icomplex2
                elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                  do ij=1,mcs
                    ilength=int(cs(4,ij))
                    lprev=int(cs(14,ij))
                    call nident(ics(lprev+1),inode2,
     &                   ilength,id)
                    if(id.gt.0) then
                      if(ics(lprev+id).eq.inode2) then
                        icomplex2=ij
                        exit
                      endif
                    endif
                  enddo
                endif
                idof2=nactdof(nodempc(2,index2),inode2)
                if((idof1.gt.0).and.(idof2.gt.0)) then
                  value=coefmpc(index1)*coefmpc(index2)*
     &                 val/coefmpc(ist)/coefmpc(ist)
                  if((icomplex1.eq.0).and.
     &                 (icomplex2.eq.0)) then
                    call add_sm_st(au,ad,jq,
     &                   irow,idof1,idof2,value,i0,i0)
                    call add_sm_st(au,ad,jq,
     &                   irow,idof1+ner,idof2+ner,value,
     &                   i0,i0)
                  elseif((icomplex1.ne.0).and.
     &                   (icomplex2.ne.0)) then
                    if(icomplex1.eq.icomplex2) then
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1,idof2,value,i0,i0)
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1+ner,idof2+ner,value,
     &                     i0,i0)
                    else
                      tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                     +cs(16,icomplex1)*cs(16,icomplex2)
                      walue=value*tr
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1,idof2,walue,i0,i0)
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1+ner,idof2+ner,walue,
     &                     i0,i0)
                      ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                     -cs(15,icomplex2)*cs(16,icomplex1)
                      walue=value*ti
                      call add_sm_st(au,ad,jq,irow
     &                     ,idof1,idof2+ner,walue,i0,i0)
                      walue=-walue
                      call add_sm_st(au,ad,jq,irow
     &                     ,idof1+ner,idof2,walue,i0,i0)
                    endif
                  elseif((icomplex1.eq.0).or.
     &                   (icomplex2.eq.0)) then
                    if(icomplex2.ne.0) then
                      walue=value*cs(15,icomplex2)
                    else
                      walue=value*cs(15,icomplex1)
                    endif
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1,idof2,walue,i0,i0)
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1+ner,idof2+ner,walue,i0,i0)
                    if(icomplex2.ne.0) then
                      walue=value*cs(16,icomplex2)
                    else
                      walue=-value*cs(16,icomplex1)
                    endif
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1,idof2+ner,walue,i0,i0)
                    walue=-walue
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1+ner,idof2,walue,i0,i0)
                  endif
                endif
                index2=nodempc(3,index2)
                if(index2.eq.0) exit
              enddo
              index1=nodempc(3,index1)
              if(index1.eq.0) exit
            enddo
          else
!     
!     MPC id1 / MPC id2
!     
            ist1=ipompc(id1)
            index1=nodempc(3,ist1)
            if(index1.eq.0) return
            do
              inode1=nodempc(1,index1)
              icomplex1=0
              if(labmpc(id1)(1:6).eq.'CYCLIC') then
                read(labmpc(id1)(7:20),'(i14)') icomplex1
              elseif(labmpc(id1)(1:9).eq.'SUBCYCLIC') then
                do ij=1,mcs
                  ilength=int(cs(4,ij))
                  lprev=int(cs(14,ij))
                  call nident(ics(lprev+1),inode1,
     &                 ilength,id)
                  if(id.gt.0) then
                    if(ics(lprev+id).eq.inode1) then
                      icomplex1=ij
                      exit
                    endif
                  endif
                enddo
              endif
              idof1=nactdof(nodempc(2,index1),inode1)
              ist2=ipompc(id2)
              index2=nodempc(3,ist2)
              if(index2.eq.0) then
                index1=nodempc(3,index1)
                if(index1.eq.0) then
                  exit
                else
                  return
                endif
              endif
              do
                inode2=nodempc(1,index2)
                icomplex2=0
                if(labmpc(id2)(1:6).eq.'CYCLIC') then
                  read(labmpc(id2)(7:20),'(i14)') icomplex2
                elseif(labmpc(id2)(1:9).eq.'SUBCYCLIC') then
                  do ij=1,mcs
                    ilength=int(cs(4,ij))
                    lprev=int(cs(14,ij))
                    call nident(ics(lprev+1),inode2,
     &                   ilength,id)
                    if(id.gt.0) then
                      if(ics(lprev+id).eq.inode2) then
                        icomplex2=ij
                        exit
                      endif
                    endif
                  enddo
                endif
                idof2=nactdof(nodempc(2,index2),inode2)
                if((idof1.gt.0).and.(idof2.gt.0)) then
                  value=coefmpc(index1)*coefmpc(index2)*
     &                 val/coefmpc(ist1)/coefmpc(ist2)
                  if(idof1.eq.idof2) then
                    value=2.d0*value
                  endif
                  if((icomplex1.eq.0).and.
     &                 (icomplex2.eq.0)) then
                    call add_sm_st(au,ad,jq,
     &                   irow,idof1,idof2,value,i0,i0)
                    call add_sm_st(au,ad,jq,
     &                   irow,idof1+ner,idof2+ner,value,
     &                   i0,i0)
                  elseif((icomplex1.ne.0).and.
     &                   (icomplex2.ne.0)) then
                    if(icomplex1.eq.icomplex2) then
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1,idof2,value,i0,i0)
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1+ner,idof2+ner,value,
     &                     i0,i0)
                    else
                      tr=cs(15,icomplex1)*cs(15,icomplex2)
     &                     +cs(16,icomplex1)*cs(16,icomplex2)
c     write(*,*) 'tr= ',tr
                      walue=value*tr
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1,idof2,walue,i0,i0)
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1+ner,idof2+ner,walue,
     &                     i0,i0)
                      ti=cs(15,icomplex1)*cs(16,icomplex2)
     &                     -cs(15,icomplex2)*cs(16,icomplex1)
                      walue=value*ti
                      call add_sm_st(au,ad,jq,irow
     &                     ,idof1,idof2+ner,walue,i0,i0)
                      walue=-walue
                      call add_sm_st(au,ad,jq,irow
     &                     ,idof1+ner,idof2,walue,i0,i0)
                    endif
                  elseif((icomplex1.eq.0).or.
     &                   (icomplex2.eq.0)) then
                    if(icomplex2.ne.0) then
                      walue=value*cs(15,icomplex2)
                    else
                      walue=value*cs(15,icomplex1)
                    endif
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1,idof2,walue,i0,i0)
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1+ner,idof2+ner,walue,i0,i0)
                    if(idof1.ne.idof2) then
                      if(icomplex2.ne.0) then
                        walue=value*cs(16,icomplex2)
                      else
                        walue=-value*cs(16,icomplex1)
                      endif
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1,idof2+ner,walue,
     &                     i0,i0)
                      walue=-walue
                      call add_sm_st(au,ad,jq,
     &                     irow,idof1+ner,idof2,walue,
     &                     i0,i0)
                    endif
                  endif
                endif
                index2=nodempc(3,index2)
                if(index2.eq.0) exit
              enddo
              index1=nodempc(3,index1)
              if(index1.eq.0) exit
            enddo
          endif
        endif
      endif
!     
      return
      end
