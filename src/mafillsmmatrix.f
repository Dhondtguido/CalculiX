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
      subroutine mafillsmmatrix(ipompc,nodempc,coefmpc,nmpc,
     &     ad,au,nactdof,jq,irow,neq,nmethod,mi,rhsi,
     &     k,m,node1,node2,jj,ll,val,istiff)
!     
!     filling the stiffness matrix in spare matrix format (sm)
!     for a substructure (superelement)
!     
      implicit none
!     
      integer rhsi
!     
      integer ipompc(*),nodempc(3,*),jq(*),mi(*),nactdof(0:mi(2),*),
     &     irow(*),icolumn,nmpc,neq(2),nmethod,k,m,jj,istiff,
     &     ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &     mpc1,mpc2,index1,index2,node1,node2,i0
!     
      real*8 coefmpc(*),val,ad(*),au(*),value
!     
      jdof1=nactdof(k,node1)
      jdof2=nactdof(m,node2)
!     
!     check whether one of the DOF belongs to a SPC or MPC
!     
      if((jdof1.gt.0).and.(jdof2.gt.0)) then
        call add_sm_st(au,ad,jq,irow,jdof1,jdof2,
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
        if(nmpc.gt.0) then
          if(idof2.ne.2*(idof2/2)) then
!     
!     regular DOF / MPC
!     
            id=(-idof2+1)/2
            ist=ipompc(id)
            index=nodempc(3,ist)
            if(index.eq.0) return
            do
              idof2=nactdof(nodempc(2,index),nodempc(1,index))
              value=-coefmpc(index)*val/coefmpc(ist)
              if(idof1.eq.idof2) value=2.d0*value
              if(idof2.gt.0) then
                call add_sm_st(au,ad,jq,irow,idof1,
     &               idof2,value,i0,i0)
              endif
              index=nodempc(3,index)
              if(index.eq.0) exit
            enddo
            return
          endif
        endif
!     
!     regular DOF / SPC
!
        if(istiff.eq.1) then
          if(rhsi.eq.1) then
          elseif(nmethod.eq.2) then
            value=val
            icolumn=neq(2)-idof2/2
            call add_bo_st(au,jq,irow,idof1,icolumn,value)
          endif
        endif
      else
        idof1=jdof1
        idof2=jdof2
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
              idof1=nactdof(nodempc(2,index1),
     &             nodempc(1,index1))
              index2=index1
              do
                idof2=nactdof(nodempc(2,index2),
     &               nodempc(1,index2))
                value=coefmpc(index1)*coefmpc(index2)*
     &               val/coefmpc(ist)/coefmpc(ist)
                if((idof1.gt.0).and.(idof2.gt.0)) then
                  call add_sm_st(au,ad,jq,irow,
     &                 idof1,idof2,value,i0,i0)
                endif
!     
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
              idof1=nactdof(nodempc(2,index1),
     &             nodempc(1,index1))
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
                idof2=nactdof(nodempc(2,index2),
     &               nodempc(1,index2))
                value=coefmpc(index1)*coefmpc(index2)*
     &               val/coefmpc(ist1)/coefmpc(ist2)
                if(idof1.eq.idof2) value=2.d0*value
                if((idof1.gt.0).and.(idof2.gt.0)) then
                    call add_sm_st(au,ad,jq,irow,
     &                   idof1,idof2,value,i0,i0)
                endif
!     
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
