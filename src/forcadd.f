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
      subroutine forcadd(node,i,val,nodeforc,ndirforc,xforc,
     &     nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &     ikforc,ilforc,isector,add,user,idefforc,ipompc,nodempc,
     &     nmpc,ikmpc,ilmpc,labmpc)
!     
!     adds a cload condition to the data base
!     
      implicit none
!     
      logical add,user
!     
      character*20 labmpc(*)
!     
      integer nodeforc(2,*),ndirforc(*),node,i,nforc,nforc_,j,
     &     iamforc(*),iamplitude,nam,ntrans,inotr(2,*),itr,idf(3),
     &     ikforc(*),ilforc(*),idof,id,k,isector,idefforc(*),ipompc(*),
     &     nodempc(3,*),nmpc,ikmpc(*),ilmpc(*)
!     
      real*8 xforc(*),val,trab(7,*),a(3,3),co(3,*)
!     
      if(ntrans.eq.0) then
        itr=0
      else
        itr=inotr(1,node)
      endif
!     
!     checking for boundary conditions on rotational dofs of
!     distributing couplings 
!     
      if((i.ge.4).and.(i.le.6)) then
!     
!     rotational dof
!     
        idof=8*(node-1)+i
        call nident(ikmpc,idof,nmpc,id)
        if(id.gt.0) then
          if(ikmpc(id).eq.idof) then
            if(labmpc(ilmpc(id))(1:14).eq.'ROTTRACOUPLING') then
              node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
              i=nodempc(2,nodempc(3,ipompc(ilmpc(id))))
              itr=0
            endif
          endif
        endif
      endif
!     
!     a transformation is taken care of at the time of setting up
!     the matrix system (mafillsmforc.f)
!     
      idof=8*(node-1)+i
      call nident(ikforc,idof,nforc,id)
      if(id.gt.0) then
        do
          if(ikforc(id).eq.idof) then
            k=ilforc(id)
            if(nodeforc(2,k).eq.isector) then
              if(add.or.(idefforc(k).eq.1)) then
                if(nam.gt.0) then
                  if(iamforc(k).ne.iamplitude) then
                    write(*,*) '*ERROR in forcadd:'
                    write(*,*) '       it is not allowed to '
                    write(*,*)'       define two concentrated'
                    write(*,*) '       loads/fluxes'
                    write(*,*) '       different amplitudes '
                    write(*,*) '       in one step'
                    write(*,*) 'node: ',node,' dof:',i
                    call exit(201)
                  endif
                endif
                xforc(k)=xforc(k)+val
              else
                xforc(k)=val
                if(.not.user) idefforc(k)=1
              endif
              if(nam.gt.0) iamforc(k)=iamplitude
              return
            endif
            id=id-1
            if(id.eq.0) exit
          else
            exit
          endif
        enddo
      endif
!     
      nforc=nforc+1
      if(nforc.gt.nforc_) then
        write(*,*) '*ERROR in forcadd: increase nforc_'
        call exit(201)
      endif
      nodeforc(1,nforc)=node
      nodeforc(2,nforc)=isector
      ndirforc(nforc)=i
      xforc(nforc)=val
      if(.not.user) idefforc(nforc)=1
      if(nam.gt.0) iamforc(nforc)=iamplitude
!     
!     updating ikforc and ilforc
!     
      do j=nforc,id+2,-1
        ikforc(j)=ikforc(j-1)
        ilforc(j)=ilforc(j-1)
      enddo
      ikforc(id+1)=idof
      ilforc(id+1)=nforc
!     
      return
      end
      
