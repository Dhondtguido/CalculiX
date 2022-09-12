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
      subroutine modifympc(inodestet,nnodestet,co,doubleglob,
     &     integerglob,ipompc,nodempc,coefmpc,nmpc,nmpc_,labmpc,mpcfree,
     &     ikmpc,ilmpc,jq,irow,icol,loc,irowt,jqt,itemp,au,ixcol,
     &     ikboun,nboun,nodeboun,mpcrfna,mpcrfnb,nodempcref,
     &     coefmpcref,memmpcref_,mpcfreeref,maxlenmpcref,memmpc_,
     &     maxlenmpc,istep)
!     
!     generating MPC's connecting the nodes in the old tet mesh in
!     which SPC's or point forces were defined with the nodes in the
!     refined tet mesh
!     
      implicit none
!
      logical spc
!     
      character*20 labmpc(*)
!     
      integer inodestet(*),nnodestet,integerglob(*),nktet,netet,ne,nkon,
     &     nfaces,nfield,nselect,imastset,iselect(1),nterms,idir1,
     &     nelem,ialset(1),mpcfreeold,jq(*),irow(*),icol(*),
     &     iendset(1),istartset(1),konl(20),loopa,loc(*),irowt(*),
     &     node,i,j,k,m,ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &     ikmpc(*),ilmpc(*),idof,id,jqt(*),nzs,nosol,nodeboun(*),
     &     kflag,itemp(*),maxrow,ixcol(*),isoltot,newnode,iy,ic,idummy,
     &     nlength,kopt,ixcolnew,idnew,i1,ikboun(*),nboun,mpcrfna,
     &     mpcrfnb,index,index1,node1,mpc,idof1,idof2,id1,id2,
     &     nodempcref(3,*),memmpcref_,mpcfreeref,maxlenmpcref,memmpc_,
     &     maxlenmpc,istep
!     
      real*8 co(3,*),doubleglob(*),coords(3),ratio(20),value,au(*),
     &     coefmpc(*),depcoef,dist,xopt,coef1,coefmpcref(*)
!
      if(istep.gt.1) then
!     
!     restoring the original equations
!     
        memmpc_=memmpcref_
        mpcfree=mpcfreeref
        maxlenmpc=maxlenmpcref
!     
!     mpcfreeref=-1 indicates that the MPC's have changed and that
!     a new copy is to be taken before the next call to cascade.c            
!     
        mpcfreeref=-1
!     
        do i=1,memmpc_
          do j=1,3
            nodempc(j,i)=nodempcref(j,i)
          enddo
          coefmpc(i)=coefmpcref(i)
        enddo
      endif
!
!     restore the original form of the MPC's
!     an MPC with number in between mpcrfna and mpcrfnb (i.e. MPC's
!     corresponding to the connection of the old mesh with the new
!     mesh) is in its orignal form if the coefficient of the first
!     term is exactly 1.d0
!
      do i=mpcrfna,mpcrfnb
        index=ipompc(i)
        if(coefmpc(index).eq.1.d0) cycle
!
!     MPC has been reordered and has to be brought in its original
!     form
!
        node1=nodempc(1,index)
        idir1=nodempc(2,index)
        coef1=coefmpc(index)
        idof1=8*(node1-1)+idir1
        index1=index
        call nident(ikmpc,idof1,nmpc,id1)
!
        index=nodempc(3,index)
!
        do
          if(index.eq.0) then
            write(*,*) '*ERROR in modifympc: equation ',i
            write(*,*) '       cannot be modified'
            call exit(201)
          endif
          if(coefmpc(index).eq.1.d0) then
!
!     restore the leading term of the equation
!
            nodempc(1,index1)=nodempc(1,index)
            coefmpc(index1)=coefmpc(index)
!
!     restore the independent term of the equation
!     id1>id2 since the refined mesh has higher node numbers
!     than the unrefined mesh
!
            idof2=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikmpc,idof2,nmpc,id2)
            do j=id1,id2+2,-1
              ikmpc(j)=ikmpc(j-1)
              ilmpc(j)=ilmpc(j-1)
            enddo
            ikmpc(id2+1)=idof2
            ilmpc(id2+1)=i
!
            nodempc(1,index)=node1
            coefmpc(index)=coef1
            exit
          else
            index=nodempc(3,index)
          endif
        enddo
      enddo
!        
!     compile all dependent nodes in jq, irow..
!     taking the equations for the first dof (x-direction) only      
!
      nzs=0
      nnodestet=0
      maxrow=0
!
      do i=mpcrfna,mpcrfnb
        index=ipompc(i)
!
!       x-direction
!
        if(nodempc(2,index).ne.1) cycle
!
        nnodestet=nnodestet+1
        inodestet(nnodestet)=nodempc(1,index)
        jq(nnodestet)=nzs+1
        index=nodempc(3,index)
        do
          if(index.eq.0) exit
          nzs=nzs+1
          irow(nzs)=nodempc(1,index)
          maxrow=max(maxrow,irow(nzs))
          au(nzs)=coefmpc(index)
          irowt(nzs)=nnodestet
          loc(nzs)=nzs
          index=nodempc(3,index)
        enddo
      enddo
      jq(nnodestet+1)=nzs+1
!
!     determine the nodes for which the dependent term in the MPC's
!     has to be modified: these are MPC's for which the dependent term
!     is subject to other SPC's or MPC's
!
!     SPC's
!
      do i=1,nboun
        node=nodeboun(i)
        call nident(inodestet,node,nnodestet,id)
        if(id.gt.0) then
          if(inodestet(id).eq.node) then
            icol(id)=1
          endif
        endif
      enddo
!
!     MPC's
!
      do i=1,mpcrfna-1
        index=ipompc(i)
        node=nodempc(1,index)
        call nident(inodestet,node,nnodestet,id)
        if(id.gt.0) then
          if(inodestet(id).eq.node) then
            icol(id)=1
          endif
        endif
      enddo
!
!     for the nodes i for which the equations have to be modified
!     icol(i) is nonzero
!
      isoltot=0
      do i=1,nnodestet
        if(icol(i).eq.1) then
          icol(i)=jq(i+1)-jq(i)
        else
          isoltot=isoltot+1
        endif
      enddo
!     
!     copying irow into itemp, sorting itemp and field loc and
!     irowt along
!     
      do j=1,nzs
        itemp(j)=irow(j)
      enddo
!     
      kflag=2
      call isortiii(itemp,loc,irowt,nzs,kflag)
!     
!     determining jqt (columns per row); jqt and irowt is the
!     equivalent (for row per row treatment) of jq and irow (for
!     column per column treatment); one can consider jqt and irowt
!     to be the jq and irow fields for the transpose of the matrix
!     (therefore the appended letter "t")
!     
      j=1
      do m=1,nzs
        if(itemp(m).ge.j) then
          do k=j,itemp(m)
            jqt(k)=m
          enddo
          j=itemp(m)+1
        endif
      enddo
      jqt(j)=nzs+1
!     
!     creating a unique number consisting of the column number and the
!     number of rows per column; if sorted, the result is the same as
!     for sorting icol.
!     
      do j=1,nnodestet
        if(icol(j).gt.0) then
          ixcol(j)=(nnodestet+1)*icol(j)+j
        endif
      enddo
!     
!     sorting ixcol 
!     
      kflag=1
      call isortii(ixcol,idummy,nnodestet,kflag)
!     
!     reordering the equations: the first term is to be a
!     new node instead of an old node; starting with columns
!     with only 1 nonzero value
!     
      nosol=0
      do
        if(isoltot.eq.nnodestet) exit
        isoltot=isoltot+1
!
!       local dependent node number
!
        i=ixcol(isoltot)
     &       -int(ixcol(isoltot)/(nnodestet+1))*(nnodestet+1)
        ixcol(isoltot)=0
        icol(i)=0
!     
!       global dependent node number
!     
        node=inodestet(i)
!
!       determining the optimal row for the exchange
!
        kopt=0
        xopt=0.d0
        do k=jq(i),jq(i+1)-1
          if(irow(k).gt.0) then
            if(dabs(au(k)).gt.xopt) then
              kopt=k
              xopt=dabs(au(k))
            endif
          endif
        enddo
        k=kopt
!     
        if(k.gt.0) then
!     
!     row dof
!     
          newnode=irow(k)
          depcoef=au(k)
!
!         generating MPC's
!
          do j=1,3
!
!     modify ikmpc to take the switch from id2 to id1 as
!     dependent node into account
!
            idof1=8*(node-1)+j
            call nident(ikmpc,idof1,nmpc,id1)
            do
              if(labmpc(ilmpc(id1)).eq.'RM                  ') exit
              id=id-1
            enddo
            mpc=ilmpc(id1)
            index=ipompc(mpc)
            idof2=8*(newnode-1)+j
            call nident(ikmpc,idof2,nmpc,id2)
            do m=id1,id2-1
              ikmpc(m)=ikmpc(m+1)
              ilmpc(m)=ilmpc(m+1)
            enddo
            ikmpc(id2)=idof2
            ilmpc(id2)=mpc
!
!     inserting the new dependent term
!
            nodempc(1,index)=newnode
            coefmpc(index)=depcoef
!
!     inserting the new independent term
!
            do
              index=nodempc(3,index)
              if(index.eq.0) then
                write(*,*) '*ERROR in modifympc: equation ',mpc
                write(*,*) '       cannot be modified'
                call exit(201)
              endif
              if(nodempc(1,index).eq.newnode) then
                nodempc(1,index)=node
                coefmpc(index)=1.d0
                exit
              endif
            enddo
          enddo
!
!     remove the indpendent node from the data base:          
!     tagging the row dof in all other columns
!     
          do m=jqt(newnode),jqt(newnode+1)-1
            if(irowt(m).eq.i) cycle
            ic=irowt(m)
!
!           no SPC/MPC or already treated
!
            if(icol(ic).eq.0) cycle
!
!     setting negative sign in irow to show that the row
!     is not available any more (already used as dependent dof)
!
            irow(loc(m))=-irow(loc(m))
            iy=(nnodestet+1)*icol(ic)+ic
            nlength=nnodestet-isoltot
            call nident(ixcol(isoltot+1),iy,nlength,id)
            if(ixcol(isoltot+id).ne.iy) then
              write(*,*) '*ERROR in modifympc:'
              write(*,*) '       data base corrupt'
              call exit(201)
            endif
            icol(ic)=icol(ic)-1
            ixcolnew=ixcol(isoltot+id)-(nnodestet+1)
!
!           not treated equation has no independent terms any more
!
            if(ixcolnew.lt.(nnodestet+1)) then
              nosol=nosol+1
              if(nosol.eq.1) then
                write(*,*) '*WARNING in modifympc; failed to connect'
                write(*,*) '         node ',inodestet(ic),
     &               ' of unrefined mesh'
                write(*,*)
     &'         to the refined mesh; loads and boundary conditions'
                write(*,*)
     &'         in this node are not taken into account; a list of'
                write(*,*)
     &'         not connected nodes is stored in'
                write(*,*)
     &'         WarnNodeMissRefineConnection.nam'
              else
                write(*,*) '*WARNING in modifympc; failed to connect'
                write(*,*) '         node ',inodestet(ic),
     &               ' of unrefined mesh'
              endif
              write(23,*)
     &             '*NSET,NSET=WarnNodeMissRefineConnection'
              write(23,*) inodestet(ic)
            endif
            call nident(ixcol(isoltot+1),ixcolnew,id-1,idnew)
            do i1=isoltot+id,isoltot+idnew+2,-1
              ixcol(i1)=ixcol(i1-1)
            enddo
            ixcol(isoltot+idnew+1)=ixcolnew
          enddo
        endif
      enddo
      write(*,*)
      write(*,*) '*INFO in modifympc: no solution for ',nosol,' nodes.'
      write(*,*)
      close(23)
!     
      return
      end
