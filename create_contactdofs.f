!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
!     
!     MASSLESS DYNAMIC CONTACT
!     
C     >    *MASSLESS DYNAMIC CONTACT*: create the fields kslav,lslav catalogueing the matrix degrees
C     >    of freedom corresponding to the slave degrees of freedom and
C     >    ktot,ltot catalogueing the matrix degrees of freedom corresponding
C     >    to the slave+master degrees of freedom
C     >    @param kslav     1D array of slave nodes indices SORTED of size *neqslav*
C     >    @param lslav     1D array of slave nodes labels  SORTED of size *neqslav* i.e. 261 = node 26 DOF 1
C     >    @param ktot      1D array of slave+master  nodes indices SORTED of size *neqtot*
C     >    @param ltot      1D array of slave+master  nodes labels  SORTED of size *neqtot* i.e. 261 = node 26 DOF 1
C     >    @param nslavs    total of slave nodes
C     >    @param islavnode see docu of massless.c 
C     >    @param nmasts    total of master nodes
C     >    @param imastnode see docu of massless.c 
C     >    @param nactdof   see docu of massless.c
C     >    @param mi        see docu of massless.c
C     >    @param neqslav   total of slave equations (DOF)
C     >    @param neqtot    total of master+slave equations (DOF)
C     >    @see massless.c
C     >    @author Guido Dhondt
C     
      subroutine create_contactdofs(kslav,lslav,ktot,ltot,nslavs,
     &     islavnode,nmasts,imastnode,nactdof,mi,neqtot,
     &     nslavnode,fric,tieset,tietol,ntie,elcon,ncmat_,ntmat_)
!     
      implicit none
!
      character*81 tieset(3,*)
!     
      integer kslav(*),lslav(*),ktot(*),ltot(*),i,j,k,nslavs,node,
     &     islavnode(*),nmasts,imastnode(*),mi(*),nactdof(0:mi(2),*),
     &     idof,neqslav,neqtot,kflag,nslavnode(*),ntie,imat,ncmat_,
     &     ntmat_
!
      real*8 fric(*),tietol(4,*),elcon(0:ncmat_,ntmat_,*)
!
!     assigning the appropriate friction coefficient to each slave node
!
      do i=1,ntie
        if(tieset(1,i)(81:81).ne.'C') cycle
        imat=int(tietol(2,i))
        if(ncmat_.ge.7) then
          do j=nslavnode(i)+1,nslavnode(i+1)
            fric(j)=elcon(6,1,imat)
          enddo
        else
          do j=nslavnode(i)+1,nslavnode(i+1)
            fric(j)=0.d0
          enddo
        endif
      enddo
!     
!     collecting the slave degrees of freedom in kslav
!     lslav contains the corresponding node and direction in the
!     for 10*node+dir
!     
      k=0
!     
      do i=1,nslavs
        node=islavnode(i)
        do j=1,3
          idof=nactdof(j,node)
          if(idof.le.0) then
            write(*,*) '*ERROR in create_contactdofs'
            write(*,*) '       a SPC and/or MPC was defined'
            write(*,*) '       in slave node: ',node,
     &           ' belonging to a massless'
            write(*,*) '       contact pair definition; this is'
            write(*,*) '       not allowed'
            call exit(201)
          endif
          k=k+1
          kslav(k)=idof
          lslav(k)=10*node+j
          ktot(k)=kslav(k)
          ltot(k)=lslav(k)
        enddo
      enddo
      neqslav=k
!     
      if(neqslav.ne.3*nslavs) then
        write(*,*) '*ERROR in create_contactdofs'
        write(*,*) '       neqslav.ne.3*nslavs'
        call exit(201)
      endif
!     
!     sorting kslav in ascending order and lslav along
!     
      kflag=2
      call isortii(kslav,lslav,neqslav,kflag)
!     
!     collecting the slave+master degrees of freedom in ktot
!     ltot contains the corresponding node and direction in the
!     for 10*node+dir
!     
      do i=1,nmasts
        node=imastnode(i)
        do j=1,3
          idof=nactdof(j,node)
          if(idof.le.0) then
            write(*,*) '*ERROR in create_contactdofs'
            write(*,*) '       a SPC and/or MPC was defined'
            write(*,*) '       in master node: ',node,
     &           ' belonging to a massless'
            write(*,*) '       contact pair definition; this is'
            write(*,*) '       not allowed'
            call exit(201)
          endif
          k=k+1
          ktot(k)=idof
          ltot(k)=10*node+j
        enddo
      enddo
      neqtot=k
!     
      if(neqtot.ne.3*(nslavs+nmasts)) then
        write(*,*) '*ERROR in create_contactdofs'
        write(*,*) '       neqtot.ne.3*(nslavs+nmasts)'
        call exit(201)
      endif
!     
!     sorting ktot in ascending order and ltot along
!     
      call isortii(ktot,ltot,neqtot,kflag)
!     
      return
      end
