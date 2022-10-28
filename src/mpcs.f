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
      subroutine mpcs(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nset_,nalset,nalset_,ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,lakon,ipkon,kon,nk,nk_,
     &  nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,iperturb,ne_,
     &  co,xboun,ctrl,typeboun,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *MPC
!
      implicit none
!
      character*1 typeboun(*),inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*),label
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),ipompc(*),
     &  nodempc(3,*),idirref,ier,
     &  nset,nset_,nalset,nalset_,nmpc,nmpc_,mpcfree,nk,nk_,ikmpc(*),
     &  ilmpc(*),ipkon(*),kon(*),i,node,ipos,istep,istat,n,ne_,
     &  j,k,nodeboun(*),ndirboun(*),ikboun(*),ilboun(*),ipoinpc(0:*),
     &  nboun,nboun_,key,iperturb(*),istart,inode,m,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*)
!
      real*8 coefmpc(3,*),co(3,*),xboun(*),ctrl(*)
!
      idirref=0
!
      if(istep.gt.0) then
         write(*,*) 
     &     '*ERROR reading *MPC: *MPC should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      if(iperturb(1).eq.1) then
         write(*,*) '*ERROR reading *MPC: the *MPC option'
         write(*,*) '       cannot be used in a perturbation step'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *MPC: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*MPC%")
      enddo
!
      istart=0
      inode=0
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!
         if(istart.eq.0) then
            label=textpart(1)(1:20)
            istart=2
         else
            istart=1
         endif
!
         do i=istart,n
            read(textpart(i)(1:10),'(i10)',iostat=istat) node
            if(istat.gt.0) then
               noset=textpart(i)(1:80)
               noset(81:81)=' '
               ipos=index(noset,' ')
               noset(ipos:ipos)='N'
c               do j=1,nset
c                  if(noset.eq.set(j)) then
               call cident81(set,noset,nset,j)
               if(j.gt.0) then
                 if(noset.eq.set(j)) then
                   m=iendset(j)-istartset(j)+1
                   do k=1,m
                     node=ialset(istartset(j)+k-1)
                     inode=inode+1
                     if(label(1:8).eq.'STRAIGHT') then
                       call straightmpc(ipompc,nodempc,coefmpc,
     &                      labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                      nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                      nboun,nboun_,xboun,inode,node,co,
     &                      typeboun)
                     elseif(label(1:5).eq.'PLANE') then
                       call planempc(ipompc,nodempc,coefmpc,
     &                      labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                      nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                      nboun,nboun_,xboun,inode,node,co,
     &                      typeboun)
                     elseif(label(1:4).eq.'BEAM') then
                       call beammpc(ipompc,nodempc,
     &                      labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                      nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                      nboun,nboun_,inode,node,co,
     &                      typeboun)
                     else
                       call usermpc(ipompc,nodempc,coefmpc,
     &                      labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                      nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                      nboun,nboun_,inode,node,co,label,
     &                      typeboun,iperturb,node,idirref,xboun)
                     endif
                   enddo
c     exit
                 endif
               endif
c     endif
c     enddo
               if(j.gt.nset) then
                  noset(ipos:ipos)=' '
                  write(*,*) '*ERROR in nosets: node set ',
     &                 noset
                  write(*,*) '       has not been defined yet'
                  ier=1
                  return
               endif
            else
               inode=inode+1
               if(node.eq.0) then
                  inode=inode-1
                  cycle
               endif
               if(label(1:8).eq.'STRAIGHT') then
                  call straightmpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,xboun,inode,node,co,
     &                 typeboun)
               elseif(label(1:5).eq.'PLANE') then
                  call planempc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,xboun,inode,node,co,typeboun)
               elseif(label(1:4).eq.'BEAM') then
                  call beammpc(ipompc,nodempc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,inode,node,co,typeboun)
               else
                  call usermpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,inode,node,co,label,
     &                 typeboun,iperturb,node,idirref,xboun)
               endif
            endif
         enddo
!
      enddo
!
!     nonhomogeneous term for user MPC
!
      if((label(1:8).ne.'STRAIGHT').and.(label(1:5).ne.'PLANE').and.
     &   (label(1:4).ne.'BEAM'))
     &     then
         node=0
         call usermpc(ipompc,nodempc,coefmpc,
     &        labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &        nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &        nboun,nboun_,inode,node,co,label,typeboun,
     &        iperturb,node,idirref,xboun)
      else
!
!     the *MPC option implies a nonlinear geometric 
!     calculation for all MPC's except MEANROT MPC's
!
         iperturb(2)=1
         write(*,*) '*INFO reading *MPC: nonlinear geometric'
         write(*,*) '      effects are turned on'
         write(*,*)
         if(iperturb(1).eq.0) iperturb(1)=2
      endif
!
      return
      end
