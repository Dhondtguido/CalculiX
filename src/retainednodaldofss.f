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
      subroutine retainednodaldofss(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nodeboun,ndirboun,xboun,nboun,nboun_,nk,
     &  iamboun,nam,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &  mpcfree,inotr,trab,ikboun,ilboun,ikmpc,ilmpc,nk_,
     &  co,labmpc,typeboun,istat,n,iline,ipol,
     &  inl,ipoinp,inp,nmethod,iperturb,
     &  ipoinpc,vold,mi,istep,ier)
!
!     reading the input deck: *RETAINED NODAL DOFS
!
      implicit none
!
      logical fixed
!
      character*1 typeboun(*),type,inpc(*)
      character*20 labmpc(*),label
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nodeboun(*),
     &  ndirboun(*),ntransl,istep,ier,id,
     &  nset,nboun,nboun_,istat,n,i,j,k,l,ibounstart,ibounend,
     &  key,nk,iamboun(*),nam,iamplitude,ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,mpcfree,inotr(2,*),ikboun(*),ilboun(*),ikmpc(*),
     &  ilmpc(*),nk_,ipos,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  nmethod,iperturb(*),ipoinpc(0:*),ktrue,mi(*)
!
      real*8 xboun(*),bounval,coefmpc(*),trab(7,*),co(3,*),
     &  vold(0:mi(2),*)
!
      iamplitude=0
      ntransl=0
      type='C'
      fixed=.false.
      label='                    '
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *RETAINED NODAL DOFS:'
         write(*,*) '       *RETAINED NODAL DOFS can only be used'
         write(*,*) '       within a STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:9).eq.'SORTED=NO') then
         else
            write(*,*) 
     &'*WARNING reading *RETAINED NODAL DOFS: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*RETAINED NODAL DOFS%")
         endif
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*RETAINED NODAL DOFS%",ier)
            return
         endif
!     
         if(textpart(3)(1:1).eq.' ') then
            ibounend=ibounstart
         else
            read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*RETAINED NODAL DOFS%",ier)
               return
            endif
         endif
!     
         bounval=0.d0
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if((l.gt.nk).or.(l.le.0)) then
               write(*,*) '*ERROR reading *RETAINED NODAL DOFS:'
               write(*,*) '       node ',l,' is not defined'
               ier=1
               return
            endif
            ktrue=l
            call bounadd(l,ibounstart,ibounend,bounval,
     &        nodeboun,ndirboun,xboun,nboun,nboun_,
     &        iamboun,iamplitude,nam,ipompc,nodempc,
     &        coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &        ntransl,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &        type,typeboun,nmethod,iperturb,fixed,vold,ktrue,mi,
     &        label)
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
c            do i=1,nset
c               if(set(i).eq.noset) exit
c            enddo
            call cident81(set,noset,nset,id)
            i=nset+1
            if(id.gt.0) then
              if(noset.eq.set(id)) then
                i=id
              endif
            endif
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *RETAINED NODAL DOFS:'
               write(*,*) '       node set ',noset
               write(*,*) '       has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*RETAINED NODAL DOFS%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  k=ialset(j)
                  ktrue=k
                  call bounadd(k,ibounstart,ibounend,bounval,
     &               nodeboun,ndirboun,xboun,nboun,nboun_,
     &               iamboun,iamplitude,nam,ipompc,nodempc,
     &               coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &               ntransl,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &               type,typeboun,nmethod,iperturb,fixed,vold,ktrue,
     &               mi,label)
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     ktrue=k
                     call bounadd(k,ibounstart,ibounend,bounval,
     &                 nodeboun,ndirboun,xboun,nboun,nboun_,
     &                 iamboun,iamplitude,nam,ipompc,nodempc,
     &                 coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &                 ntransl,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,
     &                 labmpc,type,typeboun,nmethod,iperturb,fixed,
     &                 vold,ktrue,mi,label)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

