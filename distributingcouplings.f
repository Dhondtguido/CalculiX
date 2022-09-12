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
      subroutine distributingcouplings(inpc,textpart,ipompc,nodempc,
     &  coefmpc,nmpc,nmpc_,mpcfree,nk,ikmpc,ilmpc,
     &  labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,lakon,
     &  kon,ipkon,set,nset,istartset,iendset,ialset,co,ier)
!
!     reading the input deck: *DISTRIBUTING COUPLING
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*)
      character*81 set(*),elset,noset
      character*132 textpart(16)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &  n,i,j,key,nk,node,ier,
     &  mpcfreeold,ikmpc(*),ilmpc(*),id,idof,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),irefnode,
     &  k,ipos,kon(*),ipkon(*),nset,idir,newmpc,
     &  istartset(*),iendset(*),ialset(*),indexm,ielem
!
      real*8 coefmpc(*),co(3,*),weight,totweight
!
      elset(1:1)=' '
      do i=2,n
         if(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) '*WARNING reading *DISTRIBUTING COUPLING:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*DISTRIBUTING COUPLING%")
         endif
      enddo
!
      if(elset(1:1).eq.' ') then
         write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
         write(*,*) '       no element set given'
         call inputerror(inpc,ipoinpc,iline,
     &        "*DISTRIBUTING COUPLING%",ier)
         return
      endif
!
!     check whether the element set exists
!
c      do i=1,nset
c         if(set(i).eq.elset) exit
c      enddo
      call cident81(set,elset,nset,id)
      i=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          i=id
        endif
      endif
      if(i.gt.nset) then
         write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
         write(*,*) '       element set ',elset(1:ipos-1),
     &         ' is not defined'
         ier=1
         return
      endif
!
!     check whether only one element belongs
!     to the set
!
      if(istartset(i).ne.iendset(i)) then
         write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
         write(*,*) '       element set ',elset(1:ipos-1),
     &        ' contains more than one element'
         ier=1
         return
      endif
!
!     check whether the element is a DCOUP3D element
!
      ielem=ialset(istartset(i))
      if(lakon(ielem)(1:7).ne.'DCOUP3D') then
         write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
         write(*,*) '       element ',ielem,' is not a'
         write(*,*) '       DCOUP3D element'
         ier=1
         return
      endif
!
!     the reference node belongs to the DCOUP3D element
!
      irefnode=kon(ipkon(ielem)+1)
      newmpc=0
      totweight=0.d0
!
!     generate a MPC for dof 1
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!     
         read(textpart(1)(1:10),'(i10)',iostat=istat) node
         if(istat.eq.0) then
            if(node.gt.nk) then
               write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
               write(*,*) '       node ',node,' is not defined'
               ier=1
               return
            endif
!     
!           if first node : new MPC
!     
            if(newmpc.eq.0) then
               idof=8*(node-1)+1
               call nident(ikmpc,idof,nmpc,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
                     write(*,*) '       dof 1 of node ',node,
     &                    ' is already used'
                     ier=1
                     return
                  endif
               endif
!
               nmpc=nmpc+1
               if(nmpc.gt.nmpc_) then
                  write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
                  write(*,*) '       increase nmpc_'
                  ier=1
                  return
               endif
               ipompc(nmpc)=mpcfree
               labmpc(nmpc)='                    '
               ipompc(nmpc)=mpcfree
!     
!              updating ikmpc and ilmpc
!     
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
!
               newmpc=1
            endif
!
!           reading the weight
!
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) weight
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTING COUPLING%",ier)
               return
            endif
            totweight=totweight+weight
!
!           new term in MPC
!            
            nodempc(1,mpcfree)=node
            nodempc(2,mpcfree)=1
            coefmpc(mpcfree)=weight
            mpcfree=nodempc(3,mpcfree)
!
         else
!
!           node set
!
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) weight
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTING COUPLING%",ier)
               return
            endif
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
               write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
               write(*,*) '       node set ',noset
               write(*,*) '       has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*DISTRIBUTING COUPLING%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  node=ialset(j)
                  totweight=totweight+weight
!     
                  if(newmpc.eq.0) then
                     idof=8*(node-1)+1
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,*) 
     &                       '*ERROR reading *DISTRIBUTING COUPLING:'
                           write(*,*) '       dof 1 of node ',node,
     &                          ' is already used'
                           ier=1
                           return
                        endif
                     endif
!
                     nmpc=nmpc+1
                     if(nmpc.gt.nmpc_) then
                        write(*,*) 
     &                    '*ERROR reading *DISTRIBUTING COUPLING:'
                        write(*,*) '       increase nmpc_'
                        ier=1
                        return
                     endif
                     ipompc(nmpc)=mpcfree
                     labmpc(nmpc)='                    '
                     ipompc(nmpc)=mpcfree
!     
!              updating ikmpc and ilmpc
!     
                     do k=nmpc,id+2,-1
                        ikmpc(k)=ikmpc(k-1)
                        ilmpc(k)=ilmpc(k-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
!
                     newmpc=1
                  endif
!
!           new term in MPC
!            
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=1
                  coefmpc(mpcfree)=weight
                  mpcfree=nodempc(3,mpcfree)
!
               else
                  node=ialset(j-2)
                  do
                     node=node-ialset(j)
                     if(node.ge.ialset(j-1)) exit
                     totweight=totweight+weight
!
!                    new term in MPC
!            
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=1
                     coefmpc(mpcfree)=weight
                     mpcfree=nodempc(3,mpcfree)
                  enddo
               endif
            enddo
         endif
      enddo
!
!     reference node
!
      nodempc(1,mpcfree)=irefnode
      nodempc(2,mpcfree)=1
      coefmpc(mpcfree)=-totweight
      mpcfreeold=mpcfree
      mpcfree=nodempc(3,mpcfree)
      nodempc(3,mpcfreeold)=0
!
!     dofs 2 and 3
!
      do idir=2,3
!
         indexm=ipompc(nmpc)
         node=nodempc(1,indexm)
!
         idof=8*(node-1)+idir
         call nident(ikmpc,idof,nmpc,id)
         if(id.gt.0) then
            if(ikmpc(id).eq.idof) then
               write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
               write(*,*) '       dof',idir,' of node ',node,
     &              ' is already used'
               ier=1
               return
            endif
         endif
!     
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR reading *DISTRIBUTING COUPLING:'
            write(*,*) '       increase nmpc_'
            ier=1
            return
         endif
         ipompc(nmpc)=mpcfree
         labmpc(nmpc)='                    '
         ipompc(nmpc)=mpcfree
!     
!     updating ikmpc and ilmpc
!     
         do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
         enddo
         ikmpc(id+1)=idof
         ilmpc(id+1)=nmpc
!
         do
            nodempc(1,mpcfree)=nodempc(1,indexm)
            nodempc(2,mpcfree)=idir
            coefmpc(mpcfree)=coefmpc(indexm)
            if(nodempc(3,indexm).eq.0) then
               mpcfreeold=mpcfree
               mpcfree=nodempc(3,mpcfree)
               nodempc(3,mpcfreeold)=0
               exit
            else
               mpcfree=nodempc(3,mpcfree)
               indexm=nodempc(3,indexm)
            endif
         enddo
      enddo
!
      return
      end

