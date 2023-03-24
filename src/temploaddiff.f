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
      subroutine temploaddiff(xforcold,xforc,xforcact,iamforc,nforc,
     &     xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
     &     xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,
     &     amta,namta,nam,ampli,time,reltime,ttime,dtime,ithermal,
     &     nmethod,xbounold,xboun,xbounact,iamboun,nboun,
     &     nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
     &     co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
     &     xforcdiff,xloaddiff,xbodydiff,t1diff,xboundiff,iabsload,
     &     iprescribedboundary,ntrans,trab,inotr,veold,nactdof,bcont,fn,
     &     ipobody,iponoel,inoel,ipkon,kon,lakon,ielprop,prop,ielmat,
     &     shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
!     
!     calculates the loading at a given time and the difference with
!     the last call of temploaddiff: is needed in the modal dynamic
!     procedure (dyna.c, dynacont.c; speeds up the calculation)
!     
      implicit none
!     
      logical gasnode
!     
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 amname(*)
!     
      integer iamforc(*),iamload(2,*),iamt1(*),nelemload(2,*),
     &     nam,i,istart,iend,id,nforc,nload,nk,namta(3,*),ithermal(*),
     &     nmethod,iamt1i,iamboun(*),nboun,iamforci,iambouni,
     &     iamloadi1,iamloadi2,ibody(3,*),itg(*),ntg,idof,inotr(2,*),
     &     nbody,iambodyi,nodeboun(*),ndirboun(*),nodeforc(2,*),
     &     ndirforc(*),istep,iinc,msecpt,node,j,ikboun(*),ilboun(*),
     &     ipresboun,mi(*),iabsload,iprescribedboundary,ntrans,
     &     nactdof(0:mi(2),*),ipobody(2,*),iponoel(*),inoel(2,*),
     &     ipkon(*),kon(*),ielprop(*),ielmat(mi(3),*),nshcon(*),
     &     nrhcon(*),ncocon(2,*),ntmat_
!     
      real*8 xforc(*),xforcact(*),xload(2,*),xloadact(2,*),
     &     t1(*),t1act(*),amta(2,*),ampli(*),time,xforcdiff(*),
     &     xforcold(*),xloadold(2,*),t1old(*),reltime,xloaddiff(2,*),
     &     xbounold(*),xboun(*),xbounact(*),ttime,dtime,reftime,
     &     xbody(7,*),xbodyold(7,*),xbodydiff(7,*),t1diff(*),
     &     xbodyact(7,*),co(3,*),vold(0:mi(2),*),abqtime(2),coords(3),
     &     xboundiff(*),trab(7,*), veold(0:mi(2),*),bcont(*),
     &     prop(*),shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),
     &     cocon(0:6,ntmat_,*),fn(0:mi(2),*)
!     
!     
!     
      data msecpt /1/
!     
!     if an amplitude is active, the loading is scaled according to
!     the actual time. If no amplitude is active, then the load is
!     - scaled according to the relative time for a static step
!     - applied as a step loading for a dynamic step
!     
!     calculating all amplitude values for the current time
!     
      do i=1,nam
        if(namta(3,i).lt.0) then
          reftime=ttime+time
        else
          reftime=time
        endif
        if(abs(namta(3,i)).ne.i) then
          reftime=reftime-amta(1,namta(1,i))
          istart=namta(1,abs(namta(3,i)))
          iend=namta(2,abs(namta(3,i)))
          if(istart.eq.0) then
            call uamplitude(reftime,amname(namta(3,i)),ampli(i))
            cycle
          endif
        else
          istart=namta(1,i)
          iend=namta(2,i)
          if(istart.eq.0) then
            call uamplitude(reftime,amname(i),ampli(i))
            cycle
          endif
        endif
        call identamta(amta,reftime,istart,iend,id)
        if(id.lt.istart) then
          ampli(i)=amta(2,istart)
        elseif(id.eq.iend) then
          ampli(i)=amta(2,iend)
        else
          ampli(i)=amta(2,id)+(amta(2,id+1)-amta(2,id))
     &         *(reftime-amta(1,id))/(amta(1,id+1)-amta(1,id))
        endif
      enddo
!     
!     scaling the boundary conditions
!     
      if(iprescribedboundary.eq.1) then
        do i=1,nboun
          if((xboun(i).lt.1.2357111318d0).and.
     &         (xboun(i).gt.1.2357111316d0)) then
!     
!     user subroutine for boundary conditions
!     
            node=nodeboun(i)
!     
!     check whether node is a gasnode
!     
            gasnode=.false.
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                gasnode=.true.
              endif
            endif
!     
            abqtime(1)=time
            abqtime(2)=ttime+time
!     
!     a gasnode cannot move (displacement DOFs are used
!     for other purposes, e.g. mass flow and pressure)
!     
            if(gasnode) then
              do j=1,3
                coords(j)=co(j,node)
              enddo
            else
              do j=1,3
                coords(j)=co(j,node)+vold(j,node)
              enddo
            endif
!     
            if(iabsload.eq.0) then
              xboundiff(i)=xbounact(i)
            else
              xboundiff(i)=xbounact(i)-xboundiff(i)
            endif
            if(ndirboun(i).eq.0) then
              call utemp(xbounact(i),msecpt,istep,iinc,abqtime,node,
     &             coords,vold,mi,iponoel,inoel,
     &             ipobody,xbodyact,ibody,ipkon,kon,
     &             lakon,ielprop,prop,ielmat,
     &             shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
            else
              call uboun(xbounact(i),istep,iinc,abqtime,node,
     &             ndirboun(i),coords,vold,mi,iponoel,inoel,
     &             ipobody,xbodyact,ibody,ipkon,kon,
     &             lakon,ielprop,prop,ielmat,
     &             shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon)
            endif
            xboundiff(i)=xbounact(i)-xboundiff(i)
            cycle
          endif
!     
          if(nam.gt.0) then
            iambouni=iamboun(i)
          else
            iambouni=0
          endif
!     
          if(iabsload.eq.0) then
            xboundiff(i)=xbounact(i)
          else
            xboundiff(i)=xbounact(i)-xboundiff(i)
          endif
          if(iambouni.gt.0) then
            if(amname(iambouni).eq.'RAMP12357111317') then
              xbounact(i)=xbounold(i)+
     &             (xboun(i)-xbounold(i))*reltime
            else
              xbounact(i)=xboun(i)*ampli(iambouni)
            endif
          elseif(nmethod.eq.1) then
            xbounact(i)=xbounold(i)+
     &           (xboun(i)-xbounold(i))*reltime
          else
            xbounact(i)=xboun(i)
          endif
          xboundiff(i)=xbounact(i)-xboundiff(i)
        enddo
      endif
!     
!     scaling the loading
!     
      do i=1,nforc
        if(ndirforc(i).eq.0) then
          if((xforc(i).lt.1.2357111318d0).and.
     &         (xforc(i).gt.1.2357111316d0)) then
            iabsload=2
!     
!     user subroutine for the concentrated heat flux
!     
            node=nodeforc(1,i)
!     
!     check whether node is a gasnode
!     
            gasnode=.false.
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
              if(itg(id).eq.node) then
                gasnode=.true.
              endif
            endif
!     
            abqtime(1)=time
            abqtime(2)=ttime+time
!     
!     a gasnode cannot move (displacement DOFs are used
!     for other purposes, e.g. mass flow and pressure)
!     
            if(gasnode) then
              do j=1,3
                coords(j)=co(j,node)
              enddo
            else
              do j=1,3
                coords(j)=co(j,node)+vold(j,node)
              enddo
            endif
!     
            if(iabsload.eq.0) then
              xforcdiff(i)=xforcact(i)
            else
              xforcdiff(i)=xforcact(i)-xforcdiff(i)
            endif
            call cflux(xforcact(i),msecpt,istep,iinc,abqtime,node,
     &           coords,vold,mi,ipkon,kon,lakon,iponoel,inoel,
     &           ielprop,prop,ielmat,shcon,nshcon,rhcon,nrhcon,
     &           ntmat_,cocon,ncocon)
            xforcdiff(i)=xforcact(i)-xforcdiff(i)
            cycle
          endif
        else
          if((xforc(i).lt.1.2357111318d0).and.
     &         (xforc(i).gt.1.2357111316d0)) then
            iabsload=2
!     
!     user subroutine for the concentrated force
!     
            node=nodeforc(1,i)
!     
            abqtime(1)=time
            abqtime(2)=ttime+time
!     
            do j=1,3
              coords(j)=co(j,node)+vold(j,node)
            enddo
!     
            if(iabsload.eq.0) then
              xforcdiff(i)=xforcact(i)
            else
              xforcdiff(i)=xforcact(i)-xforcdiff(i)
            endif
            call cload(xforcact(i),istep,iinc,abqtime,node,
     &           ndirforc(i),coords,vold,mi,ntrans,trab,inotr,veold)
            xforcdiff(i)=xforcact(i)-xforcdiff(i)
            cycle
          endif
        endif
        if(nam.gt.0) then
          iamforci=iamforc(i)
        else
          iamforci=0
        endif
!     
        if(iabsload.eq.0) then
          xforcdiff(i)=xforcact(i)
        else
          xforcdiff(i)=xforcact(i)-xforcdiff(i)
        endif
        if(iamforci.gt.0) then
          if(amname(iamforci).eq.'RAMP12357111317') then
            xforcact(i)=xforcold(i)+
     &           (xforc(i)-xforcold(i))*reltime
          else
            xforcact(i)=xforc(i)*ampli(iamforci)
          endif
        elseif(nmethod.eq.1) then
          xforcact(i)=xforcold(i)+
     &         (xforc(i)-xforcold(i))*reltime
        else
          xforcact(i)=xforc(i)
        endif
        xforcdiff(i)=xforcact(i)-xforcdiff(i)
      enddo
!     
      do i=1,nload
!     
!     check for dload subroutine
!     
        if(sideload(i)(3:4).eq.'NU') then
          iabsload=2
          cycle
        endif
!     
        ipresboun=0
!     
!     check for pressure boundary conditions
!     
        if(sideload(i)(3:4).eq.'NP') then
          node=nelemload(2,i)
          idof=8*(node-1)+2
          call nident(ikboun,idof,nboun,id)
          if(id.gt.0) then
            if(ikboun(id).eq.idof) then
              ipresboun=1
              if(iabsload.eq.0) then
                xloaddiff(1,i)=xloadact(1,i)
              else
                xloaddiff(1,i)=xloadact(1,i)-xloaddiff(1,i)
              endif
              xloadact(1,i)=xbounact(ilboun(id))
              xloaddiff(1,i)=xloadact(1,i)-xloaddiff(1,i)
            endif
          endif
        endif
!     
        if(ipresboun.eq.0) then
          if(nam.gt.0) then
            iamloadi1=iamload(1,i)
            iamloadi2=iamload(2,i)
          else
            iamloadi1=0
            iamloadi2=0
          endif
!     
          if(iabsload.eq.0) then
            xloaddiff(1,i)=xloadact(1,i)
          else
            xloaddiff(1,i)=xloadact(1,i)-xloaddiff(1,i)
          endif
          if(iamloadi1.gt.0) then
            if(amname(iamloadi1).eq.'RAMP12357111317') then
              xloadact(1,i)=xloadold(1,i)+
     &             (xload(1,i)-xloadold(1,i))*reltime
            else
              xloadact(1,i)=xload(1,i)*ampli(iamloadi1)
            endif
          elseif(nmethod.eq.1) then
            xloadact(1,i)=xloadold(1,i)+
     &           (xload(1,i)-xloadold(1,i))*reltime
          else
            xloadact(1,i)=xload(1,i)
          endif
          xloaddiff(1,i)=xloadact(1,i)-xloaddiff(1,i)
!     
          if(iabsload.eq.0) then
            xloaddiff(2,i)=xloadact(1,i)
          else
            xloaddiff(2,i)=xloadact(2,i)-xloaddiff(2,i)
          endif
          if(iamloadi2.gt.0) then
            if(amname(iamloadi2).eq.'RAMP12357111317') then
              xloadact(2,i)=xloadold(2,i)+
     &             (xload(2,i)-xloadold(2,i))*reltime
            else
              xloadact(2,i)=xload(2,i)*ampli(iamloadi2)
            endif
          elseif(nmethod.eq.1) then
            xloadact(2,i)=xload(2,i)
          else
            xloadact(2,i)=xload(2,i)
          endif
          xloaddiff(2,i)=xloadact(2,i)-xloaddiff(2,i)
        endif
      enddo
!     
      do i=1,nbody
        if(nam.gt.0) then
          iambodyi=ibody(2,i)
        else
          iambodyi=0
        endif
!     
        if(iabsload.eq.0) then
          xbodydiff(1,i)=xbodyact(1,i)
        else
          xbodydiff(1,i)=xbodyact(1,i)-xbodydiff(1,i)
        endif
        if(iambodyi.gt.0) then
          if(amname(iambodyi).eq.'RAMP12357111317') then
            xbodyact(1,i)=xbodyold(1,i)+
     &           (xbody(1,i)-xbodyold(1,i))*reltime
          else
            xbodyact(1,i)=xbody(1,i)*ampli(iambodyi)
          endif
        elseif(nmethod.eq.1) then
          xbodyact(1,i)=xbodyold(1,i)+
     &         (xbody(1,i)-xbodyold(1,i))*reltime
        else
          xbodyact(1,i)=xbody(1,i)
        endif
        xbodydiff(1,i)=xbodyact(1,i)-xbodydiff(1,i)
      enddo
!     
      return
      end
