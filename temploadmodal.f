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
      subroutine temploadmodal(amta,namta,nam,ampli,time,ttime,dtime,
     &  xbounold,xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,
     &  amname,reltime)
!
!     calculates the SPC boundary conditions at a given time for
!     a modal dynamic procedure; used to calculate the velocity and
!     acceleration by use of finite differences
!
      implicit none
!
      character*80 amname(*)
!
      integer nam,i,istart,iend,id,namta(3,*),
     &  iamboun(*),nboun,iambouni,nodeboun(*),ndirboun(*)
!
      real*8 amta(2,*),ampli(*),time,reltime,
     &  xbounold(*),xboun(*),xbounact(*),ttime,dtime,reftime
!
!     if an amplitude is active, the loading is scaled according to
!     the actual time. If no amplitude is active, then the load is
!     applied as a step loading
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
     &           *(reftime-amta(1,id))/(amta(1,id+1)-amta(1,id))
         endif
      enddo
!
!     scaling the boundary conditions
!
      do i=1,nboun
         if(nam.gt.0) then
            iambouni=iamboun(i)
         else
            iambouni=0
         endif
         if(iambouni.gt.0) then
            if(amname(iambouni).eq.'RAMP12357111317') then
               xbounact(i)=xbounold(i)+
     &              (xboun(i)-xbounold(i))*reltime
            else
               xbounact(i)=xboun(i)*ampli(iambouni)
            endif
         else
            xbounact(i)=xboun(i)
         endif
      enddo
!
      return
      end
