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
      subroutine checktime(itpamp,namta,tinc,ttime,amta,tmin,inext,itp,
     &  istep,tper)
!
!     checks whether tmin does not exceed the first time point,
!     in case a time points amplitude is active
!
      implicit none
!
      integer namta(3,*),itpamp,id,inew,inext,istart,iend,itp,istep
!
      real*8 amta(2,*),tinc,ttime,tmin,reftime,tper
!
      if(itpamp.gt.0) then
!
!        identifying the location in the time points amplitude
!        of the starting time of the step
!
!        for time points amplitudes based on total time the inext
!        value from the previous step should be used starting with the
!        second step
!
         if((namta(3,itpamp).ge.0).or.(inext.eq.0)) then
            if(namta(3,itpamp).lt.0) then
               reftime=ttime
            else
               reftime=0
            endif
            istart=namta(1,itpamp)
            iend=namta(2,itpamp)
            call identamta(amta,reftime,istart,iend,id)
            if(id.lt.istart) then
               inext=istart
            else
               inext=id+1
            endif
         endif
!
!        identifying the location in the time points amplitude
!        of the starting point increased by tinc
!
         if(namta(3,itpamp).lt.0) then
            reftime=ttime+tinc
         else
            reftime=tinc
         endif
         istart=namta(1,itpamp)
         iend=namta(2,itpamp)
         call identamta(amta,reftime,istart,iend,id)
         if(id.lt.istart) then
            inew=istart
         else
            inew=id+1
         endif
!
!        if the end of the new increment is less than a time
!        point by less than 1.e-6 (theta-value) dtheta is
!        enlarged up to this time point
!
         if((inext.eq.inew).and.(inew.le.iend)) then
            if(amta(1,inew)-reftime.lt.1.d-6*tper) inew=inew+1
         endif
!
!        if the next time point precedes tinc or tmin
!        appropriate action must be taken
!
         if(inew.gt.inext) then
            if(namta(3,itpamp).lt.0) then
               tinc=amta(1,inext)-ttime
            else
               tinc=amta(1,inext)
            endif
            inext=inext+1
            itp=1
            if(tinc.lt.tmin) then
               write(*,*) '*ERROR in checktime: a time point'
               write(*,*) '       precedes the minimum time tmin'
               call exit(201)
            else
               write(*,*) '*WARNING in checktime: a time point'
               write(*,*) '         precedes the initial time'
               write(*,*) '         increment tinc; tinc is'
               write(*,*) '         decreased to ',tinc
            endif
         endif
      endif
!
      return
      end
