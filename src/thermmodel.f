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
      subroutine thermmodel(amat,iel,iint,kode,coconloc,vkl,
     &  dtime,time,ttime,mi,nstate_,xstateini,xstate,qflux,xstiff,
     &  iorien,pgauss,orab,t1l,t1lold,vold,co,lakonl,konl,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,nmethod,iperturb)
!
      implicit none
!
      character*8 lakonl
      character*80 amat
!
      integer iel,iint,kode,mi(*),nstate_,iorien,ntgrd,ncoconst,
     &  layer,kspt,kstep,kinc,kal(2,6),konl(20),ipompc(*),nmethod,
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),j3,j2,j4,jj,j,i,j1,
     &  iperturb(*)
!
      real*8 coconloc(*),vkl(0:3,3),dtime,time,ttime,cond(6),
     &  xstateini(nstate_,mi(1),*),xstate(nstate_,mi(1),*),qflux(3),
     &  pgauss(3),orab(7,*),abqtime(2),u,dudt,dudg(3),dfdt(3),
     &  dfdg(3,3),dtemp,dtemdx(3),predef(1),dpred(1),pnewdt,
     &  skl(3,3),t1lold,xstiff(27,mi(1),*),xa(3,3),vold(0:mi(2),*),
     &  co(3,*),coefmpc(*),t1l
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      if(kode.eq.1) then
!
!         linear isotropic
!
         do i=1,3
            cond(i)=coconloc(1)
         enddo
         do i=4,6
            cond(i)=0.d0
         enddo
!
         do i=1,3
            qflux(i)=-coconloc(1)*vkl(0,i)
         enddo
!
      elseif((kode.eq.3).or.(kode.eq.6)) then
         if((kode.eq.3).and.(iorien.eq.0)) then
!
!           orthotropic
!
            do i=1,3
               cond(i)=coconloc(i)
            enddo
            do i=4,6
               cond(i)=0.d0
            enddo
!     
            do i=1,3
               qflux(i)=-coconloc(i)*vkl(0,i)
            enddo
!
         else
            if(iorien.ne.0) then
!
!              transformation due to special orientation
!
!              calculating the transformation matrix
!
               call transformatrix(orab(1,iorien),pgauss,skl)
!
!              modifying the conductivity constants
!
               if(kode.eq.3) then
                  do j=4,6
                     coconloc(j)=0.d0
                  enddo
               endif
!     
               xa(1,1)=coconloc(1)
               xa(1,2)=coconloc(4)
               xa(1,3)=coconloc(5)
               xa(2,1)=coconloc(4)
               xa(2,2)=coconloc(2)
               xa(2,3)=coconloc(6)
               xa(3,1)=coconloc(5)
               xa(3,2)=coconloc(6)
               xa(3,3)=coconloc(3)
!
               do jj=1,6
                  coconloc(jj)=0.d0
                  j1=kal(1,jj)
                  j2=kal(2,jj)
                  do j3=1,3
                     do j4=1,3
                        coconloc(jj)=coconloc(jj)+
     &                       xa(j3,j4)*skl(j1,j3)*skl(j2,j4)
                     enddo
                  enddo
               enddo
            endif
!
!           anisotropy
!
            do i=1,6
               cond(i)=coconloc(i)
            enddo
!     
            qflux(1)=-coconloc(1)*vkl(0,1)-coconloc(4)*vkl(0,2)-
     &           coconloc(5)*vkl(0,3)
            qflux(2)=-coconloc(4)*vkl(0,1)-coconloc(2)*vkl(0,2)-
     &           coconloc(6)*vkl(0,3)
            qflux(3)=-coconloc(5)*vkl(0,1)-coconloc(6)*vkl(0,2)-
     &           coconloc(3)*vkl(0,3)
!
         endif
      else
!
!        user material
!
         ncoconst=-kode-100
!
         do i=1,nstate_
            xstate(i,iint,iel)=xstateini(i,iint,iel)
         enddo
!
         abqtime(1)=time-dtime
         abqtime(2)=ttime+time-dtime
!
         ntgrd=3
         dtemp=t1l-t1lold
         do i=1,3
            dtemdx(i)=vkl(0,i)
         enddo
!
         call umatht(u,dudt,dudg,qflux,dfdt,dfdg,xstate(1,iint,iel),
     &     t1lold,dtemp,
     &     dtemdx,abqtime,dtime,predef,dpred,amat,ntgrd,nstate_,
     &     coconloc,ncoconst,pgauss,pnewdt,iel,iint,layer,kspt,
     &     kstep,kinc,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi)
!
         cond(1)=dfdg(1,1)
         cond(2)=dfdg(2,2)
         cond(3)=dfdg(3,3)
         cond(4)=dfdg(1,2)
         cond(5)=dfdg(1,3)
         cond(6)=dfdg(2,3)
!
      endif
!
      if(((nmethod.ne.4).or.(iperturb(1).ne.0)).and.
     &    (nmethod.ne.5)) then
         do i=1,6
            xstiff(21+i,iint,iel)=cond(i)
         enddo
      endif
!
      return
      end
