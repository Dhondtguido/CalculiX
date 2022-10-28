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
      subroutine resultsforc_se(nk,dfn,nactdofinv,ipompc,nodempc,
     &  coefmpc,nmpc,mi,fmpc,calcul_fn,calcul_f,idesvar,df,
     &  jqs,irows,distmin)
!
!     calculating the equation system internal force vector
!     (one entry for each active degree of freedom)
!
      implicit none
!
      integer mi(*),nactdofinv(*),ipompc(*),nodempc(3,*),nk,i,j,
     &  nmpc,ist,ndir,node,index,calcul_fn,calcul_f,idesvar,jqs(*),
     &  irows(*),idir,inode,idof,mt
!
      real*8 dfn(0:mi(2),*),coefmpc(*),fmpc(*),forcempc,distmin,
     &  df(*),val
!
!
!
      mt=mi(2)+1
!
!     subtracting the mpc force (for each linear mpc there is one
!     force; the actual force in a node belonging to the mpc is
!     obtained by multiplying this force with the nodal coefficient.
!     The force has to be subtracted from f, since it does not
!     appear on the rhs of the equation system
!
      if(calcul_fn.eq.1)then
        do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.gt.3) cycle
            forcempc=dfn(ndir,node)/coefmpc(ist)
            fmpc(i)=forcempc
            dfn(ndir,node)=0.d0
            index=nodempc(3,ist)
            if(index.eq.0) cycle
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               dfn(ndir,node)=dfn(ndir,node)-
     &coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         enddo
      endif
!
!     calculating the system force vector
!
      if(calcul_f.eq.1) then
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            inode=nactdofinv(idof)
            node=inode/mt+1
            idir=inode-mt*(inode/mt)
            val=-dfn(idir,node)/distmin
            call add_bo_st(df,jqs,irows,irows(j),
     &                           idesvar,val)
         enddo
      endif
!
      return
      end
