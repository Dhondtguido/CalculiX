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
      subroutine filterforward_imp(ad,au,adb,aub,feasdir,gradproj,rhs,
     &   ndesi,nodedesi,iflag,jq,irow,objectset)           
!     
!     calculates the values of the filter matrix
!     
      implicit none
!
      character*81 objectset(5,*)
!
      integer ndesi,nodedesi(*),idesvar,jdesvar,inode,jnode,iflag,j,
     &   jq(*),irow(*),istat
!     
      real*8 adb(*),aub(*),feasdir(2,*),gradproj(3,*),rhs(*),ad(*),
     &   au(*),filterrad
! 
! 
      if(iflag.eq.0) then
!
!        Multiply the design update vector with the mass matrix
!        and copy to rhs
!        loop over all columns for main- and sub-diagonal
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
!
!           multipliying the main diagonal
            rhs(idesvar)=rhs(idesvar)+gradproj(3,inode)*adb(idesvar)                      
!
!           loop over all rows in sub-diagonal
            do j=jq(idesvar),jq(idesvar+1)-1
               jdesvar=irow(j)
               jnode=nodedesi(jdesvar)
!
!              loop over all rows in sub-diagonal
               rhs(jdesvar)=rhs(jdesvar)+gradproj(3,inode)*aub(j)
!
!              loop over all columns in top-diagonal
               rhs(idesvar)=rhs(idesvar)+gradproj(3,jnode)*aub(j)

            enddo
         enddo
!   
!        Assembly of filter matrix
         read(objectset(2,1)(21:40),'(f20.0)',iostat=istat) filterrad     
!
         do idesvar=1,ndesi
            ad(idesvar)=filterrad**2*ad(idesvar)+adb(idesvar)
            do j=jq(idesvar),jq(idesvar+1)-1
               au(j)=filterrad**2*au(j)+aub(j)
            enddo
         enddo
!         
      elseif(iflag.eq.1) then
!
!        copy the filtered and unfiltered sensitivities from rhs 
!        to feasdir
         do idesvar=1,ndesi
            inode=nodedesi(idesvar)
            feasdir(1,inode)=gradproj(3,inode)
            feasdir(2,inode)=rhs(idesvar)
         enddo
      endif
!
      return        
      end
