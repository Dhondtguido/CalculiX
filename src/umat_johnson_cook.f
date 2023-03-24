!
!
!     Based on a proposal by Lukas Mazza, Matthias Mueller and
!                            Otto-Ernst Bernhardi
!
!     See also "Writing User Subroutines with ABAQUS", HKS.
!
!     Permission for use with the GPL license granted by Prof. Otto
!     Bernhardi in March 2023
!
      subroutine umat_johnson_cook(stre,statev,ddsdde,spd,rpl,
     &     stran,dtime,temp,props,nprops,noel,npt,kstep,kinc,icmd) 
!     
      implicit none
!     
      integer newton,nprops,kinc,i,j,k,noel,npt,kstep,icmd
!     
      real*8 stre(6),statev(*),ddsdde(6,6),stran(6),
     &     props(nprops),el(6),epl(6),flow(6),toler,spd,rpl,dtime,
     &     temp,a,b,c,dep,ebulk3,eff_um,eff_um2,
     &     eff_um3,eff_hard,eff_al,um,um2,um3,al,xm,ee,xn,un,ep,
     &     hard,svminv,rhs,p,svm,syield,syield0,constant,
     &     epsdot0,tmelt,ttransition
!     
      data newton,toler/40,1.d-6/ 
!     
!     local arrays 
!-----------------------------------------------------------------------
!     el -- elastic strains 
!     epl -- plastic strains 
!     flow  -- direction of plastic flow 
!-----------------------------------------------------------------------
!     
!-----------------------------------------------------------------------
!     umat for johnson-cook model 
!-----------------------------------------------------------------------
!     props(1) - young's modules 
!     props(2) - poisson ratio 
!     props(3) - inelastic heat fraction 
!     props(4) - A 
!     props(5) - B 
!     props(6) - n 
!     props(7) - C 
!     props(8) - epsilon_0^dot 
!     props(9) - theta_melt 
!     props(10) - theta_transition
!     props(11) - m 
!-----------------------------------------------------------------------
!     depvars 
!-----------------------------------------------------------------------
!     statev(1)   - equivalent plastic strain 
!     statev(2-7) - plastic strain tensor 
!     statev(8)   - plastic strain rate
!     statev(9)   - number of interations k
!-----------------------------------------------------------------------
!     
!     elastic properties 
!     
      ee=props(1) 
      un=props(2) 
      ebulk3=ee/(1.d0-2.d0*un) 
      um2=ee/(1.d0+un) 
      um=um2/2.d0 
      um3=um*3.d0 
      al=(ebulk3-um2)/3.d0 
!     
!     reload plastic strain from depvar 
!     
      ep=statev(1) 
      do i=1,6 
        epl(i)=statev(i+1)
        el(i)=stran(i)-epl(i)
      end do 
      dep=statev(8)*dtime
!     
!     calculate predictor stress from elastic strains 
!     
      stre(1)=(al+um2)*el(1)+al*el(2)+al*el(3)
      stre(2)=al*el(1)+(al+um2)*el(2)+al*el(3)
      stre(3)=al*el(1)+al*el(2)+(al+um2)*el(3)
      stre(4)=um*el(4)
      stre(5)=um*el(5)
      stre(6)=um*el(6)
!     
!     calculate mises stress 
!     
      svm=(stre(1)-stre(2))*(stre(1)-stre(2))+ 
     &     (stre(2)-stre(3))*(stre(2)-stre(3))+ 
     &     (stre(3)-stre(1))*(stre(3)-stre(1)) 
      do i=4,6 
        svm=svm+6.d0*stre(i)*stre(i) 
      end do 
      svm=dsqrt(svm/2.d0) 
!     
!     read johnson-cook parameter 
!     
      a=props(4) 
      b=props(5) 
      xn=props(6) 
      c=props(7)
      epsdot0=props(8)
      tmelt=props(9)
      ttransition=props(10)
      xm=props(11)
!
      if((tmelt.le.ttransition).or.(temp.lt.ttransition)) then
        constant=1.d0
      elseif(temp.lt.tmelt) then
        constant=1.d0-((temp-ttransition)/(tmelt-ttransition))**xm
      else
        constant=0.d0
      endif
!     
!     call userhard subroutine to get hardening rate and yield stress 
!     
      call userhard(syield0,hard,ep,dep,dtime,a,b,xn,c,epsdot0,
     &     constant)
!     
!     determine if yielding happened 
!     
      if(svm.gt.syield0) then 
!     
!     plastic response,determine flow direction 
!     
        p=(stre(1)+stre(2)+stre(3))/3.d0 
        svminv=1.d0/svm 
        do i=1,3 
          flow(i)=svminv*(stre(i)-p) 
        end do 
        do i=4,6 
          flow(i)=svminv*stre(i) 
        end do 
!     
!     newton iteration 
!     
        dep=(svm-syield0)/um3
        do k=1,newton 
          call userhard(syield,hard,ep+dep,dep,dtime,a,b,xn,c,
     &         epsdot0,constant)
          rhs=svm-um3*dep-syield
          dep=dep+rhs/(um3+hard)  
          if(dabs(rhs).le.dabs(toler*syield0)+toler) exit
        end do
!
!       check whether the loop converged
!
        if(k.gt.newton) then
          write(*,*) '*ERROR in umat_johnson_cook.f:'
          write(*,100) newton 
 100      format('      plasticity algorithm did not ',
     &         'converge after ',i3,' iterations') 
          write(*,*) '      k=',k
          write(*,*) '      a=',a,' b=',b
          write(*,*) '      xn=',xn,' c=',c
          write(*,*) '      smises=',svm
          write(*,*) '      um3=',um3
          write(*,*) '      ep=',ep
          write(*,*) '      dep=',dep
          write(*,*) '      dtime=',dtime
          write(*,*) '      epsdot=',dep/dtime
          write(*,*) '      syield=',syield,'syield0=',syield0
          write(*,*) '      hard=',hard
          write(*,*) '      rhs=',rhs
          call exit(201)
        endif
!     
!       check for very large equivalent plastic strain rate     
!     
        if(dep.gt.1000.d0) then
          write(*,*) '*WARNING in umat_johnson_cook.f: very large' 
          write(*,*) '         effective plastic strain, ep=',ep 
        endif
!        
        eff_hard=um3*hard/(um3+hard) 
!     
!     calculate stress and update strains 
!     
        do i=1,3 
          stre(i)=flow(i)*syield+p 
          epl(i)=epl(i)+3.d0*flow(i)*dep/2.d0 
        end do 
        do i=4,6 
          stre(i)=flow(i)*syield 
          epl(i)=epl(i)+3.d0*flow(i)*dep 
        end do 
        ep=ep+dep 
        spd=dep*(syield0+syield)/2.d0 
        rpl=props(3)*spd/dtime 
!     
!     jacobian 
!     
        if(icmd.ne.3) then
          eff_um=um*syield/svm 
          eff_um2=eff_um*2.d0 
          eff_um3=eff_um*3.d0
          eff_al=(ebulk3-eff_um2)/3.d0
!     
          do i=4,6 
            do j=1,i-1 
              ddsdde(i,j)=0.d0 
            end do 
          end do 
          do i=1,3 
            do j=1,i-1 
              ddsdde(j,i)=eff_al 
            end do 
            ddsdde(i,i)=eff_um2+eff_al
          end do 
          do i=4,6 
            ddsdde(i,i)=eff_um 
          end do 
          do i=1,6 
            do j=1,i 
              ddsdde(j,i)=ddsdde(j,i)+flow(j)*flow(i)*(eff_hard-eff_um3)
            end do 
          end do
        endif
      else
        dep=0.d0
!     
!     elastic jocobian 
!
        if(icmd.ne.3) then
          do i=4,6 
            do j=1,i-1 
              ddsdde(i,j)=0.d0 
            end do 
          end do 
!     
          do i=1,3 
            do j=1,i-1 
              ddsdde(j,i)=al 
            end do 
            ddsdde(i,i)=al+um2
          end do 
          do i=4,6 
            ddsdde(i,i)=um 
          end do
        endif
      end if 
!     
!     store strains in state variable array 
!     
      statev(1)=ep 
      do i=1,6 
        statev(i+1)=epl(i) 
      end do 
      statev(8)=dep/dtime
      statev(9)=1.d0*k
!     
      return 
      end 
!     
      subroutine userhard(syield,hard,ep,dep,dtime,a,b,xn,c,
     &     epsdot0,constant)
!     
      implicit none
!     
      real*8 syield,hard,ep,dep,dtime,a,b,xn,c,epsdot,delta,
     &     epsdot0,constant,xstrain,xrate,dxstrain
!     
!     calculate current yield stress and hardening rate 
!     
      epsdot=dep/(dtime*epsdot0)
      delta=0.0001d0
!     
      if(epsdot.ge.1.d0)then
        if(ep.gt.delta) then
          dxstrain=b*ep**(xn-1.d0)
          xstrain=a+dxstrain*ep
          xrate=1.d0+c*dlog(epsdot)
!     
          syield=xstrain*xrate
          hard=xstrain*c/(epsdot*dtime)
     &         +xrate*xn*dxstrain
        else
          dxstrain=b*delta**(xn-1.d0)
          xstrain=a+ep*dxstrain
          xrate=1.d0+c*dlog(epsdot)
!     
          syield=xstrain*xrate
          hard=xstrain*c/(epsdot*dtime)
     &         +xrate*dxstrain
        endif
!     
      else
        if(ep.gt.delta) then
          dxstrain=b*ep**(xn-1.d0)
          xstrain=a+dxstrain*ep
          xrate=1.d0+c*(epsdot-1.d0)
!     
          syield=xstrain*xrate
          hard=xstrain*c/(epsdot0*dtime)
     &         +xrate*xn*dxstrain
        else
          dxstrain=b*delta**(xn-1.d0)
          xstrain=a+ep*dxstrain
          xrate=1.d0+c*(epsdot-1.d0)
!     
          syield=xstrain*xrate
          hard=xstrain*c/(epsdot0*dtime)
     &         +xrate*dxstrain
        endif
      endif
!
      syield=syield*constant
      hard=hard*constant
!     
      return 
      end        
