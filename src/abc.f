        d1=dabs(xsj(1))
        d2=dabs(xsj(2))
        d3=dabs(xsj(3))
!
        if((d3.gt.d2).and.(d3.gt.d1)) then
          xsi(1,1)=xs(2,2)/xsj(3)
          xsi(2,2)=xs(1,1)/xsj(3)
          xsi(1,2)=-xs(1,2)/xsj(3)
          xsi(2,1)=-xs(2,1)/xsj(3)
          if(d2.gt.d1) then
            if(d2.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(2))
            else
              xsjmod=xsj(2)
            endif
            xsi(2,3)=xs(1,1)/(-xsjmod)
            xsi(1,3)=-xs(1,2)/(-xsjmod)
          else
            if(d1.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(1))
            else
              xsjmod=xsj(1)
            endif
            xsi(2,3)=xs(2,1)/xsjmod
            xsi(1,3)=-xs(2,2)/xsjmod
          endif
        elseif((d2.gt.d1).and.(d2.gt.d3)) then
          xsi(1,1)=xs(3,2)/(-xsj(2))
          xsi(2,3)=xs(1,1)/(-xsj(2))
          xsi(1,3)=-xs(1,2)/(-xsj(2))
          xsi(2,1)=-xs(3,1)/(-xsj(2))
          if(d1.gt.d3) then
            if(d1.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(1))
            else
              xsjmod=xsj(1)
            endif
            xsi(1,2)=xs(3,2)/xsjmod
            xsi(2,2)=-xs(3,1)/xsjmod
          else
            if(d3.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(3))
            else
              xsjmod=xsj(3)
            endif
            xsi(1,2)=-xs(1,2)/xsjmod
            xsi(2,2)=xs(1,1)/xsjmod
          endif
        else
          xsi(1,2)=xs(3,2)/xsj(1)
          xsi(2,3)=xs(2,1)/xsj(1)
          xsi(1,3)=-xs(2,2)/xsj(1)
          xsi(2,2)=-xs(3,1)/xsj(1)
          if(d3.gt.d2) then
            if(d3.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(3))
            else
              xsjmod=xsj(3)
            endif
            xsi(1,1)=xs(2,2)/xsjmod
            xsi(2,1)=-xs(2,1)/xsjmod
          else
            if(d2.lt.1.d-30) then
              xsjmod=dsign(1.d-30,xsj(2))
            else
              xsjmod=xsj(2)
            endif
            xsi(1,1)=xs(3,2)/(-xsjmod)
            xsi(2,1)=-xs(3,1)/(-xsjmod)
          endif
        endif
