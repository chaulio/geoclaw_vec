subroutine riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR, &
        bL,bR,uL,uR,vL,vR,phiL,phiR,s1,s2,drytol,g,& !sw,fw)
sw1, sw2, sw3,fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33)
    !dir$ attributes forceinline :: riemann_fwave
    !dir$ attributes vector: uniform(meqn,mwaves,drytol,g) :: riemann_fwave
 
    ! solve shallow water equations given single left and right states
    ! solution has two waves.
    ! flux - source is decomposed.
    implicit none

    integer, parameter :: DP = kind(1.d0)
    !input
    integer meqn,mwaves

    real(kind=DP) :: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,s1,s2
    real(kind=DP) :: hvL,hvR,vL,vR
    real(kind=DP) :: drytol,g

    !real(kind=DP) :: sw(mwaves)
    !real(kind=DP) :: fw(meqn,mwaves)
    real(kind=8), intent(inout) :: sw1, sw2, sw3
    real(kind=8), intent(inout) :: fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33

    !local
    real(kind=DP) :: delh,delhu,delphi,delb,delhdecomp,delphidecomp
    real(kind=DP) :: deldelh,deldelphi
    real(kind=DP) :: beta1,beta2

    !determine del vectors
    delh = hR-hL
    delhu = huR-huL
    delphi = phiR-phiL
    delb = bR-bL

    deldelphi = -g*0.5d0*(hR+hL)*delb
    delphidecomp = delphi - deldelphi

    if (s2 /= s1) then
        !flux decomposition
        beta1 = (s2*delhu - delphidecomp)/(s2-s1)
        beta2 = (delphidecomp - s1*delhu)/(s2-s1)
    else
        beta1 = 0.d0
        beta2 = 0.d0
    endif        

!    sw(1) = s1
!    sw(2) = 0.5d0*(s1+s2)
!    sw(3) = s2
!    ! 1st nonlinear wave
!    fw(1,1) = beta1
!    fw(2,1) = beta1*s1
!    fw(3,1) = beta1*vL
!    ! 2nd nonlinear wave
!    fw(1,3) = beta2
!    fw(2,3) = beta2*s2
!    fw(3,3) = beta2*vR
!    ! advection of transverse wave
!    fw(1,2) = 0.d0
!    fw(2,2) = 0.d0
!    fw(3,2) = hR*uR*vR - hL*uL*vL -fw(3,1)-fw(3,3)
    sw1 = s1
    sw2 = 0.5d0*(s1+s2)
    sw3 = s2
    ! 1st nonlinear wave
    fw11 = beta1
    fw21 = beta1*s1
    fw31 = beta1*vL
    ! 2nd nonlinear wave
    fw13 = beta2
    fw23 = beta2*s2
    fw33 = beta2*vR
    ! advection of transverse wave
    fw12 = 0.d0
    fw22 = 0.d0
    fw32 = hR*uR*vR - hL*uL*vL -fw31-fw33
    return
end subroutine

subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2, &
        maxiter,drytol,g)
    !determine the Riemann structure (wave-type in each family)
    implicit none

    integer, parameter :: DP = kind(1.d0)
    !input
    real(kind=DP) :: hL,hR,uL,uR,drytol,g
    integer maxiter

    !output
    real(kind=DP) :: s1m,s2m
    logical rare1,rare2

    !local
    real(kind=DP) :: hm,u1m,u2m,um,delu
    real(kind=DP) :: h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
    integer iter

    !Test for Riemann structure
    h_min = min(hR,hL)
    h_max = max(hR,hL)
    delu = uR-uL

    if (h_min.le.drytol) then
        hm=0.d0
        um=0.d0
        s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
        s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
        if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
        else
            rare1=.true.
            rare2=.false.
        endif

    else
        F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
        F_max= delu + &
            (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

        if (F_min.gt.0.d0) then !2-rarefactions
            hm=(1.d0/(16.d0*g))* &
                max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
            um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

            s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
            s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
            rare1=.true.
            rare2=.true.
        elseif (F_max.le.0.d0) then !2 shocks
            ! root finding using a Newton iteration on sqrt(h)===
            h0=h_max
            do iter=1,maxiter
                gL=sqrt(.5d0*g*(1/h0 + 1/hL))
                gR=sqrt(.5d0*g*(1/h0 + 1/hR))
                F0=delu+(h0-hL)*gL + (h0-hR)*gR
                dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+ &
                    gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
                slope=2.d0*sqrt(h0)*dfdh
                h0=(sqrt(h0)-F0/slope)**2
            enddo
            hm=h0
            u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
            u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
            um=.5d0*(u1m+u2m)
            s1m=u1m-sqrt(g*hm)
            s2m=u2m+sqrt(g*hm)
            rare1=.false.
            rare2=.false.
        else !one shock one rarefaction
            h0=h_min
            do iter=1,maxiter
                F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
                    + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
                slope=(F_max-F0)/(h_max-h_min)
                h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
                um    = uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
                s1m   = uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
                s2m   = uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
                rare1 = .true.
                rare2 = .false.
            else
                s2m   = uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
                s1m   = uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
                um    = uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
                rare2 = .true.
                rare1 = .false.
            endif
        endif
    endif

    return
end subroutine

