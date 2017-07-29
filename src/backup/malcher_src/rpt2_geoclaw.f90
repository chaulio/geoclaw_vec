! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx, &
        ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
    ! =====================================================
    use geoclaw_module, only: g => grav, tol => dry_tolerance
    use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

    implicit none
    !
    !     # Riemann solver in the transverse direction using an einfeldt
    !     Jacobian.

    !-----------------------last modified 1/10/05----------------------

    integer, parameter :: DP = kind(1.d0)

    real(kind=DP) :: ql(1-mbc:maxm+mbc, meqn)
    real(kind=DP) :: qr(1-mbc:maxm+mbc, meqn)
    real(kind=DP) :: asdq(meqn,1-mbc:maxm+mbc)
    real(kind=DP) :: bmasdq(meqn,1-mbc:maxm+mbc)
    real(kind=DP) :: bpasdq(meqn,1-mbc:maxm+mbc)
    real(kind=DP) :: aux1(1-mbc:maxm+mbc,maux)
    real(kind=DP) :: aux2(1-mbc:maxm+mbc,maux)
    real(kind=DP) :: aux3(1-mbc:maxm+mbc,maux)
    real(kind=DP) :: s(3)
    real(kind=DP) :: r(3,3)
    real(kind=DP) :: beta(3)
    real(kind=DP) :: abs_tol
    real(kind=DP) :: hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
    real(kind=DP) :: uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
    real(kind=DP) :: delf1,delf2,delf3
    real(kind=DP) :: dxdcm,dxdcp,topo1,topo3,eta

    integer :: ixy,maxm,meqn,maux,mwaves,mbc,mx,imp
    integer :: i,m,mw,mu,mv

    abs_tol=tol

    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    !DIR$ VECTOR ALIGNED
    !$OMP SIMD PRIVATE(hl,hr,hul,hur,hvl,hvr,s1,s2,s3, &
    !$OMP& dxdcm,dxdcp,delf1,delf2,delf3,topo1,topo3,beta,s,r)
    do i=2-mbc,mx+mbc
        hl=qr(i-1,1) 
        hr=ql(i,1) 
        hul=qr(i-1,mu) 
        hur=ql(i,mu) 
        hvl=qr(i-1,mv) 
        hvr=ql(i,mv)

        do mw=1,mwaves
            s(mw)=0.d0
            beta(mw)=0.d0
            do m=1,meqn
                r(m,mw)=0.d0
            enddo
        enddo
        
        ! check and see if cell that transverse waves are going in is high and dry
        if (imp.eq.1) then
            eta = qr(i-1,1)  + aux2(i-1,1)
            topo1 = aux1(i-1,1)
            topo3 = aux3(i-1,1)
        else
            eta = ql(i,1) + aux2(i,1)
            topo1 = aux1(i,1)
            topo3 = aux3(i,1)
        endif

        if (eta<max(topo1,topo3)) then
            hl = 0.d0
            hr = 0.d0
        endif 

        call solve_single_rpt(g, tol, hl, hul, hvl, hr, hur, hvr, s1, s2, s3)

        s(1)=s1
        s(2)=s2
        s(3)=s3
        !=======================Determine asdq decomposition (beta)============
        delf1=asdq(1,i)
        delf2=asdq(mu,i)
        delf3=asdq(mv, i)

        beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
        beta(2) = -s2*delf1 + delf2
        beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
        !======================End =================================================

        !=====================Set-up eigenvectors===================================
        r(1,1) = 1.d0
        r(2,1) = s2
        r(3,1) = s1

        r(1,2) = 0.d0
        r(2,2) = 1.d0
        r(3,2) = 0.d0

        r(1,3) = 1.d0
        r(2,3) = s2
        r(3,3) = s3
        !============================================================================
        90      continue
        !============= compute fluctuations==========================================

!        bmasdq(1,i)=0.0d0
!        bmasdq(2,i)=0.0d0
!        bmasdq(3,i)=0.0d0
!        bpasdq(1,i)=0.0d0
!        bpasdq(2,i)=0.0d0
!        bpasdq(3,i)=0.0d0

        dxdcp = 1.d0
        dxdcm = 1.d0

        if (coordinate_system.eq.2) then
            if (ixy.eq.2) then
                dxdcp=(earth_radius*deg2rad)
                dxdcm = dxdcp
            else
                if (imp.eq.1) then
                    dxdcp = earth_radius*cos(aux3(i-1,3))*deg2rad
                    dxdcm = earth_radius*cos(aux1(i-1,3))*deg2rad
                else
                    dxdcp = earth_radius*cos(aux3(i,3))*deg2rad
                    dxdcm = earth_radius*cos(aux1(i,3))*deg2rad
                endif
            endif
        endif

        bmasdq(:,i)=0.0d0
        bpasdq(:,i)=0.0d0
        do mw=1,3
            if (s(mw)<0.d0) then
                bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
            elseif (s(mw)>0.d0) then
                bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
            endif
        enddo
        !========================================================================
    enddo
    return
end subroutine 


subroutine solve_single_rpt(g, tol, hl, hul, hvl, hr, hur, hvr, &
        s1, s2, s3)
    implicit none
    integer, parameter :: DP = kind(1.d0)

    !real(kind=DP), intent(inout) :: dxdcm, dxdcp
    real(kind=DP), intent(in) :: g, tol, hul, hvl, hur, hvr
    real(kind=DP), intent(out) :: s1, s2, s3
    real(kind=DP), intent(inout) :: hl, hr 
    real(kind=DP) :: ul, vl, ur, vr, roe1, roe3, s1l, s3r,uhat,vhat,hhat, &
        sqhl, sqhr
    !===========determine velocity from momentum===========================
    if (hl<tol) then
        hl=0.d0
        ul=0.d0
        vl=0.d0
    else
        ul=hul/hl
        vl=hvl/hl
    endif

    if (hr<tol) then
        hr=0.d0
        ur=0.d0
        vr=0.d0
    else
        ur=hur/hr
        vr=hvr/hr
    endif


    if (hl <= tol .and. hr <= tol) then
        s1 = 0.d0
        s2 = 0.d0
        s3 = 0.d0
    else
        !=====Determine some speeds necessary for the Jacobian=================
        sqhr = sqrt(hr)
        sqhl = sqrt(hl)
    !---------------------------------------
        vhat=(vr*sqhr)/(sqhr+sqhl) + (vl*sqhl)/(sqhr+sqhl)
        uhat=(ur*sqhr)/(sqhr+sqhl) + (ul*sqhl)/(sqhr+sqhl)
    !---------------------------------------
        hhat=(hr+hl)*0.5d0
    !---------------------------------------
        roe1=vhat-sqrt(g*hhat)
        roe3=vhat+sqrt(g*hhat)
    !---------------------------------------
        s1l=vl-sqrt(g*hl)
        s3r=vr+sqrt(g*hr)

        s1=min(roe1,s1l)
        s3=max(roe3,s3r)
        s2=0.5d0*(s1+s3)
    endif
end subroutine
