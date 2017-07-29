!======================================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,&
        ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!======================================================================
!
! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!     with topography:
!     #        h_t + (hu)_x + (hv)_y = 0                           #
!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.

    !  Note that the i'th Riemann problem has left state qr(:,i-1)
!     and right state ql(:,i)
    !  From the basic clawpack routines, this routine is called with
    !     ql = qr
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                           !
    !      # This Riemann solver is for the shallow water equations.            !
    !                                                                           !
    !       It allows the user to easily select a Riemann solver in             !
    !       riemannsolvers_geo.f. this routine initializes all the variables    !
    !       for the shallow water equations, accounting for wet dry boundary    !
    !       dry cells, wave speeds etc.                                         !
    !                                                                           !
    !           David George, Vancouver WA, Feb. 2009                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use geoclaw_module, only: g => grav, drytol => dry_tolerance
    use geoclaw_module, only: earth_radius, deg2rad
    use amr_module, only: mcapa
#ifdef USEPAPI
    use papi_module
#endif

    implicit none
    integer, parameter :: DP = kind(1.d0)

    !input
    integer maxm,meqn,maux,mwaves,mbc,mx,ixy

    real(kind=DP) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=DP) :: s(mwaves,1-mbc:maxm+mbc)
    real(kind=DP) :: ql(1-mbc:maxm+mbc, meqn)
    real(kind=DP) :: qr(1-mbc:maxm+mbc, meqn)
    real(kind=DP) :: auxl(1-mbc:maxm+mbc,maux)
    real(kind=DP) :: auxr(1-mbc:maxm+mbc,maux)
    real(kind=DP) :: apdq(meqn,1-mbc:maxm+mbc)
    real(kind=DP) :: amdq(meqn,1-mbc:maxm+mbc)

    !local only
    integer m,i,mw,maxiter,mu,nv
    real(kind=DP) :: wall(3)
    !real(kind=DP), dimension(3,3) :: fw
    !real(kind=DP) :: sw(3)
    real(kind=DP) :: sw1, sw2, sw3
    real(kind=DP) :: fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33

    real(kind=DP) :: hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
    real(kind=DP) :: bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    real(kind=DP) :: s1m,s2m
    real(kind=DP) :: hstar,hstartest,hstarHLL,sLtest,sRtest
    real(kind=DP) :: tw, dxdc
    real(kind=DP) :: sqghl, sqghr

    logical :: rare1,rare2
    ! Status variable for negative input
    logical :: negative_input = .false.

    interface
        subroutine solve_single_layer_rp(drytol, hL, hR, huL, huR, hvL, hvR, bL, bR, &
            fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33, sw1, sw2, sw3)
        !!! - dir$ attributes vector: uniform(drytol) :: solve_single_layer_rp
        !$OMP DECLARE SIMD(solve_single_layer_rp) UNIFORM(drytol)
            real(kind=8), intent(in) :: drytol
            real(kind=8), intent(inout) :: sw1, sw2, sw3
            real(kind=8), intent(inout) :: fw11, fw12, fw13, fw21, fw22, fw23,fw31, fw32, fw33
            real(kind=8), intent(inout) :: hL, hR, huL, huR, hvL, hvR, bL, bR
        end subroutine
    end interface

    !-----------------------Initializing-----------------------------------
    !set normal direction
    if (ixy.eq.1) then
        mu=2
        nv=3
    else
        mu=3
        nv=2
    endif

!    !zero (small) negative values if they exist
!    !DIR$ VECTOR ALIGNED
!    do i=2-mbc,mx+mbc
!        if (ql(i,1).lt.0.d0) then
!            ql(i,1)=0.d0
!            ql(i,2)=0.d0
!            ql(i,3)=0.d0
!            negative_input = .true.
!        endif
!    enddo
!
!    ! Inform of a bad riemann problem from the start
!    if (negative_input) then
!        write (*,*) 'Negative input for hl,hr!'
!    endif

    !----------------------------------------------------------------------
    !loop through Riemann problems at each grid cell
    !DIR$ VECTOR ALIGNED 
    !$OMP SIMD PRIVATE(hL,hR,huL,huR,hvL,hvR,bL,bR, &
    !$OMP fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33,sw1,sw2,sw3)
    do i=2-mbc,mx+mbc
        ! Riemann problem variables. Since the cells of the i-1th and ith
        ! Rimann problem overlap, the corresponding variables need to be copied
        ! before changing them (this is done when
        ! adding the f-waves and speeds.
        hL  = qr(i-1,1) 
        hR  = ql(i,1) 
        huL = qr(i-1,mu) 
        huR = ql(i,mu) 
        hvL = qr(i-1,nv) 
        hvR = ql(i,nv)
        bL  = auxr(i-1,1)
        bR  = auxl(i,1)

        !dir$ forceinline 
        call solve_single_layer_rp(drytol, hL, hR, huL, huR, hvL, hvR, bL, bR, &
            fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33, sw1, sw2, sw3)

        s(1,i) = sw1
        s(2,i) = sw2
        s(3,i) = sw3
        ! mw=1
        fwave(1, 1,i) = fw11
        fwave(mu,1,i) = fw21
        fwave(nv,1,i) = fw31
        ! mw=2
        fwave(1, 2,i) = fw12
        fwave(mu,2,i) = fw22
        fwave(nv,2,i) = fw32
        ! mw=3
        fwave(1, 3,i) = fw13
        fwave(mu,3,i) = fw23
        fwave(nv,3,i) = fw33
    enddo

    !==========Capacity for mapping from latitude longitude to physical space====
    if (mcapa > 0) then
        if (ixy == 1) then
            dxdc = earth_radius*deg2rad
            ! OMP SIMD GIVES SPEEDUP!!!
            !DIR$ VECTOR ALIGNED
            !$OMP SIMD COLLAPSE(1)
            do i=2-mbc,mx+mbc
                do mw=1,3
                    s(mw,i) = dxdc * s(mw,i)
                    fwave(1,mw,i) = dxdc*fwave(1,mw,i)
                    fwave(2,mw,i) = dxdc*fwave(2,mw,i)
                    fwave(3,mw,i) = dxdc*fwave(3,mw,i)
                enddo
            enddo
        else
            !DIR$ VECTOR ALIGNED 
            !$OMP SIMD PRIVATE(dxdc) COLLAPSE(1)
            do i=2-mbc,mx+mbc
                dxdc=earth_radius*cos(auxl(i,3))*deg2rad
                do mw=1,3
                    s(mw,i) = dxdc * s(mw,i)
                    fwave(1,mw,i)=dxdc*fwave(1,mw,i)
                    fwave(2,mw,i)=dxdc*fwave(2,mw,i)
                    fwave(3,mw,i)=dxdc*fwave(3,mw,i)
                enddo
            enddo
        endif
    endif
    !===============================================================================

    !DIR$ VECTOR ALIGNED 
    !$OMP SIMD COLLAPSE(1)
    do i=1-mbc,mx+mbc
        do m=1,3
            amdq(m,i) = 0.d0
        enddo
    enddo

    !DIR$ VECTOR ALIGNED 
    !$OMP SIMD COLLAPSE(1)
    do i=1-mbc,mx+mbc
        do m=1,3
            apdq(m,i) = 0.d0
        enddo
    enddo

    !============= compute fluctuations=============================================

    !DIR$ VECTOR ALIGNED 
    do i=2-mbc,mx+mbc
        do  mw=1,3
            if (s(mw,i) < 0.d0) then
                amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
            else if (s(mw,i) > 0.d0) then
                apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
            else
                amdq(1:3,i) = amdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
                apdq(1:3,i) = apdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
            endif
        enddo
    enddo

    return
end subroutine
