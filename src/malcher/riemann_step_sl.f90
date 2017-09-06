subroutine solve_single_layer_rp(drytol, hL, hR, huL, huR, hvL, hvR, bL, bR, & !fw, sw)
    fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33, sw1, sw2, sw3)
    !$OMP DECLARE SIMD(solve_single_layer_rp) UNIFORM(drytol)
    use geoclaw_module, only: g => grav
    implicit none

    ! Input
    real(kind=8), intent(in) :: drytol

    ! Output
    !real(kind=8), intent(in out) :: fw(3, 3), sw(3)
    real(kind=8), intent(inout) :: sw1, sw2, sw3
    real(kind=8), intent(inout) :: fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33

    ! Locals
    integer :: mw
    real(kind=8), intent(inout) :: hL, hR, huL, huR, hvL, hvR, bL, bR
    real(kind=8) :: uL, uR, vL, vR
    real(kind=8) :: phiL, phiR
    real(kind=8) :: hstar, hstartest, s1m, s2m, rare1, rare2, uhat, chat, sE1, sE2
    real(kind=8) :: sqhL, sqhR
    real(kind=8) :: wall1, wall2, wall3

    ! ========================================
    !  Begin Snipped Code From rpn2_geoclaw.f
    ! ========================================
    ! For completely dry states, do not skip problem (hinders
    ! vectorization), but rather solve artificial 0-valued problem.
    if (hL <= drytol .and. hR <= drytol) then
        hL = drytol
        hR = drytol
    endif
    !check for wet/dry boundary
    if (hR>drytol) then
        uR=huR/hR
        vR=hvR/hR
        !phiR = 0.5d0*g*hR**2 + huR**2/hR
        phiR = 0.5d0*g*(hR*hR) + (huR*huR)/hR
    else
        uR = 0.d0
        vR = 0.d0
        phiR = 0.d0
        hR = 0.d0
        huR = 0.d0
        hvR = 0.d0
    endif

    if (hL>drytol) then
        uL=huL/hL
        vL=hvL/hL
        !phiL = 0.5d0*g*hL**2 + huL**2/hL
        phiL = 0.5d0*g*(hL*hL) + (huL*huL)/hL
    else
        uL=0.d0
        vL=0.d0
        phiL = 0.d0
        hL=0.d0
        huL=0.d0
        hvL=0.d0
    endif

    wall1 = 1.d0
    wall2 = 1.d0
    wall3 = 1.d0
#if 1
     !if (hR<=drytol) then
    if (hR<drytol) then
    !dir$ forceinline
        call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,&
                 rare1,rare2,1,drytol,g)
        hstartest=max(hL,hstar)
        if (hstartest+bL<bR) then !right state should become ghost values that mirror left for wall problem
            ! bR=hstartest+bL
            wall2=0.d0
            wall3=0.d0
            hR=hL
            huR=-huL
            bR=bL
            phiR=phiL
            uR=-uL
            vR=vL
        elseif (hL+bL<bR) then
            bR=hL+bL
        endif
    !elseif (hL<=drytol) then ! right surface is lower than left topo
    elseif (hL<drytol) then ! right surface is lower than left topo
        !dir$ forceinline
        call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,&
                rare1,rare2,1,drytol,g)
        hstartest=max(hR,hstar)
        if (hstartest+bR<bL) then  !left state should become ghost values that mirror right
            ! bL=hstartest+bR
            wall1=0.d0
            wall2=0.d0
            hL=hR
            huL=-huR
            bL=bR
            phiL=phiR
            uL=-uR
            vL=vR
        elseif (hR+bR<bL) then
            bL=hR+bR
        endif
    endif
#endif
    sqhL = sqrt(g*hL)
    sqhR = sqrt(g*hR)
    !determine wave speeds

    if (sqhR+sqhL > drytol) then
        uhat= (sqhL*uL + sqhR*uR)/(sqhR+sqhL) ! Roe average
        chat= sqrt(g*0.5d0*(hR+hL)) ! Roe average
        sE1 = min(uL-sqhL,uhat-chat) ! Eindfeldt speed 1 wave
        sE2 = max(uR+sqhR,uhat+chat) ! Eindfeldt speed 2 wave
    else
        !uhat = 0.d0
        !chat = 0.d0
        sE1 = 0.d0
        sE2 = 0.d0
    endif

    ! Computations for min/max arguments in-place
    !sRoe1=uhat-chat ! Roe wave speed 1 wave
    !sRoe2=uhat+chat ! Roe wave speed 2 wave
    !sL=uL-sqhL ! 1 wave speed of left state
    !sR=uR+sqhR ! 2 wave speed of right state
    !sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
    !sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

    !--------------------end initializing...finally----------
    !solve Riemann problem.

    ! dir$ forceinline
    call riemann_fwave(3,3,hL,hR,huL,huR,hvL,hvR, &
            bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,& !sw,fw)
    sw1, sw2, sw3,fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33)

    ! eliminate ghost fluxes for wall
    sw1 = sw1 * wall1
    sw2 = sw2 * wall2
    sw3 = sw3 * wall3

    ! mw = 1
    fw11=fw11*wall1 
    fw21=fw21*wall1
    fw31=fw31*wall1
    ! mw = 2
    fw12=fw12*wall2 
    fw22=fw22*wall2
    fw32=fw32*wall2
    ! mw = 3
    fw13=fw13*wall3 
    fw23=fw23*wall3
    fw33=fw33*wall3
    !    do mw=1,3
    !        sw(mw)  =sw(mw)*wall(mw)
    !        fw(1,mw)=fw(1,mw)*wall(mw) 
    !        fw(2,mw)=fw(2,mw)*wall(mw)
!        fw(3,mw)=fw(3,mw)*wall(mw)
    !    enddo

    !    do mw=1,mwaves
    !        s(i,mw)=sw(mw)
    !        fwave(1,mw,i)=fw(1,mw)
    !        fwave(mu,mw,i)=fw(2,mw)
!        fwave(nv,mw,i)=fw(3,mw)
    !    enddo

end subroutine solve_single_layer_rp
