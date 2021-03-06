!======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx, &
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

!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
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

      use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      use storm_module, only: pressure_forcing, pressure_index
      use statistics

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(1-mbc:maxm+mbc,maux)
      double precision  auxr(1-mbc:maxm+mbc,maux)

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      !double precision fw(3,3)
      !double precision sw(3)
      double precision fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33
      double precision sw1, sw2, sw3

      double precision hR,hL,huR,huL,hvR,hvL,pL,pR
      double precision bR,bL

      double precision tw,dxdc
      
      double precision sqrt_ghL, sqrt_ghR

      logical rare1,rare2
      
      interface
        subroutine solve_rpn2(hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR,sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33) 
        !$OMP DECLARE SIMD(solve_rpn2) 
            real(kind=8), intent(inout) :: sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33
            real(kind=8), intent(inout) :: hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR
        end subroutine
      end interface
      
      !DIR$ ASSUME_ALIGNED ql:64, qr:64, auxl:64, auxr:64, fwave:64, s:64
    
    
      
      call rpn2_start_timer()

      ! In case there is no pressure forcing
      pL = 0.d0
      pR = 0.d0
      
      !set normal direction
      if (ixy.eq.1) then
         mu=2
         nv=3
      else
         mu=3
         nv=2
      endif

      !loop through Riemann problems at each grid cell
#if !defined(__MIC__)
      ! for some reason this causes runtime error on MICs,
      ! even though the arrays are really aligned
      !DIR$ VECTOR ALIGNED 
#endif
#if defined(_OPENMP)
      !$OMP SIMD PRIVATE(hL,hR,huL,huR,hvL,hvR,bL,bR,pL,pR, &
      !$OMP& fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33,sw1,sw2,sw3)
#else
      !DIR$ SIMD
#endif
      do i=2-mbc,mx+mbc

         !Riemann problem variables
         hL = qr(i-1,1) 
         hR = ql(i,1) 
         huL = qr(i-1,mu) 
         huR = ql(i,mu) 
         bL = auxr(i-1,1)
         bR = auxl(i,1)
         pL = 0.d0
         pR = 0.d0
         if (pressure_forcing) then
             pL = auxr(i-1,pressure_index)
             pR = auxl(i, pressure_index)
         end if

         hvL=qr(i-1,nv) 
         hvR=ql(i,nv)
         
         sw1 = 0.0d0
         sw2 = 0.0d0
         sw3 = 0.0d0
         fw11 = 0.0d0
         fw21 = 0.0d0
         fw31 = 0.0d0
         fw12 = 0.0d0
         fw22 = 0.0d0
         fw32 = 0.0d0
         fw13 = 0.0d0
         fw23 = 0.0d0
         fw33 = 0.0d0   
         
         !DIR$ FORCEINLINE
         call solve_rpn(hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR, sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33) 
        

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

              do i=2-mbc,mx+mbc
                  do mw=1,3
                      s(mw,i) = dxdc * s(mw,i)
                      fwave(1,mw,i) = dxdc*fwave(1,mw,i)
                      fwave(2,mw,i) = dxdc*fwave(2,mw,i)
                      fwave(3,mw,i) = dxdc*fwave(3,mw,i)
                  enddo
              enddo
          else

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

      do i=1-mbc,mx+mbc
          do m=1,3
              amdq(m,i) = 0.d0
          enddo
      enddo


      do i=1-mbc,mx+mbc
          do m=1,3
              apdq(m,i) = 0.d0
         enddo
      enddo 

      !============= compute fluctuations=============================================


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
     
      call rpn2_stop_timer()
      return
      end subroutine
      
      
      
subroutine solve_rpn(hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR, sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33)      
      !$OMP DECLARE SIMD(solve_rpn) 
      use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      use storm_module, only: pressure_forcing, pressure_index
      use statistics

      implicit none

      !input
      double precision hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR
      
      !output
      !double precision sw(3), fw(3,3)
      double precision sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33

      !local only
      integer m,i,mw,maxiter
      double precision wall1, wall2, wall3
      double precision uR,uL,vR,vL,phiR,phiL
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc      
      double precision sqrt_ghL, sqrt_ghR
      logical rare1,rare2
         
         ! For completely dry states, do not skip problem (hinders
         ! vectorization), but rather solve artificial 0-valued problem.
         if (hL < drytol .and. hR < drytol) then
            hL = 0
            hR = 0
         endif

         !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         endif

         wall1 = 1.d0
         wall2 = 1.d0
         wall3 = 1.d0
         if (hR.le.drytol) then
            !DIR$ FORCEINLINE
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                                       rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
               !bR=hstartest+bL
               wall2=0.d0
               wall3=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            !DIR$ FORCEINLINE
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                                       rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
              !bL=hstartest+bR
               wall1=0.d0
               wall2=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif
         
         ! pre-compute square roots
         sqrt_ghL = sqrt(g*hL)
         sqrt_ghR = sqrt(g*hR)

         !determine wave speeds
         sL=uL-sqrt_ghL ! 1 wave speed of left state
         sR=uR+sqrt_ghR ! 2 wave speed of right state

         uhat=(sqrt_ghL*uL + sqrt_ghR*uR)/(sqrt_ghR+sqrt_ghL) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         !DIR$ FORCEINLINE
         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
             huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2, &
             drytol,g,rho,sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33)

!          call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
!      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
!      &     rho,sw,fw)

!          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
!      &      bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,
!      &      fw)



          ! For completely dry states, waves should be always zero!
          ! (these problems could have been skipped, but to allow 
          ! vectorization they weren't)
          if (hL < drytol .and. hR < drytol) then
             sw1 = 0.0d0
             sw2 = 0.0d0
             sw3 = 0.0d0
             fw11 = 0.0d0
             fw21 = 0.0d0
             fw31 = 0.0d0
             fw12 = 0.0d0
             fw22 = 0.0d0
             fw32 = 0.0d0
             fw13 = 0.0d0
             fw23 = 0.0d0
             fw33 = 0.0d0            
          endif


          !eliminate ghost fluxes for wall
!           do mw=1,3
!              sw(mw)=sw(mw)*wall(mw)
!  
!                 fw(1,mw)=fw(1,mw)*wall(mw) 
!                 fw(2,mw)=fw(2,mw)*wall(mw)
!                 fw(3,mw)=fw(3,mw)*wall(mw)
!           enddo

            sw1 = sw1 * wall1
            sw2 = sw2 * wall2
            sw3 = sw3 * wall3
            !mw=1
            fw11 = fw11 * wall1
            fw21 = fw21 * wall1
            fw31 = fw31 * wall1
            !mw=2
            fw12 = fw12 * wall2
            fw22 = fw22 * wall2
            fw32 = fw32 * wall2
            !mw=3
            fw13 = fw13 * wall3
            fw23 = fw23 * wall3
            fw33 = fw33 * wall3

end subroutine
