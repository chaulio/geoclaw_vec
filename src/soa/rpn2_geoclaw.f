! choose here how many Riemann Problems will
! be solved at once by the vectorized solver
#define _SOLVER_CHUNK_SIZE 2

c======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                 ql,qr,auxl,auxr,fwave,s,amdq,apdq)
c======================================================================
c
c Solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
c
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
      integer i,j,mw,mu,nv
      double precision fw(_SOLVER_CHUNK_SIZE,3,3)
      double precision sw(_SOLVER_CHUNK_SIZE,3)

      double precision hR(_SOLVER_CHUNK_SIZE),hL(_SOLVER_CHUNK_SIZE)
      double precision huR(_SOLVER_CHUNK_SIZE),huL(_SOLVER_CHUNK_SIZE)
      double precision hvR(_SOLVER_CHUNK_SIZE),hvL(_SOLVER_CHUNK_SIZE)
      double precision pL(_SOLVER_CHUNK_SIZE),pR(_SOLVER_CHUNK_SIZE)
      double precision bR(_SOLVER_CHUNK_SIZE),bL(_SOLVER_CHUNK_SIZE)
      double precision dxdc
      
      double precision sqrt_ghL(_SOLVER_CHUNK_SIZE), sqrt_ghR(_SOLVER_CHUNK_SIZE)

      logical rare1(_SOLVER_CHUNK_SIZE),rare2(_SOLVER_CHUNK_SIZE)
      
      
      call rpn2_start_timer()
      
      !set normal direction
      if (ixy.eq.1) then
         mu=2
         nv=3
      else
         mu=3
         nv=2
      endif
      
      !zero (small) negative values if they exist
      do i = 2-mbc, mx+mbc
         if (qr(i-1,1).lt.0.d0) then
               qr(i-1,1)=0.d0
               qr(i-1,2)=0.d0
               qr(i-1,3)=0.d0
         endif

         if (ql(i,1).lt.0.d0) then
               ql(i,1)=0.d0
               ql(i,2)=0.d0
               ql(i,3)=0.d0
         endif
      end do
      
      do i = 2-mbc, mx+mbc, _SOLVER_CHUNK_SIZE
        ! copy Riemann problem variables to arrays
        do j = 0, _SOLVER_CHUNK_SIZE - 1 ! j starts on 0 to make things easier, but actual arrays start on 1!
            if (i + j <= mx+mbc) then
                hL(j+1) = qr(i+j-1,1)
                hR(j+1) = ql(i+j,1) 
                huL(j+1) = qr(i+j-1,mu) 
                huR(j+1) = ql(i+j,mu) 
                bL(j+1) = auxr(i+j-1,1)
                bR(j+1) = auxl(i+j,1)
                
                ! In case there is no pressure forcing
                pL(j+1) = 0.d0
                pR(j+1) = 0.d0
                if (pressure_forcing) then
                    pL(j+1) = auxr(i+j-1,pressure_index)
                    pR(j+1) = auxl(i+j, pressure_index)
                end if

                hvL(j+1)=qr(i+j-1,nv) 
                hvR(j+1)=ql(i+j,nv)
!         print *, "========================="
!         print *, "Solving RP number", i+j
!         print *, "HL", hl(j+1)
!         print *, "HR", hR(j+1)
!         print *, "HUL", HUL(j+1)
!         print *, "HUR", HUR(j+1)
!         print *, "HVL", HVL(j+1)
!         print *, "HVR", HVR(j+1)
!         print *, "BL", BL(j+1)
!         print *, "BR", BR(j+1)

            else
                ! skip!
                !print *, "skipping", j+1
                hL(j+1) = 0.d0
                hR(j+1) = 0.d0
                huL(j+1) = 0.d0
                huR(j+1) = 0.d0
                hvL(j+1) = 0.d0
                hvR(j+1) = 0.d0
            end if

        end do
        
       ! initialize output arrays
       fw(:,:,:) = 0.d0
       sw(:,:) = 0.d0
       
       ! solve problems in the arrays
       call rpn2_vec(hL, hR, huL, huR, hvL, hvR, bL, bR, fw, sw)
      
       ! copy waves computed by solver
       do j = 0, _SOLVER_CHUNK_SIZE - 1 ! j starts on 0 to make things easier, but actual arrays start on 1!
         if (i + j <= mx+mbc) then !avoid inexistent indices
            !print *, 2-mbc, mx+mbc, "Solving: " , i + j, j
!          print *, "========================="
!          print *, "Solution for RP number", i+j 
!          print *, "fw", fw(j+1,:,:)
!          print *, "sw", sw(j+1,:)
!           print *, "working on " , i+j, "of", mx+mbc
           do mw=1,mwaves 
             s(mw,i+j)=sw(j+1,mw)
             fwave(1,mw,i+j)=fw(j+1,1,mw)
             fwave(mu,mw,i+j)=fw(j+1,2,mw)
             fwave(nv,mw,i+j)=fw(j+1,3,mw)
           enddo
         end if
       end do
      end do
      




     
     
c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxl(i,3))*deg2rad
          endif

          do mw=1,mwaves
c             if (s(mw,i) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
c                endif
               s(mw,i)=dxdc*s(mw,i)
               fwave(1,mw,i)=dxdc*fwave(1,mw,i)
               fwave(2,mw,i)=dxdc*fwave(2,mw,i)
               fwave(3,mw,i)=dxdc*fwave(3,mw,i)
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         amdq(1:3,:) = 0.d0
         apdq(1:3,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
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
!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo

      call rpn2_stop_timer()
      
      end subroutine
      
      
      
      subroutine rpn2_vec(hL, hR, huL, huR, hvL, hvR, bL, bR, fw, sw)

      use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      
      implicit none

      !input
      double precision hL(_SOLVER_CHUNK_SIZE), hR(_SOLVER_CHUNK_SIZE)
      double precision huL(_SOLVER_CHUNK_SIZE), huR(_SOLVER_CHUNK_SIZE)
      double precision hvL(_SOLVER_CHUNK_SIZE), hvR(_SOLVER_CHUNK_SIZE)
      double precision bL(_SOLVER_CHUNK_SIZE), bR(_SOLVER_CHUNK_SIZE)
      double precision fw(_SOLVER_CHUNK_SIZE,3,3), sw(_SOLVER_CHUNK_SIZE,3)
      
      !local
      integer i, maxiter, mw
      double precision uL(_SOLVER_CHUNK_SIZE), uR(_SOLVER_CHUNK_SIZE)
      double precision vL(_SOLVER_CHUNK_SIZE), vR(_SOLVER_CHUNK_SIZE) 
      double precision phiL(_SOLVER_CHUNK_SIZE), phiR(_SOLVER_CHUNK_SIZE)
      double precision pL(_SOLVER_CHUNK_SIZE), pR(_SOLVER_CHUNK_SIZE)

      double precision wall(_SOLVER_CHUNK_SIZE,3)
      double precision hstar(_SOLVER_CHUNK_SIZE), hstartest(_SOLVER_CHUNK_SIZE)
      double precision s1m(_SOLVER_CHUNK_SIZE), s2m(_SOLVER_CHUNK_SIZE)
      logical rare1(_SOLVER_CHUNK_SIZE), rare2(_SOLVER_CHUNK_SIZE)
      
      double precision sqrt_ghL(_SOLVER_CHUNK_SIZE), sqrt_ghR(_SOLVER_CHUNK_SIZE)
      double precision sL(_SOLVER_CHUNK_SIZE), sR(_SOLVER_CHUNK_SIZE)
      double precision uhat(_SOLVER_CHUNK_SIZE), chat(_SOLVER_CHUNK_SIZE)
      double precision sRoe1(_SOLVER_CHUNK_SIZE), sRoe2(_SOLVER_CHUNK_SIZE)
      double precision sE1(_SOLVER_CHUNK_SIZE), sE2(_SOLVER_CHUNK_SIZE)
      
      
         !check for wet/dry boundary
         where (hR.gt.drytol)
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         elsewhere
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         end where

         where (hL.gt.drytol) 
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         elsewhere
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         end where
         
         wall(:,1) = 1.d0
         wall(:,2) = 1.d0
         wall(:,3) = 1.d0
         
         do i=1, _SOLVER_CHUNK_SIZE !TODO: vectorize this part with riemanntype subroutine...
            if (hR(i).le.drytol) then
                call riemanntype(hL(i),hL(i),uL(i),-uL(i),hstar(i),
     &                  s1m(i),s2m(i),rare1(i),rare2(i),1,drytol,g)
                hstartest(i)=max(hL(i),hstar(i))
                if (hstartest(i)+bL(i).lt.bR(i)) then !right state should become ghost values that mirror left for wall problem
                    !bR=hstartest+bL
                    wall(i,2)=0.d0
                    wall(i,3)=0.d0
                    hR(i)=hL(i)
                    huR(i)=-huL(i)
                    bR(i)=bL(i)
                    phiR(i)=phiL(i)
                    uR(i)=-uL(i)
                    vR(i)=vL(i)
                elseif (hL(i)+bL(i).lt.bR(i)) then
                    bR(i)=hL(i)+bL(i)
                endif
            elseif (hL(i).le.drytol) then ! right surface is lower than left topo
                call riemanntype(hR(i),hR(i),-uR(i),uR(i),hstar(i),
     &                    s1m(i),s2m(i),rare1(i),rare2(i),1,drytol,g)
                hstartest(i)=max(hR(i),hstar(i))
                if (hstartest(i)+bR(i).lt.bL(i)) then  !left state should become ghost values that mirror right
                    !bL(i)=hstartest(i)+bR(i)
                    wall(i,1)=0.d0
                    wall(i,2)=0.d0
                    hL(i)=hR(i)
                    huL(i)=-huR(i)
                    bL(i)=bR(i)
                    phiL(i)=phiR(i)
                    uL(i)=-uR(i)
                    vL(i)=vR(i)
                elseif (hR(i)+bR(i).lt.bL(i)) then
                    bL(i)=hR(i)+bR(i)
                endif
            endif
         end do
         
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
         
         ! using vectorized version of this solver!
           call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,
     &        huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,
     &               pL,pR,sE1,sE2,drytol,g,rho,sw,fw)

C          call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
C      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
C      &     rho,sw,fw)

C          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
C      &      bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,
C      &      fw)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(:,mw)=sw(:,mw)*wall(:,mw)

               fw(:,1,mw)=fw(:,1,mw)*wall(:,mw) 
               fw(:,2,mw)=fw(:,2,mw)*wall(:,mw)
               fw(:,3,mw)=fw(:,3,mw)*wall(:,mw)
         enddo
      end subroutine
      
      
      subroutine rpn2_orig(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                 ql,qr,auxl,auxr,fwave,s,amdq,apdq)

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
      double precision fw(3,3)
      double precision sw(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL,pL,pR
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL
      double precision dxdc
      
      double precision sqrt_ghL, sqrt_ghR

      logical rare1,rare2
     

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
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
              s(mw,i)=0.d0
                 fwave(1,mw,i)=0.d0
                 fwave(2,mw,i)=0.d0
                 fwave(3,mw,i)=0.d0
         enddo

         !zero (small) negative values if they exist
         if (qr(i-1,1).lt.0.d0) then
               qr(i-1,1)=0.d0
               qr(i-1,2)=0.d0
               qr(i-1,3)=0.d0
         endif

         if (ql(i,1).lt.0.d0) then
               ql(i,1)=0.d0
               ql(i,2)=0.d0
               ql(i,3)=0.d0
         endif

         !skip problem if in a completely dry area
         if (qr(i-1,1) <= drytol .and. ql(i,1) <= drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(i-1,1) 
         hR = ql(i,1) 
         huL = qr(i-1,mu) 
         huR = ql(i,mu) 
         bL = auxr(i-1,1)
         bR = auxl(i,1)
         if (pressure_forcing) then
             pL = auxr(i-1,pressure_index)
             pR = auxl(i, pressure_index)
         end if

         hvL=qr(i-1,nv) 
         hvR=ql(i,nv)
         
#define solver_starts_here

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
         

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
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
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
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

         call riemann_aug_JCP_orig(maxiter,3,3,hL,hR,huL,
     &        huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,
     &                        sE1,sE2,drytol,g,rho,sw,fw)

C          call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
C      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
C      &     rho,sw,fw)

C          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
C      &      bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,
C      &      fw)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)

               fw(1,mw)=fw(1,mw)*wall(mw) 
               fw(2,mw)=fw(2,mw)*wall(mw)
               fw(3,mw)=fw(3,mw)*wall(mw)
         enddo
         
#define solver_ends_here


         do mw=1,mwaves
            s(mw,i)=sw(mw)
            fwave(1,mw,i)=fw(1,mw)
            fwave(mu,mw,i)=fw(2,mw)
            fwave(nv,mw,i)=fw(3,mw)
!            write(51,515) sw(mw),fw(1,mw),fw(2,mw),fw(3,mw)
!515         format("++sw",4e25.16)
         enddo

 30      continue
      enddo

c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxl(i,3))*deg2rad
          endif

          do mw=1,mwaves
c             if (s(mw,i) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
c                endif
               s(mw,i)=dxdc*s(mw,i)
               fwave(1,mw,i)=dxdc*fwave(1,mw,i)
               fwave(2,mw,i)=dxdc*fwave(2,mw,i)
               fwave(3,mw,i)=dxdc*fwave(3,mw,i)
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         amdq(1:3,:) = 0.d0
         apdq(1:3,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
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
!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo


      return
      end subroutine


      
! this subroutine has been moved from geoclaw_riemann_utils.f to here
c-----------------------------------------------------------------------
      subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,
     &   hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,
     &   g,rho,sw,fw)

      ! solve shallow water equations given single left and right states
      ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

      ! To use the original solver call with maxiter=1.

      ! This solver allows iteration when maxiter > 1. The iteration seems to help with
      ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
      ! due to loss of hyperbolicity.

      implicit none

      !input
      integer meqn,mwaves,maxiter
      double precision fw(_SOLVER_CHUNK_SIZE,meqn,mwaves)
      double precision sw(_SOLVER_CHUNK_SIZE,mwaves)
      double precision hL(_SOLVER_CHUNK_SIZE),hR(_SOLVER_CHUNK_SIZE)
      double precision huL(_SOLVER_CHUNK_SIZE), huR(_SOLVER_CHUNK_SIZE)
      double precision bL(_SOLVER_CHUNK_SIZE), bR(_SOLVER_CHUNK_SIZE)
      double precision uL(_SOLVER_CHUNK_SIZE), uR(_SOLVER_CHUNK_SIZE)
      double precision phiL(_SOLVER_CHUNK_SIZE), phiR(_SOLVER_CHUNK_SIZE) 
      double precision sE1(_SOLVER_CHUNK_SIZE), sE2(_SOLVER_CHUNK_SIZE)
      double precision hvL(_SOLVER_CHUNK_SIZE), hvR(_SOLVER_CHUNK_SIZE)
      double precision vL(_SOLVER_CHUNK_SIZE), vR(_SOLVER_CHUNK_SIZE)
      double precision pL(_SOLVER_CHUNK_SIZE), pR(_SOLVER_CHUNK_SIZE)
      double precision drytol,g,rho


      !local
      integer m,mw,k,iter,i
      double precision A(_SOLVER_CHUNK_SIZE,3,3)
      double precision r(_SOLVER_CHUNK_SIZE,3,3)
      double precision lambda(_SOLVER_CHUNK_SIZE,3)
      double precision del(_SOLVER_CHUNK_SIZE,3)
      double precision beta(_SOLVER_CHUNK_SIZE,3)

      double precision delh(_SOLVER_CHUNK_SIZE),delhu(_SOLVER_CHUNK_SIZE)
      double precision delphi(_SOLVER_CHUNK_SIZE),delb(_SOLVER_CHUNK_SIZE)
      double precision delnorm(_SOLVER_CHUNK_SIZE)
      double precision rare1st(_SOLVER_CHUNK_SIZE),rare2st(_SOLVER_CHUNK_SIZE)
      double precision sdelta(_SOLVER_CHUNK_SIZE),raremin(_SOLVER_CHUNK_SIZE)
      double precision raremax(_SOLVER_CHUNK_SIZE)
      double precision criticaltol(_SOLVER_CHUNK_SIZE),convergencetol(_SOLVER_CHUNK_SIZE)
      double precision raretol(_SOLVER_CHUNK_SIZE)
      double precision criticaltol_2(_SOLVER_CHUNK_SIZE),hustar_interface(_SOLVER_CHUNK_SIZE)
      double precision s1s2bar(_SOLVER_CHUNK_SIZE),s1s2tilde(_SOLVER_CHUNK_SIZE),hbar(_SOLVER_CHUNK_SIZE)
      double precision hLstar(_SOLVER_CHUNK_SIZE),hRstar(_SOLVER_CHUNK_SIZE),hustar(_SOLVER_CHUNK_SIZE)
      double precision huRstar(_SOLVER_CHUNK_SIZE),huLstar(_SOLVER_CHUNK_SIZE),uRstar(_SOLVER_CHUNK_SIZE)
      double precision uLstar(_SOLVER_CHUNK_SIZE),hstarHLL(_SOLVER_CHUNK_SIZE)
      double precision deldelh(_SOLVER_CHUNK_SIZE),deldelphi(_SOLVER_CHUNK_SIZE),delP(_SOLVER_CHUNK_SIZE)
      double precision s1m(_SOLVER_CHUNK_SIZE),s2m(_SOLVER_CHUNK_SIZE),hm(_SOLVER_CHUNK_SIZE)
      double precision det1(_SOLVER_CHUNK_SIZE),det2(_SOLVER_CHUNK_SIZE),det3(_SOLVER_CHUNK_SIZE)
      double precision determinant(_SOLVER_CHUNK_SIZE)

      logical rare1(_SOLVER_CHUNK_SIZE),rare2(_SOLVER_CHUNK_SIZE),rarecorrector(_SOLVER_CHUNK_SIZE)
      logical rarecorrectortest(_SOLVER_CHUNK_SIZE), sonic(_SOLVER_CHUNK_SIZE)

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL
      delP = pR - pL
      delnorm = delh**2 + delphi**2

      do i=1,_SOLVER_CHUNK_SIZE
        call riemanntype(hL(i),hR(i),uL(i),uR(i),
     &          hm(i),s1m(i),s2m(i),rare1(i),rare2(i),1,drytol,g)
      end do

      lambda(:,1)= min(sE1,s2m) !Modified Einfeldt speed
      lambda(:,3)= max(sE2,s1m) !Modified Eindfeldt speed
      sE1=lambda(:,1)
      sE2=lambda(:,3)
      lambda(:,2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

      
      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

c     !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.false.
      rarecorrector=.false.
      where (rarecorrectortest) 
         sdelta=lambda(:,3)-lambda(:,1)
         raremin = 0.5d0
         raremax = 0.9d0
         where (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         where (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         where (rare1.or.rare2) 
            !see which rarefaction is larger
            rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
            rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
            where (max(rare1st,rare2st).gt.raremin*sdelta.and.
     &         max(rare1st,rare2st).lt.raremax*sdelta) 
                  rarecorrector=.true.
               where (rare1st.gt.rare2st) 
                  lambda(:,2)=s1m
               elsewhere (rare2st.gt.rare1st) 
                  lambda(:,2)=s2m
               elsewhere
                  lambda(:,2)=0.5d0*(s1m+s2m)
               end where
            end where
         end where
         where (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      end where

c     ## Is this correct 2-wave when rarecorrector == .true. ??
      do mw=1,mwaves
         r(:,1,mw)=1.d0
         r(:,2,mw)=lambda(:,mw)
         r(:,3,mw)=(lambda(:,mw))**2
      enddo
      where (.not.rarecorrector)
         lambda(:,2) = 0.5d0*(lambda(:,1)+lambda(:,3))
c         lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
         r(:,1,2)=0.d0
         r(:,2,2)=0.d0
         r(:,3,2)=1.d0
      end where
c     !---------------------------------------------------


c     !determine the steady state wave -------------------
      !criticaltol = 1.d-6
      ! MODIFIED:
      criticaltol = max(drytol*g, 1d-6)
      criticaltol_2 = sqrt(criticaltol)
      deldelh = -delb
      deldelphi = -0.5d0 * (hR + hL) * (g * delb + delp / rho)

c     !determine a few quanitites needed for steady state wave if iterated
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,maxiter
         !determine steady state wave (this will be subtracted from the delta vectors)
         where (min(hLstar,hRstar).lt.drytol.and.rarecorrector)
            rarecorrector=.false.
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            lambda(:,2) = 0.5d0*(lambda(:,1)+lambda(:,3))
c           lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
            r(:,1,2)=0.d0
            r(:,2,2)=0.d0
            r(:,3,2)=1.d0
         end where

         hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
         s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
         s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

c        !find if sonic problem
         ! MODIFIED from 5.3.1 version
         ! TODO: try the one-line-solution for this - see Samoa/SWE/SWE_Patch_solvers.f90 lines 427 and 428
         sonic = .false.
         where (abs(s1s2bar) <= criticaltol) 
            sonic = .true.
         elsewhere (s1s2bar*s1s2tilde <= criticaltol**2)
            sonic = .true.
         elsewhere (s1s2bar*sE1*sE2 <= criticaltol**2) 
            sonic = .true.
         elsewhere (min(abs(sE1),abs(sE2)) < criticaltol_2) 
            sonic = .true.
         elsewhere (sE1 <  criticaltol_2 .and. s1m > -criticaltol_2) 
            sonic = .true.
         elsewhere (sE2 > -criticaltol_2 .and. s2m <  criticaltol_2) 
            sonic = .true.
         elsewhere ((uL+dsqrt(g*hL))*(uR+dsqrt(g*hR)) < 0.d0) 
            sonic = .true.
         elsewhere ((uL- dsqrt(g*hL))*(uR- dsqrt(g*hR)) < 0.d0) 
            sonic = .true.
         end where

c        !find jump in h, deldelh
         where (sonic)
            deldelh =  -delb
         elsewhere
            deldelh = delb*g*hbar/s1s2bar
         end where
c        !find bounds in case of critical state resonance, or negative states
         where (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) 
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
         elsewhere (sE1.ge.criticaltol)
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
            deldelh = max(deldelh,-hL)
         elsewhere (sE2.le.-criticaltol) 
            deldelh = min(deldelh,hR)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
         end where

c        !find jump in phi, deldelphi
         where (sonic) 
            deldelphi = -g*hbar*delb
         elsewhere
            deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
         end where
c        !find bounds in case of critical state resonance, or negative states
         deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
         deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))
         deldelphi = deldelphi - hbar * delp / rho

         del(:,1)=delh-deldelh
         del(:,2)=delhu
         del(:,3)=delphi-deldelphi

c        !Determine determinant of eigenvector matrix========
         det1=r(:,1,1)*(r(:,2,2)*r(:,3,3)-r(:,2,3)*r(:,3,2))
         det2=r(:,1,2)*(r(:,2,1)*r(:,3,3)-r(:,2,3)*r(:,3,1))
         det3=r(:,1,3)*(r(:,2,1)*r(:,3,2)-r(:,2,2)*r(:,3,1))
         determinant=det1-det2+det3

c        !solve for beta(k) using Cramers Rule=================
         do k=1,3
            do mw=1,3
                  A(:,1,mw)=r(:,1,mw)
                  A(:,2,mw)=r(:,2,mw)
                  A(:,3,mw)=r(:,3,mw)
            enddo
            A(:,1,k)=del(:,1)
            A(:,2,k)=del(:,2)
            A(:,3,k)=del(:,3)
            det1=A(:,1,1)*(A(:,2,2)*A(:,3,3)-A(:,2,3)*A(:,3,2))
            det2=A(:,1,2)*(A(:,2,1)*A(:,3,3)-A(:,2,3)*A(:,3,1))
            det3=A(:,1,3)*(A(:,2,1)*A(:,3,2)-A(:,2,2)*A(:,3,1))
            beta(:,k)=(det1-det2+det3)/determinant
         enddo

         !exit if things aren't changing --> not anymore, with vectorization we can't leave early!
         !if (abs(del(:,1)**2+del(:,3)**2-delnorm).lt.convergencetol) exit
         
         delnorm = del(:,1)**2+del(:,3)**2
         !find new states qLstar and qRstar on either side of interface
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar
         do mw=1,mwaves
            where (lambda(:,mw).lt.0.d0) 
               hLstar= hLstar + beta(:,mw)*r(:,1,mw)
               huLstar= huLstar + beta(:,mw)*r(:,2,mw)
            end where
         enddo
         do mw=mwaves,1,-1
            where (lambda(:,mw).gt.0.d0)
               hRstar= hRstar - beta(:,mw)*r(:,1,mw)
               huRstar= huRstar - beta(:,mw)*r(:,2,mw)
            end where
         enddo

         where (hLstar.gt.drytol)
            uLstar=huLstar/hLstar
         elsewhere
            hLstar=max(hLstar,0.d0)
            uLstar=0.d0
         end where
         where (hRstar.gt.drytol)
            uRstar=huRstar/hRstar
         elsewhere
            hRstar=max(hRstar,0.d0)
            uRstar=0.d0
         end where

      enddo ! end iteration on Riemann problem

      do mw=1,mwaves
         sw(:,mw)=lambda(:,mw)
         fw(:,1,mw)=beta(:,mw)*r(:,2,mw)
         fw(:,2,mw)=beta(:,mw)*r(:,3,mw)
         fw(:,3,mw)=beta(:,mw)*r(:,2,mw)
      enddo
      !find transverse components (ie huv jumps).
      ! MODIFIED from 5.3.1 version
      fw(:,3,1)=fw(:,3,1)*vL
      fw(:,3,3)=fw(:,3,3)*vR
      fw(:,3,2)= 0.d0
 
      hustar_interface = huL + fw(:,1,1)   ! = huR - fw(:,1,3)
      where (hustar_interface <= 0.0d0)
          fw(:,3,1) = fw(:,3,1) + (hR*uR*vR - hL*uL*vL - 
     &                        fw(:,3,1)- fw(:,3,3))
      else where
          fw(:,3,3) = fw(:,3,3) + (hR*uR*vR - hL*uL*vL - 
     &                        fw(:,3,1)- fw(:,3,3))
      end where

      
      ! only necessary in vectorized version, because early exit was not possible:
      ! zero the waves for DryDry cases!
      where (hL < drytol .and. hR < drytol)
        fw(:,1,1) = 0.0d0
        fw(:,1,2) = 0.0d0
        fw(:,1,3) = 0.0d0
        fw(:,2,1) = 0.0d0
        fw(:,2,2) = 0.0d0
        fw(:,2,3) = 0.0d0
        fw(:,3,1) = 0.0d0
        fw(:,3,2) = 0.0d0
        fw(:,3,3) = 0.0d0
        sw(:,1) = 0.0d0
        sw(:,2) = 0.0d0
        sw(:,3) = 0.0d0
      end where
      

      return

      end !subroutine riemann_aug_JCP-------------------------------------------------