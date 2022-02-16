# include "Flash.h"
# include "constants.h"
# include "DGFV.h"

module dgfv_kernels_mod

use Hydro_data, only: N_NODES
use Hydro_data, only: NXE,NYE,NZE
use Hydro_data, only: NXN,NYN,NZN

use Hydro_data, only: NXE_LO,NXE_HI
use Hydro_data, only: NYE_LO,NYE_HI
use Hydro_data, only: NZE_LO,NZE_HI

use Hydro_data, only: NXN_LO,NXN_HI
use Hydro_data, only: NYN_LO,NYN_HI
use Hydro_data, only: NZN_LO,NZN_HI

use Hydro_data, only: NXB_LO,NXB_HI
use Hydro_data, only: NYB_LO,NYB_HI
use Hydro_data, only: NZB_LO,NZB_HI

use Hydro_data, only: NOR,SOU,WES,EAS,FRO,BAC

use Hydro_data, only: NPRIM_VARS
use Hydro_data, only: NCONS_VARS

use Hydro_data, only: DENS_PRIM
use Hydro_data, only: VELX_PRIM
use Hydro_data, only: VELY_PRIM
use Hydro_data, only: VELZ_PRIM
use Hydro_data, only: ENER_PRIM
use Hydro_data, only: MAGX_PRIM
use Hydro_data, only: MAGY_PRIM
use Hydro_data, only: MAGZ_PRIM
use Hydro_data, only: GLMP_PRIM
use Hydro_data, only: PRES_PRIM
use Hydro_data, only: GAMC_PRIM

use Hydro_data, only: DENS_CONS
use Hydro_data, only: MOMX_CONS
use Hydro_data, only: MOMY_CONS
use Hydro_data, only: MOMZ_CONS
use Hydro_data, only: ENER_CONS
use Hydro_data, only: MAGX_CONS
use Hydro_data, only: MAGY_CONS
use Hydro_data, only: MAGZ_CONS
use Hydro_data, only: GLMP_CONS
use Hydro_data, only: PRES_CONS
use Hydro_data, only: GAMC_CONS

use Hydro_data, only: SPEC_PRIM_BEGIN
use Hydro_data, only: SPEC_PRIM_END

use Hydro_data, only: MSCL_PRIM_BEGIN
use Hydro_data, only: MSCL_PRIM_END

use Hydro_data, only: SPEC_CONS_BEGIN
use Hydro_data, only: SPEC_CONS_END

use Hydro_data, only: MSCL_CONS_BEGIN
use Hydro_data, only: MSCL_CONS_END

contains

!! ================================================================== !!

pure function check_positivity(state) result(ok)

    use dgfv_fluxes_mod, only: isvalid_cons

    implicit none

    real, intent(in) :: state(NCONS_VARS,NXB,NYB,NZB)
    logical          :: ok

    integer :: i,j,k

    ok = .true.
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        ok = ok .and. isvalid_cons(state(:,i,j,k))
    end do; end do; end do;

end function

subroutine enforce_positivity(state)

    use dgfv_fluxes_mod, only: isvalid_cons
    use dgfv_fluxes_mod, only: cons2prim
    use Driver_interface, only: Driver_abortFlash
    !use Grid_data, ONLY: gr_smallrho

    implicit none

    real, intent(inout) :: state(NCONS_VARS,NXB,NYB,NZB)

    logical :: ok

    real :: squeeze
    real :: means(NCONS_VARS)
    real :: highs(NCONS_VARS,NXB,NYB,NZB)

    real, parameter :: ws = 1.0/REAL(NXB*NYB*NZB)

    integer :: i,j,k,h

    ok = .true.
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        ok = ok .and. isvalid_cons(state(:,i,j,k))
    end do; end do; end do;

    if (ok) return

    do h = 1,NCONS_VARS
        means(h) = ws*sum(state(h,:,:,:))
    end do

    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        do h = 1,NCONS_VARS
            highs(h,i,j,k) = state(h,i,j,k) - means(h)
        end do
    end do; end do; end do;

    squeeze = 1.0
    do while (squeeze > 0.0)

        squeeze = max(0.0,squeeze - 0.1)

        do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
            do h = 1,NCONS_VARS
                state(h,i,j,k) = means(h) + squeeze*highs(h,i,j,k)
            end do
        end do; end do; end do;

        ok = .true.
        do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
            ok = ok .and. isvalid_cons(state(:,i,j,k))
        end do; end do; end do;

        if (ok) return

    end do

    if (.not.ok) then
        write (*,'(a30,99(ES12.4))') 'Unphysical block averages:', squeeze, cons2prim(state(:,1,1,1))
        call Driver_abortFlash('## Game over! Could not repair cells due to unphysical state averages of the whole block! ##')
    end if

end subroutine

!! subroutine enforce_positivity(state)
!!
!!     use dgfv_fluxes_mod, only: isvalid_cons
!!     use Driver_interface, only: Driver_abortFlash
!!
!!     implicit none
!!
!!     real, intent(inout) :: state(NCONS_VARS,NXB,NYB,NZB)
!!
!!     integer, parameter :: nlevels = NINT(LOG(REAL(NXB))/LOG(2.0))
!!
!!     integer :: level, r,p,q,h
!!     integer :: i,j,k, ii,jj,kk
!!     integer :: mx,my,mz
!!     integer :: nx,ny,nz
!!
!!     logical :: ok
!!
!!     real :: ws, squeeze
!!     real :: means(NCONS_VARS)
!!     real :: tiles(NCONS_VARS,1+K1D,1+K2D,1+K3D)
!!     real :: hiles(NCONS_VARS,1+K1D,1+K2D,1+K3D)
!!
!!     real :: massfracs(NSPECIES_FLUX)
!!
!!     ok = .true.
!!     do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
!!         ok = ok .and. isvalid_cons(state(:,i,j,k))
!!
!!         ! if (.not.ok) then
!!         !     write (*,'(a30,99(ES12.4))') 'Unphysical block averages:', state(:,i,j,k)
!!         !     call Driver_abortFlash('## Game over! Could not repair cells due to unphysical state averages of the whole block! ##')
!!         ! end if
!!     end do; end do; end do;
!!
!!     if (ok) return
!!
!!     ! if (.not.ok) then
!!     !     write (*,'(a30,99(ES12.4))') 'Unphysical block averages:', squeeze, state(:,1,1,1)
!!     !     call Driver_abortFlash('## Game over! Could not repair cells due to unphysical state averages of the whole block! ##')
!!     ! end if
!!
!!     !write (*,*) '# Correcting posivity! #'
!!
!!     do level = 1,nlevels
!!
!! # if NDIM > 0
!!         mx = 2**(level-1)
!!         nx = 2**(nlevels-level)
!! # else
!!         mx = 1
!!         nx = 1
!! # endif
!! # if NDIM > 1
!!         my = mx
!!         ny = nx
!! # else
!!         my = 1
!!         ny = 1
!! # endif
!! # if NDIM > 2
!!         mz = mx
!!         nz = nx
!! # else
!!         mz = 1
!!         nz = 1
!! # endif
!!
!!         ws = 1.0/REAL(nx**NDIM)
!!
!!         do r = 1,mz; do p = 1,my; do q = 1,mx
!!
!!             tiles = 0.0
!!             do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                 do k = 1 + (2*(r-1) + kk-1)*nz, (2*(r-1) + kk)*nz
!!                 do j = 1 + (2*(p-1) + jj-1)*ny, (2*(p-1) + jj)*ny
!!                 do i = 1 + (2*(q-1) + ii-1)*nx, (2*(q-1) + ii)*nx
!!                     do h = 1,NCONS_VARS
!!                         tiles(h,ii,jj,kk) = tiles(h,ii,jj,kk) + ws*state(h,i,j,k)
!!                     end do
!!                 end do; end do; end do;
!!             end do; end do; end do;
!!
!!             ok = .true.
!!             do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                 ok = ok .and. isvalid_cons(tiles(:,ii,jj,kk))
!!             end do; end do; end do;
!!             if (ok) cycle
!!
!!             means = 0.0
!!             do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                 do h = 1,NCONS_VARS
!!                     means(h) = means(h) + 1.0/REAL(2**NDIM)*tiles(h,ii,jj,kk)
!!                 end do
!!             end do; end do; end do;
!!
!!             do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                 do h = 1,NCONS_VARS
!!                     hiles(h,ii,jj,kk) = tiles(h,ii,jj,kk) - means(h)
!!                 end do
!!             end do; end do; end do;
!!
!!             squeeze = 1.0
!!             do while (squeeze > 0.0)
!!                 squeeze = max(0.0,squeeze - 0.1)
!!
!!                 do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                     do h = 1,NCONS_VARS
!!                         tiles(h,ii,jj,kk) = means(h) + squeeze*hiles(h,ii,jj,kk)
!!                     end do
!!                 end do; end do; end do;
!!
!!                 ok = .true.
!!                 do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                     ok = ok .and. isvalid_cons(tiles(:,ii,jj,kk))
!!                 end do; end do; end do;
!!                 if (ok) exit
!!             end do
!!
!!             do kk = 1,1+K3D; do jj = 1,1+K2D; do ii = 1,1+K1D;
!!                 do k = 1 + (2*(r-1) + kk-1)*nz, (2*(r-1) + kk)*nz
!!                 do j = 1 + (2*(p-1) + jj-1)*ny, (2*(p-1) + jj)*ny
!!                 do i = 1 + (2*(q-1) + ii-1)*nx, (2*(q-1) + ii)*nx
!!                     do h = 1,NCONS_VARS
!!                         state(h,i,j,k) = means(h) + squeeze * (state(h,i,j,k) - means(h))
!!                     end do
!!                 end do; end do; end do;
!!             end do; end do; end do;
!!
!!         end do; end do; end do !! tiles
!!
!!     end do !! levels
!!
!! # if NSPECIES > 0
!!     do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
!!         massfracs = state(SPEC_CONS_BEGIN:SPEC_CONS_END,i,j,k)/state(DENS_CONS,i,j,k)
!!         massfracs = merge(1e-12,massfracs,massfracs < 1e-12)
!!         massfracs = massfracs/sum(massfracs)
!!         state(SPEC_CONS_BEGIN:SPEC_CONS_END,i,j,k) = massfracs*state(DENS_CONS,i,j,k)
!!     end do; end do; end do;
!! # endif
!!
!!     ok = .true.
!!     do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
!!         ok = ok .and. isvalid_cons(state(:,i,j,k))
!!     end do; end do; end do;
!!
!!     if (.not.ok) then
!!         write (*,'(a30,99(ES12.4))') 'Unphysical block averages:', squeeze, state(:,1,1,1)
!!         call Driver_abortFlash('## Game over! Could not repair cells due to unphysical state averages of the whole block! ##')
!!     end if
!!
!! end subroutine

subroutine calc_fv_surface_flux(state,flux)

    use dgfv_fluxes_mod, only: xriemann => xriem_fv_prim
    use dgfv_fluxes_mod, only: yriemann => yriem_fv_prim
    use dgfv_fluxes_mod, only: zriemann => zriem_fv_prim

    use Hydro_data, only: dgfv_fv_active

    implicit none

    real, intent(in)    :: state(NPRIM_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real, intent(out)   :: flux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)

    !! Interpolated left/minux and right/plus face variables.
    real :: xuP(NPRIM_VARS,0:NXB,NYB,NZB)
    real :: xuM(NPRIM_VARS,0:NXB,NYB,NZB)

    real :: yuP(NPRIM_VARS,NXB,0:NYB,NZB)
    real :: yuM(NPRIM_VARS,NXB,0:NYB,NZB)

    real :: zuP(NPRIM_VARS,NXB,NYB,0:NZB)
    real :: zuM(NPRIM_VARS,NXB,NYB,0:NZB)

    integer :: i,j,k,h, ii,jj,kk

    !! Needed for mono-central interpolation scheme.
    real, parameter :: posvars(*) = (/DENS_PRIM,PRES_PRIM,ENER_PRIM/)
    real, parameter :: lowfactor = 1e-3

    !! ------------------------------------------------------------------ !!
    !! Reconstruct interface values.

# if HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_FIRST_ORDER
    !! First order interpolation, i.e. just copying.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xuP(:,i,j,k) = state(:,i  ,j,k)
        xuM(:,i,j,k) = state(:,i+1,j,k)
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(:,i,j,k) = state(:,i,j  ,k)
        yuM(:,i,j,k) = state(:,i,j+1,k)
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(:,i,j,k) = state(:,i,j,k  )
        zuM(:,i,j,k) = state(:,i,j,k+1)
    end do; end do; end do;
# endif

# elif HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_MINMOD
!! Second order minmod-limited reconstruction.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xuP(:,i,j,k) = state(:,i  ,j,k) + 0.5*minmod2(&
                  state(:,i  ,j,k) - state(:,i-1,j,k), &
                  state(:,i+1,j,k) - state(:,i  ,j,k))

        xuM(:,i,j,k) = state(:,i+1,j,k) - 0.5*minmod2(&
                  state(:,i+1,j,k) - state(:,i  ,j,k), &
                  state(:,i+2,j,k) - state(:,i+1,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(:,i,j,k) = state(:,i,j  ,k) + 0.5*minmod2(&
                  state(:,i,j  ,k) - state(:,i,j-1,k), &
                  state(:,i,j+1,k) - state(:,i,j  ,k))

        yuM(:,i,j,k) = state(:,i,j+1,k) - 0.5*minmod2(&
                  state(:,i,j+1,k) - state(:,i,j  ,k), &
                  state(:,i,j+2,k) - state(:,i,j+1,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(:,i,j,k) = state(:,i,j,k  ) + 0.5*minmod2(&
                  state(:,i,j,k  ) - state(:,i,j,k-1), &
                  state(:,i,j,k+1) - state(:,i,j,k  ))

        zuM(:,i,j,k) = state(:,i,j,k+1) - 0.5*minmod2(&
                  state(:,i,j,k+1) - state(:,i,j,k  ), &
                  state(:,i,j,k+2) - state(:,i,j,k+1))
    end do; end do; end do;
# endif

!! Deprecated.
!! # if NSPECIES > 0
!!     !! First order interpolation for species.
!! # if NDIM > 0
!!     do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
!!         xuP(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i  ,j,k)
!!         xuM(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i+1,j,k)
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 1
!!     do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
!!         yuP(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j  ,k)
!!         yuM(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j+1,k)
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 2
!!     do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
!!         zuP(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k  )
!!         zuM(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k+1)
!!     end do; end do; end do;
!! # endif
!! # endif
!! 
!! # if NMASS_SCALARS > 0
!! # if NDIM > 0
!!     do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
!!         xuP(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i  ,j,k)
!!         xuM(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i+1,j,k)
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 1
!!     do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
!!         yuP(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j  ,k)
!!         yuM(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j+1,k)
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 2
!!     do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
!!         zuP(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k  )
!!         zuM(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k+1)
!!     end do; end do; end do;
!! # endif
!! # endif

# elif HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_MONO_CENTRAL
!! Second order monotonized-central reconstruction.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xuP(:,i,j,k) = state(:,i  ,j,k) + 0.25*minmod3(&
                  4.0*(state(:,i  ,j,k) - state(:,i-1,j,k)), &
                      (state(:,i+1,j,k) - state(:,i-1,j,k)), &
                  4.0*(state(:,i+1,j,k) - state(:,i  ,j,k)))

        xuM(:,i,j,k) = state(:,i+1,j,k) - 0.25*minmod3(&
                  4.0*(state(:,i+1,j,k) - state(:,i  ,j,k)), &
                      (state(:,i+2,j,k) - state(:,i  ,j,k)), &
                  4.0*(state(:,i+2,j,k) - state(:,i+1,j,k)))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(:,i,j,k) = state(:,i,j  ,k) + 0.25*minmod3(&
                  4.0*(state(:,i,j  ,k) - state(:,i,j-1,k)), &
                      (state(:,i,j+1,k) - state(:,i,j-1,k)), &
                  4.0*(state(:,i,j+1,k) - state(:,i,j  ,k)))

        yuM(:,i,j,k) = state(:,i,j+1,k) - 0.25*minmod3(&
                  4.0*(state(:,i,j+1,k) - state(:,i,j  ,k)), &
                      (state(:,i,j+2,k) - state(:,i,j  ,k)), &
                  4.0*(state(:,i,j+2,k) - state(:,i,j+1,k)))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(:,i,j,k) = state(:,i,j,k  ) + 0.25*minmod3(&
                  4.0*(state(:,i,j,k  ) - state(:,i,j,k-1)), &
                      (state(:,i,j,k+1) - state(:,i,j,k-1)), &
                  4.0*(state(:,i,j,k+1) - state(:,i,j,k  )))

        zuM(:,i,j,k) = state(:,i,j,k+1) - 0.25*minmod3(&
                  4.0*(state(:,i,j,k+1) - state(:,i,j,k  )), &
                      (state(:,i,j,k+2) - state(:,i,j,k  )), &
                  4.0*(state(:,i,j,k+2) - state(:,i,j,k+1)))
    end do; end do; end do;
# endif

    !! To prevent undershoots, limit interface states to lower thresholds.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xuP(posvars,i,j,k) = max(lowfactor*state(posvars,i  ,j,k),xuP(posvars,i,j,k))
        xuM(posvars,i,j,k) = max(lowfactor*state(posvars,i+1,j,k),xuM(posvars,i,j,k))
    end do; end do; end do;
# if NSPECIES > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB; do h = SPEC_PRIM_BEGIN,SPEC_PRIM_END;
        xuP(h,i,j,k) = max(lowfactor*state(h,i  ,j,k),xuP(h,i,j,k))
        xuM(h,i,j,k) = max(lowfactor*state(h,i+1,j,k),xuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB; do h = MSCL_PRIM_BEGIN,MSCL_PRIM_END;
        xuP(h,i,j,k) = max(lowfactor*state(h,i  ,j,k),xuP(h,i,j,k))
        xuM(h,i,j,k) = max(lowfactor*state(h,i+1,j,k),xuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 0 */

# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(posvars,i,j,k) = max(lowfactor*state(posvars,i,j  ,k),yuP(posvars,i,j,k))
        yuM(posvars,i,j,k) = max(lowfactor*state(posvars,i,j+1,k),yuM(posvars,i,j,k))
    end do; end do; end do;
# if NSPECIES > 0
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB; do h = SPEC_PRIM_BEGIN,SPEC_PRIM_END;
        yuP(h,i,j,k) = max(lowfactor*state(h,i,j  ,k),yuP(h,i,j,k))
        yuM(h,i,j,k) = max(lowfactor*state(h,i,j+1,k),yuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB; do h = MSCL_PRIM_BEGIN,MSCL_PRIM_END;
        yuP(h,i,j,k) = max(lowfactor*state(h,i,j  ,k),yuP(h,i,j,k))
        yuM(h,i,j,k) = max(lowfactor*state(h,i,j+1,k),yuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 1 */

# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(posvars,i,j,k) = max(lowfactor*state(posvars,i,j,k  ),zuP(posvars,i,j,k))
        zuM(posvars,i,j,k) = max(lowfactor*state(posvars,i,j,k+1),zuM(posvars,i,j,k))
    end do; end do; end do;
# if NSPECIES > 0
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB; do h = SPEC_PRIM_BEGIN,SPEC_PRIM_END;
        zuP(h,i,j,k) = max(lowfactor*state(h,i,j,k  ),zuP(h,i,j,k))
        zuM(h,i,j,k) = max(lowfactor*state(h,i,j,k+1),zuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB; do h = MSCL_PRIM_BEGIN,MSCL_PRIM_END;
        zuP(h,i,j,k) = max(lowfactor*state(h,i,j,k  ),zuP(h,i,j,k))
        zuM(h,i,j,k) = max(lowfactor*state(h,i,j,k+1),zuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 2 */
# else
# error: Unknown reconstruction mode for FV scheme.
# endif /* HYDRO_FV_RECONSTRUCTION */

    !! ------------------------------------------------------------------ !!
    !! Compute interface fluxes

# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        flux(:,i,j,k,IAXIS) = xriemann(xuP(:,i,j,k),xuM(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        flux(:,i,j,k,JAXIS) = yriemann(yuP(:,i,j,k),yuM(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        flux(:,i,j,k,KAXIS) = zriemann(zuP(:,i,j,k),zuM(:,i,j,k))
    end do; end do; end do;
# endif

    contains

    pure elemental function minmod2(a1,a2) result(mm)

        implicit none

        real, intent(in) :: a1,a2
        real             :: mm

        if ((sign(1.0,a1)*sign(1.0,a2)) > 0.0) then
            mm = sign(min(abs(a1),abs(a2)),a1)
        else
            mm = 0.0
        end if

    end function

    pure elemental function minmod3(a1,a2,a3) result(mm)

        implicit none

        real, intent(in) :: a1, a2, a3
        real             :: mm

        if (((sign(1.0,a1)*sign(1.0,a2)) > 0.0) .and. ((sign(1.0,a2)*sign(1.0,a3)) > 0.0)) then
            mm = sign(min(abs(a1),abs(a2),abs(a3)),a1)
        else
            mm = 0.0
        end if

    end function

end subroutine

!! ================================================================== !!

subroutine calc_fv_rhs(state,flux,sdx,rhs,erhs)

    use dgfv_fluxes_mod, only: Non_Conservative_xFlux_prim
    use dgfv_fluxes_mod, only: Non_Conservative_yFlux_prim
    use dgfv_fluxes_mod, only: Non_Conservative_zFlux_prim

    use Hydro_data, only: DGFV_MHD_NON_CONS
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION

    use dgfv_fluxes_mod, only: cons2evec
    use dgfv_fluxes_mod, only: evec2cons 
    use dgfv_fluxes_mod, only: prim2evec

    use dgfv_fluxes_mod, only: prim2xpot
    use dgfv_fluxes_mod, only: prim2ypot
    use dgfv_fluxes_mod, only: prim2zpot

    use dgfv_fluxes_mod, only: xefluxjump_prim
    use dgfv_fluxes_mod, only: yefluxjump_prim
    use dgfv_fluxes_mod, only: zefluxjump_prim

    use dgfv_fluxes_mod, only: xefluxmean_prim
    use dgfv_fluxes_mod, only: yefluxmean_prim
    use dgfv_fluxes_mod, only: zefluxmean_prim

    use dgfv_fluxes_mod, only: xentrflux_prim
    use dgfv_fluxes_mod, only: yentrflux_prim
    use dgfv_fluxes_mod, only: zentrflux_prim

    ! use dgfv_fluxes_mod, only: xriemann_species_surface
    ! use dgfv_fluxes_mod, only: yriemann_species_surface
    ! use dgfv_fluxes_mod, only: zriemann_species_surface

    implicit none

    real, intent(in)    :: state(NPRIM_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real, intent(in)    :: flux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
    real, intent(in)    :: sdx(NDIM)
    real, intent(out)   :: rhs(NFLUXES,NXB,NYB,NZB)
    real, intent(out)   :: erhs(NXB,NYB,NZB)

    ! real :: corrected_flux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)

    !! Interpolated left/minux and right/plus face variables.
    real :: xuP(NPRIM_VARS,0:NXB,NYB,NZB)
    real :: xuM(NPRIM_VARS,0:NXB,NYB,NZB)

    real :: yuP(NPRIM_VARS,NXB,0:NYB,NZB)
    real :: yuM(NPRIM_VARS,NXB,0:NYB,NZB)

    real :: zuP(NPRIM_VARS,NXB,NYB,0:NZB)
    real :: zuM(NPRIM_VARS,NXB,NYB,0:NZB)

    integer :: i,j,k,h, ii,jj,kk

    !! Reconstructed variables for non-conservative terms.
    ! integer, parameter :: recovars(*) = (/DENS_PRIM,PRES_PRIM,ENER_PRIM,VELX_PRIM:VELZ_PRIM,MAGX_PRIM:MAGZ_PRIM,GLMP_PRIM/)
    ! integer, parameter :: recovars(*) = (/MAGX_PRIM:MAGZ_PRIM,GLMP_PRIM/)

    !! Needed for mono-central interpolation scheme.
    real, parameter :: posvars(*) = (/DENS_PRIM,PRES_PRIM,ENER_PRIM/)
    real, parameter :: lowfactor = 1e-3

    !! ------------------------------------------------------------------ !!
    !! Compute right-hand-side (RHS).

    ! corrected_flux = flux

! # if NSPECIES > 0
! # if NDIM > 0
!     do k = 1,NZB; do j = 1,NYB;
!         corrected_flux(:,  0,j,k,IAXIS) = xriemann_species_surface(flux(:,  0,j,k,IAXIS),state(:,  0,j,k),state(:,    1,j,k))
!         corrected_flux(:,NXB,j,k,IAXIS) = xriemann_species_surface(flux(:,NXB,j,k,IAXIS),state(:,NXB,j,k),state(:,NXB+1,j,k))
!     end do; end do;
! # endif
! # if NDIM > 1
!     do k = 1,NZB; do i = 1,NXB;
!         corrected_flux(:,i,  0,k,JAXIS) = yriemann_species_surface(flux(:,i,  0,k,JAXIS),state(:,i,  0,k),state(:,i,    1,k))
!         corrected_flux(:,i,NYB,k,JAXIS) = yriemann_species_surface(flux(:,i,NYB,k,JAXIS),state(:,i,NYB,k),state(:,i,NYB+1,k))
!     end do; end do;
! # endif
! # if NDIM > 2
!     do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
!         corrected_flux(:,i,j,  0,KAXIS) = zriemann_species_surface(flux(:,i,j,  0,KAXIS),state(:,i,j,  0),state(:,i,j,    1))
!         corrected_flux(:,i,j,NZB,KAXIS) = zriemann_species_surface(flux(:,i,j,NZB,KAXIS),state(:,i,j,NZB),state(:,i,j,NZB+1))
!     end do; end do; end do;
! # endif
! # endif

    rhs = 0.0

    !! Standard Finite-Volume method.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB
        rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(IAXIS)*(flux(:,i-1,j,k,IAXIS) - flux(:,i,j,k,IAXIS))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(JAXIS)*(flux(:,i,j-1,k,JAXIS) - flux(:,i,j,k,JAXIS))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(KAXIS)*(flux(:,i,j,k-1,KAXIS) - flux(:,i,j,k,KAXIS))
    end do; end do; end do;
# endif

    if (DGFV_MHD_NON_CONS .or. DGFV_ENTROPY_CORRECTION) then
        !! ------------------------------------------------------------------ !!
        !! Reconstruct interface values.

# if HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_FIRST_ORDER
        !! First order interpolation, i.e. just copying.
# if NDIM > 0
        do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
            xuP(:,i,j,k) = state(:,i  ,j,k)
            xuM(:,i,j,k) = state(:,i+1,j,k)
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
            yuP(:,i,j,k) = state(:,i,j  ,k)
            yuM(:,i,j,k) = state(:,i,j+1,k)
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
            zuP(:,i,j,k) = state(:,i,j,k  )
            zuM(:,i,j,k) = state(:,i,j,k+1)
        end do; end do; end do;
# endif

# elif HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_MINMOD
!! Second order minmod-limited reconstruction.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xuP(:,i,j,k) = state(:,i  ,j,k) + 0.5*minmod2(&
                  state(:,i  ,j,k) - state(:,i-1,j,k), &
                  state(:,i+1,j,k) - state(:,i  ,j,k))

        xuM(:,i,j,k) = state(:,i+1,j,k) - 0.5*minmod2(&
                  state(:,i+1,j,k) - state(:,i  ,j,k), &
                  state(:,i+2,j,k) - state(:,i+1,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(:,i,j,k) = state(:,i,j  ,k) + 0.5*minmod2(&
                  state(:,i,j  ,k) - state(:,i,j-1,k), &
                  state(:,i,j+1,k) - state(:,i,j  ,k))

        yuM(:,i,j,k) = state(:,i,j+1,k) - 0.5*minmod2(&
                  state(:,i,j+1,k) - state(:,i,j  ,k), &
                  state(:,i,j+2,k) - state(:,i,j+1,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(:,i,j,k) = state(:,i,j,k  ) + 0.5*minmod2(&
                  state(:,i,j,k  ) - state(:,i,j,k-1), &
                  state(:,i,j,k+1) - state(:,i,j,k  ))

        zuM(:,i,j,k) = state(:,i,j,k+1) - 0.5*minmod2(&
                  state(:,i,j,k+1) - state(:,i,j,k  ), &
                  state(:,i,j,k+2) - state(:,i,j,k+1))
    end do; end do; end do;
# endif

# elif HYDRO_FV_RECONSTRUCTION == HYDRO_FV_RECONSTRUCTION_MONO_CENTRAL
    !! Second order monotonized-central reconstruction.
# if NDIM > 0
        do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
            xuP(:,i,j,k) = state(:,i  ,j,k) + 0.25*minmod3(&
                      4.0*(state(:,i  ,j,k) - state(:,i-1,j,k)), &
                          (state(:,i+1,j,k) - state(:,i-1,j,k)), &
                      4.0*(state(:,i+1,j,k) - state(:,i  ,j,k)))

            xuM(:,i,j,k) = state(:,i+1,j,k) - 0.25*minmod3(&
                      4.0*(state(:,i+1,j,k) - state(:,i  ,j,k)), &
                          (state(:,i+2,j,k) - state(:,i  ,j,k)), &
                      4.0*(state(:,i+2,j,k) - state(:,i+1,j,k)))
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
            yuP(:,i,j,k) = state(:,i,j  ,k) + 0.25*minmod3(&
                      4.0*(state(:,i,j  ,k) - state(:,i,j-1,k)), &
                          (state(:,i,j+1,k) - state(:,i,j-1,k)), &
                      4.0*(state(:,i,j+1,k) - state(:,i,j  ,k)))

            yuM(:,i,j,k) = state(:,i,j+1,k) - 0.25*minmod3(&
                      4.0*(state(:,i,j+1,k) - state(:,i,j  ,k)), &
                          (state(:,i,j+2,k) - state(:,i,j  ,k)), &
                      4.0*(state(:,i,j+2,k) - state(:,i,j+1,k)))
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
            zuP(:,i,j,k) = state(:,i,j,k  ) + 0.25*minmod3(&
                      4.0*(state(:,i,j,k  ) - state(:,i,j,k-1)), &
                          (state(:,i,j,k+1) - state(:,i,j,k-1)), &
                      4.0*(state(:,i,j,k+1) - state(:,i,j,k  )))

            zuM(:,i,j,k) = state(:,i,j,k+1) - 0.25*minmod3(&
                      4.0*(state(:,i,j,k+1) - state(:,i,j,k  )), &
                          (state(:,i,j,k+2) - state(:,i,j,k  )), &
                      4.0*(state(:,i,j,k+2) - state(:,i,j,k+1)))
        end do; end do; end do;
# endif

        !! To prevent undershoots, limit interface states to lower thresholds.
# if NDIM > 0
        do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
            xuP(posvars,i,j,k) = max(lowfactor*state(posvars,i  ,j,k),xuP(posvars,i,j,k))
            xuM(posvars,i,j,k) = max(lowfactor*state(posvars,i+1,j,k),xuM(posvars,i,j,k))
        end do; end do; end do;
# if NSPECIES > 0
        do k = 1,NZB; do j = 1,NYB; do i = 0,NXB; do h = SPECIES_BEGIN,SPECIES_END;
            xuP(h,i,j,k) = max(lowfactor*state(h,i  ,j,k),xuP(h,i,j,k))
            xuM(h,i,j,k) = max(lowfactor*state(h,i+1,j,k),xuM(h,i,j,k))
        end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB; do h = MASS_SCALARS_BEGIN,MASS_SCALARS_END;
        xuP(h,i,j,k) = max(lowfactor*state(h,i  ,j,k),xuP(h,i,j,k))
        xuM(h,i,j,k) = max(lowfactor*state(h,i+1,j,k),xuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 0 */

# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yuP(posvars,i,j,k) = max(lowfactor*state(posvars,i,j  ,k),yuP(posvars,i,j,k))
        yuM(posvars,i,j,k) = max(lowfactor*state(posvars,i,j+1,k),yuM(posvars,i,j,k))
    end do; end do; end do;
# if NSPECIES > 0
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB; do h = SPECIES_BEGIN,SPECIES_END;
        yuP(h,i,j,k) = max(lowfactor*state(h,i,j  ,k),yuP(h,i,j,k))
        yuM(h,i,j,k) = max(lowfactor*state(h,i,j+1,k),yuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB; do h = MASS_SCALARS_BEGIN,MASS_SCALARS_END;
        yuP(h,i,j,k) = max(lowfactor*state(h,i,j  ,k),yuP(h,i,j,k))
        yuM(h,i,j,k) = max(lowfactor*state(h,i,j+1,k),yuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 1 */

# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zuP(posvars,i,j,k) = max(lowfactor*state(posvars,i,j,k  ),zuP(posvars,i,j,k))
        zuM(posvars,i,j,k) = max(lowfactor*state(posvars,i,j,k+1),zuM(posvars,i,j,k))
    end do; end do; end do;
# if NSPECIES > 0
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB; do h = SPECIES_BEGIN,SPECIES_END;
        zuP(h,i,j,k) = max(lowfactor*state(h,i,j,k  ),zuP(h,i,j,k))
        zuM(h,i,j,k) = max(lowfactor*state(h,i,j,k+1),zuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# if NMASS_SCALARS > 0
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB; do h = MASS_SCALARS_BEGIN,MASS_SCALARS_END;
        zuP(h,i,j,k) = max(lowfactor*state(h,i,j,k  ),zuP(h,i,j,k))
        zuM(h,i,j,k) = max(lowfactor*state(h,i,j,k+1),zuM(h,i,j,k))
    end do; end do; end do; end do;
# endif
# endif /* NDIM > 2 */
# else
# error: Unknown reconstruction mode for FV scheme.
# endif /* HYDRO_FV_RECONSTRUCTION */
    end if !! DGFV_MHD_NON_CONS

    !! ------------------------------------------------------------------ !!
    !! Add non-conservative source terms.

    if (DGFV_MHD_NON_CONS) then
        associate(u => state)
        do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
# if NDIM > 0
            rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(IAXIS)*Non_Conservative_xFlux_prim(u(:,i,j,k),xuP(:,i-1,j,k),xuM(:,i-1,j,k))
            rhs(:,i,j,k) = rhs(:,i,j,k) - sdx(IAXIS)*Non_Conservative_xFlux_prim(u(:,i,j,k),xuP(:,i  ,j,k),xuM(:,i  ,j,k))
# endif
# if NDIM > 1
            rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(JAXIS)*Non_Conservative_yFlux_prim(u(:,i,j,k),yuP(:,i,j-1,k),yuM(:,i,j-1,k))
            rhs(:,i,j,k) = rhs(:,i,j,k) - sdx(JAXIS)*Non_Conservative_yFlux_prim(u(:,i,j,k),yuP(:,i,j  ,k),yuM(:,i,j  ,k))
# endif
# if NDIM > 2
            rhs(:,i,j,k) = rhs(:,i,j,k) + sdx(KAXIS)*Non_Conservative_zFlux_prim(u(:,i,j,k),zuP(:,i,j,k-1),zuM(:,i,j,k-1))
            rhs(:,i,j,k) = rhs(:,i,j,k) - sdx(KAXIS)*Non_Conservative_zFlux_prim(u(:,i,j,k),zuP(:,i,j,k  ),zuM(:,i,j,k  ))
# endif
        end do; end do; end do;
        end associate
    end if

    erhs = 0.0
# if 0
    !! ------------------------------------------------------------------ !!
    !! Calculate entropy production rate per cell.

    ! if (DGFV_ENTROPY_CORRECTION) then

    block
    real :: xejumps(0:NXB,1:NYB,1:NZB)
    real :: yejumps(1:NXB,0:NYB,1:NZB)
    real :: zejumps(1:NXB,1:NYB,0:NZB)
    real :: fvrhs(1:NXB,1:NYB,1:NZB)

# if 0
    associate(u => state)
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xejumps(i,j,k) = 0.5*sdx(IAXIS)*xefluxjump_prim(xuP(:,i,j,k),xuM(:,i,j,k),flux(:,i,j,k,IAXIS))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        yejumps(i,j,k) = 0.5*sdx(JAXIS)*yefluxjump_prim(yuP(:,i,j,k),yuM(:,i,j,k),flux(:,i,j,k,JAXIS))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zejumps(i,j,k) = 0.5*sdx(KAXIS)*zefluxjump_prim(zuP(:,i,j,k),zuM(:,i,j,k),flux(:,i,j,k,KAXIS))
    end do; end do; end do;
# endif
    end associate


    erhs = 0.0

    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
# if NDIM > 0
        erhs(i,j,k) = erhs(i,j,k) + xejumps(i-1,j,k) + xejumps(i,j,k)
# endif
# if NDIM > 1
        erhs(i,j,k) = erhs(i,j,k) + yejumps(i,j-1,k) + yejumps(i,j,k)
# endif
# if NDIM > 2
        erhs(i,j,k) = erhs(i,j,k) + zejumps(i,j,k-1) + zejumps(i,j,k)
# endif
    end do; end do; end do;

# else 

    associate(u => state)
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        !xejumps(i,j,k) = xefluxmean_prim(xuP(:,i,j,k),xuM(:,i,j,k),flux(:,i,j,k,IAXIS))
        xejumps(i,j,k) = xentrflux_prim(flux(:,i,j,k,IAXIS),xuP(:,i,j,k),xuM(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        !yejumps(i,j,k) = yefluxmean_prim(yuP(:,i,j,k),yuM(:,i,j,k),flux(:,i,j,k,JAXIS))
        yejumps(i,j,k) = yentrflux_prim(flux(:,i,j,k,JAXIS),yuP(:,i,j,k),yuM(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        !zejumps(i,j,k) = zefluxmean_prim(zuP(:,i,j,k),zuM(:,i,j,k),flux(:,i,j,k,KAXIS))
        zejumps(i,j,k) = zentrflux_prim(flux(:,i,j,k,KAXIS),zuP(:,i,j,k),zuM(:,i,j,k))
    end do; end do; end do;
# endif
    end associate

    erhs = 0.0

    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
# if NDIM > 0
        erhs(i,j,k) = erhs(i,j,k) + sdx(IAXIS)*(xejumps(i-1,j,k) - xejumps(i,j,k))
# endif
# if NDIM > 1
        erhs(i,j,k) = erhs(i,j,k) + sdx(JAXIS)*(yejumps(i,j-1,k) - yejumps(i,j,k))
# endif
# if NDIM > 2
        erhs(i,j,k) = erhs(i,j,k) + sdx(KAXIS)*(zejumps(i,j,k-1) - zejumps(i,j,k))
# endif
    end do; end do; end do;

    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        erhs(i,j,k) = dot_product(prim2evec(state(:,i,j,k)),rhs(:,i,j,k)) - erhs(i,j,k)
        !erhs(i,j,k) = dot_product(prim2evec(state(:,i,j,k)),rhs(:,i,j,k))
        !erhs(i,j,k) = fvrhs(i,j,k) - erhs(i,j,k)
        !erhs(i,j,k) = erhs(i,j,k)
    end do; end do; end do;
# endif

    end block
# endif

    contains

    pure elemental function minmod2(a1,a2) result(mm)

        implicit none

        real, intent(in) :: a1,a2
        real             :: mm

        if ((sign(1.0,a1)*sign(1.0,a2)) > 0.0) then
            mm = sign(min(abs(a1),abs(a2)),a1)
        else
            mm = 0.0
        end if

    end function

    pure elemental function minmod3(a1,a2,a3) result(mm)

        implicit none

        real, intent(in) :: a1, a2, a3
        real             :: mm

        if (((sign(1.0,a1)*sign(1.0,a2)) > 0.0) .and. ((sign(1.0,a2)*sign(1.0,a3)) > 0.0)) then
            mm = sign(min(abs(a1),abs(a2),abs(a3)),a1)
        else
            mm = 0.0
        end if

    end function

end subroutine

subroutine calc_species_rhs(state,flux,sdx,rhs)

    ! use dgfv_fluxes_mod, only: xriemann_species_surface
    ! use dgfv_fluxes_mod, only: yriemann_species_surface
    ! use dgfv_fluxes_mod, only: zriemann_species_surface

    implicit none

    real, intent(in)    :: state(NPRIM_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real, intent(in)    :: flux(0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
    real, intent(in)    :: sdx(NDIM)
    real, intent(inout) :: rhs(NFLUXES,NXB,NYB,NZB)

# if NSPECIES > 0
    integer, parameter :: sb = SPECIES_FLUX_BEGIN
    integer, parameter :: se = SPECIES_FLUX_END

    real :: xsP(NSPECIES,0:NXB,NYB,NZB)
    real :: xsM(NSPECIES,0:NXB,NYB,NZB)

    real :: ysP(NSPECIES,NXB,0:NYB,NZB)
    real :: ysM(NSPECIES,NXB,0:NYB,NZB)

    real :: zsP(NSPECIES,NXB,NYB,0:NZB)
    real :: zsM(NSPECIES,NXB,NYB,0:NZB)

    real :: xsf(NSPECIES,0:NXB,NYB,NZB)
    real :: ysf(NSPECIES,NXB,0:NYB,NZB)
    real :: zsf(NSPECIES,NXB,NYB,0:NZB)
# endif

# if NMASS_SCALARS > 0
    integer, parameter :: mb = MASS_SCALARS_FLUX_BEGIN
    integer, parameter :: me = MASS_SCALARS_FLUX_END

    real :: xmP(NMASS_SCALARS,0:NXB,NYB,NZB)
    real :: xmM(NMASS_SCALARS,0:NXB,NYB,NZB)

    real :: ymP(NMASS_SCALARS,NXB,0:NYB,NZB)
    real :: ymM(NMASS_SCALARS,NXB,0:NYB,NZB)

    real :: zmP(NMASS_SCALARS,NXB,NYB,0:NZB)
    real :: zmM(NMASS_SCALARS,NXB,NYB,0:NZB)

    real :: xmf(NMASS_SCALARS,0:NXB,NYB,NZB)
    real :: ymf(NMASS_SCALARS,NXB,0:NYB,NZB)
    real :: zmf(NMASS_SCALARS,NXB,NYB,0:NZB)
# endif

    integer :: i,j,k

    !! ------------------------------------------------------------------ !!

# if NSPECIES > 0
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xsP(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i  ,j,k)
        xsM(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i+1,j,k)
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        ysP(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j  ,k)
        ysM(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j+1,k)
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zsP(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k  )
        zsM(:,i,j,k) = state(SPEC_PRIM_BEGIN:SPEC_PRIM_END,i,j,k+1)
    end do; end do; end do;
# endif

    !! Robust upwind scheme.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xsf(:,i,j,k) = flux(i,j,k,IAXIS)*merge(xsP(:,i,j,k),xsM(:,i,j,k),flux(i,j,k,IAXIS) > 0.0)
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        ysf(:,i,j,k) = flux(i,j,k,JAXIS)*merge(ysP(:,i,j,k),ysM(:,i,j,k),flux(i,j,k,JAXIS) > 0.0)
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zsf(:,i,j,k) = flux(i,j,k,KAXIS)*merge(zsP(:,i,j,k),zsM(:,i,j,k),flux(i,j,k,KAXIS) > 0.0)
    end do; end do; end do;
# endif

    rhs(sb:se,:,:,:) = 0.0

    !! Standard Finite-Volume method.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB
        rhs(sb:se,i,j,k) = rhs(sb:se,i,j,k) + sdx(IAXIS)*(xsf(:,i-1,j,k) - xsf(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(sb:se,i,j,k) = rhs(sb:se,i,j,k) + sdx(JAXIS)*(ysf(:,i,j-1,k) - ysf(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(sb:se,i,j,k) = rhs(sb:se,i,j,k) + sdx(KAXIS)*(zsf(:,i,j,k-1) - zsf(:,i,j,k))
    end do; end do; end do;
# endif
# endif !! NSPECIES > 0

# if NMASS_SCALARS > 0
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xmP(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i  ,j,k)
        xmM(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i+1,j,k)
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        ymP(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j  ,k)
        ymM(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j+1,k)
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zmP(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k  )
        zmM(:,i,j,k) = state(MSCL_PRIM_BEGIN:MSCL_PRIM_END,i,j,k+1)
    end do; end do; end do;
# endif

    !! Robust upwind scheme.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 0,NXB;
        xmf(:,i,j,k) = flux(i,j,k,IAXIS)*merge(xmP(:,i,j,k),xmM(:,i,j,k),flux(i,j,k,IAXIS) > 0.0)
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 0,NYB; do i = 1,NXB;
        ymf(:,i,j,k) = flux(i,j,k,JAXIS)*merge(ymP(:,i,j,k),ymM(:,i,j,k),flux(i,j,k,JAXIS) > 0.0)
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZB; do j = 1,NYB; do i = 1,NXB;
        zmf(:,i,j,k) = flux(i,j,k,KAXIS)*merge(zmP(:,i,j,k),zmM(:,i,j,k),flux(i,j,k,KAXIS) > 0.0)
    end do; end do; end do;
# endif

    rhs(mb:me,:,:,:) = 0.0

    !! Standard Finite-Volume method.
# if NDIM > 0
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB
        rhs(mb:me,i,j,k) = rhs(mb:me,i,j,k) + sdx(IAXIS)*(xmf(:,i-1,j,k) - xmf(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(mb:me,i,j,k) = rhs(mb:me,i,j,k) + sdx(JAXIS)*(ymf(:,i,j-1,k) - ymf(:,i,j,k))
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
        rhs(mb:me,i,j,k) = rhs(mb:me,i,j,k) + sdx(KAXIS)*(zmf(:,i,j,k-1) - zmf(:,i,j,k))
    end do; end do; end do;
# endif
# endif !! NMASS_SCALARS > 0

end subroutine

!! ================================================================== !!

subroutine calc_dg_surface_flux(state,flux,sfblend,vlblend)

    use dgfv_fluxes_mod, only: xriemann => xriem_dg_prim
    use dgfv_fluxes_mod, only: yriemann => yriem_dg_prim
    use dgfv_fluxes_mod, only: zriemann => zriem_dg_prim

    use Hydro_data, only: recoMat   !! mean values -> node values reconstruction matrix
    use Hydro_data, only: projMat   !! node values -> mean values projection matrix
    use Hydro_data, only: projTam   !! node values -> mean values projection matrix

    use Hydro_data, only: diffMat

    use Hydro_data, only: toBoundaryVecM  !! interpolation operator (left/minus face)
    use Hydro_data, only: toBoundaryVecP  !! interpolation operator (right/plus face)

    use Hydro_data, only: toMidpointsMat  !! interpolation operator from nodes to midpoints

    use Hydro_data, only: NXE,NYE,NZE
    use Hydro_data, only: NXN,NYN,NZN

    use dgfv_fluxes_mod, only: prim2cons
    use dgfv_fluxes_mod, only: cons2prim
    use dgfv_fluxes_mod, only: cons2evec
    use dgfv_fluxes_mod, only: evec2prim

    use dgfv_fluxes_mod, only: isvalid_prim
    use dgfv_fluxes_mod, only: isvalid_cons

    use Hydro_data, only: DGFV_FV_ACTIVE
    use Hydro_data, only: DGFV_FLUX_DIFFERENCING
    use Hydro_data, only: DGFV_ENTROPY_BOUNDARY_PROJECTION
    !use Hydro_data, only: DGFV_ENTROPY_CORRECTION
    !use Hydro_data, only: DGFV_ENTROPY_CORRECTION_ACTIVE

    use Hydro_data, only: DGFV_NODE_TYPE
    use Hydro_data, only: DGFV_NODE_TYPE_GAUSS

    implicit none

    real, intent(in)    :: state(NCONS_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)

    real, intent(out)   :: flux (NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
    real, intent(out)   :: sfblend(                NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
    real, intent(out)   :: vlblend(                           NXE,NYE,NZE)

    real :: blend(NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: ufv (NPRIM_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: cfv (max(NPRIM_VARS,NCONS_VARS), NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: udg (max(NPRIM_VARS,NCONS_VARS), NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: edg (max(NPRIM_VARS,NCONS_VARS), NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: atmp(max(NPRIM_VARS,NCONS_VARS), NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: btmp(max(NPRIM_VARS,NCONS_VARS), NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    !integer, parameter :: jumpvars(*) = (/DENS_PRIM,PRES_PRIM/)
    !integer, parameter :: jumpvars(*) = (/DENS_PRIM,PRES_PRIM,MAGX_PRIM,MAGY_PRIM/)
    !integer, parameter :: jumpvars(*) = (/PRES_PRIM,MAGX_PRIM,MAGY_PRIM/)
    !integer, parameter :: jumpvars(*) = (/DENS_PRIM/)
    integer, parameter :: jumpvars(*) = (/PRES_PRIM/)

    real :: itf (max(NPRIM_VARS,NCONS_VARS), NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)
    real :: imp (size(jumpvars), NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)
    !real :: etf (max(NPRIM_VARS,NCONS_VARS), NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)

    real :: jumpFV (size(jumpvars), NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI, IAXIS:KAXIS)
    real :: jumpDG (size(jumpvars), NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI, IAXIS:KAXIS)

    real :: rim (NFLUXES,      NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

    integer :: i,j,k,h, ii,jj,kk,ll

    logical :: ok

    !! ------------------------------------------------------------------ !!
    !! Split block into DG elements.

    do k = NZE_LO,NZE_HI; do j = NYE_LO,NYE_HI; do i = NXE_LO,NXE_HI;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
            ufv(:,ii,jj,kk,i,j,k) = state(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*(k-1)+kk)
        end do; end do; end do
    end do; end do; end do

    do k = NZE_LO,NZE_HI; do j = NYE_LO,NYE_HI; do i = NXE_LO,NXE_HI;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
            cfv(:,ii,jj,kk,i,j,k) = prim2cons(ufv(:,ii,jj,kk,i,j,k))
        end do; end do; end do
    end do; end do; end do

    !! ------------------------------------------------------------------ !!
    !! Transform to nodal space.

# if NDIM == 1
# define atmp udg
# elif NDIM == 2
# define btmp udg
# elif NDIM == 3
# define ctmp udg
# endif

# if NDIM > 0
    atmp = 0.0
# endif
# if NDIM > 1
    btmp = 0.0
# endif
# if NDIM > 2
    ctmp = 0.0
# endif

!! Inner elements.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,j,k) = atmp(h,ii,jj,kk,i,j,k) + recoMat(ii,ll)*cfv(h,ll,jj,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,j,k) = btmp(h,ii,jj,kk,i,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,j,k) = ctmp(h,ii,jj,kk,i,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM > 0
!! NORTH side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,NXE_LO,j,k) = atmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(ii,ll)*cfv(h,ll,jj,kk,NXE_LO,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,NXE_LO,j,k) = btmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_LO,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,NXE_LO,j,k) = ctmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_LO,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif

!! SOUTH side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,NXE_HI,j,k) = atmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(ii,ll)*cfv(h,ll,jj,kk,NXE_HI,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,NXE_HI,j,k) = btmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_HI,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,NXE_HI,j,k) = ctmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_HI,j,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# endif /* NDIM > 0 */

# if NDIM > 1
!! WEST side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,NYE_LO,k) = atmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(ii,ll)*cfv(h,ll,jj,kk,i,NYE_LO,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,NYE_LO,k) = btmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_LO,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,NYE_LO,k) = ctmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_LO,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif

!! EAST side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,NYE_HI,k) = atmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(ii,ll)*cfv(h,ll,jj,kk,i,NYE_HI,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,NYE_HI,k) = btmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_HI,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,NYE_HI,k) = ctmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_HI,k)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# endif /* NDIM > 1 */

# if NDIM > 2
!! FRONT side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,j,NZE_LO) = atmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(ii,ll)*cfv(h,ll,jj,kk,i,j,NZE_LO)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,j,NZE_LO) = btmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_LO)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,j,NZE_LO) = ctmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_LO)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif

!! BACK side.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,j,NZE_HI) = atmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(ii,ll)*cfv(h,ll,jj,kk,i,j,NZE_HI)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,j,NZE_HI) = btmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_HI)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,j,NZE_HI) = ctmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_HI)
        end do; end do; end do; end do; end do;
    end do; end do;
# endif
# endif /* NDIM > 2 */

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

    vlblend = 1.0

    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        ok = .true.
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
            ok = ok .and. isvalid_cons(udg(:,ii,jj,kk,i,j,k))
        end do; end do; end do;
        if (.not. ok) vlblend(i,j,k) = 0.0
    end do; end do; end do;

    !! ------------------------------------------------------------------ !!

    !if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS .and. DGFV_ENTROPY_BOUNDARY_PROJECTION) then
    if (DGFV_ENTROPY_BOUNDARY_PROJECTION) then

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,j,k) = cons2evec(udg(:,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,NXE_LO,j,k) = cons2evec(udg(:,ii,jj,kk,NXE_LO,j,k))
                edg(:,ii,jj,kk,NXE_HI,j,k) = cons2evec(udg(:,ii,jj,kk,NXE_HI,j,k))
            end do; end do; end do;
        end do; end do
# endif
# if NDIM > 1
        do k = 1,NZE; do i = 1,NXE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,NYE_LO,k) = cons2evec(udg(:,ii,jj,kk,i,NYE_LO,k))
                edg(:,ii,jj,kk,i,NYE_HI,k) = cons2evec(udg(:,ii,jj,kk,i,NYE_HI,k))
            end do; end do; end do;
        end do; end do
# endif
# if NDIM > 2
        do j = 1,NYE; do i = 1,NXE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,j,NZE_LO) = cons2evec(udg(:,ii,jj,kk,i,j,NZE_LO))
                edg(:,ii,jj,kk,i,j,NZE_HI) = cons2evec(udg(:,ii,jj,kk,i,j,NZE_HI))
            end do; end do; end do;
        end do; end do
# endif

        !! ------------------------------------------------------------------ !!
        !! Interpolate to element boundaries.

        itf = 0.0

# if NDIM > 0
        !! Interpolate to boundaries in x-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*edg(h,ll,jj,kk,i+1,j,k)
                itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*edg(h,ll,jj,kk,i  ,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        !! Interpolate to boundaries in y-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*edg(h,ii,ll,kk,i,j+1,k)
                itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*edg(h,ii,ll,kk,i,j  ,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        !! Interpolate to boundaries in z-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*edg(h,ii,jj,ll,i,j,k+1)
                itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*edg(h,ii,jj,ll,i,j,k  )
            end do; end do; end do; end do;
        end do; end do; end do;
# endif

        !! ------------------------------------------------------------------ !!
        !! Transform back to primitive variables.

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                itf(:,jj,kk,i,j,k,NOR) = evec2prim(itf(:,jj,kk,i,j,k,NOR))
                itf(:,jj,kk,i,j,k,SOU) = evec2prim(itf(:,jj,kk,i,j,k,SOU))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                itf(:,ii,kk,i,j,k,WES) = evec2prim(itf(:,ii,kk,i,j,k,WES))
                itf(:,ii,kk,i,j,k,EAS) = evec2prim(itf(:,ii,kk,i,j,k,EAS))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                itf(:,ii,jj,i,j,k,FRO) = evec2prim(itf(:,ii,jj,i,j,k,FRO))
                itf(:,ii,jj,i,j,k,BAC) = evec2prim(itf(:,ii,jj,i,j,k,BAC))
            end do; end do;
        end do; end do; end do;
# endif

    else

        !! ------------------------------------------------------------------ !!
        !! Interpolate to element boundaries.

        itf = 0.0

# if NDIM > 0
        !! Interpolate to boundaries in x-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*udg(h,ll,jj,kk,i+1,j,k)
                itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*udg(h,ll,jj,kk,i  ,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        !! Interpolate to boundaries in y-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*udg(h,ii,ll,kk,i,j+1,k)
                itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*udg(h,ii,ll,kk,i,j  ,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        !! Interpolate to boundaries in z-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*udg(h,ii,jj,ll,i,j,k+1)
                itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*udg(h,ii,jj,ll,i,j,k  )
            end do; end do; end do; end do;
        end do; end do; end do;
# endif

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                itf(:,jj,kk,i,j,k,NOR) = cons2prim(itf(:,jj,kk,i,j,k,NOR))
                itf(:,jj,kk,i,j,k,SOU) = cons2prim(itf(:,jj,kk,i,j,k,SOU))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                itf(:,ii,kk,i,j,k,WES) = cons2prim(itf(:,ii,kk,i,j,k,WES))
                itf(:,ii,kk,i,j,k,EAS) = cons2prim(itf(:,ii,kk,i,j,k,EAS))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                itf(:,ii,jj,i,j,k,FRO) = cons2prim(itf(:,ii,jj,i,j,k,FRO))
                itf(:,ii,jj,i,j,k,BAC) = cons2prim(itf(:,ii,jj,i,j,k,BAC))
            end do; end do;
        end do; end do; end do;
# endif
    end if !! (DGFV_FLUX_DIFFERENCING .and. DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS)

    !! ------------------------------------------------------------------ !!
    !! Interpolate to midpoints for indicator calculation.

    if (DGFV_FV_ACTIVE) then

# if NDIM == 1

        imp = itf(jumpvars,:,:,:,:,:,:)

# elif NDIM == 2
        imp = 0.0

        !! Interpolate to midpoints.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                imp(h,jj,kk,i,j,k,NOR) = imp(h,jj,kk,i,j,k,NOR) + toMidpointsMat(jj,ll)*itf(jumpvars(h),ll,kk,i,j,k,NOR)
                imp(h,jj,kk,i,j,k,SOU) = imp(h,jj,kk,i,j,k,SOU) + toMidpointsMat(jj,ll)*itf(jumpvars(h),ll,kk,i,j,k,SOU)
            end do; end do; end do; end do;
        end do; end do; end do;

        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                imp(h,ii,kk,i,j,k,WES) = imp(h,ii,kk,i,j,k,WES) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,kk,i,j,k,WES)
                imp(h,ii,kk,i,j,k,EAS) = imp(h,ii,kk,i,j,k,EAS) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,kk,i,j,k,EAS)
            end do; end do; end do; end do;
        end do; end do; end do;

# elif NDIM == 3
        imp = 0.0

        !! Reconstruction in y-direction.
        atmp = 0.0; btmp = 0.0;
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                atmp(h,1,jj,kk,i,j,k) = atmp(h,1,jj,kk,i,j,k) + toMidpointsMat(jj,ll)*itf(jumpvars(h),ll,kk,i,j,k,NOR)
                btmp(h,1,jj,kk,i,j,k) = btmp(h,1,jj,kk,i,j,k) + toMidpointsMat(jj,ll)*itf(jumpvars(h),ll,kk,i,j,k,SOU)
            end do; end do; end do; end do;
        end do; end do; end do;

        !! Reconstruction in z-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                imp(h,jj,kk,i,j,k,NOR) = imp(h,jj,kk,i,j,k,NOR) + toMidpointsMat(kk,ll)*atmp(h,1,jj,ll,i,j,k)
                imp(h,jj,kk,i,j,k,SOU) = imp(h,jj,kk,i,j,k,SOU) + toMidpointsMat(kk,ll)*btmp(h,1,jj,ll,i,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;

        !! Reconstruction in x-direction.
        atmp = 0.0; btmp = 0.0;
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                atmp(h,ii,1,kk,i,j,k) = atmp(h,ii,1,kk,i,j,k) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,kk,i,j,k,WES)
                btmp(h,ii,1,kk,i,j,k) = btmp(h,ii,1,kk,i,j,k) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,kk,i,j,k,EAS)
            end do; end do; end do; end do;
        end do; end do; end do;

        !! Reconstruction in z-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                imp(h,ii,kk,i,j,k,WES) = imp(h,ii,kk,i,j,k,WES) + toMidpointsMat(kk,ll)*atmp(h,ii,1,ll,i,j,k)
                imp(h,ii,kk,i,j,k,EAS) = imp(h,ii,kk,i,j,k,EAS) + toMidpointsMat(kk,ll)*btmp(h,ii,1,ll,i,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;

        !! Reconstruction in x-direction.
        atmp = 0.0; btmp = 0.0;
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                atmp(h,ii,jj,1,i,j,k) = atmp(h,ii,jj,1,i,j,k) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,jj,i,j,k,FRO)
                btmp(h,ii,jj,1,i,j,k) = btmp(h,ii,jj,1,i,j,k) + toMidpointsMat(ii,ll)*itf(jumpvars(h),ll,jj,i,j,k,BAC)
            end do; end do; end do; end do;
        end do; end do; end do;

        !! Reconstruction in y-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,size(jumpvars);
                imp(h,ii,jj,i,j,k,FRO) = imp(h,ii,jj,i,j,k,FRO) + toMidpointsMat(jj,ll)*atmp(h,ii,ll,1,i,j,k)
                imp(h,ii,jj,i,j,k,BAC) = imp(h,ii,jj,i,j,k,BAC) + toMidpointsMat(jj,ll)*btmp(h,ii,ll,1,i,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif

        !! Calculate jumps.
# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN;
                jumpDG(:,jj,kk,i,j,k,IAXIS) = (imp(:,jj,kk,i,j,k,NOR) - imp(:,jj,kk,i,j,k,SOU)) / &
                                                max(1.0,abs(imp(:,jj,kk,i,j,k,NOR)),abs(imp(:,jj,kk,i,j,k,SOU)))
                jumpFV(:,jj,kk,i,j,k,IAXIS) = (ufv(jumpvars,1,jj,kk,i+1,j,k) - ufv(jumpvars,NXN,jj,kk,i,j,k)) / &
                                                max(1.0,abs(imp(:,jj,kk,i,j,k,NOR)),abs(imp(:,jj,kk,i,j,k,SOU)))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                jumpDG(:,ii,kk,i,j,k,JAXIS) = (imp(:,ii,kk,i,j,k,WES) - imp(:,ii,kk,i,j,k,EAS)) / &
                                                max(1.0,abs(imp(:,ii,kk,i,j,k,WES)),abs(imp(:,ii,kk,i,j,k,EAS)))
                jumpFV(:,ii,kk,i,j,k,JAXIS) = (ufv(jumpvars,ii,1,kk,i,j+1,k) - ufv(jumpvars,ii,NYN,kk,i,j,k)) / &
                                                max(1.0,abs(imp(:,ii,kk,i,j,k,WES)),abs(imp(:,ii,kk,i,j,k,EAS)))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                jumpDG(:,ii,jj,i,j,k,KAXIS) = (imp(:,ii,jj,i,j,k,FRO) - imp(:,ii,jj,i,j,k,BAC)) / &
                                                max(1.0,abs(imp(:,ii,jj,i,j,k,FRO)),abs(imp(:,ii,jj,i,j,k,BAC)))
                jumpFV(:,ii,jj,i,j,k,KAXIS) = (ufv(jumpvars,ii,jj,1,i,j,k+1) - ufv(jumpvars,ii,jj,NZN,i,j,k)) / &
                                                max(1.0,abs(imp(:,ii,jj,i,j,k,FRO)),abs(imp(:,ii,jj,i,j,k,BAC)))
            end do; end do;
        end do; end do; end do;
# endif

        sfblend = 1.0

        !! Calculate blending.
# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            sfblend(i,j,k,IAXIS) = calc_sfblend(jumpFV(:,:,:,i,j,k,IAXIS),jumpDG(:,:,:,i,j,k,IAXIS))
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            sfblend(i,j,k,JAXIS) = calc_sfblend(jumpFV(:,:,:,i,j,k,JAXIS),jumpDG(:,:,:,i,j,k,JAXIS))
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            sfblend(i,j,k,KAXIS) = calc_sfblend(jumpFV(:,:,:,i,j,k,KAXIS),jumpDG(:,:,:,i,j,k,KAXIS))
        end do; end do; end do;
# endif

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            ok = .true.
            do kk = 1,NZN; do jj = 1,NYN;
                ok = ok .and. isvalid_prim(itf(:,jj,kk,i,j,k,SOU)) .and. isvalid_prim(itf(:,jj,kk,i,j,k,NOR))
            end do; end do;
            if (.not. ok) sfblend(i,j,k,IAXIS) = 0.0
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            ok = .true.
            do kk = 1,NZN; do ii = 1,NXN;
                ok = ok .and. isvalid_prim(itf(:,ii,kk,i,j,k,EAS)) .and. isvalid_prim(itf(:,ii,kk,i,j,k,WES))
            end do; end do;
            if (.not. ok) sfblend(i,j,k,JAXIS) = 0.0
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            ok = .true.
            do jj = 1,NYN; do ii = 1,NXN;
                ok = ok .and. isvalid_prim(itf(:,ii,jj,i,j,k,BAC)) .and. isvalid_prim(itf(:,ii,jj,i,j,k,FRO))
            end do; end do;
            if (.not. ok) sfblend(i,j,k,KAXIS) = 0.0
        end do; end do; end do;
# endif

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            vlblend(i,j,k) = min(vlblend(i,j,k), &
# if NDIM > 0
                 0.5*(sfblend(i-1,j,k,IAXIS) + sfblend(i,j,k,IAXIS)) &
# endif
# if NDIM > 1
                ,0.5*(sfblend(i,j-1,k,JAXIS) + sfblend(i,j,k,JAXIS)) &
# endif
# if NDIM > 2
                ,0.5*(sfblend(i,j,k-1,KAXIS) + sfblend(i,j,k,KAXIS)) &
# endif
            )
        ;end do; end do; end do;

    else

        sfblend = 1.0
        vlblend = 1.0

    end if !! HYDRO_FV_ACTIVE

    !! ------------------------------------------------------------------ !!
    !! Calculate surface fluxes.

    rim = 0.0

# if NDIM > 0
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        if (sfblend(i,j,k,IAXIS) > 0.0) then
            do kk = 1,NZN; do jj = 1,NYN;
                rim(:,jj,kk,i,j,k,IAXIS) = xriemann(itf(:,jj,kk,i,j,k,SOU),itf(:,jj,kk,i,j,k,NOR))
            end do; end do;
        end if
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        if (sfblend(i,j,k,JAXIS) > 0.0) then
            do kk = 1,NZN; do ii = 1,NXN;
                rim(:,ii,kk,i,j,k,JAXIS) = yriemann(itf(:,ii,kk,i,j,k,EAS),itf(:,ii,kk,i,j,k,WES))
            end do; end do;
        end if
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        if (sfblend(i,j,k,KAXIS) > 0.0) then
            do jj = 1,NYN; do ii = 1,NXN;
                rim(:,ii,jj,i,j,k,KAXIS) = zriemann(itf(:,ii,jj,i,j,k,BAC),itf(:,ii,jj,i,j,k,FRO))
            end do; end do;
        end if
    end do; end do; end do;
# endif

    !! ------------------------------------------------------------------ !!
    !! Transform to median space.

    flux = 0.0

# if NDIM == 1
    flux(:, :,:, :,:,:, IAXIS) = rim(:, :,:, :,:,:, IAXIS)

# elif NDIM == 2
    !! Project in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            flux(h,jj,kk,i,j,k,IAXIS) = flux(h,jj,kk,i,j,k,IAXIS) + projMat(jj,ll)*rim(h,ll,kk,i,j,k,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Project in x-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            flux(h,ii,kk,i,j,k,JAXIS) = flux(h,ii,kk,i,j,k,JAXIS) + projMat(ii,ll)*rim(h,ll,kk,i,j,k,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

# elif NDIM == 3
    !! Reconstruction in y-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,1,jj,kk,i,j,k) = atmp(h,1,jj,kk,i,j,k) + projMat(jj,ll)*rim(h,ll,kk,i,j,k,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruction in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            flux(h,jj,kk,i,j,k,IAXIS) = flux(h,jj,kk,i,j,k,IAXIS) + projMat(kk,ll)*atmp(h,1,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruction in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,1,kk,i,j,k) = atmp(h,ii,1,kk,i,j,k) + projMat(ii,ll)*rim(h,ll,kk,i,j,k,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruction in z-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            flux(h,ii,kk,i,j,k,JAXIS) = flux(h,ii,kk,i,j,k,JAXIS) + projMat(kk,ll)*atmp(h,ii,1,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruction in x-direction.
    atmp = 0.0
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,jj,1,i,j,k) = atmp(h,ii,jj,1,i,j,k) + projMat(ii,ll)*rim(h,ll,jj,i,j,k,KAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruction in y-direction.
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            flux(h,ii,jj,i,j,k,KAXIS) = flux(h,ii,jj,i,j,k,KAXIS) + projMat(jj,ll)*atmp(h,ii,ll,1,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

    contains

# if 1
    pure function calc_sfblend(jumpFV,jumpDG) result(blend)

        implicit none

        real, intent(in)    :: jumpFV(size(jumpvars),NYN,NZN)
        real, intent(in)    :: jumpDG(size(jumpvars),NYN,NZN)
        real                :: blend

        !! DGFV4
        ! real, parameter :: tauA = 100.0
        ! real, parameter :: tauB = 1.0

        real, parameter :: tauA = 20.0
        real, parameter :: tauB = 1.0

        real :: jvFV(NYN,NZN)
        real :: jvDG(NYN,NZN)
        real :: jabs(NYN,NZN)

        integer :: j,k,h

        blend = 1.0

        do h = 1,size(jumpvars)
            do k = 1,NZN; do j = 1,NYN;
                jvFV(j,k) = jumpFV(h,j,k)
                jvDG(j,k) = jumpDG(h,j,k)
            end do; end do;

            jabs = abs(jvFV)

            blend = min(blend,minval(1.0 - max(0.0,min(tauA*(abs(jvFV - jvDG) - tauB*jabs)/max(jabs,1.0),1.0))))
            if (blend > 0.999) blend = 1.0
        end do

    end function
# else
    function calc_sfblend(jumpFV,jumpDG) result(blend)

        implicit none

        real, intent(in)    :: jumpFV(size(jumpvars),NYN,NZN)
        real, intent(in)    :: jumpDG(size(jumpvars),NYN,NZN)
        real                :: blend

        !! XXX
        !real, parameter :: rmin = 0.0
        !real, parameter :: rmax = 2.0
        !real, parameter :: gamb = 1.0

        ! This is somewhat blunt. There must be a better solution.
        ! (Maybe use std of var?)
        ! The standard sfblend also does this, though.
        !real, parameter :: jmin = 0.1

        ! PARAMETER dgfv_blend_asymptotic     BOOLEAN    TRUE
        ! PARAMETER dgfv_blend_clampmode      INTEGER    0
        ! PARAMETER dgfv_blend_as_jmin        REAL       0.1
        ! PARAMETER dgfv_blend_as_rmin        REAL       0.0
        ! PARAMETER dgfv_blend_as_rmax        REAL       2.0
        ! PARAMETER dgfv_blend_as_gamb        REAL       1.0

        ! real, PARAMETER dgfv_blend_asymptotic     BOOLEAN    TRUE
        integer, PARAMETER :: dgfv_blend_clampmode      =   0
        real, PARAMETER :: dgfv_blend_as_jmin        =       0.1
        real, PARAMETER :: dgfv_blend_as_rmin        =       0.0
        real, PARAMETER :: dgfv_blend_as_rmax        =       2.0
        real, PARAMETER :: dgfv_blend_as_gamb        =       1.0

        integer, parameter  :: clampmode = dgfv_blend_clampmode

        real    :: jvFV(NYN,NZN)
        real    :: jvDG(NYN,NZN)
        real    :: jK(NYN,NZN)
        real    :: alpha_rmin(NYN,NZN)
        real    :: alpha_rmax(NYN,NZN)
        real    :: alpha_side(NYN,NZN)
        logical :: alpha_mask(NYN,NZN)
        real :: alpha(NYN,NZN)

        integer :: j,k,h

        blend = 1.0

        do h = 1,size(jumpvars)
            do k = 1,NZN; do j = 1,NYN;
                jvFV(j,k) = jumpFV(h,j,k)
                jvDG(j,k) = jumpDG(h,j,k)
            end do; end do;

            jK = sign(max(abs(jvFV),dgfv_blend_as_jmin)/max(abs(jvDG-jvFV),1e-10), jvFV*(jvDG-jvFV))
            alpha_rmin = jK*(dgfv_blend_as_rmin-1.0)
            alpha_rmax = jK*(dgfv_blend_as_rmax-1.0)
            alpha_mask = jK<0.0
            alpha_side = merge(alpha_rmin, alpha_rmax, alpha_mask)
            alpha = min(1.0, max(0.0, alpha_side))
            blend = min(blend, minval(alpha)**dgfv_blend_as_gamb)

            if ((blend > 0.9).AND.BTEST(clampmode,0)) blend = 1.0
            if ((blend < 0.1).AND.BTEST(clampmode,1)) blend = 0.1

        end do

    end function
# endif

end subroutine

!! ================================================================== !!

subroutine calc_dg_rhs_weak_form(state,fvflux,blend,sdx,rhs,dgerhs,densflux)

    use dgfv_fluxes_mod, only: xflux => xflux_prim
    use dgfv_fluxes_mod, only: yflux => yflux_prim
    use dgfv_fluxes_mod, only: zflux => zflux_prim

    use Hydro_data, only: recoMat   !! mean values -> node values reconstruction matrix
    use Hydro_data, only: projMat   !! node values -> mean values projection matrix

    use Hydro_data, only: toBoundaryVecM  !! interpolation operator (left/minus face)
    use Hydro_data, only: toBoundaryVecP  !! interpolation operator (right/plus face)

    use Hydro_data, only: weakdiffMat   !! weak-form differentation matrix
    use Hydro_data, only: surfVecM      !! surface (vector) operator (left/minus face)
    use Hydro_data, only: surfVecP      !! surface (vector) operator (right/plus face)

    use Hydro_data, only: dgvolumes
    use Hydro_data, only: DGFV_FV_ACTIVE
    use Hydro_data, only: DGFV_MHD_NON_CONS
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION
    use Hydro_data, only: DGFV_ENTROPY_BOUNDARY_PROJECTION
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION_ACTIVE
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION_CUTOFF
    use Hydro_data, only: DGFV_NODE_TYPE
    use Hydro_data, only: DGFV_NODE_TYPE_GAUSS
    use dgfv_fluxes_mod, only: cons2prim

    use dgfv_fluxes_mod, only: cons2evec
    use dgfv_fluxes_mod, only: evec2cons 
    ! use dgfv_fluxes_mod, only: prim2evec

    ! use dgfv_fluxes_mod, only: cons2xpot
    ! use dgfv_fluxes_mod, only: cons2ypot
    ! use dgfv_fluxes_mod, only: cons2zpot

    use dgfv_fluxes_mod, only: xentrflux_cons
    use dgfv_fluxes_mod, only: yentrflux_cons
    use dgfv_fluxes_mod, only: zentrflux_cons

    use dgfv_fluxes_mod, only: isvalid_cons

    implicit none

    !!                     Conservative state variables.
    real, intent(in)    :: state(NCONS_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real, intent(in)    :: fvflux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
    real, intent(inout) :: blend(NXE,NYE,NZE)
    real, intent(in)    :: sdx(NDIM)

    real, intent(out)   :: rhs(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real, intent(out)   :: dgerhs(NXE,NYE,NZE)
    real, intent(out)   :: densflux(NXN,NYN,NZN, NXE,NYE,NZE,IAXIS:KAXIS)

    real :: ufv (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: udg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: pdg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: atmp(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: btmp(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: rim(NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

    real :: xfx(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: yfy(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: zfz(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real :: vfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: sfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real:: dgrhs(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    integer :: i,j,k,h, ii,jj,kk,ll

    !! ------------------------------------------------------------------ !!
    !! Split block into DG elements.

    do k = NZE_LO,NZE_HI; do j = NYE_LO,NYE_HI; do i = NXE_LO,NXE_HI;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
            ufv(:,ii,jj,kk,i,j,k) = state(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*(k-1)+kk)
        end do; end do; end do
    end do; end do; end do

    !! ------------------------------------------------------------------ !!
    !! Transform to nodal space.

# if NDIM == 1
# define atmp udg
# elif NDIM == 2
# define btmp udg
# elif NDIM == 3
# define ctmp udg
# endif

# if NDIM > 0
    atmp = 0.0
# endif
# if NDIM > 1
    btmp = 0.0
# endif
# if NDIM > 2
    ctmp = 0.0
# endif

!! Inner elements.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,j,k) = atmp(h,ii,jj,kk,i,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,j,k) = btmp(h,ii,jj,kk,i,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,j,k) = ctmp(h,ii,jj,kk,i,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM > 0
!! NORTH side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,NXE_LO,j,k) = atmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,NXE_LO,j,k) = btmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,NXE_LO,j,k) = ctmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! SOUTH side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,NXE_HI,j,k) = atmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,NXE_HI,j,k) = btmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,NXE_HI,j,k) = ctmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 0 */

# if NDIM > 1
!! WEST side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,NYE_LO,k) = atmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,NYE_LO,k) = btmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,NYE_LO,k) = ctmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! EAST side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,NYE_HI,k) = atmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,NYE_HI,k) = btmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,NYE_HI,k) = ctmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 1 */

# if NDIM > 2
!! FRONT side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,j,NZE_LO) = atmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,j,NZE_LO) = btmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,j,NZE_LO) = ctmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! BACK side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,j,NZE_HI) = atmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,j,NZE_HI) = btmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,j,NZE_HI) = ctmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 2 */

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

    !! ------------------------------------------------------------------ !!
    !! Translate to primitive variables.

    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN
            pdg(:,ii,jj,kk,i,j,k) = cons2prim(udg(:,ii,jj,kk,i,j,k))
        end do; end do; end do;
    end do; end do; end do;

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

    !! ------------------------------------------------------------------ !!
    !! Calculate volume flux.

    !! volume flux
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
            xfx(:,ii,jj,kk,i,j,k) = xflux(pdg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 1
            yfy(:,ii,jj,kk,i,j,k) = yflux(pdg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 2
            zfz(:,ii,jj,kk,i,j,k) = zflux(pdg(:,ii,jj,kk,i,j,k))
# endif
        end do; end do; end do;
    end do; end do; end do;

    !! ------------------------------------------------------------------ !!
    !! Reconstruct surface flux.

    rim = 0.0

# if NDIM == 1
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = fvflux(h,NXN*i,j,k,IAXIS)
        end do; end do; end do;
    end do; end do; end do;

# elif NDIM == 2
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = rim(h,jj,kk,i,j,k,IAXIS) + recoMat(jj,ll)*fvflux(h,NXN*i,NYN*(j-1)+ll,NZN*(k-1)+kk,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,kk,i,j,k,JAXIS) = rim(h,ii,kk,i,j,k,JAXIS) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*j,NZN*(k-1)+kk,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

# elif NDIM == 3
    !! Reconstruct in y-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,1,jj,kk,i,j,k) = atmp(h,1,jj,kk,i,j,k) + recoMat(jj,ll)*fvflux(h,NXN*i,NYN*(j-1)+ll,NZN*(k-1)+kk,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = rim(h,jj,kk,i,j,k,IAXIS) + recoMat(kk,ll)*atmp(h,1,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruct in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,1,kk,i,j,k) = atmp(h,ii,1,kk,i,j,k) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*j,NZN*(k-1)+kk,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,kk,i,j,k,JAXIS) = rim(h,ii,kk,i,j,k,JAXIS) + recoMat(kk,ll)*atmp(h,ii,1,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruct in x-direction.
    atmp = 0.0
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,jj,1,i,j,k) = atmp(h,ii,jj,1,i,j,k) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*(j-1)+jj,NZN*k,KAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in y-direction.
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,jj,i,j,k,KAXIS) = rim(h,ii,jj,i,j,k,KAXIS) + recoMat(jj,ll)*atmp(h,ii,ll,1,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

    !! ------------------------------------------------------------------ !!
    !! Compute RHS, i.e. apply the weak-form DG operator.

    !! Apply volume operator.
    vfl = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            vfl(h,ii,jj,kk,i,j,k) = vfl(h,ii,jj,kk,i,j,k) &
# if NDIM > 0
                                        + sdx(IAXIS)*weakdiffMat(ii,ll)*xfx(h,ll,jj,kk,i,j,k) &
# endif
# if NDIM > 1
                                        + sdx(JAXIS)*weakdiffMat(jj,ll)*yfy(h,ii,ll,kk,i,j,k) &
# endif
# if NDIM > 2
                                        + sdx(KAXIS)*weakdiffMat(kk,ll)*zfz(h,ii,jj,ll,i,j,k) &
# endif
        ;end do; end do; end do; end do; end do;
    end do; end do; end do;

    !! Apply surface operator.
    sfl = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do h = 1,NFLUXES;
            sfl(h,ii,jj,kk,i,j,k) = sfl(h,ii,jj,kk,i,j,k) &
# if NDIM > 0
                                        + sdx(IAXIS)*surfVecM(ii)*rim(h,jj,kk,i-1,j,k,IAXIS) &
                                        + sdx(IAXIS)*surfVecP(ii)*rim(h,jj,kk,i  ,j,k,IAXIS) &
# endif
# if NDIM > 1
                                        + sdx(JAXIS)*surfVecM(jj)*rim(h,ii,kk,i,j-1,k,JAXIS) &
                                        + sdx(JAXIS)*surfVecP(jj)*rim(h,ii,kk,i,j  ,k,JAXIS) &
# endif
# if NDIM > 2
                                        + sdx(KAXIS)*surfVecM(kk)*rim(h,ii,jj,i,j,k-1,KAXIS) &
                                        + sdx(KAXIS)*surfVecP(kk)*rim(h,ii,jj,i,j,k  ,KAXIS) &
# endif
        ;end do; end do; end do; end do;
    end do; end do; end do;

    !! ------------------------------------------------------------------ !!
    !! Collect all flux contributions.

    dgrhs = vfl + sfl

    !! ------------------------------------------------------------------ !!
    !! Interpolate to element boundaries.

    if (DGFV_MHD_NON_CONS) then

        block
        real :: means(NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
        real :: itf(NCONS_VARS, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)

        real :: divB(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: divG(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: ncs (NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

        itf = 0.0

# if NDIM > 0
        !! Interpolate to boundaries in x-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES;
                itf(MAGX_PRIM,jj,kk,i,j,k,NOR) = itf(MAGX_PRIM,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*udg(MAGX_CONS,ll,jj,kk,i+1,j,k)
                itf(MAGX_PRIM,jj,kk,i,j,k,SOU) = itf(MAGX_PRIM,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*udg(MAGX_CONS,ll,jj,kk,i  ,j,k)
                itf(GLMP_PRIM,jj,kk,i,j,k,NOR) = itf(GLMP_PRIM,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*udg(GLMP_CONS,ll,jj,kk,i+1,j,k)
                itf(GLMP_PRIM,jj,kk,i,j,k,SOU) = itf(GLMP_PRIM,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*udg(GLMP_CONS,ll,jj,kk,i  ,j,k)
            end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        !! Interpolate to boundaries in y-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES;
                itf(MAGY_PRIM,ii,kk,i,j,k,WES) = itf(MAGY_PRIM,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*udg(MAGY_CONS,ii,ll,kk,i,j+1,k)
                itf(MAGY_PRIM,ii,kk,i,j,k,EAS) = itf(MAGY_PRIM,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*udg(MAGY_CONS,ii,ll,kk,i,j  ,k)
                itf(GLMP_PRIM,ii,kk,i,j,k,WES) = itf(GLMP_PRIM,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*udg(GLMP_CONS,ii,ll,kk,i,j+1,k)
                itf(GLMP_PRIM,ii,kk,i,j,k,EAS) = itf(GLMP_PRIM,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*udg(GLMP_CONS,ii,ll,kk,i,j  ,k)
            end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        !! Interpolate to boundaries in z-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
                itf(MAGZ_PRIM,ii,jj,i,j,k,FRO) = itf(MAGZ_PRIM,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*udg(MAGZ_CONS,ii,jj,ll,i,j,k+1)
                itf(MAGZ_PRIM,ii,jj,i,j,k,BAC) = itf(MAGZ_PRIM,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*udg(MAGZ_CONS,ii,jj,ll,i,j,k  )
                itf(GLMP_PRIM,ii,jj,i,j,k,FRO) = itf(GLMP_PRIM,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*udg(GLMP_CONS,ii,jj,ll,i,j,k+1)
                itf(GLMP_PRIM,ii,jj,i,j,k,BAC) = itf(GLMP_PRIM,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*udg(GLMP_CONS,ii,jj,ll,i,j,k  )
            end do; end do; end do;
        end do; end do; end do;
# endif

        !! ------------------------------------------------------------------ !!
        !! Calculate surface averages.

        means = 0.0
# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                means(MAGX_FLUX,jj,kk,i,j,k,IAXIS) = 0.5*(itf(MAGX_PRIM,jj,kk,i,j,k,SOU) + itf(MAGX_PRIM,jj,kk,i,j,k,NOR))
                means(GLMP_FLUX,jj,kk,i,j,k,IAXIS) = 0.5*(itf(GLMP_PRIM,jj,kk,i,j,k,SOU) + itf(GLMP_PRIM,jj,kk,i,j,k,NOR))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                means(MAGY_FLUX,ii,kk,i,j,k,JAXIS) = 0.5*(itf(MAGY_PRIM,ii,kk,i,j,k,EAS) + itf(MAGY_PRIM,ii,kk,i,j,k,WES))
                means(GLMP_FLUX,ii,kk,i,j,k,JAXIS) = 0.5*(itf(GLMP_PRIM,ii,kk,i,j,k,EAS) + itf(GLMP_PRIM,ii,kk,i,j,k,WES))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                means(MAGZ_FLUX,ii,jj,i,j,k,KAXIS) = 0.5*(itf(MAGZ_PRIM,ii,jj,i,j,k,BAC) + itf(MAGZ_PRIM,ii,jj,i,j,k,FRO))
                means(GLMP_FLUX,ii,jj,i,j,k,KAXIS) = 0.5*(itf(GLMP_PRIM,ii,jj,i,j,k,BAC) + itf(GLMP_PRIM,ii,jj,i,j,k,FRO))
            end do; end do;
        end do; end do; end do;
# endif

        !! Apply volume operator: magnetic field divergence
        divB = 0.0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
                divB(ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k) &
# if NDIM > 0
                + sdx(IAXIS)*weakdiffMat(ii,ll)*pdg(MAGX_PRIM,ll,jj,kk,i,j,k) &
# endif
# if NDIM > 1
                + sdx(JAXIS)*weakdiffMat(jj,ll)*pdg(MAGY_PRIM,ii,ll,kk,i,j,k) &
# endif
# if NDIM > 2
                + sdx(KAXIS)*weakdiffMat(kk,ll)*pdg(MAGZ_PRIM,ii,jj,ll,i,j,k) &
# endif
            ;end do; end do; end do; end do;
        end do; end do; end do;

        !! Apply surface operator: magnetic field divergence
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                divB(ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k) &
# if NDIM > 0
                + sdx(IAXIS)*surfVecM(ii)*means(MAGX_FLUX,jj,kk,i-1,j,k,IAXIS) &
                + sdx(IAXIS)*surfVecP(ii)*means(MAGX_FLUX,jj,kk,i  ,j,k,IAXIS) &
# endif
# if NDIM > 1
                + sdx(JAXIS)*surfVecM(jj)*means(MAGY_FLUX,ii,kk,i,j-1,k,JAXIS) &
                + sdx(JAXIS)*surfVecP(jj)*means(MAGY_FLUX,ii,kk,i,j  ,k,JAXIS) &
# endif
# if NDIM > 2
                + sdx(KAXIS)*surfVecM(kk)*means(MAGZ_FLUX,ii,jj,i,j,k-1,KAXIS) &
                + sdx(KAXIS)*surfVecP(kk)*means(MAGZ_FLUX,ii,jj,i,j,k  ,KAXIS) &
# endif
            ;end do; end do; end do;
        end do; end do; end do;

        !! ------------------------------------------------------------------ !!

        !! Apply volume operator: GLM gradient.
        divG = 0.0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
                divG(ii,jj,kk,i,j,k) = divG(ii,jj,kk,i,j,k) &
# if NDIM > 0
                + sdx(IAXIS)*weakdiffMat(ii,ll)*pdg(GLMP_PRIM,ll,jj,kk,i,j,k)*pdg(VELX_PRIM,ii,jj,kk,i,j,k) &
# endif
# if NDIM > 1
                + sdx(JAXIS)*weakdiffMat(jj,ll)*pdg(GLMP_PRIM,ii,ll,kk,i,j,k)*pdg(VELY_PRIM,ii,jj,kk,i,j,k) &
# endif
# if NDIM > 2
                + sdx(KAXIS)*weakdiffMat(kk,ll)*pdg(GLMP_PRIM,ii,jj,ll,i,j,k)*pdg(VELZ_PRIM,ii,jj,kk,i,j,k) &
# endif
            ;end do; end do; end do; end do;
        end do; end do; end do;

        !! Apply surface operator: GLM gradient.
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                divG(ii,jj,kk,i,j,k) = divG(ii,jj,kk,i,j,k) &
# if NDIM > 0
                + sdx(IAXIS)*surfVecM(ii)*means(GLMP_FLUX,jj,kk,i-1,j,k,IAXIS)*pdg(VELX_PRIM,ii,jj,kk,i,j,k) &
                + sdx(IAXIS)*surfVecP(ii)*means(GLMP_FLUX,jj,kk,i  ,j,k,IAXIS)*pdg(VELX_PRIM,ii,jj,kk,i,j,k) &
# endif
# if NDIM > 1
                + sdx(JAXIS)*surfVecM(jj)*means(GLMP_FLUX,ii,kk,i,j-1,k,JAXIS)*pdg(VELY_PRIM,ii,jj,kk,i,j,k) &
                + sdx(JAXIS)*surfVecP(jj)*means(GLMP_FLUX,ii,kk,i,j  ,k,JAXIS)*pdg(VELY_PRIM,ii,jj,kk,i,j,k) &
# endif
# if NDIM > 2
                + sdx(KAXIS)*surfVecM(kk)*means(GLMP_FLUX,ii,jj,i,j,k-1,KAXIS)*pdg(VELZ_PRIM,ii,jj,kk,i,j,k) &
                + sdx(KAXIS)*surfVecP(kk)*means(GLMP_FLUX,ii,jj,i,j,k  ,KAXIS)*pdg(VELZ_PRIM,ii,jj,kk,i,j,k) &
# endif
            ;end do; end do; end do;
        end do; end do; end do;

        !! Calculate non-conservative terms.
        ncs = 0.0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                ncs(XMOM_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(MAGX_PRIM,ii,jj,kk,i,j,k)
                ncs(YMOM_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(MAGY_PRIM,ii,jj,kk,i,j,k)
                ncs(ZMOM_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(MAGZ_PRIM,ii,jj,kk,i,j,k)

                ncs(ENER_FLUX,ii,jj,kk,i,j,k) &
                    = divB(ii,jj,kk,i,j,k)*SUM(pdg(VELX_PRIM:VELZ_PRIM,ii,jj,kk,i,j,k)*pdg(MAGX_PRIM:MAGZ_PRIM,ii,jj,kk,i,j,k)) &
                    + divG(ii,jj,kk,i,j,k)*pdg(GLMP_PRIM,ii,jj,kk,i,j,k)

                ncs(MAGX_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(VELX_PRIM,ii,jj,kk,i,j,k)
                ncs(MAGY_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(VELY_PRIM,ii,jj,kk,i,j,k)
                ncs(MAGZ_FLUX,ii,jj,kk,i,j,k) = divB(ii,jj,kk,i,j,k)*pdg(VELZ_PRIM,ii,jj,kk,i,j,k)

                ncs(GLMP_FLUX,ii,jj,kk,i,j,k) = divG(ii,jj,kk,i,j,k)
            end do; end do; end do;
        end do; end do; end do;

        dgrhs = dgrhs + ncs

        end block
    end if !! DGFV_MHD_NON_CONS

    if (DGFV_ENTROPY_CORRECTION) then

        !! ------------------------------------------------------------------ !!
        !! Calculate entropy production rate per DG (sub-)element.

        block
        real  :: erhs(NXE,NYE,NZE)

        real :: edg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
        real :: sdg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

        real :: highs (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
        real :: means (NCONS_VARS,              NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

        real :: evfl(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: esfl(         NXN,NYN,NZN, NXE,NYE,NZE)

        ! real :: ievec(NCONS_VARS, NYN,NZN,     NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
        real :: vevec(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        !real :: iepot(           NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
        real :: ieflx(           NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

        real :: devec(NPROP_FLUX,  NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        real :: aevec(NPROP_FLUX,  NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        real :: mevec(NPROP_FLUX,               NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        real :: tevec(          NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        real :: tesum(                       NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        !! used for validation
        real :: dgevfl(         NXN,NYN,NZN, NXE,NYE,NZE)

        real :: itf(NCONS_VARS, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)

        real :: squeeze(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        logical :: ok(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE), first

        !! ------------------------------------------------------------------ !!
        !! Interpolate to element boundaries.


        ! if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS.and. DGFV_ENTROPY_BOUNDARY_PROJECTION) then
        if (DGFV_ENTROPY_BOUNDARY_PROJECTION) then
 
            squeeze = 1.0
            first = .true.
            sdg = udg

            do while (all(squeeze > 0.0))

                ok = .true.
                itf = 0.0

                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        edg(:,ii,jj,kk,i,j,k) = cons2evec(sdg(:,ii,jj,kk,i,j,k))
                    end do; end do; end do;
                end do; end do; end do;

# if NDIM > 0
                do k = 1,NZE; do j = 1,NYE
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        edg(:,ii,jj,kk,NXE_LO,j,k) = cons2evec(sdg(:,ii,jj,kk,NXE_LO,j,k))
                        edg(:,ii,jj,kk,NXE_HI,j,k) = cons2evec(sdg(:,ii,jj,kk,NXE_HI,j,k))
                    end do; end do; end do;
                end do; end do
# endif
# if NDIM > 1
                do k = 1,NZE; do i = 1,NXE
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        edg(:,ii,jj,kk,i,NYE_LO,k) = cons2evec(sdg(:,ii,jj,kk,i,NYE_LO,k))
                        edg(:,ii,jj,kk,i,NYE_HI,k) = cons2evec(sdg(:,ii,jj,kk,i,NYE_HI,k))
                    end do; end do; end do;
                end do; end do
# endif
# if NDIM > 2
                do j = 1,NYE; do i = 1,NXE
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        edg(:,ii,jj,kk,i,j,NZE_LO) = cons2evec(sdg(:,ii,jj,kk,i,j,NZE_LO))
                        edg(:,ii,jj,kk,i,j,NZE_HI) = cons2evec(sdg(:,ii,jj,kk,i,j,NZE_HI))
                    end do; end do; end do;
                end do; end do
# endif

# if NDIM > 0
                !! Interpolate to boundaries in x-direction.
                do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                    do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                        itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*edg(h,ll,jj,kk,i+1,j,k)
                        itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*edg(h,ll,jj,kk,i  ,j,k)
                    end do; end do; end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 1
                !! Interpolate to boundaries in y-direction.
                do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                        itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*edg(h,ii,ll,kk,i,j+1,k)
                        itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*edg(h,ii,ll,kk,i,j  ,k)
                    end do; end do; end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 2
                !! Interpolate to boundaries in z-direction.
                do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                    do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                        itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*edg(h,ii,jj,ll,i,j,k+1)
                        itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*edg(h,ii,jj,ll,i,j,k  )
                    end do; end do; end do; end do;
                end do; end do; end do;
# endif

# if NDIM > 0
                do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                    do kk = 1,NZN; do jj = 1,NYN
                        itf(:,jj,kk,i,j,k,NOR) = evec2cons(itf(:,jj,kk,i,j,k,NOR))
                        itf(:,jj,kk,i,j,k,SOU) = evec2cons(itf(:,jj,kk,i,j,k,SOU))
                        if (.not.isvalid_cons(itf(:,jj,kk,i,j,k,NOR))) ok(i,j,k) = .false.
                        if (.not.isvalid_cons(itf(:,jj,kk,i,j,k,SOU))) ok(i,j,k) = .false.
                    end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 1
                do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do ii = 1,NXN;
                        itf(:,ii,kk,i,j,k,WES) = evec2cons(itf(:,ii,kk,i,j,k,WES))
                        itf(:,ii,kk,i,j,k,EAS) = evec2cons(itf(:,ii,kk,i,j,k,EAS))
                        if (.not.isvalid_cons(itf(:,ii,kk,i,j,k,WES))) ok(i,j,k) = .false.
                        if (.not.isvalid_cons(itf(:,ii,kk,i,j,k,EAS))) ok(i,j,k) = .false.
                    end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 2
                do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                    do jj = 1,NYN; do ii = 1,NXN;
                        itf(:,ii,jj,i,j,k,FRO) = evec2cons(itf(:,ii,jj,i,j,k,FRO))
                        itf(:,ii,jj,i,j,k,BAC) = evec2cons(itf(:,ii,jj,i,j,k,BAC))
                        if (.not.isvalid_cons(itf(:,ii,jj,i,j,k,FRO))) ok(i,j,k) = .false.
                        if (.not.isvalid_cons(itf(:,ii,jj,i,j,k,BAC))) ok(i,j,k) = .false.
                    end do; end do;
                end do; end do; end do;
# endif

                exit
                !! if (all(ok)) exit

                !! where (.not.ok)
                !!     squeeze = max(0.0,squeeze - 0.1)
                !! end where

                !! ! write (*,*) 'squeeze', minval(squeeze)

                !! if (first) then
                !!     means = 0.0
                !!     do k = 0,NZE; do j = 0,NYE; do i = 0,NXE;
                !!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                !!             means(:,i,j,k) = means(:,i,j,k) + dgvolumes(ii,jj,kk)*sdg(:,ii,jj,kk,i,j,k)
                !!         end do; end do; end do;
                !!     end do; end do; end do;

                !!     do k = 0,NZE; do j = 0,NYE; do i = 0,NXE;
                !!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                !!             highs(:,ii,jj,kk,i,j,k) = sdg(:,ii,jj,kk,i,j,k) - means(:,i,j,k)
                !!         end do; end do; end do;
                !!     end do; end do; end do;

                !!     first = .false.
                !! end if

                !! do k = 0,NZE; do j = 0,NYE; do i = 0,NXE;
                !!     do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                !!         sdg(:,ii,jj,kk,i,j,k) = means(:,i,j,k) + squeeze(i,j,k)*highs(:,ii,jj,kk,i,j,k)
                !!     end do; end do; end do;
                !! end do; end do; end do;

            end do

        else

            itf = 0.0
# if NDIM > 0
            !! Interpolate to boundaries in x-direction.
            do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                    itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*udg(h,ll,jj,kk,i+1,j,k)
                    itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*udg(h,ll,jj,kk,i  ,j,k)
                end do; end do; end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 1
            !! Interpolate to boundaries in y-direction.
            do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                    itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*udg(h,ii,ll,kk,i,j+1,k)
                    itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*udg(h,ii,ll,kk,i,j  ,k)
                end do; end do; end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 2
            !! Interpolate to boundaries in z-direction.
            do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                    itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*udg(h,ii,jj,ll,i,j,k+1)
                    itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*udg(h,ii,jj,ll,i,j,k  )
                end do; end do; end do; end do;
            end do; end do; end do;
# endif
        end if !! DGFV_NODE_TYPE

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                ieflx(jj,kk,i,j,k,IAXIS) = sdx(IAXIS)*xentrflux_cons(-1,rim(:,jj,kk,i,j,k,IAXIS),itf(:,jj,kk,i,j,k,SOU),itf(:,jj,kk,i,j,k,NOR))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                ieflx(ii,kk,i,j,k,JAXIS) = sdx(JAXIS)*yentrflux_cons( 1,rim(:,ii,kk,i,j,k,JAXIS),itf(:,ii,kk,i,j,k,EAS),itf(:,ii,kk,i,j,k,WES))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                ieflx(ii,jj,i,j,k,KAXIS) = sdx(KAXIS)*zentrflux_cons( 1,rim(:,ii,jj,i,j,k,KAXIS),itf(:,ii,jj,i,j,k,BAC),itf(:,ii,jj,i,j,k,FRO))
            end do; end do;
        end do; end do; end do;
# endif

        !! Apply surface operator.
        esfl = 0.0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                esfl(ii,jj,kk,i,j,k) = esfl(ii,jj,kk,i,j,k) &
# if NDIM > 0
                    + surfVecM(ii)*ieflx(jj,kk,i-1,j,k,IAXIS) &
                    + surfVecP(ii)*ieflx(jj,kk,i  ,j,k,IAXIS) &
# endif
# if NDIM > 1
                    + surfVecM(jj)*ieflx(ii,kk,i,j-1,k,JAXIS) &
                    + surfVecP(jj)*ieflx(ii,kk,i,j  ,k,JAXIS) &
# endif
# if NDIM > 2
                    + surfVecM(kk)*ieflx(ii,jj,i,j,k-1,KAXIS) &
                    + surfVecP(kk)*ieflx(ii,jj,i,j,k  ,KAXIS) &
# endif
            ;end do; end do; end do;
        end do; end do; end do;


        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                vevec(:,ii,jj,kk,i,j,k) = cons2evec(udg(:,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

        if (DGFV_ENTROPY_CORRECTION_ACTIVE) then

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN
                    evfl(ii,jj,kk,i,j,k) = dot_product(vevec(1:NPROP_FLUX,ii,jj,kk,i,j,k),dgrhs(1:NPROP_FLUX,ii,jj,kk,i,j,k))
                end do; end do; end do;
            end do; end do; end do;

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                erhs(i,j,k) = sum(dgvolumes*(evfl(:,:,:,i,j,k) - esfl(:,:,:,i,j,k)))
            end do; end do; end do;

            mevec = 0.0
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                    mevec(:,i,j,k) = mevec(:,i,j,k) + dgvolumes(ii,jj,kk)*vevec(1:NPROP_FLUX,ii,jj,kk,i,j,k)
                end do; end do; end do;
            end do; end do; end do;

            devec = 0.0
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                    devec(:,ii,jj,kk,i,j,k) = vevec(1:NPROP_FLUX,ii,jj,kk,i,j,k) - mevec(:,i,j,k)
                end do; end do; end do;
            end do; end do; end do;

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                    tevec(ii,jj,kk,i,j,k) = dot_product(devec(:,ii,jj,kk,i,j,k),devec(:,ii,jj,kk,i,j,k)) + 1e-20
                end do; end do; end do;
            end do; end do; end do;

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                tesum(i,j,k) = sum(dgvolumes*tevec(:,:,:,i,j,k))
            end do; end do; end do;

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                    aevec(:,ii,jj,kk,i,j,k) = devec(:,ii,jj,kk,i,j,k)/tesum(i,j,k)
                end do; end do; end do;
            end do; end do; end do;

            if (DGFV_FV_ACTIVE) then
                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    if (blend(i,j,k) == 0.0 .or. isnan(erhs(i,j,k)) .or. any(isnan(aevec(:,:,:,:,i,j,k)))) then
                        erhs(i,j,k) = 0.0
                        aevec(:,:,:,:,i,j,k) = 0.0
                    end if
                end do; end do; end do;
            end if

            !! cutoff
            if (DGFV_ENTROPY_CORRECTION_CUTOFF) then
                erhs = max(0.0,erhs)
            end if

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN
                    dgrhs(1:NPROP_FLUX,ii,jj,kk,i,j,k) = dgrhs(1:NPROP_FLUX,ii,jj,kk,i,j,k) - erhs(i,j,k) * aevec(:,ii,jj,kk,i,j,k)
                end do; end do; end do;
            end do; end do; end do;

        end if !! DGFV_ENTROPY_CORRECTION_ACTIVE

        !! For validation.
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN
                dgevfl(ii,jj,kk,i,j,k) = dot_product(vevec(1:NPROP_FLUX,ii,jj,kk,i,j,k),dgrhs(1:NPROP_FLUX,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            dgerhs(i,j,k) = sum(dgvolumes*(dgevfl(:,:,:,i,j,k) - esfl(:,:,:,i,j,k)))
        end do; end do; end do;

        if (DGFV_FV_ACTIVE) then
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                if (blend(i,j,k) == 0.0 .or. isnan(dgerhs(i,j,k))) then
                    dgerhs(i,j,k) = 0.0
                end if
            end do; end do; end do;
        end if
        end block
    end if !! DGFV_ENTROPY_CORRECTION


    if (DGFV_FV_ACTIVE) then
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            if (blend(i,j,k) == 0.0 .or. any(isnan(dgrhs(:,:,:,:,i,j,k))) &
                    .or. any(dgrhs(:,:,:,:,i,j,k)-1.0 == dgrhs(:,:,:,:,i,j,k))) then

                blend(i,j,k) = 0.0
                dgrhs(:,:,:,:,i,j,k) = 0.0
            end if
        end do; end do; end do;
    end if

    !! ------------------------------------------------------------------ !!
    !! Project to median space.

# if NDIM == 1
# define atmp rhs
# elif NDIM == 2
# define btmp rhs
# elif NDIM == 3
# define ctmp rhs
# endif

# if NDIM > 0
    !! Projection in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,jj,kk,i,j,k) = atmp(h,ii,jj,kk,i,j,k) + projMat(ii,ll)*dgrhs(h,ll,jj,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Projection in y-direction.
    btmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            btmp(h,ii,jj,kk,i,j,k) = btmp(h,ii,jj,kk,i,j,k) + projMat(jj,ll)*atmp(h,ii,ll,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Projection in z-direction.
    ctmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            ctmp(h,ii,jj,kk,i,j,k) = ctmp(h,ii,jj,kk,i,j,k) + projMat(kk,ll)*btmp(h,ii,jj,ll,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

    !! ------------------------------------------------------------------ !!
    !! ------------------------------------------------------------------ !!

!     !! volume flux
!     do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
! # if NDIM > 0
!             xfx(:,ii,jj,kk,i,j,k) = xflux(pdg(:,ii,jj,kk,i,j,k))
! # endif
! # if NDIM > 1
!             yfy(:,ii,jj,kk,i,j,k) = yflux(pdg(:,ii,jj,kk,i,j,k))
! # endif
! # if NDIM > 2
!             zfz(:,ii,jj,kk,i,j,k) = zflux(pdg(:,ii,jj,kk,i,j,k))
! # endif
!         end do; end do; end do;
!     end do; end do; end do;




# if NSPECIES > 0 || NMASS_SCALARS > 0
    !! ------------------------------------------------------------------ !!
    !! Project to median space.

    block
    real :: atmp(NXN,NYN,NZN, NXE,NYE,NZE)
    real :: btmp(NXN,NYN,NZN, NXE,NYE,NZE)
    real :: temp(NXN,NYN,NZN, NXE,NYE,NZE)

# if NDIM > 0

# if NDIM == 1
# define atmp temp
# elif NDIM == 2
# define btmp temp
# elif NDIM == 3
# define ctmp temp
# endif

# if NDIM > 0
    !! Projection in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            atmp(ii,jj,kk,i,j,k) = atmp(ii,jj,kk,i,j,k) + projMat(ii,ll)*xfx(DENS_FLUX,ll,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Projection in y-direction.
    btmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            btmp(ii,jj,kk,i,j,k) = btmp(ii,jj,kk,i,j,k) + projMat(jj,ll)*atmp(ii,ll,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Projection in z-direction.
    ctmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            ctmp(ii,jj,kk,i,j,k) = ctmp(ii,jj,kk,i,j,k) + projMat(kk,ll)*btmp(ii,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif
    
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN-1;
            densflux(ii,jj,kk,i,j,k,IAXIS) = 0.5*(temp(ii,jj,kk,i,j,k) + temp(ii+1,jj,kk,i,j,k))
        end do; end do; end do;
    end do; end do; end do;
# endif !! NDIM > 0

# if NDIM > 1

# if NDIM == 1
# define atmp temp
# elif NDIM == 2
# define btmp temp
# elif NDIM == 3
# define ctmp temp
# endif

# if NDIM > 0
    !! Projection in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            atmp(ii,jj,kk,i,j,k) = atmp(ii,jj,kk,i,j,k) + projMat(ii,ll)*yfy(DENS_FLUX,ll,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Projection in y-direction.
    btmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            btmp(ii,jj,kk,i,j,k) = btmp(ii,jj,kk,i,j,k) + projMat(jj,ll)*atmp(ii,ll,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Projection in z-direction.
    ctmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            ctmp(ii,jj,kk,i,j,k) = ctmp(ii,jj,kk,i,j,k) + projMat(kk,ll)*btmp(ii,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif
    
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN-1; do ii = 1,NXN;
            densflux(ii,jj,kk,i,j,k,JAXIS) = 0.5*(temp(ii,jj,kk,i,j,k) + temp(ii,jj+1,kk,i,j,k))
        end do; end do; end do;
    end do; end do; end do;
# endif !! NDIM > 0

# if NDIM > 2

# if NDIM == 1
# define atmp temp
# elif NDIM == 2
# define btmp temp
# elif NDIM == 3
# define ctmp temp
# endif

# if NDIM > 0
    !! Projection in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            atmp(ii,jj,kk,i,j,k) = atmp(ii,jj,kk,i,j,k) + projMat(ii,ll)*zfz(DENS_FLUX,ll,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Projection in y-direction.
    btmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            btmp(ii,jj,kk,i,j,k) = btmp(ii,jj,kk,i,j,k) + projMat(jj,ll)*atmp(ii,ll,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Projection in z-direction.
    ctmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
            ctmp(ii,jj,kk,i,j,k) = ctmp(ii,jj,kk,i,j,k) + projMat(kk,ll)*btmp(ii,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif
    
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN-1; do jj = 1,NYN; do ii = 1,NXN;
            densflux(ii,jj,kk,i,j,k,KAXIS) = 0.5*(temp(ii,jj,kk,i,j,k) + temp(ii,jj,kk+1,i,j,k))
        end do; end do; end do;
    end do; end do; end do;
# endif !! NDIM > 0

    end block
# else

    densflux = 0.0

# endif !! NSPECIES > 0 || NMASS_SCALARS > 0

end subroutine

!! ================================================================== !!

subroutine calc_dg_rhs_flux_differencing(state,fvflux,blend,sdx,rhs,dgerhs)

    use dgfv_fluxes_mod, only: xflux => xflux_prim
    use dgfv_fluxes_mod, only: yflux => yflux_prim
    use dgfv_fluxes_mod, only: zflux => zflux_prim

    use dgfv_fluxes_mod, only: TWO_POINT_XFLUX_cons
    use dgfv_fluxes_mod, only: TWO_POINT_YFLUX_cons
    use dgfv_fluxes_mod, only: TWO_POINT_ZFLUX_cons

    use dgfv_fluxes_mod, only: non_Conservative_xFlux_cons
    use dgfv_fluxes_mod, only: non_Conservative_yFlux_cons
    use dgfv_fluxes_mod, only: non_Conservative_zFlux_cons

    use Hydro_data, only: recoMat   !! mean values -> node values reconstruction matrix
    use Hydro_data, only: projMat   !! node values -> mean values projection matrix

    use Hydro_data, only: toBoundaryVecM  !! interpolation operator (left/minus face)
    use Hydro_data, only: toBoundaryVecP  !! interpolation operator (right/plus face)

    use Hydro_data, only: weakdiffMat   !! weak differencing matrix
    use Hydro_data, only: fluxdiffMat   !! flux differencing matrix
    use Hydro_data, only: surfVecM      !! surface (vector) operator (left/minus face)
    use Hydro_data, only: surfVecP      !! surface (vector) operator (right/plus face)

    use Hydro_data, only: DGFV_FV_ACTIVE
    use Hydro_data, only: DGFV_MHD_NON_CONS

    use Hydro_data, only: dgvolumes
    use Hydro_data, only: xweights
    use Hydro_data, only: yweights
    use Hydro_data, only: zweights

    use dgfv_fluxes_mod, only: isvalid_cons
    use dgfv_fluxes_mod, only: cons2prim
    use dgfv_fluxes_mod, only: cons2evec
    use dgfv_fluxes_mod, only: evec2cons

    use dgfv_fluxes_mod, only: cons2xpot
    use dgfv_fluxes_mod, only: cons2ypot
    use dgfv_fluxes_mod, only: cons2zpot

    use dgfv_fluxes_mod, only: xefluxmean_cons
    use dgfv_fluxes_mod, only: yefluxmean_cons
    use dgfv_fluxes_mod, only: zefluxmean_cons

    use dgfv_fluxes_mod, only: xentrflux_cons
    use dgfv_fluxes_mod, only: yentrflux_cons
    use dgfv_fluxes_mod, only: zentrflux_cons

    use Hydro_data, only: DGFV_ENTROPY_BOUNDARY_PROJECTION
    use Hydro_data, only: DGFV_NODE_TYPE
    use Hydro_data, only: DGFV_NODE_TYPE_GAUSS

    implicit none

    !!                     Conservative state variables.
    real, intent(in)    :: state(NCONS_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real, intent(in)    :: fvflux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
    real, intent(inout) :: blend(NXE,NYE,NZE)
    real, intent(in)    :: sdx(NDIM)

    real, intent(out)   :: rhs(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real, intent(out)   :: dgerhs(NXE,NYE,NZE)

    real :: ufv (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: udg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: pdg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: atmp(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
    real :: btmp(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

    real :: itf(NCONS_VARS, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)
    real :: rim(NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

    real :: xfx(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: yfy(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: zfz(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real :: xfxm(NFLUXES, N_NODES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: yfym(NFLUXES, N_NODES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: zfzm(NFLUXES, N_NODES, NXN,NYN,NZN, NXE,NYE,NZE)

    real :: xvfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: yvfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: zvfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real :: xsfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: ysfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)
    real :: zsfl(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real :: divB(         NXN,NYN,NZN, NXE,NYE,NZE)
    real :: divG(         NXN,NYN,NZN, NXE,NYE,NZE)
    real :: ncs (NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    real:: dgrhs(NFLUXES, NXN,NYN,NZN, NXE,NYE,NZE)

    integer :: i,j,k,h, ii,jj,kk,ll

    logical :: ok

    !! ------------------------------------------------------------------ !!
    !! Split block into DG elements.

    do k = NZE_LO,NZE_HI; do j = NYE_LO,NYE_HI; do i = NXE_LO,NXE_HI;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
            ufv(:,ii,jj,kk,i,j,k) = state(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*(k-1)+kk)
        end do; end do; end do
    end do; end do; end do

    !! ------------------------------------------------------------------ !!
    !! Reconstruct from mean values to node values.

# if NDIM == 1
# define atmp udg
# elif NDIM == 2
# define btmp udg
# elif NDIM == 3
# define ctmp udg
# endif

# if NDIM > 0
    atmp = 0.0
# endif
# if NDIM > 1
    btmp = 0.0
# endif
# if NDIM > 2
    ctmp = 0.0
# endif

!! Inner elements.
# if NDIM > 0
    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            atmp(h,ii,jj,kk,i,j,k) = atmp(h,ii,jj,kk,i,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            btmp(h,ii,jj,kk,i,j,k) = btmp(h,ii,jj,kk,i,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
            ctmp(h,ii,jj,kk,i,j,k) = ctmp(h,ii,jj,kk,i,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM > 0
!! NORTH side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,NXE_LO,j,k) = atmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,NXE_LO,j,k) = btmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,NXE_LO,j,k) = ctmp(h,ii,jj,kk,NXE_LO,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_LO,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! SOUTH side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,NXE_HI,j,k) = atmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,NXE_HI,j,k) = btmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do j = 1,NYE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,NXE_HI,j,k) = ctmp(h,ii,jj,kk,NXE_HI,j,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,NXE_HI,j,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 0 */

# if NDIM > 1
!! WEST side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,NYE_LO,k) = atmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,NYE_LO,k) = btmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,NYE_LO,k) = ctmp(h,ii,jj,kk,i,NYE_LO,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_LO,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! EAST side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,NYE_HI,k) = atmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,NYE_HI,k) = btmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do k = 1,NZE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,NYE_HI,k) = ctmp(h,ii,jj,kk,i,NYE_HI,k) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,NYE_HI,k)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 1 */

# if NDIM > 2
!! FRONT side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,j,NZE_LO) = atmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,j,NZE_LO) = btmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,j,NZE_LO) = ctmp(h,ii,jj,kk,i,j,NZE_LO) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_LO)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif

!! BACK side.
# if NDIM > 0
        !! Reconstruct in x-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                atmp(h,ii,jj,kk,i,j,NZE_HI) = atmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(ii,ll)*ufv(h,ll,jj,kk,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 1
        !! Reconstruct in y-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                btmp(h,ii,jj,kk,i,j,NZE_HI) = btmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(jj,ll)*atmp(h,ii,ll,kk,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# if NDIM > 2
        !! Reconstruct in z-direction.
        do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                ctmp(h,ii,jj,kk,i,j,NZE_HI) = ctmp(h,ii,jj,kk,i,j,NZE_HI) + recoMat(kk,ll)*btmp(h,ii,jj,ll,i,j,NZE_HI)
            end do; end do; end do; end do; end do;
        end do; end do;
# endif
# endif /* NDIM > 2 */

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

    !! ------------------------------------------------------------------ !!
    !! Interpolate to element boundaries.

    if (DGFV_ENTROPY_BOUNDARY_PROJECTION) then
    ! if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS) then

        block
        real :: edg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,j,k) = cons2evec(udg(:,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,NXE_LO,j,k) = cons2evec(udg(:,ii,jj,kk,NXE_LO,j,k))
                edg(:,ii,jj,kk,NXE_HI,j,k) = cons2evec(udg(:,ii,jj,kk,NXE_HI,j,k))
            end do; end do; end do;
        end do; end do
# endif
# if NDIM > 1
        do k = 1,NZE; do i = 1,NXE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,NYE_LO,k) = cons2evec(udg(:,ii,jj,kk,i,NYE_LO,k))
                edg(:,ii,jj,kk,i,NYE_HI,k) = cons2evec(udg(:,ii,jj,kk,i,NYE_HI,k))
            end do; end do; end do;
        end do; end do
# endif
# if NDIM > 2
        do j = 1,NYE; do i = 1,NXE
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                edg(:,ii,jj,kk,i,j,NZE_LO) = cons2evec(udg(:,ii,jj,kk,i,j,NZE_LO))
                edg(:,ii,jj,kk,i,j,NZE_HI) = cons2evec(udg(:,ii,jj,kk,i,j,NZE_HI))
            end do; end do; end do;
        end do; end do
# endif

        itf = 0.0
# if NDIM > 0
        !! Interpolate to boundaries in x-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*edg(h,ll,jj,kk,i+1,j,k)
                itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*edg(h,ll,jj,kk,i  ,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        !! Interpolate to boundaries in y-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*edg(h,ii,ll,kk,i,j+1,k)
                itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*edg(h,ii,ll,kk,i,j  ,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        !! Interpolate to boundaries in z-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*edg(h,ii,jj,ll,i,j,k+1)
                itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*edg(h,ii,jj,ll,i,j,k  )
            end do; end do; end do; end do;
        end do; end do; end do;
# endif

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                itf(:,jj,kk,i,j,k,NOR) = evec2cons(itf(:,jj,kk,i,j,k,NOR))
                itf(:,jj,kk,i,j,k,SOU) = evec2cons(itf(:,jj,kk,i,j,k,SOU))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                itf(:,ii,kk,i,j,k,WES) = evec2cons(itf(:,ii,kk,i,j,k,WES))
                itf(:,ii,kk,i,j,k,EAS) = evec2cons(itf(:,ii,kk,i,j,k,EAS))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                itf(:,ii,jj,i,j,k,FRO) = evec2cons(itf(:,ii,jj,i,j,k,FRO))
                itf(:,ii,jj,i,j,k,BAC) = evec2cons(itf(:,ii,jj,i,j,k,BAC))
            end do; end do;
        end do; end do; end do;
# endif
        end block

    else

        itf = 0.0
# if NDIM > 0
        !! Interpolate to boundaries in x-direction.
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,jj,kk,i,j,k,NOR) = itf(h,jj,kk,i,j,k,NOR) + toBoundaryVecM(ll)*udg(h,ll,jj,kk,i+1,j,k)
                itf(h,jj,kk,i,j,k,SOU) = itf(h,jj,kk,i,j,k,SOU) + toBoundaryVecP(ll)*udg(h,ll,jj,kk,i  ,j,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        !! Interpolate to boundaries in y-direction.
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,kk,i,j,k,WES) = itf(h,ii,kk,i,j,k,WES) + toBoundaryVecM(ll)*udg(h,ii,ll,kk,i,j+1,k)
                itf(h,ii,kk,i,j,k,EAS) = itf(h,ii,kk,i,j,k,EAS) + toBoundaryVecP(ll)*udg(h,ii,ll,kk,i,j  ,k)
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        !! Interpolate to boundaries in z-direction.
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NCONS_VARS;
                itf(h,ii,jj,i,j,k,FRO) = itf(h,ii,jj,i,j,k,FRO) + toBoundaryVecM(ll)*udg(h,ii,jj,ll,i,j,k+1)
                itf(h,ii,jj,i,j,k,BAC) = itf(h,ii,jj,i,j,k,BAC) + toBoundaryVecP(ll)*udg(h,ii,jj,ll,i,j,k  )
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
    end if !! DGFV_NODE_TYPE

    !! ------------------------------------------------------------------ !!
    !! Calculate volume flux.

    xfxm = 0.0
    yfym = 0.0
    zfzm = 0.0

    !! symmetry-exploiting form
# if NDIM > 0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = ii+1,NXN;
            xfxm(:,ll,ii,jj,kk,i,j,k) = TWO_POINT_XFLUX_cons(udg(:,ll,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k))
            xfxm(:,ii,ll,jj,kk,i,j,k) = xfxm(:,ll,ii,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = jj+1,NYN;
            yfym(:,ll,ii,jj,kk,i,j,k) = TWO_POINT_YFLUX_cons(udg(:,ii,ll,kk,i,j,k),udg(:,ii,jj,kk,i,j,k))
            yfym(:,jj,ii,ll,kk,i,j,k) = yfym(:,ll,ii,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = kk+1,NZN;
            zfzm(:,ll,ii,jj,kk,i,j,k) = TWO_POINT_ZFLUX_cons(udg(:,ii,jj,ll,i,j,k),udg(:,ii,jj,kk,i,j,k))
            zfzm(:,kk,ii,jj,ll,i,j,k) = zfzm(:,ll,ii,jj,kk,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

    if (DGFV_MHD_NON_CONS) then
# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,NXN;
                xfxm(:,ll,ii,jj,kk,i,j,k) = xfxm(:,ll,ii,jj,kk,i,j,k) &
                    + non_Conservative_xFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),udg(:,ll,jj,kk,i,j,k))
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,NYN;
                yfym(:,ll,ii,jj,kk,i,j,k) = yfym(:,ll,ii,jj,kk,i,j,k) &
                    + non_Conservative_yFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),udg(:,ii,ll,kk,i,j,k))
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,NZN;
                zfzm(:,ll,ii,jj,kk,i,j,k) = zfzm(:,ll,ii,jj,kk,i,j,k) &
                    + non_Conservative_zFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,ll,i,j,k))
            end do; end do; end do; end do;
        end do; end do; end do;
# endif
    end if !! DGFV_MHD_NON_CONS

    xvfl = 0.0
    yvfl = 0.0
    zvfl = 0.0

    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
# if NDIM > 0
           xvfl(h,ii,jj,kk,i,j,k) = xvfl(h,ii,jj,kk,i,j,k) + fluxDiffMat(ii,ll) * xfxm(h,ll,ii,jj,kk,i,j,k)
# endif
# if NDIM > 1
           yvfl(h,ii,jj,kk,i,j,k) = yvfl(h,ii,jj,kk,i,j,k) + fluxDiffMat(jj,ll) * yfym(h,ll,ii,jj,kk,i,j,k)
# endif
# if NDIM > 2
           zvfl(h,ii,jj,kk,i,j,k) = zvfl(h,ii,jj,kk,i,j,k) + fluxDiffMat(kk,ll) * zfzm(h,ll,ii,jj,kk,i,j,k)
# endif
        end do; end do; end do; end do; end do;
    end do; end do; end do;

    !! ------------------------------------------------------------------ !!
    !! ------------------------------------------------------------------ !!
    !! Surface fluxes.

    !! ------------------------------------------------------------------ !!
    !! Reconstruct block surface flux.

    rim = 0.0

# if NDIM == 1
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = fvflux(h,NXN*i,j,k,IAXIS)
        end do; end do; end do;
    end do; end do; end do;

# elif NDIM == 2
    !! Reconstruct in y-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = rim(h,jj,kk,i,j,k,IAXIS) + recoMat(jj,ll)*fvflux(h,NXN*i,NYN*(j-1)+ll,NZN*(k-1)+kk,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in x-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,kk,i,j,k,JAXIS) = rim(h,ii,kk,i,j,k,JAXIS) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*j,NZN*(k-1)+kk,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

# elif NDIM == 3
    !! Reconstruct in y-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,1,jj,kk,i,j,k) = atmp(h,1,jj,kk,i,j,k) + recoMat(jj,ll)*fvflux(h,NXN*i,NYN*(j-1)+ll,NZN*(k-1)+kk,IAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,jj,kk,i,j,k,IAXIS) = rim(h,jj,kk,i,j,k,IAXIS) + recoMat(kk,ll)*atmp(h,1,jj,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruct in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,1,kk,i,j,k) = atmp(h,ii,1,kk,i,j,k) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*j,NZN*(k-1)+kk,JAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in z-direction.
    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
        do kk = 1,NZN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,kk,i,j,k,JAXIS) = rim(h,ii,kk,i,j,k,JAXIS) + recoMat(kk,ll)*atmp(h,ii,1,ll,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;


    !! Reconstruct in x-direction.
    atmp = 0.0
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,jj,1,i,j,k) = atmp(h,ii,jj,1,i,j,k) + recoMat(ii,ll)*fvflux(h,NXN*(i-1)+ll,NYN*(j-1)+jj,NZN*k,KAXIS)
        end do; end do; end do; end do;
    end do; end do; end do;

    !! Reconstruct in y-direction.
    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
        do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            rim(h,ii,jj,i,j,k,KAXIS) = rim(h,ii,jj,i,j,k,KAXIS) + recoMat(jj,ll)*atmp(h,ii,ll,1,i,j,k)
        end do; end do; end do; end do;
    end do; end do; end do;
# endif

    !! ------------------------------------------------------------------ !!
    !! Surface flux differencing term.

    if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS) then

        block
        real :: xsflVecM(NFLUXES,NYN,NZN,NXE,NYE,NZE)
        real :: xsflVecP(NFLUXES,NYN,NZN,NXE,NYE,NZE)

        real :: ysflVecM(NFLUXES,NXN,NZN,NXE,NYE,NZE)
        real :: ysflVecP(NFLUXES,NXN,NZN,NXE,NYE,NZE)

        real :: zsflVecM(NFLUXES,NXN,NYN,NXE,NYE,NZE)
        real :: zsflVecP(NFLUXES,NXN,NYN,NXE,NYE,NZE)

        real :: xsflM(NFLUXES,NXN,NYN,NZN,NXE,NYE,NZE)
        real :: xsflP(NFLUXES,NYN,NYN,NZN,NXE,NYE,NZE)

        real :: ysflM(NFLUXES,NXN,NYN,NZN,NXE,NYE,NZE)
        real :: ysflP(NFLUXES,NYN,NYN,NZN,NXE,NYE,NZE)

        real :: zsflM(NFLUXES,NXN,NYN,NZN,NXE,NYE,NZE)
        real :: zsflP(NFLUXES,NYN,NYN,NZN,NXE,NYE,NZE)

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
                xsflM(:,ii,jj,kk,i,j,k) = TWO_POINT_XFLUX_cons(itf(:,jj,kk,i-1,j,k,NOR),udg(:,ii,jj,kk,i,j,k))
                xsflP(:,ii,jj,kk,i,j,k) = TWO_POINT_XFLUX_cons(itf(:,jj,kk,i  ,j,k,SOU),udg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 1
                ysflM(:,ii,jj,kk,i,j,k) = TWO_POINT_YFLUX_cons(itf(:,ii,kk,i,j-1,k,WES),udg(:,ii,jj,kk,i,j,k))
                ysflP(:,ii,jj,kk,i,j,k) = TWO_POINT_YFLUX_cons(itf(:,ii,kk,i,j  ,k,EAS),udg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 2
                zsflM(:,ii,jj,kk,i,j,k) = TWO_POINT_ZFLUX_cons(itf(:,ii,jj,i,j,k-1,FRO),udg(:,ii,jj,kk,i,j,k))
                zsflP(:,ii,jj,kk,i,j,k) = TWO_POINT_ZFLUX_cons(itf(:,ii,jj,i,j,k  ,BAC),udg(:,ii,jj,kk,i,j,k))
# endif
            end do; end do; end do;
        end do; end do; end do;

        xsflVecM = 0.0
        xsflVecP = 0.0

        ysflVecM = 0.0
        ysflVecP = 0.0

        zsflVecM = 0.0
        zsflVecP = 0.0

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
                xsflVecM(:,jj,kk,i,j,k) = xsflVecM(:,jj,kk,i,j,k) - toBoundaryVecM(ii)*xsflM(:,ii,jj,kk,i,j,k)
                xsflVecP(:,jj,kk,i,j,k) = xsflVecP(:,jj,kk,i,j,k) - toBoundaryVecP(ii)*xsflP(:,ii,jj,kk,i,j,k)
# endif
# if NDIM > 1
                ysflVecM(:,ii,kk,i,j,k) = ysflVecM(:,ii,kk,i,j,k) - toBoundaryVecM(jj)*ysflM(:,ii,jj,kk,i,j,k)
                ysflVecP(:,ii,kk,i,j,k) = ysflVecP(:,ii,kk,i,j,k) - toBoundaryVecP(jj)*ysflP(:,ii,jj,kk,i,j,k)
# endif
# if NDIM > 2
                zsflVecM(:,ii,jj,i,j,k) = zsflVecM(:,ii,jj,i,j,k) - toBoundaryVecM(kk)*zsflM(:,ii,jj,kk,i,j,k)
                zsflVecP(:,ii,jj,i,j,k) = zsflVecP(:,ii,jj,i,j,k) - toBoundaryVecP(kk)*zsflP(:,ii,jj,kk,i,j,k)
# endif
            end do; end do; end do;
        end do; end do; end do;

        if (DGFV_MHD_NON_CONS) then
# if NDIM > 0
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN;
                    xsflVecM(:,jj,kk,i,j,k) = xsflVecM(:,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(itf(:,jj,kk,i-1,j,k,NOR),itf(:,jj,kk,i-1,j,k,NOR),itf(:,jj,kk,i-1,j,k,SOU))
                    xsflVecP(:,jj,kk,i,j,k) = xsflVecP(:,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,NOR))
                end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 1
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do ii = 1,NXN;
                    ysflVecM(:,ii,kk,i,j,k) = ysflVecM(:,ii,kk,i,j,k) + Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j-1,k,WES),itf(:,ii,kk,i,j-1,k,WES),itf(:,ii,kk,i,j-1,k,EAS))
                    ysflVecP(:,ii,kk,i,j,k) = ysflVecP(:,ii,kk,i,j,k) + Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,WES))
                end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 2
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do jj = 1,NYN; do ii = 1,NXN;
                    zsflVecM(:,ii,jj,i,j,k) = zsflVecM(:,ii,jj,i,j,k) + Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k-1,FRO),itf(:,ii,jj,i,j,k-1,FRO),itf(:,jj,kk,i,j,k-1,BAC))
                    zsflVecP(:,ii,jj,i,j,k) = zsflVecP(:,ii,jj,i,j,k) + Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k  ,BAC),itf(:,ii,jj,i,j,k  ,BAC),itf(:,jj,kk,i,j,k  ,FRO))
                end do; end do;
            end do; end do; end do;
# endif


            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
                    xsflVecM(:,jj,kk,i,j,k) = xsflVecM(:,jj,kk,i,j,k) - toBoundaryVecM(ii)*Non_Conservative_xFlux_cons(itf(:,jj,kk,i-1,j,k,NOR),itf(:,jj,kk,i-1,j,k,NOR),udg(:,ii,jj,kk,i,j,k))
                    xsflVecP(:,jj,kk,i,j,k) = xsflVecP(:,jj,kk,i,j,k) - toBoundaryVecP(ii)*Non_Conservative_xFlux_cons(itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,SOU),udg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 1
                    ysflVecM(:,ii,kk,i,j,k) = ysflVecM(:,ii,kk,i,j,k) - toBoundaryVecM(jj)*Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j-1,k,WES),itf(:,ii,kk,i,j-1,k,WES),udg(:,ii,jj,kk,i,j,k))
                    ysflVecP(:,ii,kk,i,j,k) = ysflVecP(:,ii,kk,i,j,k) - toBoundaryVecP(jj)*Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,EAS),udg(:,ii,jj,kk,i,j,k))
# endif
# if NDIM > 2
                    zsflVecM(:,ii,jj,i,j,k) = zsflVecM(:,ii,jj,i,j,k) - toBoundaryVecM(kk)*Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k-1,FRO),itf(:,ii,jj,i,j,k-1,FRO),udg(:,ii,jj,kk,i,j,k))
                    zsflVecP(:,ii,jj,i,j,k) = zsflVecP(:,ii,jj,i,j,k) - toBoundaryVecP(kk)*Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k  ,BAC),itf(:,ii,jj,i,j,k  ,BAC),udg(:,ii,jj,kk,i,j,k))
# endif
                end do; end do; end do;
            end do; end do; end do;

        end if !! DGFV_MHD_NON_CONS

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
                xsflM(:,ii,jj,kk,i,j,k) = xsflM(:,ii,jj,kk,i,j,k) + rim(:,jj,kk,i-1,j,k,IAXIS) + xsflVecM(:,jj,kk,i,j,k)
                xsflP(:,ii,jj,kk,i,j,k) = xsflP(:,ii,jj,kk,i,j,k) + rim(:,jj,kk,i  ,j,k,IAXIS) + xsflVecP(:,jj,kk,i,j,k)
# endif
# if NDIM > 1
                ysflM(:,ii,jj,kk,i,j,k) = ysflM(:,ii,jj,kk,i,j,k) + rim(:,ii,kk,i,j-1,k,JAXIS) + ysflVecM(:,ii,kk,i,j,k)
                ysflP(:,ii,jj,kk,i,j,k) = ysflP(:,ii,jj,kk,i,j,k) + rim(:,ii,kk,i,j  ,k,JAXIS) + ysflVecP(:,ii,kk,i,j,k)
# endif
# if NDIM > 2
                zsflM(:,ii,jj,kk,i,j,k) = zsflM(:,ii,jj,kk,i,j,k) + rim(:,ii,jj,i,j,k-1,KAXIS) + zsflVecM(:,ii,jj,i,j,k)
                zsflP(:,ii,jj,kk,i,j,k) = zsflP(:,ii,jj,kk,i,j,k) + rim(:,ii,jj,i,j,k  ,KAXIS) + zsflVecP(:,ii,jj,i,j,k)
# endif
            end do; end do; end do;
        end do; end do; end do;

        if (DGFV_MHD_NON_CONS) then
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
# if NDIM > 0
                    xsflM(:,ii,jj,kk,i,j,k) = xsflM(:,ii,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,jj,kk,i-1,j,k,NOR))
                    xsflP(:,ii,jj,kk,i,j,k) = xsflP(:,ii,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,jj,kk,i  ,j,k,SOU))
# endif
# if NDIM > 1
                    ysflM(:,ii,jj,kk,i,j,k) = ysflM(:,ii,jj,kk,i,j,k) + Non_Conservative_yFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,ii,kk,i,j-1,k,WES))
                    ysflP(:,ii,jj,kk,i,j,k) = ysflP(:,ii,jj,kk,i,j,k) + Non_Conservative_yFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,ii,kk,i,j  ,k,EAS))
# endif
# if NDIM > 2
                    zsflM(:,ii,jj,kk,i,j,k) = zsflM(:,ii,jj,kk,i,j,k) + Non_Conservative_zFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,ii,jj,i,j,k-1,FRO))
                    zsflP(:,ii,jj,kk,i,j,k) = zsflP(:,ii,jj,kk,i,j,k) + Non_Conservative_zFlux_cons(udg(:,ii,jj,kk,i,j,k),udg(:,ii,jj,kk,i,j,k),itf(:,ii,jj,i,j,k  ,BAC))
# endif
                end do; end do; end do;
            end do; end do; end do;
        end if !! DGFV_MHD_NON_CONS

        !! Apply surface operator.
        xsfl = 0.0
        ysfl = 0.0
        zsfl = 0.0

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do h = 1,NFLUXES;
# if NDIM > 0
                xsfl(h,ii,jj,kk,i,j,k) = xsfl(h,ii,jj,kk,i,j,k) + surfVecM(ii)*xsflM(h,ii,jj,kk,i,j,k) + surfVecP(ii)*xsflP(h,ii,jj,kk,i,j,k)
# endif
# if NDIM > 1
                ysfl(h,ii,jj,kk,i,j,k) = ysfl(h,ii,jj,kk,i,j,k) + surfVecM(jj)*ysflM(h,ii,jj,kk,i,j,k) + surfVecP(jj)*ysflP(h,ii,jj,kk,i,j,k)
# endif
# if NDIM > 2
                zsfl(h,ii,jj,kk,i,j,k) = zsfl(h,ii,jj,kk,i,j,k) + surfVecM(kk)*zsflM(h,ii,jj,kk,i,j,k) + surfVecP(kk)*zsflP(h,ii,jj,kk,i,j,k)
# endif
            end do; end do; end do; end do;
        end do; end do; end do;
        end block

        !! ------------------------------------------------------------------ !!
        !! Collect all flux contributions.

# if NDIM > 0
        dgrhs =         sdx(IAXIS)*(xvfl + xsfl)
# endif
# if NDIM > 1
        dgrhs = dgrhs + sdx(JAXIS)*(yvfl + ysfl)
# endif
# if NDIM > 2
        dgrhs = dgrhs + sdx(KAXIS)*(zvfl + zsfl)
# endif

        if (DGFV_FV_ACTIVE) then
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                if (blend(i,j,k) == 0.0 .or. any(isnan(dgrhs(:,:,:,:,i,j,k))) .or. any(dgrhs(:,:,:,:,i,j,k)-1.0 == dgrhs(:,:,:,:,i,j,k))) then
                    blend(i,j,k) = 0.0
                    dgrhs(:,:,:,:,i,j,k) = 0.0
                end if
            end do; end do; end do;

            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                ok = .true.
# if NDIM > 0
                do kk = 1,NZN; do jj = 1,NYN;
                    ok = ok .and. isvalid_cons(itf(:,jj,kk,i,j,k,SOU)) .and. isvalid_cons(itf(:,jj,kk,i-1,j,k,NOR))
                end do; end do;
# endif
# if NDIM > 1
                do kk = 1,NZN; do ii = 1,NXN;
                    ok = ok .and. isvalid_cons(itf(:,ii,kk,i,j,k,EAS)) .and. isvalid_cons(itf(:,ii,kk,i,j-1,k,WES))
                end do; end do;
# endif
# if NDIM > 2
                do jj = 1,NYN; do ii = 1,NXN;
                    ok = ok .and. isvalid_cons(itf(:,ii,jj,i,j,k,BAC)) .and. isvalid_cons(itf(:,ii,jj,i,j,k-1,FRO))
                end do; end do;
# endif
                if (.not. ok) then
                    blend(i,j,k) = 0.0
                    dgrhs(:,:,:,:,i,j,k) = 0.0
                end if
            end do; end do; end do;
        end if

    else ! if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_LOBATTO) then

        block
        real :: xsflM(NFLUXES,NYN,NZN,NXE,NYE,NZE)
        real :: xsflP(NFLUXES,NYN,NZN,NXE,NYE,NZE)

        real :: ysflM(NFLUXES,NXN,NZN,NXE,NYE,NZE)
        real :: ysflP(NFLUXES,NXN,NZN,NXE,NYE,NZE)

        real :: zsflM(NFLUXES,NXN,NYN,NXE,NYE,NZE)
        real :: zsflP(NFLUXES,NXN,NYN,NXE,NYE,NZE)

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN;
                xsflM(:,jj,kk,i,j,k) = rim(:,jj,kk,i-1,j,k,IAXIS)
                xsflP(:,jj,kk,i,j,k) = rim(:,jj,kk,i  ,j,k,IAXIS)
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                ysflM(:,ii,kk,i,j,k) = rim(:,ii,kk,i,j-1,k,JAXIS)
                ysflP(:,ii,kk,i,j,k) = rim(:,ii,kk,i,j  ,k,JAXIS)
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                zsflM(:,ii,jj,i,j,k) = rim(:,ii,jj,i,j,k-1,KAXIS)
                zsflP(:,ii,jj,i,j,k) = rim(:,ii,jj,i,j,k  ,KAXIS)
            end do; end do;
        end do; end do; end do;
# endif

        if (DGFV_MHD_NON_CONS) then
# if NDIM > 0
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN;
                    xsflM(:,jj,kk,i,j,k) = xsflM(:,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(itf(:,jj,kk,i-1,j,k,NOR),itf(:,jj,kk,i-1,j,k,NOR),itf(:,jj,kk,i-1,j,k,SOU))
                    xsflP(:,jj,kk,i,j,k) = xsflP(:,jj,kk,i,j,k) + Non_Conservative_xFlux_cons(itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,NOR))
                end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 1
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do ii = 1,NXN;
                    ysflM(:,ii,kk,i,j,k) = ysflM(:,ii,kk,i,j,k) + Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j-1,k,WES),itf(:,ii,kk,i,j-1,k,WES),itf(:,ii,kk,i,j-1,k,EAS))
                    ysflP(:,ii,kk,i,j,k) = ysflP(:,ii,kk,i,j,k) + Non_Conservative_yFlux_cons(itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,WES))
                end do; end do;
            end do; end do; end do;
# endif
# if NDIM > 2
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do jj = 1,NYN; do ii = 1,NXN;
                    zsflM(:,ii,jj,i,j,k) = zsflM(:,ii,jj,i,j,k) + Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k-1,FRO),itf(:,ii,jj,i,j,k-1,FRO),itf(:,jj,kk,i,j,k-1,BAC))
                    zsflP(:,ii,jj,i,j,k) = zsflP(:,ii,jj,i,j,k) + Non_Conservative_zFlux_cons(itf(:,ii,jj,i,j,k  ,BAC),itf(:,ii,jj,i,j,k  ,BAC),itf(:,jj,kk,i,j,k  ,FRO))
                end do; end do;
            end do; end do; end do;
# endif
        end if

        !! Apply surface operator.
        xsfl = 0.0
        ysfl = 0.0
        zsfl = 0.0

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN;
                xsfl(:,  1,jj,kk,i,j,k) = surfVecM(  1)*xsflM(:,jj,kk,i,j,k)
                xsfl(:,NXN,jj,kk,i,j,k) = surfVecP(NXN)*xsflP(:,jj,kk,i,j,k)
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                ysfl(:,ii,  1,kk,i,j,k) = surfVecM(  1)*ysflM(:,ii,kk,i,j,k)
                ysfl(:,ii,NYN,kk,i,j,k) = surfVecP(NYN)*ysflP(:,ii,kk,i,j,k)
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                zsfl(:,ii,jj,  1,i,j,k) = surfVecM(  1)*zsflM(:,ii,jj,i,j,k)
                zsfl(:,ii,jj,NZN,i,j,k) = surfVecP(NZN)*zsflP(:,ii,jj,i,j,k)
            end do; end do;
        end do; end do; end do;
# endif
        end block

# if NDIM > 0
        dgrhs =         sdx(IAXIS)*(xvfl + xsfl)
# endif
# if NDIM > 1
        dgrhs = dgrhs + sdx(JAXIS)*(yvfl + ysfl)
# endif
# if NDIM > 2
        dgrhs = dgrhs + sdx(KAXIS)*(zvfl + zsfl)
# endif

        if (DGFV_FV_ACTIVE) then
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                if (blend(i,j,k) == 0.0 .or. any(isnan(dgrhs(:,:,:,:,i,j,k))) &
                        .or. any(dgrhs(:,:,:,:,i,j,k)-1.0 == dgrhs(:,:,:,:,i,j,k))) then
                    blend(i,j,k) = 0.0
                    dgrhs(:,:,:,:,i,j,k) = 0.0
                end if
            end do; end do; end do;
        end if

    end if !! DGFV_NODE_TYPE

    ! if (DGFV_ENTROPY_CORRECTION) then

        !! ------------------------------------------------------------------ !!
        !! Calculate entropy production rate per DG (sub-)element.

        block
        ! real  :: erhs(NXE,NYE,NZE)

        ! real :: edg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
        ! real :: sdg (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)

        ! real :: highs (NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
        ! real :: means (NCONS_VARS,              NXE_LO:NXE_HI, NYE_LO:NYE_HI, NZE_LO:NZE_HI)
        real :: means(NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

        real :: eflux(         NXE,NYE,NZE)
        real :: erhs(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: theta(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: evfl(         NXN,NYN,NZN, NXE,NYE,NZE)
        real :: esfl(         NXN,NYN,NZN, NXE,NYE,NZE)

        real :: ievec(NCONS_VARS, NYN,NZN,     NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
        real :: vevec(NCONS_VARS, NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        real :: iepot(           NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
        real :: ieflx(           NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

        ! real :: devec(NPROP_FLUX,  NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        ! real :: aevec(NPROP_FLUX,  NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        ! real :: mevec(NPROP_FLUX,               NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        ! real :: tevec(          NXN,NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        ! real :: tesum(                       NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)

        !! used for validation
        real :: dgevfl(         NXN,NYN,NZN, NXE,NYE,NZE)

        ! real :: itf(NCONS_VARS, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, NOR:BAC)
        ! real :: squeeze(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE)
        ! logical :: ok(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE), first

        !! Calculate means.

!!         ievec = 0
!!         iepot = 0
!! 
!! # if NDIM > 0
!!         do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
!!             do kk = 1,NZN; do jj = 1,NYN
!!                 ievec(:,jj,kk,i,j,k,IAXIS) = 0.5*(cons2evec(itf(:,jj,kk,i,j,k,SOU)) + cons2evec(itf(:,jj,kk,i,j,k,NOR)))
!!                 iepot(  jj,kk,i,j,k,IAXIS) = 0.5*(cons2xpot(itf(:,jj,kk,i,j,k,SOU)) + cons2xpot(itf(:,jj,kk,i,j,k,NOR)))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif
!! # if NDIM > 1
!!         do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
!!             do kk = 1,NZN; do ii = 1,NXN;
!!                 ievec(:,ii,kk,i,j,k,JAXIS) = 0.5*(cons2evec(itf(:,ii,kk,i,j,k,EAS)) + cons2evec(itf(:,ii,kk,i,j,k,WES)))
!!                 iepot(  ii,kk,i,j,k,JAXIS) = 0.5*(cons2ypot(itf(:,ii,kk,i,j,k,EAS)) + cons2ypot(itf(:,ii,kk,i,j,k,WES)))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif
!! # if NDIM > 2
!!         do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
!!             do jj = 1,NYN; do ii = 1,NXN;
!!                 ievec(:,ii,jj,i,j,k,KAXIS) = 0.5*(cons2evec(itf(:,ii,jj,i,j,k,BAC)) + cons2evec(itf(:,ii,jj,i,j,k,FRO)))
!!                 iepot(  ii,jj,i,j,k,KAXIS) = 0.5*(cons2zpot(itf(:,ii,jj,i,j,k,BAC)) + cons2zpot(itf(:,ii,jj,i,j,k,FRO)))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif

# if NDIM > 0
        do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
            do kk = 1,NZN; do jj = 1,NYN
                !ieflx(jj,kk,i,j,k,IAXIS) = sdx(IAXIS)*(dot_product(ievec(1:NPROP_FLUX,jj,kk,i,j,k,IAXIS),rim(1:NPROP_FLUX,jj,kk,i,j,k,IAXIS)) - iepot(jj,kk,i,j,k,IAXIS))
                ieflx(jj,kk,i,j,k,IAXIS) = sdx(IAXIS)*xentrflux_cons(-1,rim(:,jj,kk,i,j,k,IAXIS),itf(:,jj,kk,i,j,k,SOU),itf(:,jj,kk,i,j,k,NOR))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 1
        do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
            do kk = 1,NZN; do ii = 1,NXN;
                !ieflx(ii,kk,i,j,k,JAXIS) = sdx(JAXIS)*(dot_product(ievec(1:NPROP_FLUX,ii,kk,i,j,k,JAXIS),rim(1:NPROP_FLUX,ii,kk,i,j,k,JAXIS)) - iepot(ii,kk,i,j,k,JAXIS))
                ieflx(ii,kk,i,j,k,JAXIS) = sdx(JAXIS)*yentrflux_cons( 1,rim(:,ii,kk,i,j,k,JAXIS),itf(:,ii,kk,i,j,k,EAS),itf(:,ii,kk,i,j,k,WES))
            end do; end do;
        end do; end do; end do;
# endif
# if NDIM > 2
        do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
            do jj = 1,NYN; do ii = 1,NXN;
                !ieflx(ii,jj,i,j,k,KAXIS) = sdx(KAXIS)*(dot_product(ievec(1:NPROP_FLUX,ii,kk,i,j,k,KAXIS),rim(1:NPROP_FLUX,ii,jj,i,j,k,KAXIS)) - iepot(ii,jj,i,j,k,KAXIS))
                ieflx(ii,jj,i,j,k,KAXIS) = sdx(KAXIS)*zentrflux_cons( 1,rim(:,ii,jj,i,j,k,KAXIS),itf(:,ii,jj,i,j,k,BAC),itf(:,ii,jj,i,j,k,FRO))
            end do; end do;
        end do; end do; end do;
# endif

        !! Apply surface operator.
        esfl = 0.0
        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                esfl(ii,jj,kk,i,j,k) = esfl(ii,jj,kk,i,j,k) &
# if NDIM > 0
                    + surfVecM(ii)*ieflx(jj,kk,i-1,j,k,IAXIS) &
                    + surfVecP(ii)*ieflx(jj,kk,i  ,j,k,IAXIS) &
# endif
# if NDIM > 1
                    + surfVecM(jj)*ieflx(ii,kk,i,j-1,k,JAXIS) &
                    + surfVecP(jj)*ieflx(ii,kk,i,j  ,k,JAXIS) &
# endif
# if NDIM > 2
                    + surfVecM(kk)*ieflx(ii,jj,i,j,k-1,KAXIS) &
                    + surfVecP(kk)*ieflx(ii,jj,i,j,k  ,KAXIS) &
# endif
            ;end do; end do; end do;
        end do; end do; end do;

!!         eflux = 0.0
!! 
!! # if NDIM > 0
!!         do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!!             do kk = 1,NZN; do jj = 1,NYN
!!                 eflux(i,j,k) = eflux(i,j,k) + yweights(jj)*zweights(kk)*sdx(IAXIS)*xentrflux_cons(-1,rim(:,jj,kk,i-1,j,k,IAXIS),itf(:,jj,kk,i-1,j,k,SOU),itf(:,jj,kk,i-1,j,k,NOR))
!!                 eflux(i,j,k) = eflux(i,j,k) - yweights(jj)*zweights(kk)*sdx(IAXIS)*xentrflux_cons( 1,rim(:,jj,kk,i  ,j,k,IAXIS),itf(:,jj,kk,i  ,j,k,SOU),itf(:,jj,kk,i  ,j,k,NOR))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif
!! # if NDIM > 1
!!         do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!!             do kk = 1,NZN; do ii = 1,NXN;
!!                 eflux(i,j,k) = eflux(i,j,k) + xweights(ii)*zweights(kk)*sdx(JAXIS)*yentrflux_cons(-1,rim(:,ii,kk,i,j-1,k,JAXIS),itf(:,ii,kk,i,j-1,k,EAS),itf(:,ii,kk,i,j-1,k,WES))
!!                 eflux(i,j,k) = eflux(i,j,k) - xweights(ii)*zweights(kk)*sdx(JAXIS)*yentrflux_cons( 1,rim(:,ii,kk,i,j  ,k,JAXIS),itf(:,ii,kk,i,j  ,k,EAS),itf(:,ii,kk,i,j  ,k,WES))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif
!! # if NDIM > 2
!!         do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!!             do jj = 1,NYN; do ii = 1,NXN;
!!                 eflux(i,j,k) = eflux(i,j,k) + xweights(ii)*yweights(jj)*sdx(KAXIS)*zentrflux_cons(-1,rim(:,ii,jj,i,j,k-1,KAXIS),itf(:,ii,jj,i,j,k-1,BAC),itf(:,ii,jj,i,j,k-1,FRO))
!!                 eflux(i,j,k) = eflux(i,j,k) - xweights(ii)*yweights(jj)*sdx(KAXIS)*zentrflux_cons( 1,rim(:,ii,jj,i,j,k  ,KAXIS),itf(:,ii,jj,i,j,k  ,BAC),itf(:,ii,jj,i,j,k  ,FRO))
!!             end do; end do;
!!         end do; end do; end do;
!! # endif

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                vevec(:,ii,jj,kk,i,j,k) = cons2evec(udg(:,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN
                dgevfl(ii,jj,kk,i,j,k) = dot_product(vevec(1:NPROP_FLUX,ii,jj,kk,i,j,k),dgrhs(1:NPROP_FLUX,ii,jj,kk,i,j,k))
            end do; end do; end do;
        end do; end do; end do;

        do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
            dgerhs(i,j,k) = sum(dgvolumes*(dgevfl(:,:,:,i,j,k) - esfl(:,:,:,i,j,k)))
            !dgerhs(i,j,k) = sum(dgvolumes*(dgevfl(:,:,:,i,j,k) - evfl(:,:,:,i,j,k) - esfl(:,:,:,i,j,k)))
            !dgerhs(i,j,k) = sum(dgvolumes*(dgevfl(:,:,:,i,j,k))) - eflux(i,j,k)
            !dgerhs(i,j,k) = eflux(i,j,k)
        end do; end do; end do;
        end block

    ! end if !! DGFV_ENTROPY_CORRECTION


    !! ------------------------------------------------------------------ !!
    !! Project to median space.

# if NDIM == 1
# define atmp rhs
# elif NDIM == 2
# define btmp rhs
# elif NDIM == 3
# define ctmp rhs
# endif

# if NDIM > 0
    !! Projection in x-direction.
    atmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            atmp(h,ii,jj,kk,i,j,k) = atmp(h,ii,jj,kk,i,j,k) + projMat(ii,ll)*dgrhs(h,ll,jj,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 1
    !! Projection in y-direction.
    btmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            btmp(h,ii,jj,kk,i,j,k) = btmp(h,ii,jj,kk,i,j,k) + projMat(jj,ll)*atmp(h,ii,ll,kk,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif
# if NDIM > 2
    !! Projection in z-direction.
    ctmp = 0.0
    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES; do h = 1,NFLUXES;
            ctmp(h,ii,jj,kk,i,j,k) = ctmp(h,ii,jj,kk,i,j,k) + projMat(kk,ll)*btmp(h,ii,jj,ll,i,j,k)
        end do; end do; end do; end do; end do;
    end do; end do; end do;
# endif

# if NDIM == 1
# undef atmp
# elif NDIM == 2
# undef btmp
# elif NDIM == 3
# undef ctmp
# endif

end subroutine

end module


!! # if NDIM > 0
!!     do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN
!!             itf(:,jj,kk,i,j,k,NOR) = cons2prim(itf(:,jj,kk,i,j,k,NOR))
!!             itf(:,jj,kk,i,j,k,SOU) = cons2prim(itf(:,jj,kk,i,j,k,SOU))
!!         end do; end do;
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 1
!!     do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do ii = 1,NXN;
!!             itf(:,ii,kk,i,j,k,WES) = cons2prim(itf(:,ii,kk,i,j,k,WES))
!!             itf(:,ii,kk,i,j,k,EAS) = cons2prim(itf(:,ii,kk,i,j,k,EAS))
!!         end do; end do;
!!     end do; end do; end do;
!! # endif
!! # if NDIM > 2
!!     do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
!!         do jj = 1,NYN; do ii = 1,NXN;
!!             itf(:,ii,jj,i,j,k,FRO) = cons2prim(itf(:,ii,jj,i,j,k,FRO))
!!             itf(:,ii,jj,i,j,k,BAC) = cons2prim(itf(:,ii,jj,i,j,k,BAC))
!!         end do; end do;
!!     end do; end do; end do;
!! # endif

!! # if HYDRO_MHD_NON_CONS
!! # if NDIM > 0
!! !! NORTH side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) = atmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(ii,ll)*ufv(MAGX_CONS,ll,jj,kk,NXE_LO,j,k)
!!             atmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) = atmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,NXE_LO,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) = btmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(jj,ll)*atmp(MAGX_CONS,ii,ll,kk,NXE_LO,j,k)
!!             btmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) = btmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,NXE_LO,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) = ctmp(MAGX_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(kk,ll)*btmp(MAGX_CONS,ii,jj,ll,NXE_LO,j,k)
!!             ctmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) = ctmp(GLMP_CONS,ii,jj,kk,NXE_LO,j,k) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,NXE_LO,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!!
!! !! SOUTH side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) = atmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(ii,ll)*ufv(MAGX_CONS,ll,jj,kk,NXE_HI,j,k)
!!             atmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) = atmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,NXE_HI,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) = btmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(jj,ll)*atmp(MAGX_CONS,ii,ll,kk,NXE_HI,j,k)
!!             btmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) = btmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,NXE_HI,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do k = 1,NZE; do j = 1,NYE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) = ctmp(MAGX_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(kk,ll)*btmp(MAGX_CONS,ii,jj,ll,NXE_HI,j,k)
!!             ctmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) = ctmp(GLMP_CONS,ii,jj,kk,NXE_HI,j,k) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,NXE_HI,j,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # endif /* NDIM > 0 */
!!
!! # if NDIM > 1
!! !! WEST side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) = atmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(ii,ll)*ufv(MAGY_CONS,ll,jj,kk,i,NYE_LO,k)
!!             atmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) = atmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,i,NYE_LO,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) = btmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(jj,ll)*atmp(MAGY_CONS,ii,ll,kk,i,NYE_LO,k)
!!             btmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) = btmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,i,NYE_LO,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) = ctmp(MAGY_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(kk,ll)*btmp(MAGY_CONS,ii,jj,ll,i,NYE_LO,k)
!!             ctmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) = ctmp(GLMP_CONS,ii,jj,kk,i,NYE_LO,k) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,i,NYE_LO,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!!
!! !! EAST side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) = atmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(ii,ll)*ufv(MAGY_CONS,ll,jj,kk,i,NYE_HI,k)
!!             atmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) = atmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,i,NYE_HI,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) = btmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(jj,ll)*atmp(MAGY_CONS,ii,ll,kk,i,NYE_HI,k)
!!             btmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) = btmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,i,NYE_HI,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do k = 1,NZE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) = ctmp(MAGY_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(kk,ll)*btmp(MAGY_CONS,ii,jj,ll,i,NYE_HI,k)
!!             ctmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) = ctmp(GLMP_CONS,ii,jj,kk,i,NYE_HI,k) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,i,NYE_HI,k)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # endif /* NDIM > 1 */
!!
!! # if NDIM > 2
!! !! FRONT side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) = atmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(ii,ll)*ufv(MAGZ_CONS,ll,jj,kk,i,j,NZE_LO)
!!             atmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) = atmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,i,j,NZE_LO)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) = btmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(jj,ll)*atmp(MAGZ_CONS,ii,ll,kk,i,j,NZE_LO)
!!             btmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) = btmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,i,j,NZE_LO)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) = ctmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(kk,ll)*btmp(MAGZ_CONS,ii,jj,ll,i,j,NZE_LO)
!!             ctmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) = ctmp(GLMP_CONS,ii,jj,kk,i,j,NZE_LO) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,i,j,NZE_LO)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!!
!! !! BACK side.
!! # if NDIM > 0
!!     !! Reconstruct in x-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             atmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) = atmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(ii,ll)*ufv(MAGZ_CONS,ll,jj,kk,i,j,NZE_HI)
!!             atmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) = atmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(ii,ll)*ufv(GLMP_CONS,ll,jj,kk,i,j,NZE_HI)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 1
!!     !! Reconstruct in y-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             btmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) = btmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(jj,ll)*atmp(MAGZ_CONS,ii,ll,kk,i,j,NZE_HI)
!!             btmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) = btmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(jj,ll)*atmp(GLMP_CONS,ii,ll,kk,i,j,NZE_HI)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # if NDIM > 2
!!     !! Reconstruct in z-direction.
!!     do j = 1,NYE; do i = 1,NXE;
!!         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN; do ll = 1,N_NODES;
!!             ctmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) = ctmp(MAGZ_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(kk,ll)*btmp(MAGZ_CONS,ii,jj,ll,i,j,NZE_HI)
!!             ctmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) = ctmp(GLMP_CONS,ii,jj,kk,i,j,NZE_HI) + recoMat(kk,ll)*btmp(GLMP_CONS,ii,jj,ll,i,j,NZE_HI)
!!         end do; end do; end do; end do;
!!     end do; end do;
!! # endif
!! # endif /* NDIM > 2 */
!! # endif /* HYDRO_MHD_NON_CONS */
