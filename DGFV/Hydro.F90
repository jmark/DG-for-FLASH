!!****if* source/physics/Hydro/HydroMain/split/DGFV/Hydro
!!
!! NAME
!!
!!  Hydro
!!
!!
!! SYNOPSIS
!!
!!  Hydro(integer(IN) :: blockCount,
!!        integer(IN) :: blockList(blockCount)
!!        real(IN)    :: timeEndAdv,
!!        real(IN)    :: dt,
!!        real(IN)    :: dtOld,
!!        integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!!
!!  Top level function of the DGFV solver.
!!
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - end time (ignored)
!!  dt_inout   - timestep
!!  dtOld      - old timestep (ignored)
!!  sweepOrder - sweep order (ignored)
!!
!!
!! NOTE
!!
!!  Most of the arguments are dummy in this routine and they are present
!!  to be consistent with the very top layer stub function.
!!
!!  This file copied from ES solver (originally authored by D Dedrig(?))
!!  and modified by Tim Crundall in 2020 to perform evolution in an
!!  unsplit manner.
!!
!!***

! # define LABEL_RECYLCE_DGFV 4242

subroutine Hydro(blockCount, blockList,timeEndAdv, dt_inout, dtOld, sweepOrderIn)

# include "constants.h"
# include "Flash.h"
# include "Eos.h"
# include "DGFV.h"

# ifdef INDEXREORDER
# error in DGFV: Index reordering not supported!
# endif

# ifndef FIXEDBLOCKSIZE
# error in DGFV: Variable block sizes not supported!
# endif

    use Hydro_data, ONLY : hy_ch, hy_ch_local, hy_damp, &
                           hy_chscalingfactor, hy_tref, hy_dtmin
    use Timers_interface, ONLY : Timers_start, Timers_stop
    use Driver_data, only: dr_globalComm, dr_simGeneration, dr_globalMe

    use Hydro_data, ONLY : hy_cfl,  hy_xref, hy_eswitch,            &
                           hy_tref, hy_dref, hy_vref, hy_pref,      &
                           hy_eref, hy_qref, hy_bref, hy_gref,      &
                           hy_mref, hy_nref, hy_kref, hy_useGravity,&
                           hy_useTreeRay,                           &
                           hy_fluxCorrect, hy_irenorm,hy_gcMaskSize,&
                           hy_meshMe, hy_smalldens, hy_dtmin,       &
                           hy_debug, hy_dtminloc

    use Hydro_data, only: nstages
    use Hydro_data, only: nrkcaches

    use Hydro_data, only: dgfv_dg_active
    use Hydro_data, only: dgfv_fv_active
    use Hydro_data, only: dgfv_mhd_non_cons
    use Hydro_data, only: dgfv_mhd_div_cleaning
    use Hydro_data, only: dgfv_enforce_positivity
    use Hydro_data, only: dgfv_fv_at_boundary
    use Hydro_data, only: DGFV_FLUX_DIFFERENCING

    use Hydro_data, only: DGFV_RUNGE_KUTTA
    use Hydro_data, only: DGFV_RUNGE_KUTTA_EULER
    use Hydro_data, only: DGFV_RUNGE_KUTTA_RALSTON
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_22
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_33
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_34
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_45
    use Hydro_data, only: DGFV_RUNGE_KUTTA_LOW_STORAGE_45

    use Hydro_data, only: RK_LS_45_A, RK_LS_45_B, RK_LS_45_C

    use Hydro_data, only: RK_SSP_45_A1, RK_SSP_45_A2, RK_SSP_45_A3
    use Hydro_data, only: RK_SSP_45_B1, RK_SSP_45_B2, RK_SSP_45_B3
    use Hydro_data, only: RK_SSP_45_C1, RK_SSP_45_C2, RK_SSP_45_C3
    use Hydro_data, only: RK_SSP_45_D1, RK_SSP_45_D2, RK_SSP_45_D3
    use Hydro_data, only: RK_SSP_45_E1, RK_SSP_45_E2, RK_SSP_45_E3
    use Hydro_data, only: RK_SSP_45_E4, RK_SSP_45_E5, RK_SSP_45_E6

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

    use Grid_interface, only: Grid_getBlkBC

    use Grid_interface, ONLY : Grid_fillGuardCells,   &
                               Grid_getDeltas,        &
                               Grid_getBlkIndexLimits,&
                               Grid_getCellCoords,    &
                               Grid_getFluxData,      &
                               Grid_putFluxData,      &
                               Grid_conserveFluxes,   &
                               Grid_getBlkPtr,        &
                               Grid_releaseBlkPtr,    &
                               Grid_renormAbundance,  &
                               Grid_limitAbundance,   &
                               Grid_notifySolnDataUpdate

    use Grid_interface, only: Grid_getBlkRefineLevel

    use Eos_interface,     ONLY : Eos_wrapped

# ifdef GRAVITY
    use Gravity_interface, ONLY : Gravity_accelOneRow
# endif
# ifdef TREERAY
    use TreeRay_interface, ONLY : TreeRay_accelOneRow
# endif

    use Driver_interface, only: Driver_abortFlash

    use Hydro_data, only: hy_limit_abundance

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

    use dgfv_kernels_mod, only: calc_dg_surface_flux
    use dgfv_kernels_mod, only: calc_dg_rhs_weak_form
    use dgfv_kernels_mod, only: calc_dg_rhs_flux_differencing

    use dgfv_kernels_mod, only: calc_fv_surface_flux
    use dgfv_kernels_mod, only: calc_fv_rhs

    use dgfv_fluxes_mod, only: non_directional_Fastest_Signal_Speed

    use dgfv_fluxes_mod, only: prim2cons
    use dgfv_fluxes_mod, only: cons2prim

    use dgfv_kernels_mod, only: check_positivity
    use dgfv_kernels_mod, only: enforce_positivity

    use dgfv_utils_mod, only: pp

    use dgfv_fluxes_mod, only: cons2evec
    use dgfv_fluxes_mod, only: prim2evec

    use dgfv_kernels_mod, only: calc_species_rhs

    implicit none

# include "Flash_mpi.h"

    integer, INTENT(in)     :: blockCount, sweepOrderIn
    integer, INTENT(in)     :: blockList(blockCount)
    real,    INTENT(in)     :: timeEndAdv, dtOld
    real,    INTENT(inout)  :: dt_inout
    
    integer :: rkstage, error

    integer :: facesBC(2,NDIM)
    integer :: onBoundary(2,NDIM)

    integer :: blkLimits(LOW:HIGH,MDIM), blkLimitsGC(LOW:HIGH,MDIM)

    integer :: i, j, k, h, blockID, lb
    integer :: ii,jj,kk

    real :: dtmin,dt

    ! logical :: fv_recycle
    logical :: ok

    real, pointer   :: U_ptr(:,:,:,:)
    real            :: U(NUNK_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)

    real, save, allocatable :: blends      (:,:,:,:)
    real, save, allocatable :: flux    (:,:,:,:,:,:)
    real, save, allocatable :: rkcache (:,:,:,:,:,:)

    real :: Uprim(NPRIM_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)
    real :: Ucons(NCONS_VARS,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI)

    real :: Unew ( NCONS_VARS,NXB,NYB,NZB)
    real :: rhs  ( NFLUXES,NXB,NYB,NZB)
    real :: fvrhs( NFLUXES,NXB,NYB,NZB)
    real :: fverhs(        NXB,NYB,NZB)

    real :: sfflux(NFLUXES,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI,IAXIS:KAXIS)
    real :: crflux(NFLUXES)

    real :: fvflux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)

    real :: dgflux(NFLUXES, NYN,NZN, NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
    real :: dgrhs( NFLUXES, NXN,NYN,NZN, 1:NXE,1:NYE,1:NZE)
    real :: dgerhs(                      1:NXE,1:NYE,1:NZE)

    !! surface blend
    real :: sfblend(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)
    real :: sfdnelb(NXE_LO:NXE, NYE_LO:NYE, NZE_LO:NZE, IAXIS:KAXIS)

    !! volume blend
    real :: vlblend(1:NXE,1:NYE,1:NZE)
    real :: vldnelb(1:NXE,1:NYE,1:NZE)

    integer :: axis

    real :: del(NDIM)
    real :: sdx(NDIM)
    real :: deltamin

    !! max. velocity, max. lambda (i.e. signal speed)
    real :: vmax, lmax

# ifdef GRAVITY
    real :: gacx(1-NGUARD*K1D:NXB+NGUARD*K1D,NYB,NZB)
    real :: gacy(1-NGUARD*K2D:NYB+NGUARD*K2D,NXB,NZB)
    real :: gacz(1-NGUARD*K3D:NZB+NGUARD*K3D,NXB,NYB)
# endif
# ifdef TREERAY
    real :: tacx(1-NGUARD*K1D:NXB+NGUARD*K1D,NYB,NZB)
    real :: tacy(1-NGUARD*K2D:NYB+NGUARD*K2D,NXB,NZB)
    real :: tacz(1-NGUARD*K3D:NZB+NGUARD*K3D,NXB,NYB)
# endif
    real :: sacc(NXB,NYB,NZB,IAXIS:KAXIS)

!# if NSPECIES > 0 || NMASS_SCALARS > 0
    real :: fvdensflux(0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
    real :: dgdensflux(NXN,NYN,NZN, NXE,NYE,NZE,IAXIS:KAXIS)
!# endif

    integer :: allocstatus
    !character(len=1024) :: errmsg

    ! if (dr_globalMe == MASTER_PE) then
    !     write (*,*) ' ## ', timeEndAdv, dt_inout
    ! end if

    call Timers_start("DGFV")

    allocate(flux(NFLUXES,NXB_LO:NXB_HI,NYB_LO:NYB_HI,NZB_LO:NZB_HI,IAXIS:KAXIS,blockCount),stat=allocstatus)
    if (allocstatus /= 0) then
        call Driver_abortFlash(' [Hydro] Could not allocate memory for flux array.')
    end if

    if (dgfv_dg_active .and. dgfv_fv_active) then
        allocate(blends(NXE,NYE,NZE,blockCount),stat=allocstatus)
        if (allocstatus /= 0) then
            call Driver_abortFlash(' [Hydro] Could not allocate memory for blend array.')
        end if
    end if

    if (.not.(DGFV_RUNGE_KUTTA == DGFV_RUNGE_KUTTA_EULER)) then
        allocate(rkcache(NUNK_VARS,NXB,NYB,NZB,blockCount,nrkcaches),stat=allocstatus)
        if (allocstatus /= 0) then
            call Driver_abortFlash(' [Hydro] Could not allocate memory for rkcache array.')
        end if
    end if

    !! Normalize timestep.
    dt = dt_inout / hy_tref

    do rkstage = 1,nstages

        call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,eosMode=HYDRO_EOS_MODE,selectBlockType=LEAF)

        flux = 0.0
        if (dgfv_dg_active .and. dgfv_fv_active) then
            blends = 0.0
        end if

        call Timers_start("DGFV_leafs")

        !! Main loop over leaf blocks.
        do lb = 1,blockCount

            blockID = blockList(lb)

            !! ------------------------------------------------------------------ !!
            !! Get a local copy.

            call Grid_getBlkPtr(blockID,U_ptr,CENTER)
            U = U_ptr
            call Grid_releaseBlkPtr(blockID,U_ptr,CENTER)

            !! Apply unitsystem conversion (e.g. CGS, SI, ...).
            Uprim(DENS_PRIM,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
            Uprim(VELX_PRIM,:,:,:) = U(VELX_VAR,:,:,:)/hy_vref
            Uprim(VELY_PRIM,:,:,:) = U(VELY_VAR,:,:,:)/hy_vref
            Uprim(VELZ_PRIM,:,:,:) = U(VELZ_VAR,:,:,:)/hy_vref
            Uprim(PRES_PRIM,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
            Uprim(MAGX_PRIM,:,:,:) = U(MAGX_VAR,:,:,:)/hy_bref
            Uprim(MAGY_PRIM,:,:,:) = U(MAGY_VAR,:,:,:)/hy_bref
            Uprim(MAGZ_PRIM,:,:,:) = U(MAGZ_VAR,:,:,:)/hy_bref
            Uprim(GLMP_PRIM,:,:,:) = U(GLMP_VAR,:,:,:)/hy_bref
            Uprim(ENER_PRIM,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
            Uprim(GAMC_PRIM,:,:,:) = U(GAMC_VAR,:,:,:)

# if NSPECIES > 0
            Uprim(SPEC_PRIM_BEGIN:SPEC_PRIM_END,:,:,:) = U(SPECIES_BEGIN:SPECIES_END,:,:,:)
# endif
# if NMASS_SCALARS > 0
            Uprim(MSCL_PRIM_BEGIN:MSCL_PRIM_END,:,:,:) = U(MASS_SCALARS_BEGIN:MASS_SCALARS_END,:,:,:)
# endif

            !! ------------------------------------------------------------------ !!

            if (dgfv_dg_active .and. dgfv_fv_active) then
                call calc_dg_surface_flux(Uprim,dgflux,sfblend,blends(:,:,:,lb))

                if (dgfv_fv_at_boundary) then
                    call Grid_getBlkBC(blockID,facesBC)
                    if (.not.all(facesBC == NOT_BOUNDARY)) then
                        blends(:,:,:,lb) = 0.0
                    end if
                end if
            else
                call calc_dg_surface_flux(Uprim,dgflux,sfblend,vlblend)
            end if

            !! ------------------------------------------------------------------ !!

            if (dgfv_dg_active .and. dgfv_fv_active) then

                if (any(sfblend < 1.0)) then

                    call calc_fv_surface_flux(Uprim,fvflux)
                    sfdnelb = 1.0 - sfblend

# if NDIM > 0
                    flux(:,1:NXB+1,1:NYB,1:NZB,IAXIS,lb) = fvflux(:,0:NXB,1:NYB,1:NZB,IAXIS)
                    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                        do kk = 1,NZN; do jj = 1,NYN;
                            associate(baseflux => flux(:,NXN*i+1,NYN*(j-1)+jj,NZN*(k-1)+kk,IAXIS,lb), superflux => dgflux(:,jj,kk,i,j,k,IAXIS))
                            baseflux = sfdnelb(i,j,k,IAXIS)*baseflux + sfblend(i,j,k,IAXIS)*superflux
                            end associate
                        end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 1
                    flux(:,1:NXB,1:NYB+1,1:NZB,JAXIS,lb) = fvflux(:,1:NXB,0:NYB,1:NZB,JAXIS)
                    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do ii = 1,NXN;
                            associate(baseflux => flux(:,NXN*(i-1)+ii,NYN*j+1,NZN*(k-1)+kk,JAXIS,lb), superflux => dgflux(:,ii,kk,i,j,k,JAXIS))
                            baseflux = sfdnelb(i,j,k,JAXIS)*baseflux + sfblend(i,j,k,JAXIS)*superflux
                            end associate
                        end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 2
                    flux(:,1:NXB,1:NYB,1:NZB+1,KAXIS,lb) = fvflux(:,1:NXB,1:NYB,0:NZB,KAXIS)
                    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                        do jj = 1,NYN; do ii = 1,NXN;
                            associate(baseflux => flux(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*k+1,KAXIS,lb), superflux => dgflux(:,ii,jj,i,j,k,KAXIS))
                            baseflux = sfdnelb(i,j,k,KAXIS)*baseflux + sfblend(i,j,k,KAXIS)*superflux
                            end associate
                        end do; end do;
                    end do; end do; end do;
# endif
                else
# if NDIM > 0
                    do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                        do kk = 1,NZN; do jj = 1,NYN;
                            flux(:,NXN*i+1,NYN*(j-1)+jj,NZN*(k-1)+kk,IAXIS,lb) = dgflux(:,jj,kk,i,j,k,IAXIS)
                        end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 1
                    do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do ii = 1,NXN;
                            flux(:,NXN*(i-1)+ii,NYN*j+1,NZN*(k-1)+kk,JAXIS,lb) = dgflux(:,ii,kk,i,j,k,JAXIS)
                        end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 2
                    do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                        do jj = 1,NYN; do ii = 1,NXN;
                            flux(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*k+1,KAXIS,lb) = dgflux(:,ii,jj,i,j,k,KAXIS)
                        end do; end do;
                    end do; end do; end do;
# endif
                end if

            else if (dgfv_dg_active) then

# if NDIM > 0
                do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
                    do kk = 1,NZN; do jj = 1,NYN;
                        flux(:,NXN*i+1,NYN*(j-1)+jj,NZN*(k-1)+kk,IAXIS,lb) = dgflux(:,jj,kk,i,j,k,IAXIS)
                    end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 1
                do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do ii = 1,NXN;
                        flux(:,NXN*(i-1)+ii,NYN*j+1,NZN*(k-1)+kk,JAXIS,lb) = dgflux(:,ii,kk,i,j,k,JAXIS)
                    end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 2
                do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
                    do jj = 1,NYN; do ii = 1,NXN;
                        flux(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*k+1,KAXIS,lb) = dgflux(:,ii,jj,i,j,k,KAXIS)
                    end do; end do;
                end do; end do; end do;
# endif

            else if (dgfv_fv_active) then

                call calc_fv_surface_flux(Uprim,fvflux)

# if NDIM > 0
                flux(:,1:NXB+1,1:NYB,1:NZB,IAXIS,lb) = fvflux(:,0:NXB,1:NYB,1:NZB,IAXIS)
# endif
# if NDIM > 1
                flux(:,1:NXB,1:NYB+1,1:NZB,JAXIS,lb) = fvflux(:,1:NXB,0:NYB,1:NZB,JAXIS)
# endif
# if NDIM > 2
                flux(:,1:NXB,1:NYB,1:NZB+1,KAXIS,lb) = fvflux(:,1:NXB,1:NYB,0:NZB,KAXIS)
# endif
            end if
           
            if (hy_fluxCorrect) then
                do axis = 1,NDIM
                    call Grid_putFluxData(blockID,axis,flux(:,:,:,:,axis,lb),(/GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC/))
                end do
            end if

        end do !! lb = 1,blockCount

        call Timers_stop("DGFV_leafs")

        !! ====================================================================== !!
        !! ====================================================================== !!

        if (hy_fluxCorrect) then
            call Grid_conserveFluxes(ALLDIR,0)
        end if

        if (rkstage == nstages) then
            hy_dtmin = HUGE(1.0)
            hy_dtminloc = 0
            hy_ch_local = 0.0
        end if

        !! ====================================================================== !!
        !! ====================================================================== !!

        call Timers_start("DGFV_leafs")

        !! Main loop over leaf blocks.
        do lb = 1,blockCount
            blockID = blockList(lb)

            call Grid_getDeltas(blockID,del)
            del = del/hy_xref
            sdx(1:NDIM) = 1.0/del(1:NDIM)
            deltamin = minval(del(1:NDIM))

            !! ------------------------------------------------------------------ !!
            !! Extract and scale block surface fluxes.

# if NDIM > 0
            fvflux(:,0:NXB,1:NYB,1:NZB,IAXIS) = flux(:,1:NXB+1,1:NYB,1:NZB,IAXIS,lb)
# endif
# if NDIM > 1
            fvflux(:,1:NXB,0:NYB,1:NZB,JAXIS) = flux(:,1:NXB,1:NYB+1,1:NZB,JAXIS,lb)
# endif
# if NDIM > 2
            fvflux(:,1:NXB,1:NYB,0:NZB,KAXIS) = flux(:,1:NXB,1:NYB,1:NZB+1,KAXIS,lb)
# endif

            !! ------------------------------------------------------------------ !!
            !! Get a local copy.

            call Grid_getBlkPtr(blockID,U_ptr,CENTER)
            U = U_ptr
            call Grid_releaseBlkPtr(blockID,U_ptr,CENTER)

            !! Apply unitsystem conversion (e.g. CGS, SI, ...).
            Uprim(DENS_PRIM,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
            Uprim(VELX_PRIM,:,:,:) = U(VELX_VAR,:,:,:)/hy_vref
            Uprim(VELY_PRIM,:,:,:) = U(VELY_VAR,:,:,:)/hy_vref
            Uprim(VELZ_PRIM,:,:,:) = U(VELZ_VAR,:,:,:)/hy_vref
            Uprim(PRES_PRIM,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
            Uprim(MAGX_PRIM,:,:,:) = U(MAGX_VAR,:,:,:)/hy_bref
            Uprim(MAGY_PRIM,:,:,:) = U(MAGY_VAR,:,:,:)/hy_bref
            Uprim(MAGZ_PRIM,:,:,:) = U(MAGZ_VAR,:,:,:)/hy_bref
            Uprim(GLMP_PRIM,:,:,:) = U(GLMP_VAR,:,:,:)/hy_bref
            Uprim(ENER_PRIM,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
            Uprim(GAMC_PRIM,:,:,:) = U(GAMC_VAR,:,:,:)

# if NSPECIES > 0
            Uprim(SPEC_PRIM_BEGIN:SPEC_PRIM_END,:,:,:) = U(SPECIES_BEGIN:SPECIES_END,:,:,:)
# endif
# if NMASS_SCALARS > 0
            Uprim(MSCL_PRIM_BEGIN:MSCL_PRIM_END,:,:,:) = U(MASS_SCALARS_BEGIN:MASS_SCALARS_END,:,:,:)
# endif

            !! ------------------------------------------------------------------ !!
            !! Calculate conservative variables and DG surface flux.

            do k = NZB_LO,NZB_HI; do j = NYB_LO,NYB_HI; do i = NXB_LO,NXB_HI;
                Ucons(:,i,j,k) = prim2cons(Uprim(:,i,j,k))
            end do; end do; end do;

            !! ------------------------------------------------------------------ !!

            !dgerhs = HUGE(1.0)
            !fverhs = HUGE(1.0)

            dgerhs = 0.0
            fverhs = 0.0

            if (dgfv_dg_active .and. dgfv_fv_active) then
                vlblend = blends(:,:,:,lb)
                vldnelb = 1.0-vlblend
            else
                vlblend = 1.0
                vldnelb = 0.0
            end if

            if (dgfv_dg_active) then
                if (DGFV_FLUX_DIFFERENCING) then
                    call calc_dg_rhs_flux_differencing(Ucons,fvflux,vlblend,sdx,dgrhs,dgerhs)
                else
                    call calc_dg_rhs_weak_form(Ucons,fvflux,vlblend,sdx,dgrhs,dgerhs,dgdensflux)
                end if
            end if

            !! ------------------------------------------------------------------ !!

            ! fv_recycle = .false.

            ! LABEL_RECYLCE_DGFV continue !! redo computation with 100% finite-volume scheme

            ! if (fv_recycle) then
            !     vlblend = 0.0
            !     blends(:,:,:,lb) = 0.0
            !     dgrhs = 0.0
            !     rhs = 0.0
            ! end if

            !! ------------------------------------------------------------------ !!

            if (dgfv_dg_active .and. dgfv_fv_active) then

                if (any(vlblend < 1.0)) then
                    call calc_fv_rhs(Uprim,fvflux,sdx,fvrhs,fverhs)

                    ! vldnelb = 1.0 - vlblend

                    rhs = 0.0
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            rhs(:,iii,jjj,kkk) = vldnelb(i,j,k)*fvrhs(:,iii,jjj,kkk) + vlblend(i,j,k)*dgrhs(:,ii,jj,kk,i,j,k)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
                
                else

                    rhs = 0.0
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                            rhs(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*(k-1)+kk) = dgrhs(:,ii,jj,kk,i,j,k)
                        end do; end do; end do;
                    end do; end do; end do;
                end if

            else if (dgfv_dg_active) then

!                 if (.false.) then
!                 !if (any(dgerhs > 0.0)) then
!                     !! This branch is for testing the sledgehammer approach!
! 
!                     block
!                     real :: vfflux(NFLUXES,0:NXB,0:NYB,0:NZB,IAXIS:KAXIS)
!                     real :: vferhs(1:NXE,1:NYE,1:NZE)
!                     real, parameter :: eps = 1e-20
!                     real, parameter :: tau = 50.0
! 
!                     call calc_fv_surface_flux(Uprim,vfflux)
! 
!                     sfblend = 1.0
!                     sfdnelb = 1.0 - sfblend
! 
! # if NDIM > 0
!                     do k = 1,NZE; do j = 1,NYE; do i = 0,NXE;
!                         do kk = 1,NZN; do jj = 1,NYN;
!                             vfflux(:,NXN*i,NYN*(j-1)+jj,NZN*(k-1)+kk,IAXIS) = fvflux(:,NXN*i,NYN*(j-1)+jj,NZN*(k-1)+kk,IAXIS)
!                         end do; end do;
!                     end do; end do; end do;
! # endif
! # if NDIM > 1
!                     do k = 1,NZE; do j = 0,NYE; do i = 1,NXE;
!                         do kk = 1,NZN; do ii = 1,NXN;
!                             vfflux(:,NXN*(i-1)+ii,NYN*j,NZN*(k-1)+kk,JAXIS) = fvflux(:,NXN*(i-1)+ii,NYN*j,NZN*(k-1)+kk,JAXIS)
!                         end do; end do;
!                     end do; end do; end do;
! # endif
! # if NDIM > 2
!                     do k = 0,NZE; do j = 1,NYE; do i = 1,NXE;
!                         do jj = 1,NYN; do ii = 1,NXN;
!                             vfflux(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*k,KAXIS) = fvflux(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*k,KAXIS)
!                         end do; end do;
!                     end do; end do; end do;
! # endif
! 
!                     call calc_fv_rhs(Uprim,vfflux,sdx,fvrhs,fverhs)
! 
!                     vferhs = 0.0
!                     do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!                         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
!                             associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
!                             vferhs(i,j,k) = vferhs(i,j,k) + fverhs(iii,jjj,kkk)
!                             end associate
!                         end do; end do; end do;
!                     end do; end do; end do;
! 
!                     do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!                         vlblend(i,j,k) = (eps + abs(vferhs(i,j,k))) / (eps + abs(vferhs(i,j,k)) + tau*max(0.0,dgerhs(i,j,k)))
!                     end do; end do; end do;
! 
!                     vldnelb = 1.0 - vlblend
! 
!                     do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!                         dgerhs(i,j,k) = vldnelb(i,j,k)*vferhs(i,j,k) + vlblend(i,j,k)*dgerhs(i,j,k)
!                     end do; end do; end do;
! 
!                     rhs = 0.0
!                     do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
!                         do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
!                             associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
!                             rhs(:,iii,jjj,kkk) = vldnelb(i,j,k)*fvrhs(:,iii,jjj,kkk) + vlblend(i,j,k)*dgrhs(:,ii,jj,kk,i,j,k)
!                             end associate
!                         end do; end do; end do;
!                     end do; end do; end do;
!                     end block
!                     
!                 else

                    rhs = 0.0
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                            rhs(:,NXN*(i-1)+ii,NYN*(j-1)+jj,NZN*(k-1)+kk) = dgrhs(:,ii,jj,kk,i,j,k)
                        end do; end do; end do;
                    end do; end do; end do;

!                end if

            else if (dgfv_fv_active) then

                call calc_fv_rhs(Uprim,fvflux,sdx,rhs,fverhs)

            end if

            !! ------------------------------------------------------------------ !!

# if defined(GLMP_VAR)
            !! Divergence error damping.
            rhs(GLMP_FLUX,:,:,:) = rhs(GLMP_FLUX,:,:,:) - hy_damp*Uprim(GLMP_PRIM,1:NXB,1:NYB,1:NZB)
# endif

            !! ------------------------------------------------------------------ !!
            !! Calculate surface flux difference and correct via applying a source term.

            if (hy_fluxCorrect) then
                sfflux = 0.0
                crflux = 0.0

                do axis = 1,NDIM
                    call Grid_getFluxData(blockID,axis,sfflux(:,:,:,:,axis),(/GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC/))
                end do

                do h = 1,NFLUXES
# if NDIM > 0
                    crflux(h) = crflux(h) - sdx(IAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,    1,1:NYB,1:NZB,IAXIS,lb)) - sum(sfflux(h,    1,1:NYB,1:NZB,IAXIS)))
                    crflux(h) = crflux(h) + sdx(IAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,NXB+1,1:NYB,1:NZB,IAXIS,lb)) - sum(sfflux(h,NXB+1,1:NYB,1:NZB,IAXIS)))
# endif
# if NDIM > 1
                    crflux(h) = crflux(h) - sdx(JAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,1:NXB,    1,1:NZB,JAXIS,lb)) - sum(sfflux(h,1:NXB,    1,1:NZB,JAXIS)))
                    crflux(h) = crflux(h) + sdx(JAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,1:NXB,NYB+1,1:NZB,JAXIS,lb)) - sum(sfflux(h,1:NXB,NYB+1,1:NZB,JAXIS)))
# endif
# if NDIM > 2
                    crflux(h) = crflux(h) - sdx(JAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,1:NXB,1:NYB,    1,KAXIS,lb)) - sum(sfflux(h,1:NXB,1:NYB,    1,KAXIS)))
                    crflux(h) = crflux(h) + sdx(JAXIS)/REAL(NXB*NYB*NZB)*(sum(flux(h,1:NXB,1:NYB,NZB+1,KAXIS,lb)) - sum(sfflux(h,1:NXB,1:NYB,NZB+1,KAXIS)))
# endif
                end do

                do h = 1,NFLUXES
                    rhs(h,:,:,:) = rhs(h,:,:,:) + crflux(h)
                end do
            end if

# if NSPECIES > 0 || NMASS_SCALARS > 0

# if NDIM > 0
            fvdensflux(0:NXB,1:NYB,1:NZB,IAXIS) = fvflux(DENS_FLUX,0:NXB,1:NYB,1:NZB,IAXIS)
# endif
# if NDIM > 1
            fvdensflux(1:NXB,0:NYB,1:NZB,JAXIS) = fvflux(DENS_FLUX,1:NXB,0:NYB,1:NZB,JAXIS)
# endif
# if NDIM > 2
            fvdensflux(1:NXB,1:NYB,0:NZB,KAXIS) = fvflux(DENS_FLUX,1:NXB,1:NYB,0:NZB,KAXIS)
# endif

            if (hy_fluxCorrect) then
# if NDIM > 0
                fvdensflux(  0,1:NYB,1:NZB,IAXIS) = sfflux(DENS_FLUX,    1,1:NYB,1:NZB,IAXIS)
                fvdensflux(NXB,1:NYB,1:NZB,IAXIS) = sfflux(DENS_FLUX,NXB+1,1:NYB,1:NZB,IAXIS)
# endif
# if NDIM > 1
                fvdensflux(1:NXB,  0,1:NZB,JAXIS) = sfflux(DENS_FLUX,1:NXB,    1,1:NZB,JAXIS)
                fvdensflux(1:NXB,NYB,1:NZB,JAXIS) = sfflux(DENS_FLUX,1:NXB,NYB+1,1:NZB,JAXIS)
# endif
# if NDIM > 2
                fvdensflux(1:NXB,1:NYB,  0,KAXIS) = sfflux(DENS_FLUX,1:NXB,1:NYB,    1,KAXIS)
                fvdensflux(1:NXB,1:NYB,NZB,KAXIS) = sfflux(DENS_FLUX,1:NXB,1:NYB,NZB+1,KAXIS)
# endif
            end if

            if (dgfv_dg_active .and. dgfv_fv_active) then

                if (any(vlblend < 1.0)) then

# if NDIM > 0
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN-1;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,IAXIS) = vldnelb(i,j,k)*fvflux(DENS_FLUX,iii,jjj,kkk,IAXIS) + vlblend(i,j,k)*dgdensflux(ii,jj,kk,i,j,k,IAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 1
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN-1; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,JAXIS) = vldnelb(i,j,k)*fvflux(DENS_FLUX,iii,jjj,kkk,JAXIS) + vlblend(i,j,k)*dgdensflux(ii,jj,kk,i,j,k,JAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 2
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN-1; do jj = 1,NYN; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,KAXIS) = vldnelb(i,j,k)*fvflux(DENS_FLUX,iii,jjj,kkk,KAXIS) + vlblend(i,j,k)*dgdensflux(ii,jj,kk,i,j,k,KAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
                else
# if NDIM > 0
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN-1;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,IAXIS) = dgdensflux(ii,jj,kk,i,j,k,IAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 1
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN-1; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,JAXIS) = dgdensflux(ii,jj,kk,i,j,k,JAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
# if NDIM > 2
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN-1; do jj = 1,NYN; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            fvdensflux(iii,jjj,kkk,KAXIS) = dgdensflux(ii,jj,kk,i,j,k,KAXIS)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
# endif
                end if

            else if (dgfv_dg_active) then

# if NDIM > 0
                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN-1;
                        associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                        fvdensflux(iii,jjj,kkk,IAXIS) = dgdensflux(ii,jj,kk,i,j,k,IAXIS)
                        end associate
                    end do; end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 1
                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do jj = 1,NYN-1; do ii = 1,NXN;
                        associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                        fvdensflux(iii,jjj,kkk,JAXIS) = dgdensflux(ii,jj,kk,i,j,k,JAXIS)
                        end associate
                    end do; end do; end do;
                end do; end do; end do;
# endif
# if NDIM > 2
                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN-1; do jj = 1,NYN; do ii = 1,NXN;
                        associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                        fvdensflux(iii,jjj,kkk,KAXIS) = dgdensflux(ii,jj,kk,i,j,k,KAXIS)
                        end associate
                    end do; end do; end do;
                end do; end do; end do;
# endif
            end if

            call calc_species_rhs(Uprim,fvdensflux,sdx,rhs)
# endif

            !! ------------------------------------------------------------------ !!
            !! Calculate source terms.

# if defined(GRAVITY) || defined(TREERAY)
            sacc = 0.0

# ifdef GRAVITY
            if (hy_useGravity) then
# if NDIM > 0
                do k = 1,NZB; do j = 1,NYB;
                    call Gravity_accelOneRow((/j+NGUARD,k+NGUARD/),SWEEP_X,blockID,size(gacx,dim=1),gacx(:,j,k))
                end do; end do;
# endif
# if NDIM > 1
                do k = 1,NZB; do i = 1,NXB;
                    call Gravity_accelOneRow((/i+NGUARD,k+NGUARD/),SWEEP_Y,blockID,size(gacy,dim=1),gacy(:,i,k))
                end do; end do;
# endif
# if NDIM > 2
                do j = 1,NYB; do i = 1,NXB;
                    call Gravity_accelOneRow((/i+NGUARD,j+NGUARD/),SWEEP_Z,blockID,size(gacz,dim=1),gacz(:,i,j))
                end do; end do;
# endif
                do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                    sacc(i,j,k,IAXIS) = sacc(i,j,k,IAXIS) + gacx(i,j,k)/hy_gref
                    sacc(i,j,k,JAXIS) = sacc(i,j,k,JAXIS) + gacy(j,i,k)/hy_gref
                    sacc(i,j,k,KAXIS) = sacc(i,j,k,KAXIS) + gacz(k,i,j)/hy_gref
                end do; end do; end do;
            endif
# endif /* GRAVITY */

# ifdef TREERAY
            if (hy_useTreeRay) then
# if NDIM > 0
                do k = 1,NZB; do j = 1,NYB;
                    call TreeRay_accelOneRow((/j+NGUARD,k+NGUARD/),SWEEP_X,blockID,size(tacx,dim=1),tacx(:,j,k))
                end do; end do;
# endif
# if NDIM > 1
                do k = 1,NZB; do i = 1,NXB;
                    call TreeRay_accelOneRow((/i+NGUARD,k+NGUARD/),SWEEP_Y,blockID,size(tacy,dim=1),tacy(:,i,k))
                end do; end do;
# endif
# if NDIM > 2
                do j = 1,NYB; do i = 1,NXB;
                    call TreeRay_accelOneRow((/i+NGUARD,j+NGUARD/),SWEEP_Z,blockID,size(tacz,dim=1),tacz(:,i,j))
                end do; end do;
# endif
                do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                    sacc(i,j,k,IAXIS) = sacc(i,j,k,IAXIS) + tacx(i,j,k)/hy_gref
                    sacc(i,j,k,JAXIS) = sacc(i,j,k,JAXIS) + tacy(j,i,k)/hy_gref
                    sacc(i,j,k,KAXIS) = sacc(i,j,k,KAXIS) + tacz(k,i,j)/hy_gref
                end do; end do; end do;
            end if
# endif /* TREERAY */

!! Following is just for debugging.
# ifdef ACCX_VAR
            do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                U(ACCX_VAR,i,j,k) = sacc(i,j,k,IAXIS)
            end do; end do; end do;
# endif
# ifdef ACCY_VAR
            do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                U(ACCY_VAR,i,j,k) = sacc(i,j,k,JAXIS)
            end do; end do; end do;
# endif
# ifdef ACCZ_VAR
            do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                U(ACCZ_VAR,i,j,k) = sacc(i,j,k,KAXIS)
            end do; end do; end do;
# endif
       
            if (hy_useGravity .or. hy_useTreeRay) then
                do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                    rhs(XMOM_FLUX,i,j,k) = rhs(XMOM_FLUX,i,j,k) + Ucons(DENS_CONS,i,j,k)*sacc(i,j,k,IAXIS)
                    rhs(YMOM_FLUX,i,j,k) = rhs(YMOM_FLUX,i,j,k) + Ucons(DENS_CONS,i,j,k)*sacc(i,j,k,JAXIS)
                    rhs(ZMOM_FLUX,i,j,k) = rhs(ZMOM_FLUX,i,j,k) + Ucons(DENS_CONS,i,j,k)*sacc(i,j,k,KAXIS)

                    rhs(ENER_FLUX,i,j,k) = rhs(ENER_FLUX,i,j,k) + Ucons(MOMX_CONS,i,j,k)*sacc(i,j,k,IAXIS)
                    rhs(ENER_FLUX,i,j,k) = rhs(ENER_FLUX,i,j,k) + Ucons(MOMY_CONS,i,j,k)*sacc(i,j,k,JAXIS)
                    rhs(ENER_FLUX,i,j,k) = rhs(ENER_FLUX,i,j,k) + Ucons(MOMZ_CONS,i,j,k)*sacc(i,j,k,KAXIS)
                end do; end do; end do;
            end if
# endif /* defined(GRAVITY) || defined(TREERAY) */

            !! Used for convergence tests.
            ! call add_residual(blockID,timeEndAdv,dt_inout,rkstage,rhs)

            !! ------------------------------------------------------------------ !!
            !! Evolve in time.

            select case (DGFV_RUNGE_KUTTA)

                case (DGFV_RUNGE_KUTTA_EULER)
                    !! 1st order Euler timestepping.
                    Unew(1:NFLUXES,:,:,:) = Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB) + dt*rhs(1:NFLUXES,:,:,:)

                case (DGFV_RUNGE_KUTTA_RALSTON)
                    !! 2nd(/3rd) order Heun-type Runge-Kutta.
                    associate(uu => Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              un => Unew (1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              r1 => rkcache(1:NFLUXES,:,:,:,lb,1))

                    if (rkstage == 1) then
                        r1 = uu
                        un = uu + 2.0/3.0*dt*rhs
                    else
                        un = 5.0/8.0*r1 + 3.0/8.0*uu + 3.0/4.0*dt*rhs
                    end if
                    end associate

                case (DGFV_RUNGE_KUTTA_SSP_22)
                    !! 2nd order SSP Runge-Kutta.
                    associate(uu => Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              un => Unew (1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              r1 => rkcache(1:NFLUXES,:,:,:,lb,1))

                    if (rkstage == 1) then
                        r1 = uu
                        un = uu + dt*rhs
                    else
                        un = 0.5*(r1 + uu + dt*rhs)
                    end if
                    end associate

                case (DGFV_RUNGE_KUTTA_SSP_33)
                    !! 3rd order strong-stability prexerving Runge-Kutta.
                    associate(uu => Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              un => Unew (1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              r1 => rkcache(1:NFLUXES,:,:,:,lb,1))

                    select case(rkstage)
                        case(1)
                            r1 = uu
                            un = uu + dt*rhs

                        case(2)
                            un = 3.0/4.0*r1 + 1.0/4.0*uu + 1.0/4.0*dt*rhs

                        case(3)
                            un = 1.0/3.0*r1 + 2.0/3.0*uu + 2.0/3.0*dt*rhs
                    end select
                    end associate

                case (DGFV_RUNGE_KUTTA_SSP_34)
                    !! 3rd order / 4 stages strong-stability preserving Runge-Kutta.
                    associate(uu => Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              un => Unew (1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              r1 => rkcache(1:NFLUXES,:,:,:,lb,1))

                    select case(rkstage)
                        case(1)
                            r1 = uu
                            un = uu + 0.5*dt*rhs

                        case(2)
                            un = uu + 0.5*dt*rhs

                        case(3)
                            un = 2.0/3.0*r1 + 1.0/3.0*uu + 1.0/6.0*dt*rhs

                        case(4)
                            un = uu + 0.5*dt*rhs
                    end select
                    end associate

                case (DGFV_RUNGE_KUTTA_SSP_45)
                    !! 4th order strong-stability prexerving Runge-Kutta.
                    associate(uu => Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              un => Unew (1:NFLUXES,1:NXB,1:NYB,1:NZB), &
                              r1 => rkcache(1:NFLUXES,:,:,:,lb,1), &
                              r2 => rkcache(1:NFLUXES,:,:,:,lb,2))

                    select case(rkstage)
                        case(1)
                            r1 =                                    uu
                            un =                                    uu + RK_SSP_45_A3 * dt * rhs

                        case(2)
                            un = RK_SSP_45_B1 * r1 + RK_SSP_45_B2 * uu + RK_SSP_45_B3 * dt * rhs
                            r2 = RK_SSP_45_E4 * un

                        case(3)
                            un = RK_SSP_45_C1 * r1 + RK_SSP_45_C2 * uu + RK_SSP_45_C3 * dt * rhs
                            r2 = r2 + RK_SSP_45_E5 * un

                        case(4)
                            un = RK_SSP_45_D1 * r1 + RK_SSP_45_D2 * uu + RK_SSP_45_D3 * dt * rhs
                            r2 = r2 + RK_SSP_45_E6 * dt * rhs

                        case(5)
                            un = RK_SSP_45_E1 * r1 + RK_SSP_45_E2 * uu + RK_SSP_45_E3 * dt * rhs + r2
                    end select
                    end associate

                case (DGFV_RUNGE_KUTTA_LOW_STORAGE_45)
                    !! 4th order low-storage Runge-Kutta.
                    if (rkstage == 1) then
                        rkcache(1:NFLUXES,:,:,:,lb,1) = rhs(1:NFLUXES,:,:,:)
                    else
                        rkcache(1:NFLUXES,:,:,:,lb,1) = RK_LS_45_A(rkstage)*rkcache(1:NFLUXES,:,:,:,lb,1) + rhs(1:NFLUXES,:,:,:)
                    end if

                    Unew(1:NFLUXES,:,:,:) = Ucons(1:NFLUXES,1:NXB,1:NYB,1:NZB) + RK_LS_45_C(rkstage)*dt*rkcache(1:NFLUXES,:,:,:,lb,1)

            end select

            !! ------------------------------------------------------------------ !!

! # if NSPECIES > 0
!             Unew(DENS_CONS,:,:,:) = SUM(Unew(SPEC_CONS_BEGIN:SPEC_CONS_END,:,:,:))
! # endif
            Unew(GAMC_CONS,1:NXB,1:NYB,1:NZB) = Ucons(GAMC_CONS,1:NXB,1:NYB,1:NZB)

            !! ------------------------------------------------------------------ !!
            !! Enforce posivite densities and pressures.

            ! if (dgfv_dg_active .and. dgfv_fv_active .and. .not.fv_recycle) then
            !     if (.not.check_positivity(Unew)) then
            !         fv_recycle = .true.
            !         goto LABEL_RECYLCE_DGFV
            !     end if
            ! end if

            if (dgfv_enforce_positivity) then
                call enforce_positivity(Unew)
            else
                ok = check_positivity(Unew)
                if (.not.ok) then
                    call Driver_abortFlash(' [DGFV] Unphysical states detected!')
                end if
            end if

            !! ------------------------------------------------------------------ !!
            !! Convert back to primitive form.

            do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                Uprim(:,i,j,k) = cons2prim(Unew(:,i,j,k))
            end do; end do; end do;

            !! ------------------------------------------------------------------ !!
            !! Compute CFL constrained timestep.

            if (rkstage == nstages) then
                do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                    vmax = maxval(abs(Uprim(VELX_PRIM:VELZ_PRIM,i,j,k)))
                    lmax = non_directional_Fastest_Signal_Speed(Uprim(:,i,j,k))
                    dtmin = hy_CFL*hy_tref*deltamin/(vmax + lmax)

                    !! Divergence cleaning speed.
                    hy_ch_local = max(hy_ch_local,lmax)
# ifdef DTMN_VAR
                    U(DTMN_VAR,i,j,k) = dtmin
# endif
                    if (dtmin <= hy_dtmin) then
                        hy_dtmin = dtmin
                        hy_dtminloc(1) = i
                        hy_dtminloc(2) = j
                        hy_dtminloc(3) = k
                        hy_dtminloc(4) = blockID
                        hy_dtminloc(5) = hy_meshMe
                    end if
                end do; end do; end do;
            end if

            !! ------------------------------------------------------------------ !!
            !! Transform to physical domain.

            U(DENS_VAR,1:NXB,1:NYB,1:NZB) = hy_dref*Uprim(DENS_PRIM,1:NXB,1:NYB,1:NZB)
            U(VELX_VAR,1:NXB,1:NYB,1:NZB) = hy_vref*Uprim(VELX_PRIM,1:NXB,1:NYB,1:NZB)
            U(VELY_VAR,1:NXB,1:NYB,1:NZB) = hy_vref*Uprim(VELY_PRIM,1:NXB,1:NYB,1:NZB)
            U(VELZ_VAR,1:NXB,1:NYB,1:NZB) = hy_vref*Uprim(VELZ_PRIM,1:NXB,1:NYB,1:NZB)
            U(PRES_VAR,1:NXB,1:NYB,1:NZB) = hy_pref*Uprim(PRES_PRIM,1:NXB,1:NYB,1:NZB)
            U(ENER_VAR,1:NXB,1:NYB,1:NZB) = hy_eref*Uprim(ENER_PRIM,1:NXB,1:NYB,1:NZB)
            U(MAGX_VAR,1:NXB,1:NYB,1:NZB) = hy_bref*Uprim(MAGX_PRIM,1:NXB,1:NYB,1:NZB)
            U(MAGY_VAR,1:NXB,1:NYB,1:NZB) = hy_bref*Uprim(MAGY_PRIM,1:NXB,1:NYB,1:NZB)
            U(MAGZ_VAR,1:NXB,1:NYB,1:NZB) = hy_bref*Uprim(MAGZ_PRIM,1:NXB,1:NYB,1:NZB)
            U(GLMP_VAR,1:NXB,1:NYB,1:NZB) = hy_bref*Uprim(GLMP_PRIM,1:NXB,1:NYB,1:NZB)
            U(EINT_VAR,1:NXB,1:NYB,1:NZB) = hy_eref*Uprim(PRES_PRIM,1:NXB,1:NYB,1:NZB) &
                /(Uprim(GAMC_PRIM,1:NXB,1:NYB,1:NZB)-1.0)/Uprim(DENS_PRIM,1:NXB,1:NYB,1:NZB)

# if NSPECIES > 0
            U(SPECIES_BEGIN:SPECIES_END,1:NXB,1:NYB,1:NZB) = Uprim(SPEC_PRIM_BEGIN:SPEC_PRIM_END,1:NXB,1:NYB,1:NZB)
# endif
# if NMASS_SCALARS > 0
            U(MASS_SCALARS_BEGIN:MASS_SCALARS_END,1:NXB,1:NYB,1:NZB) = Uprim(MSCL_PRIM_BEGIN:MSCL_PRIM_END,1:NXB,1:NYB,1:NZB)
# endif

            !! ------------------------------------------------------------------ !!
            !! ------------------------------------------------------------------ !!

# ifdef ERHS_VAR
            if (dgfv_dg_active .and. dgfv_fv_active) then
            ! if (.false.) then

                if (any(vlblend < 0.9999)) then
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            U(ERHS_VAR,iii,jjj,kkk) = vldnelb(i,j,k)*fverhs(i,j,k) + vlblend(i,j,k)*dgerhs(i,j,k)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
                else
                    do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                        do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                            associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                            U(ERHS_VAR,iii,jjj,kkk) = dgerhs(i,j,k)
                            end associate
                        end do; end do; end do;
                    end do; end do; end do;
                end if

            ! else if (dgfv_dg_active) then
            else if (.false.) then
            ! else if (.true.) then

                block
                real :: vferhs(NXE,NYE,NZE)

                vferhs = 0.0
                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                        vferhs(i,j,k) = vferhs(i,j,k) + fverhs(iii,jjj,kkk)
                        end associate
                    end do; end do; end do;
                end do; end do; end do;


                do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                    do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                        associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                        ! U(ERHS_VAR,iii,jjj,kkk) = dgerhs(i,j,k)
                        ! U(ERHS_VAR,iii,jjj,kkk) = min(fverhs(i,j,k),dgerhs(i,j,k))
                        ! U(ERHS_VAR,iii,jjj,kkk) = vldnelb(i,j,k)*vferhs(i,j,k) + vlblend(i,j,k)*dgerhs(i,j,k)
                        U(ERHS_VAR,iii,jjj,kkk) = vldnelb(i,j,k)*fverhs(iii,jjj,kkk) + vlblend(i,j,k)*dgerhs(i,j,k)
                        end associate
                    end do; end do; end do;
                end do; end do; end do;

                end block

            ! else if (dgfv_fv_active) then
            else if (.true.) then

                do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
                    U(ERHS_VAR,i,j,k) = fverhs(i,j,k)
                end do; end do; end do;

            end if
# endif

# ifdef BLND_VAR
            do k = 1,NZE; do j = 1,NYE; do i = 1,NXE;
                do kk = 1,NZN; do jj = 1,NYN; do ii = 1,NXN;
                    associate(iii => NXN*(i-1)+ii, jjj => NYN*(j-1)+jj, kkk => NZN*(k-1)+kk)
                    U(BLND_VAR,iii,jjj,kkk) = vlblend(i,j,k)
                    end associate
                end do; end do; end do;
            end do; end do; end do;
# endif

            !! ------------------------------------------------------------------ !!
            !! Write data back to database.

            call Grid_getBlkPtr(blockID,U_ptr,CENTER)
            U_ptr = U

            !! Needed for 'Eos_wrapped' call.
            call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

            !! ------------------------------------------------------------------ !!

# if NSPECIES > 0
            if (hy_limit_abundance) then
                ! call Grid_renormAbundance(blockID,blkLimits,U_ptr)
                call Grid_limitAbundance(blkLimits,U_ptr)
            end if
# endif

            call Grid_releaseBlkPtr(blockID,U_ptr,CENTER)

            !! ------------------------------------------------------------------ !!

            call Eos_wrapped(HYDRO_EOS_MODE,blkLimits,blockID)

        end do !! lb = 1,blockCount
    
        call Timers_stop("DGFV_leafs")

        call Grid_notifySolnDataUpdate()

    end do !! rkstages

    if (allocated(flux)) deallocate(flux)
    if (allocated(blends)) deallocate(blends)
    if (allocated(rkcache)) deallocate(rkcache)

    if (DGFV_MHD_DIV_CLEANING) then
        !! Update cleaning speed for next solver round.
        call MPI_AllReduce(hy_ch_local,hy_ch,1,FLASH_REAL,MPI_MAX,dr_globalComm,error)

        hy_ch = hy_chscalingfactor*hy_ch
        hy_damp = hy_ch * 2.0/0.18

        !! Reset value for next evolution loop.
        hy_ch_local = 0.0
    end if

    call Timers_stop ("DGFV")

end subroutine Hydro

! # if 0
! subroutine add_residual(blockID,t,dt,rkstage,rhs)
! 
!     use Hydro_data, only: nstages
!     use Hydro_data, only: nrkcaches
! 
!     use Hydro_data, only: dgfv_dg_active
!     use Hydro_data, only: dgfv_fv_active
!     use Hydro_data, only: dgfv_mhd_non_cons
!     use Hydro_data, only: dgfv_mhd_div_cleaning
!     use Hydro_data, only: dgfv_enforce_positivity
! 
!     use Hydro_data, only: DGFV_RUNGE_KUTTA
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_EULER
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_RALSTON
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_22
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_33
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_45
!     use Hydro_data, only: DGFV_RUNGE_KUTTA_LOW_STORAGE_45
! 
!     use Hydro_data, only: RK_LS_45_A, RK_LS_45_B, RK_LS_45_C
! 
!     use Hydro_data, only: dt_coeffs => RK_SSP_45_dt_coeffs
! 
!     use Hydro_data, only: RK_SSP_45_A1, RK_SSP_45_A2, RK_SSP_45_A3
!     use Hydro_data, only: RK_SSP_45_B1, RK_SSP_45_B2, RK_SSP_45_B3
!     use Hydro_data, only: RK_SSP_45_C1, RK_SSP_45_C2, RK_SSP_45_C3
!     use Hydro_data, only: RK_SSP_45_D1, RK_SSP_45_D2, RK_SSP_45_D3
!     use Hydro_data, only: RK_SSP_45_E1, RK_SSP_45_E2, RK_SSP_45_E3
!     use Hydro_data, only: RK_SSP_45_E4, RK_SSP_45_E5, RK_SSP_45_E6
! 
!     use Hydro_data, only: NXE,NYE,NZE
!     use Hydro_data, only: NXN,NYN,NZN
! 
!     use Hydro_data, only: NXE_LO,NXE_HI
!     use Hydro_data, only: NYE_LO,NYE_HI
!     use Hydro_data, only: NZE_LO,NZE_HI
! 
!     use Hydro_data, only: NXN_LO,NXN_HI
!     use Hydro_data, only: NYN_LO,NYN_HI
!     use Hydro_data, only: NZN_LO,NZN_HI
! 
!     use Hydro_data, only: NXB_LO,NXB_HI
!     use Hydro_data, only: NYB_LO,NYB_HI
!     use Hydro_data, only: NZB_LO,NZB_HI
! 
!     use Hydro_data, only: NOR,SOU,WES,EAS,FRO,BAC
! 
!     use Hydro_data, only: N_NODES
!     use Hydro_data, only: dgnodes,dgvolumes
! 
!     use Simulation_manufacsol_mhd, only: calc_manufacsol
! 
!     use Grid_interface, ONLY : Grid_fillGuardCells,   &
!                                Grid_getDeltas,        &
!                                Grid_getBlkIndexLimits,&
!                                Grid_getBlkBoundBox,   &
!                                Grid_getCellCoords,    &
!                                Grid_getFluxData,      &
!                                Grid_putFluxData,      &
!                                Grid_conserveFluxes,   &
!                                Grid_getBlkPtr,        &
!                                Grid_releaseBlkPtr,    &
!                                Grid_renormAbundance,  &
!                                Grid_limitAbundance,   &
!                                Grid_notifySolnDataUpdate
! 
!     integer, intent(in) :: blockID, rkstage
!     real, intent(in)    :: t,dt
!     real, intent(inout) :: rhs(NFLUXES,NXB,NYB,NZB)
! 
!     real :: x,y,z
!   
!     real, dimension(N_NODES,N_NODES,N_NODES) :: xx,yy,zz
! 
!     integer :: i,j,k
!     integer :: ii,jj,kk
! 
!     real, dimension(N_NODES,N_NODES,N_NODES) :: h
!     real, dimension(N_NODES,N_NODES,N_NODES) :: diffh
! 
!     real :: rktime
!     real :: boundBox(2,MDIM)
!     real :: deltas(MDIM)
! 
!     real :: r(NFLUXES,N_NODES,N_NODES,N_NODES)
! 
!     ! write (*,'(a8,99(ES12.4))') 'nodes', 0.5*dgnodes
!     ! write (*,'(a8,99(ES12.4))') 'weights', s, sum(0.5*dgweights)
! 
!     ! return
! 
!     call Grid_getDeltas(blockID,deltas)
!     call Grid_getBlkBoundBox(blockId,boundBox)
! 
!     !rktime = t + dt*dt_coeffs(rkstage)
!     rktime = 0.0
! 
!     do k = 1,NZB; do j = 1,NYB; do i = 1,NXB;
! 
!         x = boundBox(1,IAXIS) + ((i-1) + 0.5)*deltas(IAXIS)
!         y = boundBox(1,JAXIS) + ((j-1) + 0.5)*deltas(JAXIS)
!         z = boundBox(1,KAXIS) + ((k-1) + 0.5)*deltas(KAXIS)
! 
!         do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
! 
!             xx(ii,jj,kk) = x + 0.5*dgnodes(ii)*deltas(IAXIS)
!             yy(ii,jj,kk) = y + 0.5*dgnodes(jj)*deltas(JAXIS)
!             zz(ii,jj,kk) = z + 0.5*dgnodes(kk)*deltas(KAXIS)
! 
!         end do; end do; end do;
! 
!         do kk = 1,N_NODES; do jj = 1,N_NODES; do ii = 1,N_NODES;
!             r(:,ii,jj,kk) = calc_manufacsol(rktime,xx(ii,jj,kk),yy(ii,jj,kk),zz(ii,jj,kk))
!         end do; end do; end do;
! 
!         ! write (*,'(99(ES12.4))') rhs(DENS_FLUX,i,j,k), sum(dgvolumes*r(DENS_FLUX,:,:,:)), rhs(DENS_FLUX,i,j,k) + sum(dgvolumes*r(DENS_FLUX,:,:,:))
! 
!         rhs(DENS_FLUX,i,j,k) = rhs(DENS_FLUX,i,j,k) + sum(dgvolumes*r(DENS_FLUX,:,:,:))
!         rhs(XMOM_FLUX,i,j,k) = rhs(XMOM_FLUX,i,j,k) + sum(dgvolumes*r(XMOM_FLUX,:,:,:))
!         rhs(YMOM_FLUX,i,j,k) = rhs(YMOM_FLUX,i,j,k) + sum(dgvolumes*r(YMOM_FLUX,:,:,:))
!         rhs(ZMOM_FLUX,i,j,k) = rhs(ZMOM_FLUX,i,j,k) + sum(dgvolumes*r(ZMOM_FLUX,:,:,:))
!         rhs(ENER_FLUX,i,j,k) = rhs(ENER_FLUX,i,j,k) + sum(dgvolumes*r(ENER_FLUX,:,:,:))
!         rhs(MAGX_FLUX,i,j,k) = rhs(MAGX_FLUX,i,j,k) + sum(dgvolumes*r(MAGX_FLUX,:,:,:))
!         rhs(MAGY_FLUX,i,j,k) = rhs(MAGY_FLUX,i,j,k) + sum(dgvolumes*r(MAGY_FLUX,:,:,:))
!         rhs(MAGZ_FLUX,i,j,k) = rhs(MAGZ_FLUX,i,j,k) + sum(dgvolumes*r(MAGZ_FLUX,:,:,:))
! 
!     end do; end do; end do;
! 
! end subroutine
! # endif
