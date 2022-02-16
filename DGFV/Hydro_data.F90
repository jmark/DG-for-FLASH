!!****if* source/physics/Hydro/HydroMain/split/ES/Hydro_data
!!
!! NAME
!!
!!  Hydro_data
!!
!!
!! SYNOPSIS
!!
!!  use Hydro_data
!!
!!
!! DESCRIPTION
!!
!!  Placeholder module containing private 8Wave MHD solver-specific data
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

module Hydro_data

# include "Flash.h"
# include "constants.h"
# include "DGFV.h"

    implicit none

    integer, parameter :: dp = selected_real_kind(15)

    integer, save :: nstages
    integer, save :: nrkcaches

    real(dp), parameter :: RK_SSP_45_A0 = 1.00000000000000_dp
    real(dp), parameter :: RK_SSP_45_A1 = 1.00000000000000_dp
    real(dp), parameter :: RK_SSP_45_A2 = 0.00000000000000_dp
    real(dp), parameter :: RK_SSP_45_A3 = 0.39175222700392_dp

    real(dp), parameter :: RK_SSP_45_B0 = 0.39175222700392_dp
    real(dp), parameter :: RK_SSP_45_B1 = 0.44437049406734_dp
    real(dp), parameter :: RK_SSP_45_B2 = 0.55562950593266_dp
    real(dp), parameter :: RK_SSP_45_B3 = 0.36841059262959_dp

    real(dp), parameter :: RK_SSP_45_C0 = 0.58607968896780_dp
    real(dp), parameter :: RK_SSP_45_C1 = 0.62010185138540_dp
    real(dp), parameter :: RK_SSP_45_C2 = 0.37989814861460_dp
    real(dp), parameter :: RK_SSP_45_C3 = 0.25189177424738_dp

    real(dp), parameter :: RK_SSP_45_D0 = 0.47454236302687_dp
    real(dp), parameter :: RK_SSP_45_D1 = 0.17807995410773_dp
    real(dp), parameter :: RK_SSP_45_D2 = 0.82192004589227_dp
    real(dp), parameter :: RK_SSP_45_D3 = 0.54497475021237_dp

    real(dp), parameter :: RK_SSP_45_E0 = 0.93501063100924_dp
    real(dp), parameter :: RK_SSP_45_E1 = 0.00683325884039_dp
    real(dp), parameter :: RK_SSP_45_E2 = 0.34833675773694_dp
    real(dp), parameter :: RK_SSP_45_E3 = 0.22600748319395_dp

    real(dp), parameter :: RK_SSP_45_E4 = 0.51723167208978_dp
    real(dp), parameter :: RK_SSP_45_E5 = 0.12759831133288_dp
    real(dp), parameter :: RK_SSP_45_E6 = 0.08460416338212_dp

    real(dp), parameter :: RK_SSP_45_dt_coeffs(5) = (/RK_SSP_45_A0,RK_SSP_45_B0,RK_SSP_45_C0,RK_SSP_45_D0,RK_SSP_45_E0/)

    real(dp), parameter :: RK_LS_45_A_1 = REAL(  0.0_dp                                , dp)
    real(dp), parameter :: RK_LS_45_A_2 = REAL( -567301805773.0_dp /  1357537059087.0_dp,dp)
    real(dp), parameter :: RK_LS_45_A_3 = REAL(-2404267990393.0_dp /  2016746695238.0_dp,dp)
    real(dp), parameter :: RK_LS_45_A_4 = REAL(-3550918686646.0_dp /  2091501179385.0_dp,dp)
    real(dp), parameter :: RK_LS_45_A_5 = REAL(-1275806237668.0_dp /   842570457699.0_dp,dp)

    real(dp), parameter :: RK_LS_45_B_1 = REAL( 0.0_dp                                  ,dp)
    real(dp), parameter :: RK_LS_45_B_2 = REAL( 1432997174477.0_dp /  9575080441755.0_dp,dp)
    real(dp), parameter :: RK_LS_45_B_3 = REAL( 2526269341429.0_dp /  6820363962896.0_dp,dp)
    real(dp), parameter :: RK_LS_45_B_4 = REAL( 2006345519317.0_dp /  3224310063776.0_dp,dp)
    real(dp), parameter :: RK_LS_45_B_5 = REAL( 2802321613138.0_dp /  2924317926251.0_dp,dp)

    real(dp), parameter :: RK_LS_45_C_1 = REAL( 1432997174477.0_dp /  9575080441755.0_dp,dp)
    real(dp), parameter :: RK_LS_45_C_2 = REAL( 5161836677717.0_dp / 13612068292357.0_dp,dp)
    real(dp), parameter :: RK_LS_45_C_3 = REAL( 1720146321549.0_dp /  2090206949498.0_dp,dp)
    real(dp), parameter :: RK_LS_45_C_4 = REAL( 3134564353537.0_dp /  4481467310338.0_dp,dp)
    real(dp), parameter :: RK_LS_45_C_5 = REAL( 2277821191437.0_dp / 14882151754819.0_dp,dp)

    real(dp), parameter :: RK_LS_45_A(5) = (/ RK_LS_45_a_1,RK_LS_45_a_2,RK_LS_45_a_3,RK_LS_45_a_4,RK_LS_45_a_5 /)
    real(dp), parameter :: RK_LS_45_B(5) = (/ RK_LS_45_b_1,RK_LS_45_b_2,RK_LS_45_b_3,RK_LS_45_b_4,RK_LS_45_b_5 /)
    real(dp), parameter :: RK_LS_45_C(5) = (/ RK_LS_45_c_1,RK_LS_45_c_2,RK_LS_45_c_3,RK_LS_45_c_4,RK_LS_45_c_5 /)

    logical, save :: dgfv_dg_active
    logical, save :: dgfv_fv_active
    logical, save :: dgfv_mhd_non_cons
    logical, save :: dgfv_mhd_div_cleaning
    logical, save :: dgfv_enforce_positivity

    logical, save :: DGFV_FLUX_DIFFERENCING
    logical, save :: DGFV_ENTROPY_CORRECTION
    logical, save :: DGFV_ENTROPY_CORRECTION_ACTIVE
    logical, save :: DGFV_ENTROPY_CORRECTION_CUTOFF
    logical, save :: DGFV_ENTROPY_BOUNDARY_PROJECTION

    logical, save :: dgfv_fv_at_boundary

    integer, save :: dgfv_runge_kutta

    integer, save      :: DGFV_NODE_TYPE
    integer, parameter :: DGFV_NODE_TYPE_GAUSS      = 1
    integer, parameter :: DGFV_NODE_TYPE_LOBATTO    = 2

    integer, parameter :: DGFV_RUNGE_KUTTA_EULER            = 1
    integer, parameter :: DGFV_RUNGE_KUTTA_RALSTON          = 2
    integer, parameter :: DGFV_RUNGE_KUTTA_SSP_22           = 3
    integer, parameter :: DGFV_RUNGE_KUTTA_SSP_33           = 4
    integer, parameter :: DGFV_RUNGE_KUTTA_SSP_34           = 5
    integer, parameter :: DGFV_RUNGE_KUTTA_SSP_45           = 6
    integer, parameter :: DGFV_RUNGE_KUTTA_LOW_STORAGE_45   = 7

    integer, save :: hy_irenorm
    real,    save :: hy_Rconst
    real,    save :: hy_cfl
    real,    save :: hy_smallpres, hy_smalldens
    logical, save :: hy_useGravity,    &
                     hy_useTreeRay,    &
                     hy_fluxCorrect,   &
                     hy_psidissipation

    logical, save :: hy_write_integral_quantities

    logical, save :: hy_limit_abundance

    ! Energy switch
    real, PARAMETER :: hy_eswitch = 1.e-6

    ! System of units used
    character(4), save :: hy_units

    ! Storage for timestep calculation
    integer,save, DIMENSION(5) :: hy_dtminloc
    real, save :: hy_dtmin

    ! Everybody should know these!
    integer, save :: hy_meshNumProcs, hy_meshMe
    integer, save :: hy_gcMaskSize=NUNK_VARS+NDIM*NFACE_VARS

    ! Constants for non-dimensionalization
    real, save :: hy_xref
    real, save :: hy_tref
    real, save :: hy_dref
    real, save :: hy_vref
    real, save :: hy_pref
    real, save :: hy_eref
    real, save :: hy_qref
    real, save :: hy_bref
    real, save :: hy_gref
    real, save :: hy_mref
    real, save :: hy_nref
    real, save :: hy_kref

    real, save :: hy_ch, hy_ch_local, hy_chscalingfactor
      real, save :: hy_damp
    !real, save :: hy_LambdaMax

    character(len=MAX_STRING_LENGTH), save :: hy_EntropyScheme_str
    integer, save :: hy_EntropyScheme

    ! logical, save :: hy_onlyDT
    logical, save :: hy_debug

    ! integer, parameter :: N_NODES = NGUARD
    integer, parameter :: N_NODES = 4

# if NDIM > 0
    !! Number of DG elements per block.
    integer, parameter :: N_ELEMS = NXB/N_NODES

    integer, parameter :: NXE = NXB/N_NODES
    integer, parameter :: NXN = N_NODES
# endif

# if NDIM > 1
    integer, parameter :: NYE = NYB/N_NODES
    integer, parameter :: NYN = N_NODES
# else
    integer, parameter :: NYE = 1
    integer, parameter :: NYN = 1 
# endif

# if NDIM > 2
    integer, parameter :: NZE = NZB/N_NODES
    integer, parameter :: NZN = N_NODES
# else
    integer, parameter :: NZE = 1
    integer, parameter :: NZN = 1 
# endif

    integer, parameter :: NXE_LO = 1-K1D
    integer, parameter :: NXE_HI = 1+NXE*K1D

    integer, parameter :: NYE_LO = 1-K2D
    integer, parameter :: NYE_HI = 1+NYE*K2D

    integer, parameter :: NZE_LO = 1-K3D
    integer, parameter :: NZE_HI = 1+NZE*K3D

    integer, parameter :: NXN_LO =   1-N_NODES*K1D
    integer, parameter :: NXN_HI = NXB+N_NODES*K1D

    integer, parameter :: NYN_LO =   1-N_NODES*K2D
    integer, parameter :: NYN_HI = NYB+N_NODES*K2D

    integer, parameter :: NZN_LO =   1-N_NODES*K3D
    integer, parameter :: NZN_HI = NZB+N_NODES*K3D

    integer, parameter :: NXB_LO =   1-NGUARD*K1D
    integer, parameter :: NXB_HI = NXB+NGUARD*K1D

    integer, parameter :: NYB_LO =   1-NGUARD*K2D
    integer, parameter :: NYB_HI = NYB+NGUARD*K2D

    integer, parameter :: NZB_LO =   1-NGUARD*K3D
    integer, parameter :: NZB_HI = NZB+NGUARD*K3D

    !! face codes
    integer, parameter :: ZEN = 0
    integer, parameter :: NOR = 1
    integer, parameter :: SOU = 2
    integer, parameter :: WES = 3
    integer, parameter :: EAS = 4
    integer, parameter :: FRO = 5
    integer, parameter :: BAC = 6

    real, save :: dgnodes(N_NODES)
    real, save :: dgweights(N_NODES)

    real, save :: dgvolumes(NXN,NYN,NZN)

    real, save :: xweights(N_NODES)
    real, save :: yweights(N_NODES)
    real, save :: zweights(N_NODES)

    real, save :: recoMat(N_NODES,N_NODES)   !! mean values -> node values reconstruction matrix
    real, save :: projMat(N_NODES,N_NODES)   !! node values -> mean values projection matrix
    real, save :: projTam(N_NODES,N_NODES)   !! node values -> mean values projection matrix

    real, save :: weakdiffMat(N_NODES,N_NODES)   !! weak-form differentation matrix
    real, save :: splitdiffMat(N_NODES,N_NODES)   !! split-form differentation matrix

    real, save :: fluxdiffMat(N_NODES,N_NODES)   !! weak-form differentation matrix

    real, save :: diffMat(N_NODES,N_NODES)   !! weak-form differentation matrix

    real, save :: surfVecM(N_NODES)  !! weak-form surface (vector) operator (left/minus face)
    real, save :: surfVecP(N_NODES)  !! weak-form surface (vector) operator (right/plus face)

    real, save :: toBoundaryVecM(N_NODES)  !! interpolation operator (left/minus face)
    real, save :: toBoundaryVecP(N_NODES)  !! interpolation operator (right/plus face)

    real, save :: toMidpointsMat(N_NODES,N_NODES)  !! interpolation operator from nodes to midpoints

    real, save :: toInnerItfsMat(N_NODES-1,N_NODES)  !! interpolation operator from nodes to midpoints

    integer, parameter :: NPRIM_VARS = NFLUXES+2

    integer, parameter :: DENS_PRIM = DENS_FLUX
    integer, parameter :: VELX_PRIM = XMOM_FLUX
    integer, parameter :: VELY_PRIM = YMOM_FLUX
    integer, parameter :: VELZ_PRIM = ZMOM_FLUX
    integer, parameter :: ENER_PRIM = ENER_FLUX
    integer, parameter :: MAGX_PRIM = MAGX_FLUX
    integer, parameter :: MAGY_PRIM = MAGY_FLUX
    integer, parameter :: MAGZ_PRIM = MAGZ_FLUX
    integer, parameter :: GLMP_PRIM = GLMP_FLUX

    integer, parameter :: SPEC_PRIM_BEGIN = SPECIES_FLUX_BEGIN
    integer, parameter :: SPEC_PRIM_END = SPECIES_FLUX_END

    integer, parameter :: MSCL_PRIM_BEGIN = MASS_SCALARS_FLUX_BEGIN
    integer, parameter :: MSCL_PRIM_END = MASS_SCALARS_FLUX_END

    integer, parameter :: PRES_PRIM = NFLUXES+1
    integer, parameter :: GAMC_PRIM = NFLUXES+2


    integer, parameter :: NCONS_VARS = NFLUXES+2

    integer, parameter :: DENS_CONS = DENS_FLUX
    integer, parameter :: MOMX_CONS = XMOM_FLUX
    integer, parameter :: MOMY_CONS = YMOM_FLUX
    integer, parameter :: MOMZ_CONS = ZMOM_FLUX
    integer, parameter :: ENER_CONS = ENER_FLUX
    integer, parameter :: MAGX_CONS = MAGX_FLUX
    integer, parameter :: MAGY_CONS = MAGY_FLUX
    integer, parameter :: MAGZ_CONS = MAGZ_FLUX
    integer, parameter :: GLMP_CONS = GLMP_FLUX

    integer, parameter :: SPEC_CONS_BEGIN = SPECIES_FLUX_BEGIN
    integer, parameter :: SPEC_CONS_END = SPECIES_FLUX_END

    integer, parameter :: MSCL_CONS_BEGIN = MASS_SCALARS_FLUX_BEGIN
    integer, parameter :: MSCL_CONS_END = MASS_SCALARS_FLUX_END

    integer, parameter :: PRES_CONS = NFLUXES+1
    integer, parameter :: GAMC_CONS = NFLUXES+2

end module Hydro_data
