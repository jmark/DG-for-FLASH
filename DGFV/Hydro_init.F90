!!****if* source/physics/Hydro/HydroMain/split/ES/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  Hydro_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the entropy-stable MHD unit
!!
!!
!! ARGUMENTS
!!
!!
!!***

subroutine Hydro_init()

    use Hydro_data, ONLY : hy_units, hy_cfl, hy_Rconst, &
                           hy_dtmin, hy_xref, hy_vref, hy_dref,      &
                           hy_mref,  hy_tref, hy_eref, hy_nref,      &
                           hy_pref,  hy_gref, hy_qref, hy_kref,      &
                           hy_bref,  hy_units,hy_useGravity,         &
                           hy_useTreeRay,                            &
                           hy_smallpres,   hy_smalldens,             &
                           hy_fluxCorrect, hy_irenorm,               &
                           hy_meshMe, hy_meshNumProcs,               &
                           hy_ch, hy_chscalingfactor, hy_ch_local, hy_damp, &
                           hy_EntropyScheme_str, hy_EntropyScheme,   &
                           hy_psidissipation, hy_debug
    ! use Hydro_data

    use Hydro_data, only: dgfv_fv_at_boundary

    use Hydro_data, only: hy_write_integral_quantities
    use Hydro_data, only: hy_limit_abundance

    use Hydro_data, only: nstages
    use Hydro_data, only: nrkcaches

    use Hydro_data, only: dgfv_dg_active
    use Hydro_data, only: dgfv_fv_active
    use Hydro_data, only: dgfv_mhd_non_cons
    use Hydro_data, only: dgfv_mhd_div_cleaning
    use Hydro_data, only: dgfv_enforce_positivity
    use Hydro_data, only: DGFV_FLUX_DIFFERENCING
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION_ACTIVE
    use Hydro_data, only: DGFV_ENTROPY_CORRECTION_CUTOFF
    use Hydro_data, only: DGFV_ENTROPY_BOUNDARY_PROJECTION

    use Hydro_data, only: dgfv_runge_kutta
    use Hydro_data, only: DGFV_RUNGE_KUTTA_EULER
    use Hydro_data, only: DGFV_RUNGE_KUTTA_RALSTON
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_22
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_33
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_34
    use Hydro_data, only: DGFV_RUNGE_KUTTA_SSP_45
    use Hydro_data, only: DGFV_RUNGE_KUTTA_LOW_STORAGE_45

    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
    use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
    use Grid_interface, ONLY:   Grid_setFluxHandling
    use Driver_data, ONLY : dr_dt, dr_restart

    use Hydro_data, only: dgnodes
    use Hydro_data, only: dgweights
    use Hydro_data, only: dgvolumes
    use Hydro_data, only: xweights
    use Hydro_data, only: yweights
    use Hydro_data, only: zweights

    use Hydro_data, only: recoMat       !! mean values -> node values reconstruction matrix
    use Hydro_data, only: projMat       !! node values -> mean values projection matrix
    use Hydro_data, only: projTam       !! node values -> mean values projection matrix

    use Hydro_data, only: weakdiffMat   !! weak-form differentation matrix
    use Hydro_data, only: surfVecM      !! weak-form surface (vector) operator (left/minus face)
    use Hydro_data, only: surfVecP      !! weak-form surface (vector) operator (right/plus face)

    use Hydro_data, only: splitdiffMat     !! split-form differentation matrix
    use Hydro_data, only: fluxdiffMat     !! split-form differentation matrix

    use Hydro_data, only: toBoundaryVecM  !! interpolation operator (left/minus face)
    use Hydro_data, only: toBoundaryVecP  !! interpolation operator (right/plus face)
    use Hydro_data, only: toMidpointsMat  !! interpolation operator from nodes to midpoints

    use Hydro_data, only: toInnerItfsMat

    use Hydro_data, only: diffMat

    use Hydro_data, only: N_NODES
    use Hydro_data, only: NXN,NYN,NZN

    use dgfv_utils_mod, only: LegendreGaussNodesAndWeights
    use dgfv_utils_mod, only: LegendreGaussLobattoNodesAndWeights
    use dgfv_utils_mod, only: fill_projection_matrix
    use dgfv_utils_mod, only: fill_Diff_Matrix
    use dgfv_utils_mod, only: fill_weak_Diff_Matrix
    use dgfv_utils_mod, only: fill_surface_vectors
    use dgfv_utils_mod, only: fill_generalized_SBP_surface_Matrix
    use dgfv_utils_mod, only: lagrange_polynomial
    use dgfv_utils_mod, only: invert_matrix
    use dgfv_utils_mod, only: pp

    ! use dgfv_fluxes_mod, only: prim2cons
    ! use dgfv_fluxes_mod, only: cons2prim
    ! use dgfv_fluxes_mod, only: cons2evec
    ! use dgfv_fluxes_mod, only: evec2prim
    ! use dgfv_fluxes_mod, only: evec2cons

    use Hydro_data, only: DGFV_NODE_TYPE
    use Hydro_data, only: DGFV_NODE_TYPE_GAUSS
    use Hydro_data, only: DGFV_NODE_TYPE_LOBATTO

    implicit none

    integer :: i,j,k

    character(len=32) :: dgfv_runge_kutta_string
    character(len=256) :: temp_string

#include "constants.h"
#include "Flash.h"

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)

  call RuntimeParameters_get('UnitSystem', hy_units)
  call RuntimeParameters_get('smlrho', hy_smalldens)
  call RuntimeParameters_get('smallP', hy_smallpres)
  call RuntimeParameters_get('cfl', hy_cfl)

  !hy_cfl = hy_cfl * 0.5 **(NDIM-1)
  hy_cfl = hy_cfl / REAL(NDIM)

  call RuntimeParameters_get('flux_correct', hy_fluxCorrect)

  if (NDIM > 1) then
     if (hy_fluxCorrect) then
        call Grid_setFluxHandling('consv_flux_densities')
     end if
  end if

# ifdef GLMP_VAR
  ! call RuntimeParameters_get("psidissipation", hy_psidissipation)
  ! if(hy_psidissipation) then
  !   if(hy_meshMe .EQ. MASTER_PE) print*,"[Hydro_init] INFO: hy_psidissipation is ENABLED"
  ! else
  !   if(hy_meshMe .EQ. MASTER_PE) print*,"[Hydro_init] INFO: hy_psidissipation is DISABLED"
  ! end if

  call RuntimeParameters_get('chscaling', hy_chscalingfactor)
1 format(" [Hydro_init] INFO: c_h scaling factor is ",F0.4)
  if (hy_meshMe .EQ. MASTER_PE) print 1, hy_chscalingfactor
# endif


  hy_useGravity = .false.
# ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
# endif
  hy_useTreeRay = .false.
# ifdef TREERAY
  call RuntimeParameters_get("useTreeRay", hy_useTreeRay)
# endif

  call RuntimeParameters_get('irenorm', hy_irenorm)
  call PhysicalConstants_get("ideal gas constant", hy_Rconst)


    if (.not.dr_restart) then
       hy_dtmin = HUGE(hy_dtmin)
    else
       hy_dtmin = dr_dt
    endif

    call RuntimeParameters_get("hy_write_integral_quantities", hy_write_integral_quantities)
    call RuntimeParameters_get("hy_limit_abundance", hy_limit_abundance)

    call RuntimeParameters_get("dgfv_fv_active", dgfv_fv_active)
    call RuntimeParameters_get("dgfv_dg_active", dgfv_dg_active)

    call RuntimeParameters_get("dgfv_fv_at_boundary", dgfv_fv_at_boundary)

    call RuntimeParameters_get("dgfv_enforce_positivity", dgfv_enforce_positivity)
    call RuntimeParameters_get("dgfv_mhd_div_cleaning", dgfv_mhd_div_cleaning)
    call RuntimeParameters_get("dgfv_mhd_non_cons", dgfv_mhd_non_cons)

    call RuntimeParameters_get("dgfv_flux_differencing",DGFV_FLUX_DIFFERENCING)
    call RuntimeParameters_get("dgfv_entropy_correction",DGFV_ENTROPY_CORRECTION)
    call RuntimeParameters_get("dgfv_entropy_correction_active",DGFV_ENTROPY_CORRECTION_ACTIVE)
    call RuntimeParameters_get("dgfv_entropy_correction_cutoff",DGFV_ENTROPY_CORRECTION_CUTOFF)
    call RuntimeParameters_get("dgfv_entropy_boundary_projection",DGFV_ENTROPY_BOUNDARY_PROJECTION)

    call RuntimeParameters_get("dgfv_runge_kutta", dgfv_runge_kutta_string)

    if (trim(dgfv_runge_kutta_string) == "EULER") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_EULER
        nstages     = 1
        nrkcaches   = 0

    else if (trim(dgfv_runge_kutta_string) == "RALSTON") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_RALSTON
        nstages     = 2
        nrkcaches   = 1

    else if (trim(dgfv_runge_kutta_string) == "SSP_22") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_SSP_22
        nstages     = 2
        nrkcaches   = 1

    else if (trim(dgfv_runge_kutta_string) == "SSP_33") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_SSP_33
        nstages     = 3
        nrkcaches   = 1

    else if (trim(dgfv_runge_kutta_string) == "SSP_34") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_SSP_34
        nstages     = 4
        nrkcaches   = 1

    else if (trim(dgfv_runge_kutta_string) == "SSP_45") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_SSP_45
        nstages     = 5
        nrkcaches   = 2

    else if (trim(dgfv_runge_kutta_string) == "LOW_STORAGE_45") then
        dgfv_runge_kutta = DGFV_RUNGE_KUTTA_LOW_STORAGE_45
        nstages     = 5
        nrkcaches   = 1

    else
        call Driver_abortFlash(" [Hydro_init]: Given Runge-Kutta scheme is of" // &
            " unknown type: Options are EULER, RALSTON, SSP_22, SSP_33, SSP_34, SSP_45 (default), LOW_STORAGE_45.")
    end if

    if (dgfv_dg_active .and. dgfv_runge_kutta < DGFV_RUNGE_KUTTA_SSP_33) then
        call Driver_abortFlash(" [Hydro_init]: You must not use the DG method with a Runge-Kutta scheme of less than third order!")
    end if

    call RuntimeParameters_get("dgfv_characteristic_time", hy_tref)
    call RuntimeParameters_get("dgfv_characteristic_density", hy_dref)
    call RuntimeParameters_get("dgfv_characteristic_dimension", hy_xref)

    hy_vref = hy_xref/hy_tref
    hy_mref = hy_xref*hy_vref
    hy_eref = hy_vref*hy_vref
    hy_nref = hy_dref*hy_vref*hy_xref
    hy_pref = hy_dref*hy_vref*hy_vref
    hy_gref = hy_vref*hy_vref/hy_xref

    hy_qref = hy_vref*hy_vref/hy_Rconst
    hy_kref = hy_dref*hy_vref*hy_xref*hy_Rconst

    if ( hy_units == "SI" .or. hy_units == "si" ) then
      hy_bref = hy_vref*sqrt(4.0*PI*hy_dref*1.e-7)
    else if ( hy_units == "CGS" .or. hy_units == "cgs" ) then
      hy_bref = hy_vref*sqrt(4.0*PI*hy_dref)
    else
      hy_bref = hy_vref*sqrt(hy_dref)
    end if

  ! write (*,*) hy_bref, hy_tref, hy_eref

  ! Initialize hyperbolic cleaning speed with zero
  ! We cannot do anything in the first step, as we don't know
  ! the maximum fast magnetoacoustic speed at the beginning of
  ! the iteration
  hy_ch = 0.0
  hy_ch_local = 0.0
  hy_damp = 0.0

  ! calculate proper hydrodynamics after start
  ! hy_onlyDT = .false.


    ! ======================================================================= !!

! # if HYDRO_DG_KERNEL == HYDRO_DG_WEAK
!     call LegendreGaussNodesAndWeights(N_NODES-1,dgnodes,dgweights)
! # elif HYDRO_DG_KERNEL == HYDRO_DG_SPLIT
!     call LegendreGaussLobattoNodesAndWeights(N_NODES-1,dgnodes,dgweights)
! # else
! # error: NoDG kernel given.
! # endif

    call RuntimeParameters_get("dgfv_node_type", temp_string)

    if (trim(temp_string) == "GAUSS") then
        DGFV_NODE_TYPE = DGFV_NODE_TYPE_GAUSS

    else if (trim(temp_string) == "LOBATTO") then
        DGFV_NODE_TYPE = DGFV_NODE_TYPE_LOBATTO

    else
        call Driver_abortFlash(" [Hydro_init]: Unknown DGFV node type: " // trim(temp_string))
    end if


    if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_GAUSS) then
        call LegendreGaussNodesAndWeights(N_NODES-1,dgnodes,dgweights)
    else if (DGFV_NODE_TYPE == DGFV_NODE_TYPE_LOBATTO) then
        call LegendreGaussLobattoNodesAndWeights(N_NODES-1,dgnodes,dgweights)
    end if

# if 0
    write (*,'(a8,99(ES12.4))') 'nodes', 0.5*dgnodes
    write (*,'(a8,99(ES12.4))') 'weights', 0.5*dgweights, sum(0.5*dgweights)
    stop
# endif

    dgvolumes = 1.0
    xweights = 1.0
    yweights = 1.0
    zweights = 1.0

# if NDIM > 0
    do i = 1,NXN; do j = 1,NYN; do k = 1,NZN;
        dgvolumes(i,j,k) = dgvolumes(i,j,k) * 0.5*dgweights(i)
    end do; end do; end do;
    
    xweights = 0.5*dgweights/REAL(N_NODES)
# endif
# if NDIM > 1
    do i = 1,NXN; do j = 1,NYN; do k = 1,NZN;
        dgvolumes(i,j,k) = dgvolumes(i,j,k) * 0.5*dgweights(j)
    end do; end do; end do;

    yweights = 0.5*dgweights/REAL(N_NODES)
# endif
# if NDIM > 2
    do i = 1,NXN; do j = 1,NYN; do k = 1,NZN;
        dgvolumes(i,j,k) = dgvolumes(i,j,k) * 0.5*dgweights(k)
    end do; end do; end do;

    zweights = 0.5*dgweights/REAL(N_NODES)
# endif

    ! write (*,*) 'x', xweights
    ! write (*,*) 'y', yweights
    ! write (*,*) 'z', zweights

    call fill_projection_matrix(N_NODES,dgnodes,dgweights,projMat)

    recoMat = invert_matrix(N_NODES,projMat)

    projTam = transpose(projMat)

    ! recoMat = 0.0
    ! do i = 1,N_NODES;
    !     recoMat(i,i) = 1.0
    ! end do

    ! write (*,*)
    ! do i = 1,N_NODES
    !     write (*,'(a8,2x,99(F10.4))') 'recoMat',recoMat(i,:)
    ! end do;
    ! write (*,*)

    call fill_diff_Matrix(N_NODES-1,dgnodes,diffMat)
    diffMat = 2.0/REAL(N_NODES)*diffMat

    ! write (*,*)
    ! do i = 1,N_NODES
    !     write (*,'(a8,2x,99(F10.4))') 'diff',diffMat(i,:)
    ! end do;
    ! write (*,*)

    call fill_weak_Diff_Matrix(N_NODES,dgnodes,dgweights,weakdiffMat)
    weakdiffMat = weakdiffMat/N_NODES

    ! write (*,*)
    ! do i = 1,N_NODES
    !     write (*,'(a8,2x,99(F10.4))') 'weakdiff',weakdiffMat(i,:)
    ! end do;
    ! write (*,*)

    call fill_diff_Matrix(N_NODES-1,dgnodes,splitdiffMat)

    splitdiffMat = 2.0*splitdiffMat
    splitDiffMat(1,1) = splitDiffMat(1,1) + 1.0/dgweights(1)
    splitDiffMat(N_NODES,N_NODES) = splitDiffMat(N_NODES,N_NODES) - 1.0/dgweights(N_NODES)
    splitdiffMat = -2.0/REAL(N_NODES) * splitdiffMat

    ! write (*,*)
    ! do i = 1,N_NODES
    !     write (*,'(a8,2x,99(F10.4))') 'splitdiff',splitdiffMat(i,:)
    ! end do;
    ! write (*,*)

    ! stop

    ! write (*,*)
    ! do i = 1,N_NODES
    !     write (*,'(a8,2x,99(F10.4))') 'splitdiff',sum(splitdiffMat(i,:))
    ! end do;
    ! write (*,*)

    call fill_surface_vectors(N_NODES,dgnodes,dgweights,surfVecM,surfVecP)

    surfVecM = -surfVecM/N_NODES
    surfVecP = -surfVecP/N_NODES

    ! write (*,'(a8,2x,99(F10.4))') 'surfvecM',surfVecM
    ! write (*,'(a8,2x,99(F10.4))') 'surfvecP',surfVecP

    do i = 1,N_NODES;
        toBoundaryVecM(i) = lagrange_polynomial(N_NODES,dgnodes,i,-1.0)
        toBoundaryVecP(i) = lagrange_polynomial(N_NODES,dgnodes,i, 1.0)
    end do

    do i = 1,N_NODES; do j = 1,N_NODES;
        ! write (*,*) i, j, -1.0 + i*2.0/REAL(N_NODES) - 1.0/REAL(N_NODES)
        toMidpointsMat(i,j) = lagrange_polynomial(N_NODES,dgnodes,j,-1.0 + i*2.0/REAL(N_NODES) - 1.0/REAL(N_NODES))
    end do; end do;

    do i = 1,N_NODES-1; do j = 1,N_NODES;
        !write (*,*) i, j, -1.0 + i*2.0/REAL(N_NODES)
        toInnerItfsMat(i,j) = lagrange_polynomial(N_NODES,dgnodes,j,-1.0 + i*2.0/REAL(N_NODES))
    end do; end do;

    ! stop

    block
    real :: diffMat(N_NODES,N_NODES)
    real :: surfMat(N_NODES,N_NODES)

    call fill_diff_Matrix(N_NODES-1,dgnodes,diffMat)
    !call pp('diffMat', diffMat)

    call fill_generalized_SBP_surface_Matrix(N_NODES,dgnodes,dgweights,surfMat)
    ! dgsurfTam = transpose(dgsurfMat)
    ! call pp('dgSurfMat', dgsurfMat)

    fluxDiffMat = -2.0*(2.0 * diffMat - surfMat)/REAL(N_NODES)
    !call pp('fluxdiffMat', fluxdiffMat*REAL(N_NODES))
    end block

    ! block
    ! real :: p(NPRIM_VARS)
    ! real :: q(NPRIM_VARS)
    ! real :: u(NCONS_VARS)
    ! real :: w(NCONS_VARS)

    ! p = 0.0

    ! p(DENS_PRIM) = 1.5
    ! p(PRES_PRIM) = 0.71428

    ! p(VELX_PRIM) = 0.11
    ! p(VELY_PRIM) = 0.22
    ! p(VELZ_PRIM) = 0.33

    ! p(MAGX_PRIM) = 0.1
    ! p(MAGY_PRIM) = 0.2
    ! p(MAGZ_PRIM) = 0.3

    ! p(GLMP_PRIM) = 0.123

    ! p(GAMC_PRIM) = 1.4

    ! p(ENER_PRIM) = p(PRES_PRIM)/(p(GAMC_PRIM)-1.0)/p(DENS_PRIM) + 0.5*sum(p(VELX_PRIM:VELZ_PRIM)**2)

    ! u = prim2cons(p)

    ! ! u = (/1.123_dp, 2.0_dp, -3.0_dp, 0.4502_dp, 10.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.42_dp/)

    ! ! u(DENS_CONS) = 1.123
    ! ! !u(PRES_CONS) = 0.71428

    ! ! u(MOMX_CONS) =  2.0
    ! ! u(MOMY_CONS) = -3.0
    ! ! u(MOMZ_CONS) =  0.4502

    ! ! u(MAGX_CONS) = 0.1
    ! ! u(MAGY_CONS) = 0.2
    ! ! u(MAGZ_CONS) = 0.3

    ! ! u(GLMP_CONS) = 0.42

    ! ! u(GAMC_CONS) = 5.0/3.0

    ! ! u(ENER_CONS) = 10.0

    ! w = cons2evec(u)
    ! !q = evec2prim(w)
    ! q = evec2cons(w)

    ! write (*,'(99(ES12.4))') p
    ! write (*,'(99(ES12.4))') cons2prim(u)
    ! write (*,*)
    ! write (*,'(99(ES12.4))') u
    ! write (*,'(99(ES12.4))') q
    ! write (*,*)
    ! write (*,'(99(ES12.4))') w

    ! end block

end subroutine Hydro_init
