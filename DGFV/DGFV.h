!! -------------------------------------------------------------------------- !!
!! -------------------------------------------------------------------------- !!
!! ACHTUNG !! ACHTUNG !! ACHTUNG !!

!! This file is considered deprecated and will be deleted in the near future. Ignore
!! this file and forget about it.

!! -------------------------------------------------------------------------- !!
!! -------------------------------------------------------------------------- !!







!! -------------------------------------------------------------------------- !!
!! Nothing to configure here.

# define HYDRO_EOS_MODE MODE_DENS_EI

# define HYDRO_FV_RECONSTRUCTION_FIRST_ORDER    1
# define HYDRO_FV_RECONSTRUCTION_MINMOD         2
# define HYDRO_FV_RECONSTRUCTION_MONO_CENTRAL   3

! # define HYDRO_RUNGE_KUTTA_EULER                1
! # define HYDRO_RUNGE_KUTTA_RALSTON              2
! # define HYDRO_RUNGE_KUTTA_SSP_22               3
! # define HYDRO_RUNGE_KUTTA_SSP_33               4
! # define HYDRO_RUNGE_KUTTA_SSP_45               5
! # define HYDRO_RUNGE_KUTTA_LOW_STORAGE_45       6
! 
! # define HYDRO_DG_WEAK  1
! # define HYDRO_DG_SPLIT 2

!! -------------------------------------------------------------------------- !!
!! Choose between three possible modes: FV-only, DG-only or blending scheme.

! # define HYDRO_FV_ACTIVE 1
! # define HYDRO_DG_ACTIVE 1

!! -------------------------------------------------------------------------- !!
!! Only use a fourth-order RK (*_45) together with DG.
!! For FV-only mode all RK-schemes are ok.

!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_EULER
!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_RALSTON
!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_SSP_22
!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_SSP_33
!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_SSP_45
!# define HYDRO_RUNGE_KUTTA_SCHEME HYDRO_RUNGE_KUTTA_LOW_STORAGE_45

!! -------------------------------------------------------------------------- !!
!! Activate MHD divergence cleaning.

!# define HYDRO_MHD_DIV_CLEANING 1

!! Activate MHD-non-conservative terms. Highly recommended for MHD simulations.
!# define HYDRO_MHD_NON_CONS 1

!! -------------------------------------------------------------------------- !!
!! Ignore those options. Experimental!

!# define HYDRO_FV_RECONSTRUCTION HYDRO_FV_RECONSTRUCTION_FIRST_ORDER
# define HYDRO_FV_RECONSTRUCTION HYDRO_FV_RECONSTRUCTION_MINMOD
!# define HYDRO_FV_RECONSTRUCTION HYDRO_FV_RECONSTRUCTION_MONO_CENTRAL

!# define HYDRO_DG_KERNEL HYDRO_DG_WEAK
# define HYDRO_DG_WEAK_ENTROPY_CORRECTION 0
!# define HYDRO_DG_KERNEL HYDRO_DG_SPLIT
