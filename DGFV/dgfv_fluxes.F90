# include "Flash.h"
# include "constants.h"
# include "DGFV.h"

module dgfv_fluxes_mod

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

use Hydro_data, only: hy_ch

contains

pure function calc_pressure_cons(u) result(p)

    !! Calculate pressure (ideal gas law) from conservative state variables.

    implicit none

    real, intent(in) :: u(NCONS_VARS)
    real             :: p,k,e,m

    associate(&
        rho => u(DENS_CONS), &
        e   => u(ENER_CONS), &
        m1  => u(MOMX_CONS), &
        m2  => u(MOMY_CONS), &
        m3  => u(MOMZ_CONS), &
        b1  => u(MAGX_CONS), &
        b2  => u(MAGY_CONS), &
        b3  => u(MAGZ_CONS), &
        gam => u(GAMC_CONS), &
        psi => u(GLMP_CONS))

    !! kinetic energy
    k = 0.5/rho*(m1*m1 + m2*m2 + m3*m3)

    !! magnetic energy (+ psi energy)
    m = 0.5*(b1*b1 + b2*b2 + b3*b3 + psi*psi)

    !! pressure
    p = (gam-1.0)*(e - k - m)
    end associate

end function

!! TODO: Calculate proper GAMC from species.

pure function prim2cons(prim) result(cons)

    !! Convert from primitive to conservative state variables.

    implicit none

    real, intent(in) :: prim(NPRIM_VARS)
    real             :: cons(NCONS_VARS)

    cons(DENS_CONS) = prim(DENS_PRIM)
    cons(MOMX_CONS) = prim(DENS_PRIM)*prim(VELX_PRIM)
    cons(MOMY_CONS) = prim(DENS_PRIM)*prim(VELY_PRIM)
    cons(MOMZ_CONS) = prim(DENS_PRIM)*prim(VELZ_PRIM)
    cons(ENER_CONS) = prim(DENS_PRIM)*prim(ENER_PRIM) &
        + 0.5*(sum(prim(MAGX_PRIM:MAGZ_PRIM)**2) + prim(GLMP_PRIM)**2)
    cons(MAGX_CONS) = prim(MAGX_PRIM)
    cons(MAGY_CONS) = prim(MAGY_PRIM)
    cons(MAGZ_CONS) = prim(MAGZ_PRIM)
    cons(GLMP_CONS) = prim(GLMP_PRIM)
    cons(GAMC_CONS) = prim(GAMC_PRIM)
    !! This is just an auxiliary variable.
    cons(PRES_CONS) = prim(PRES_PRIM)

# if NSPECIES > 0
    cons(SPEC_CONS_BEGIN:SPEC_CONS_END) &
        = prim(DENS_PRIM)*prim(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
# endif
# if NMASS_SCALARS > 0
    cons(MSCL_CONS_BEGIN:MSCL_CONS_END) &
        = prim(DENS_PRIM)*prim(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
# endif

end function

pure function cons2prim(cons) result(prim)

    !! Convert from conservative to primitive state variables.

    implicit none

    real, intent(in) :: cons(NCONS_VARS)
    real             :: prim(NPRIM_VARS)

    real :: ekin,emag,ener,dens

# if NSPECIES > 0
    !! Get total mass
    dens = SUM(cons(SPEC_CONS_BEGIN:SPEC_CONS_END))
# else
    dens = cons(DENS_CONS)
# endif

    ekin = 0.5/dens*sum(cons(MOMX_CONS:MOMZ_CONS)**2)
    emag = 0.5*(sum(cons(MAGX_CONS:MAGZ_CONS)**2) + cons(GLMP_CONS)**2)
    ener = cons(ENER_CONS) - emag

    prim(PRES_PRIM) = (ener - ekin)*(cons(GAMC_CONS)-1.0)
    prim(ENER_PRIM) = ener/dens

    prim(DENS_PRIM) = dens
    prim(VELX_PRIM) = cons(MOMX_CONS)/dens
    prim(VELY_PRIM) = cons(MOMY_CONS)/dens
    prim(VELZ_PRIM) = cons(MOMZ_CONS)/dens

    prim(MAGX_PRIM) = cons(MAGX_CONS)
    prim(MAGY_PRIM) = cons(MAGY_CONS)
    prim(MAGZ_PRIM) = cons(MAGZ_CONS)
    prim(GLMP_PRIM) = cons(GLMP_CONS)

    prim(GAMC_PRIM) = cons(GAMC_CONS)

# if NSPECIES > 0
    !! Convert to normalized mass fractions.
    prim(SPEC_PRIM_BEGIN:SPEC_PRIM_END) &
        = cons(SPEC_CONS_BEGIN:SPEC_CONS_END)/dens
# endif
# if NMASS_SCALARS > 0
    prim(MSCL_PRIM_BEGIN:MSCL_PRIM_END) &
        = cons(MSCL_CONS_BEGIN:MSCL_CONS_END)/dens
# endif

end function

pure function prim2evec(u) result(w)

    !! Computes the entropy vector for given states.

    implicit none

    real, intent(IN)  :: U(NPRIM_VARS)
    real              :: w(NFLUXES)
    real              :: s,b

    associate(&
        rho     => U(DENS_PRIM), &
        p       => U(PRES_PRIM), &
        v1      => U(VELX_PRIM), &
        v2      => U(VELY_PRIM), &
        v3      => U(VELZ_PRIM), &
        B1      => U(MAGX_PRIM), &
        B2      => U(MAGY_PRIM), &
        B3      => U(MAGZ_PRIM), &
        psi     => U(GLMP_PRIM), &
        gamc    => U(GAMC_PRIM))

    s = LOG(p) - gamc*LOG(rho)
    b = rho/p !! 2*beta

    w = 0.0

    w(DENS_FLUX) =  ((gamc-s)/(gamc - 1.0) - 0.5*b*(v1*v1+v2*v2+v3*v3))
    w(XMOM_FLUX) =  b*v1
    w(YMOM_FLUX) =  b*v2
    w(ZMOM_FLUX) =  b*v3
    w(ENER_FLUX) = -b
    w(MAGX_FLUX) =  b*B1
    w(MAGY_FLUX) =  b*B2
    w(MAGZ_FLUX) =  b*B3
    w(GLMP_FLUX) =  b*psi
    !w(GAMC_CONS) =  gamc
    end associate

    !! This is just an auxiliary variable.
    !w(PRES_CONS) = u(PRES_PRIM)

! # if NSPECIES > 0
!     w(SPEC_CONS_BEGIN:SPEC_CONS_END) &
!         = u(DENS_PRIM)*u(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
! # endif
! # if NMASS_SCALARS > 0
!     w(MSCL_CONS_BEGIN:MSCL_CONS_END) &
!         = u(DENS_PRIM)*u(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
! # endif

end function

pure function cons2evec(u) result(w)

    !! Computes the entropy vector for given states.

    implicit none

    real, intent(IN)  :: U(NCONS_VARS)
    real              :: w(NCONS_VARS)
    real              :: s,b,p

    w = u

    associate(&
        rho     => U(DENS_CONS), &
        v1      => U(MOMX_CONS)/U(DENS_CONS), &
        v2      => U(MOMY_CONS)/U(DENS_CONS), &
        v3      => U(MOMZ_CONS)/U(DENS_CONS), &
        B1      => U(MAGX_CONS), &
        B2      => U(MAGY_CONS), &
        B3      => U(MAGZ_CONS), &
        psi     => U(GLMP_CONS), &
        gamc    => U(GAMC_CONS))

    p = calc_pressure_cons(U)
    s = LOG(p) - gamc*LOG(rho)
    b = rho/p !! 2*beta

    w(DENS_CONS) =  ((gamc-s)/(gamc - 1.0) - 0.5*b*(v1*v1+v2*v2+v3*v3))
    w(MOMX_CONS) =  b*v1
    w(MOMY_CONS) =  b*v2
    w(MOMZ_CONS) =  b*v3
    w(ENER_CONS) = -b
    w(MAGX_CONS) =  b*B1
    w(MAGY_CONS) =  b*B2
    w(MAGZ_CONS) =  b*B3
    w(GLMP_CONS) =  b*psi
    end associate

end function

pure function evec2cons(w) result(u)

    real, intent(IN)  :: w(NCONS_VARS)
    real              :: u(NCONS_VARS)

    real :: rho_sp, p_srho, v1,v2,v3,v2s2, s

    real :: kappa
    real, parameter :: mu0 = 1.0

    u = w

    kappa  =  w(GAMC_CONS)
    rho_sp = -w(ENER_CONS)
    p_srho = 1.0/rho_sp

    v1 = w(MOMX_CONS)*p_srho
    v2 = w(MOMY_CONS)*p_srho
    v3 = w(MOMZ_CONS)*p_srho

    v2s2 = 0.5*(v1*v1 + v2*v2 + v3*v3)

    s = kappa - ((w(DENS_CONS) + rho_sp*v2s2) * (kappa-1.0)) 

    u(DENS_CONS) = (exp(-s)*p_srho)**(1.0/(kappa-1.0))

    u(MOMX_CONS) = v1 * u(DENS_CONS)
    u(MOMY_CONS) = v2 * u(DENS_CONS)
    u(MOMZ_CONS) = v3 * u(DENS_CONS)

    u(MAGX_CONS) = p_srho * w(MAGX_CONS)
    u(MAGY_CONS) = p_srho * w(MAGY_CONS)
    u(MAGZ_CONS) = p_srho * w(MAGZ_CONS)

    u(GLMP_CONS) = p_srho * w(GLMP_CONS)

    u(ENER_CONS) = u(DENS_CONS)*(1.0/(kappa-1.0)*p_srho + v2s2) + 0.5/mu0*(SUM(u(MAGX_CONS:MAGZ_CONS)**2) + u(GLMP_CONS)**2)

# if NSPECIES > 0
    !! Convert to normalized mass fractions.
    u(SPEC_CONS_BEGIN:SPEC_CONS_END) &
        = u(DENS_CONS)*w(SPEC_CONS_BEGIN:SPEC_CONS_END)/sum(w(SPEC_CONS_BEGIN:SPEC_CONS_END))
# endif

end function

pure function evec2prim(w) result(p)

    real, intent(IN)  :: w(NCONS_VARS)
    real              :: p(NPRIM_VARS)

    real              :: u(NCONS_VARS)

    u = evec2cons(w)
    p = cons2prim(u)

end function

pure function isvalid_prim(prim) result(ok)

    !! Check if primitive state variables are physical.

    implicit none

    real, intent(in) :: prim(NPRIM_VARS)
    logical          :: ok

    real, parameter :: tol = 1e-50

    ok = .true.

    ok = ok .and. prim(DENS_PRIM) > tol .and. (prim(DENS_PRIM)-1.0 /= prim(DENS_PRIM))
    ok = ok .and. prim(PRES_PRIM) > tol .and. (prim(PRES_PRIM)-1.0 /= prim(PRES_PRIM))

! # if NSPECIES > 0
!     ok = ok .and. all(prim(SPEC_PRIM_BEGIN:SPEC_PRIM_END) >= 0.0)
! # endif
! # if NMASS_SCALARS > 0
!     ok = ok .and. all(prim(MSCL_PRIM_BEGIN:MSCL_PRIM_END) >= 0.0)
! # endif

end function

pure function isvalid_cons(cons) result(ok)

    !! Check if conservative state variables are physical.

    implicit none

    real, intent(in) :: cons(NCONS_VARS)
    logical          :: ok

    real, parameter :: tol = 1e-50

    real :: eint

    ok = .true.

    eint = calc_eint_cons(cons)

    ok = ok .and. cons(DENS_CONS) > tol .and. (cons(DENS_CONS)-1.0 /= cons(DENS_CONS))
    ok = ok .and. eint > tol .and. (eint-1.0 /= eint)

! # if NSPECIES > 0
!     ok = ok .and. all(cons(SPEC_CONS_BEGIN:SPEC_CONS_END) >= 0.0)
! # endif
! # if NMASS_SCALARS > 0
!     ok = ok .and. all(cons(MSCL_CONS_BEGIN:MSCL_CONS_END) >= 0.0)
! # endif

    contains

    pure function calc_eint_cons(u) result(eint)

        !! Calculate internal energy from conservative variables.

        implicit none

        real, intent(in) :: u(NCONS_VARS)
        real             :: eint,k,e,m

        associate(&
            rho => u(DENS_CONS), &
            e   => u(ENER_CONS), &
            m1  => u(MOMX_CONS), &
            m2  => u(MOMY_CONS), &
            m3  => u(MOMZ_CONS), &
            b1  => u(MAGX_CONS), &
            b2  => u(MAGY_CONS), &
            b3  => u(MAGZ_CONS), &
            psi => u(GLMP_CONS))

        !! kinetic energy
        k = 0.5/rho*(m1*m1 + m2*m2 + m3*m3)

        !! magnetic energy (+ psi energy)
        m = 0.5*(b1*b1 + b2*b2 + b3*b3 + psi*psi)

        !! internal energy
        eint = e - k - m
        end associate

    end function

end function

!! Standard advective MHD flux for primitive variables.
pure function xflux_prim(u) result(f)

    real, intent(in)    :: u(NPRIM_VARS)
    real                :: f(NFLUXES)
    real                :: etot,vels(3)

    real :: bb2,vb,pt,Ep,rhov1

    integer :: s

    f = 0.0

    associate(rho   => U(DENS_PRIM), &
            v1      => U(VELX_PRIM), &
            v2      => U(VELY_PRIM), &
            v3      => U(VELZ_PRIM), &
            pres    => U(PRES_PRIM), &
            ener    => U(ENER_PRIM), &
            b1      => U(MAGX_PRIM), &
            b2      => U(MAGY_PRIM), &
            b3      => U(MAGZ_PRIM), &
            psi     => U(GLMP_PRIM))

    rhov1 = rho*v1

    vb    = b1*v1+b2*v2+b3*v3
    bb2   = b1*b1+b2*b2+b3*b3

    pt = pres + 0.5*bb2
    Ep = ener*rho + pres + bb2 

    f(DENS_FLUX) = rhov1                   

    f(XMOM_FLUX) = rhov1*v1 + pt - b1*b1
    f(YMOM_FLUX) = rhov1*v2      - b1*b2
    f(ZMOM_FLUX) = rhov1*v3      - b1*b3

    f(ENER_FLUX) = Ep*v1         - b1*vb + hy_ch*b1*psi

    f(MAGX_FLUX) = hy_ch*psi
    f(MAGY_FLUX) = v1*b2 - b1*v2
    f(MAGZ_FLUX) = v1*b3 - b1*v3

    f(GLMP_FLUX) = hy_ch*b1

! # if NSPECIES > 0
!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = f(DENS_FLUX) * U(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
! # endif
! # if NMASS_SCALARS > 0
!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = f(DENS_FLUX) * U(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
! # endif
    end associate

end function

pure function yflux_prim(u) result(f)

    real, intent(in)    :: u(NPRIM_VARS)
    real                :: f(NFLUXES)

    f = yrotate_back_flux(xflux_prim(yrotate_prim(u)))

end function

pure function zflux_prim(u) result(f)

    real, intent(in)    :: u(NPRIM_VARS)
    real                :: f(NFLUXES)

    f = zrotate_back_flux(xflux_prim(zrotate_prim(u)))

end function

pure function Fastest_Signal_Speed_prim(U) result(c)

    real, intent(in)    :: U(NPRIM_VARS) !! primitive state variables
    real                :: c

    real :: cs2,va2,ca2,pres
    real :: astar

    cs2 = U(GAMC_PRIM)*U(PRES_PRIM)/U(DENS_PRIM)   !! sound wave speed
    ca2 = U(MAGX_PRIM)**2/U(DENS_PRIM)            !! Alfen wave speed

    va2 = SUM(U(MAGY_PRIM:MAGZ_PRIM)**2)/U(DENS_PRIM) + ca2

    astar = SQRT((cs2 + va2)**2 - 4.0*cs2*ca2)

    c = SQRT(0.5*(cs2 + va2 + astar))

end function

pure function non_directional_Fastest_Signal_Speed(U) result(c)

    !! This is a little hack. Here we overestimate the
    !! max. signal speed a little bit by discarding the
    !! directional information of the magnetic field.
    !! For CLF time step estimation this is ok.

    real, intent(in)    :: U(NPRIM_VARS) !! primitive state variables
    real                :: c

    real :: cs2,va2,ca2,pres,b_2(3)
    real :: astar

    cs2 = U(GAMC_PRIM)*U(PRES_PRIM)/U(DENS_PRIM)   !! sound wave speed

    b_2 = U(MAGX_PRIM:MAGZ_PRIM)**2
    ca2 = minval(b_2)/U(DENS_PRIM)               !! Alfen wave speed

    va2 = SUM(b_2)/U(DENS_PRIM) - ca2

    astar = SQRT((cs2 + va2)**2 - 4.0*cs2*ca2)

    c = SQRT(0.5*(cs2 + va2 + astar))

end function

pure function Non_Conservative_xFlux_prim(UI,UL,UR) result(fl)

    implicit none

    real, intent(in)    :: uI(NPRIM_VARS) !! inner state
    real, intent(in)    :: uL(NPRIM_VARS) !! left state
    real, intent(in)    :: uR(NPRIM_VARS) !! right state

    real                :: FL(NFLUXES)

    real :: Bavg,Gavg
    real :: mhdF(NFLUXES) !! MHD non-cons. flux
    real :: glmF(NFLUXES) !! GLM non-cons. flux

    mhdF = 0.0
    glmF = 0.0

    Bavg = 0.5*(UL(MAGX_PRIM) + UR(MAGX_PRIM))
    Gavg = 0.5*(UL(GLMP_PRIM) + UR(GLMP_PRIM))

    !! Powell term.
    mhdF(XMOM_FLUX:ZMOM_FLUX)       = Bavg * UI(MAGX_PRIM:MAGZ_PRIM)
    mhdF(ENER_FLUX)                 = Bavg * SUM(UI(MAGX_PRIM:MAGZ_PRIM)*UI(VELX_PRIM:VELZ_PRIM))
    mhdF(MAGX_FLUX:MAGZ_FLUX)       = Bavg * UI(VELX_PRIM:VELZ_PRIM)

    !! GLM term.
    glmF((/ENER_FLUX,GLMP_FLUX/))   = Gavg * UI(VELX_PRIM)*(/UI(GLMP_PRIM),1.0/)

    FL = mhdF + glmF

end function

pure function Non_Conservative_yFlux_prim(UI,UL,UR) result(F)

    implicit none

    real, intent(in)    :: UI(NPRIM_VARS)
    real, intent(in)    :: UL(NPRIM_VARS)
    real, intent(in)    :: UR(NPRIM_VARS)
    real                :: F(NFLUXES)

    f = yrotate_back_flux(Non_Conservative_xFlux_prim(yrotate_prim(uI),yrotate_prim(uL),yrotate_prim(uR)))

end function

pure function Non_Conservative_zFlux_prim(UI,UL,UR) result(F)

    implicit none

    real, intent(in)    :: UI(NPRIM_VARS)
    real, intent(in)    :: UL(NPRIM_VARS)
    real, intent(in)    :: UR(NPRIM_VARS)
    real                :: F(NFLUXES)

    f = zrotate_back_flux(Non_Conservative_xFlux_prim(zrotate_prim(uI),zrotate_prim(uL),zrotate_prim(uR)))

end function

! pure function Non_Conservative_xpot_prim(UI,UL,UR) result(epot)
! 
!     implicit none
! 
!     real, intent(in)    :: uI(NPRIM_VARS) !! inner state
!     real, intent(in)    :: uL(NPRIM_VARS) !! left state
!     real, intent(in)    :: uR(NPRIM_VARS) !! right state
!     real                :: epot
! 
!     real :: beta
! 
!     beta = 0.5*uI(DENS_PRIM)/uI(PRES_PRIM)
! 
!     epot = -2.0*beta*sum(uI(VELX_PRIM:VELZ_PRIM)*uI(MAGX_PRIM:MAGZ_PRIM)) * 0.5*(UL(MAGX_PRIM) + UR(MAGX_PRIM))
! 
! end function
! 
! pure function Non_Conservative_ypot_prim(UI,UL,UR) result(epot)
! 
!     implicit none
! 
!     real, intent(in)    :: UI(NPRIM_VARS)
!     real, intent(in)    :: UL(NPRIM_VARS)
!     real, intent(in)    :: UR(NPRIM_VARS)
!     real                :: epot
! 
!     epot = Non_Conservative_xpot_prim(yrotate_prim(uI),yrotate_prim(uL),yrotate_prim(uR))
! 
! end function
! 
! pure function Non_Conservative_zpot_prim(UI,UL,UR) result(epot)
! 
!     implicit none
! 
!     real, intent(in)    :: UI(NPRIM_VARS)
!     real, intent(in)    :: UL(NPRIM_VARS)
!     real, intent(in)    :: UR(NPRIM_VARS)
!     real                :: epot
! 
!     epot = Non_Conservative_xpot_prim(zrotate_prim(uI),zrotate_prim(uL),zrotate_prim(uR))
! 
! end function


pure function non_Conservative_xFlux_cons(UI,UL,UR) result(fl)

    implicit none

    real, intent(in)    :: uI(NCONS_VARS) !! inner state
    real, intent(in)    :: uL(NCONS_VARS) !! left state
    real, intent(in)    :: uR(NCONS_VARS) !! right state

    real                :: FL(NFLUXES)

    real :: Bavg,Gavg
    real :: mhdF(NFLUXES) !! MHD non-cons. flux
    real :: glmF(NFLUXES) !! GLM non-cons. flux
    real :: vels(3)

    vels = UI(MOMX_CONS:MOMZ_CONS)/UI(DENS_CONS)

    mhdF = 0.0
    glmF = 0.0

    Bavg = 0.5*(UL(MAGX_CONS) + UR(MAGX_CONS))
    Gavg = 0.5*(UL(GLMP_CONS) + UR(GLMP_CONS))

    !! Powell term.
    mhdF(XMOM_FLUX:ZMOM_FLUX)       = Bavg *     UI(MAGX_CONS:MAGZ_CONS)
    mhdF(ENER_FLUX)                 = Bavg * SUM(UI(MAGX_CONS:MAGZ_CONS)*vels)
    mhdF(MAGX_FLUX:MAGZ_FLUX)       = Bavg * vels 

    !! GLM term.
    glmF((/ENER_FLUX,GLMP_FLUX/))   = Gavg * vels(1)*(/UI(GLMP_CONS),1.0/)

    FL = mhdF + glmF

end function

pure function Non_Conservative_yFlux_cons(UI,UL,UR) result(F)

    implicit none

    real, intent(in)    :: UI(NCONS_VARS)
    real, intent(in)    :: UL(NCONS_VARS)
    real, intent(in)    :: UR(NCONS_VARS)
    real                :: F(NFLUXES)

    f = yrotate_back_flux(Non_Conservative_xFlux_cons(yrotate_cons(uI),yrotate_cons(uL),yrotate_cons(uR)))

end function

pure function Non_Conservative_zFlux_cons(UI,UL,UR) result(F)

    implicit none

    real, intent(in)    :: UI(NCONS_VARS)
    real, intent(in)    :: UL(NCONS_VARS)
    real, intent(in)    :: UR(NCONS_VARS)
    real                :: F(NFLUXES)

    f = zrotate_back_flux(Non_Conservative_xFlux_cons(zrotate_cons(uI),zrotate_cons(uL),zrotate_cons(uR)))

end function

pure function cons2xpot(U) result(ep)

    implicit none

    real, intent(IN)  :: U(NCONS_VARS)
    real              :: ep

    ep = prim2xpot(cons2prim(u))

end function

pure function cons2ypot(U) result(ep)

    implicit none

    real, intent(IN)  :: U(NCONS_VARS)
    real              :: ep

    ep = cons2xpot(yrotate_cons(u))

end function

pure function cons2zpot(U) result(ep)

    implicit none

    real, intent(IN)  :: U(NCONS_VARS)
    real              :: ep

    ep = cons2xpot(zrotate_cons(u))

end function

pure function prim2xpot(U) result(ep)

    use Hydro_data, only: hy_ch

    implicit none

    real, intent(IN)  :: U(NPRIM_VARS)
    real              :: ep

    real :: beta

    beta = 0.5*u(DENS_PRIM)/u(PRES_PRIM)

    ep = u(DENS_PRIM)*u(VELX_PRIM) &
        +     beta*u(VELX_PRIM)*sum(u(MAGX_PRIM:MAGZ_PRIM)**2) &
        + 2.0*beta*u(MAGX_PRIM)*hy_ch*u(GLMP_PRIM)

end function

pure function prim2ypot(U) result(ep)

    implicit none

    real, intent(IN)  :: U(NPRIM_VARS)
    real              :: ep

    ep = prim2xpot(yrotate_prim(u))

end function

pure function prim2zpot(U) result(ep)

    implicit none

    real, intent(IN)  :: U(NPRIM_VARS)
    real              :: ep

    ep = prim2xpot(zrotate_prim(u))

end function

pure function xefluxjump_prim(U_L,U_R,F) result(jump)

    use Hydro_data, only: hy_ch

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    real :: b_L
    real :: b_R

    b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    jump = dot_product(prim2evec(U_R) - prim2evec(U_L),f) - (prim2xpot(U_R) - prim2xpot(U_L))

    jump = jump + 0.5*(u_R(MAGX_PRIM) + u_L(MAGX_PRIM)) * (&
        2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) - &
        2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function yefluxjump_prim(U_L,U_R,F) result(jump)

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    real :: b_L
    real :: b_R

    b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    jump = dot_product(prim2evec(U_R) - prim2evec(U_L),f) - (prim2ypot(U_R) - prim2ypot(U_L))

    jump = jump + 0.5*(u_R(MAGY_PRIM) + u_L(MAGY_PRIM)) * (&
        2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) - &
        2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function zefluxjump_prim(U_L,U_R,F) result(jump)

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    real :: b_L
    real :: b_R

    b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    jump = dot_product(prim2evec(U_R) - prim2evec(U_L),f) - (prim2zpot(U_R) - prim2zpot(U_L))

    jump = jump + 0.5*(u_R(MAGZ_PRIM) + u_L(MAGZ_PRIM)) * (&
        2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) - &
        2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function xefluxjump_cons(U_L,U_R,F) result(jump)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    jump = xefluxjump_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure function yefluxjump_cons(U_L,U_R,F) result(jump)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    jump = yefluxjump_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure function zefluxjump_cons(U_L,U_R,F) result(jump)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: jump

    jump = zefluxjump_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure function xefluxmean_prim(U_L,U_R,F) result(mean)

    use Hydro_data, only: hy_ch

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    real :: b_L
    real :: b_R

    b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    mean = 0.5*(dot_product(prim2evec(U_R) + prim2evec(U_L),f) - (prim2xpot(U_R) + prim2xpot(U_L)))

    mean = mean - 0.5*(u_R(MAGX_PRIM) - u_L(MAGX_PRIM)) * (&
        2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) + &
        2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function yefluxmean_prim(U_L,U_R,F) result(mean)

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    real :: b_L
    real :: b_R

    b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    mean = 0.5*(dot_product(prim2evec(U_R) + prim2evec(U_L),f) - (prim2ypot(U_R) + prim2ypot(U_L)))

    mean = mean - 0.5*(u_R(MAGY_PRIM) - u_L(MAGY_PRIM)) * (&
        2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) + &
        2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function zefluxmean_prim(U_L,U_R,F) result(mean)

    implicit none

    real, intent(IN)  :: U_L(NPRIM_VARS)
    real, intent(IN)  :: U_R(NPRIM_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    ! real :: b_L
    ! real :: b_R

    ! b_L = 0.5*u_L(DENS_PRIM)/u_L(PRES_PRIM)
    ! b_R = 0.5*u_R(DENS_PRIM)/u_R(PRES_PRIM)

    mean = 0.5*(dot_product(prim2evec(U_R) + prim2evec(U_L),f) - (prim2zpot(U_R) + prim2zpot(U_L)))

    ! mean = mean + 0.5*(u_R(MAGZ_PRIM) - u_L(MAGZ_PRIM)) * (&
    !     2.0*b_R*sum(u_R(VELX_PRIM:VELZ_PRIM)*u_R(MAGX_PRIM:MAGZ_PRIM)) + &
    !     2.0*b_L*sum(u_L(VELX_PRIM:VELZ_PRIM)*u_L(MAGX_PRIM:MAGZ_PRIM)) )

end function

pure function xefluxmean_cons(U_L,U_R,F) result(mean)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    mean = xefluxmean_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure function yefluxmean_cons(U_L,U_R,F) result(mean)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    mean = yefluxmean_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure function zefluxmean_cons(U_L,U_R,F) result(mean)

    implicit none

    real, intent(IN)  :: U_L(NCONS_VARS)
    real, intent(IN)  :: U_R(NCONS_VARS)
    real, intent(IN)  :: f(NFLUXES)
    real              :: mean

    mean = zefluxmean_prim(cons2prim(u_L),cons2prim(u_R),f)

end function

pure elemental function mean(uL,uR)

    implicit none

    real, intent(in) :: uL
    real, intent(in) :: uR
    real             :: mean

    mean = 0.5*(uL + uR)

end function

pure elemental function jump(uL,uR)

    implicit none

    real, intent(in) :: uL
    real, intent(in) :: uR
    real             :: jump

    jump = uR - uL

end function

pure function xentrflux_cons(iside,f,uL,uR) result(ef)

    implicit none

    integer, intent(in) :: iside !! -1 | + 1
    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NCONS_VARS)
    real, intent(in)    :: uR(NCONS_VARS)
    real                :: ef

    real :: ss

    real :: evL(NCONS_VARS)
    real :: evR(NCONS_VARS)

    real :: phiL(NFLUXES)
    real :: phiR(NFLUXES)

    real :: potL, potR

    evL = cons2evec(uL)
    evR = cons2evec(uR)

    potL = cons2xpot(uL)
    potR = cons2xpot(uR)

    phiL = non_Conservative_xFlux_cons(uL,uL,uR)
    phiR = non_Conservative_xFlux_cons(uR,uR,uL)

    ef = 0.0

    !ss = iside
    !ss = 1.0

    ef = ef +              (dot_product(mean(evL(1:NPROP_FLUX),evR(1:NPROP_FLUX)),f(1:NPROP_FLUX)) - mean(potL,potR))
    ef = ef +          mean(dot_product(evL(1:NPROP_FLUX),phiL(1:NPROP_FLUX)),dot_product(evR(1:NPROP_FLUX),phiR(1:NPROP_FLUX)))

    ! ef = ef + 0.5*ss *     (dot_product(jump(evL(1:NPROP_FLUX),evR(1:NPROP_FLUX)),f(1:NPROP_FLUX)) - jump(potL,potR))
    ! ef = ef + 0.5*ss * jump(dot_product(evL(1:NPROP_FLUX),phiL(1:NPROP_FLUX)),dot_product(evR(1:NPROP_FLUX),phiR(1:NPROP_FLUX)))

end function

pure function yentrflux_cons(iside,f,uL,uR) result(ef)

    implicit none

    integer, intent(in) :: iside
    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NCONS_VARS)
    real, intent(in)    :: uR(NCONS_VARS)
    real                :: ef

    ef = xentrflux_cons(iside,yrotate_flux(f),yrotate_cons(uL),yrotate_cons(uR)) 

end function

pure function zentrflux_cons(iside,f,uL,uR) result(ef)

    implicit none

    integer, intent(in) :: iside
    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NCONS_VARS)
    real, intent(in)    :: uR(NCONS_VARS)
    real                :: ef

    ef = xentrflux_cons(iside,zrotate_flux(f),zrotate_cons(uL),zrotate_cons(uR)) 

end function

pure function xentrflux_prim(f,uL,uR) result(ef)

    implicit none

    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NCONS_VARS)
    real, intent(in)    :: uR(NCONS_VARS)
    real                :: ef

    real :: ss

    real :: evL(NFLUXES)
    real :: evR(NFLUXES)

    real :: phiL(NFLUXES)
    real :: phiR(NFLUXES)

    real :: potL, potR

    evL = prim2evec(uL)
    evR = prim2evec(uR)

    potL = prim2xpot(uL)
    potR = prim2xpot(uR)

    phiL = non_Conservative_xFlux_prim(uL,uL,uR)
    phiR = non_Conservative_xFlux_prim(uR,uR,uL)

    ef = 0.0

    !ss = iside
    !ss = 1.0

    ef = ef +              (dot_product(mean(evL(1:NPROP_FLUX),evR(1:NPROP_FLUX)),f(1:NPROP_FLUX)) - mean(potL,potR))
    ef = ef +          mean(dot_product(evL(1:NPROP_FLUX),phiL(1:NPROP_FLUX)),dot_product(evR(1:NPROP_FLUX),phiR(1:NPROP_FLUX)))

    ! ef = ef + 0.5*ss *     (dot_product(jump(evL(1:NPROP_FLUX),evR(1:NPROP_FLUX)),f(1:NPROP_FLUX)) - jump(potL,potR))
    ! ef = ef + 0.5*ss * jump(dot_product(evL(1:NPROP_FLUX),phiL(1:NPROP_FLUX)),dot_product(evR(1:NPROP_FLUX),phiR(1:NPROP_FLUX)))

end function

pure function yentrflux_prim(f,uL,uR) result(ef)

    implicit none

    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NPRIM_VARS)
    real, intent(in)    :: uR(NPRIM_VARS)
    real                :: ef

    ef = xentrflux_prim(yrotate_flux(f),yrotate_prim(uL),yrotate_prim(uR)) 

end function

pure function zentrflux_prim(f,uL,uR) result(ef)

    implicit none

    real, intent(in)    :: f(NFLUXES)
    real, intent(in)    :: uL(NPRIM_VARS)
    real, intent(in)    :: uR(NPRIM_VARS)
    real                :: ef

    ef = xentrflux_prim(zrotate_flux(f),zrotate_prim(uL),zrotate_prim(uR)) 

end function


pure function Nine_wave_none_prim(U_L,U_R) result(F)

    implicit none

    real, intent(in)    :: U_L(NPRIM_VARS)
    real, intent(in)    :: U_R(NPRIM_VARS)
    real                :: F(NFLUXES)

    !! Averages for wavespeeds.
    real :: bb1A, bb2A, bb3A, bbA, caA, csA, cfA, aaA
    real :: u1A, gammaA, xx

    real :: betaL, betaR

    !! Averages for flux computations.
    real :: rhoA, rhoLN, pA, pLN
    real :: B1A,B2A,B3A
    real :: u2A,v2A,w2A,betaLN,betaA,betauA,betavA,betawA
    real :: B1B1A,B2B2A,B3B3A,B1B2A,B1B3A
    real :: v1A, w1A, uB1A, vB2A, wB3A
    real :: uB1B1A, uB2B2A, uB3B3A, vB1B2A, wB1B3A
    real :: B1psiA, psiA

    real :: F_EC(NFLUXES)

    !! MHD waves selective dissipation.
    real :: bperpA, psiFplus, psiFminus, psiSplus, psiSminus
    real :: beta2A, beta3A, alphas, alphaf, sgnb1, abeta, aA, aLN
    real :: aabbA, ca, cs, cf

    associate(&
        rhoL    => U_L(DENS_PRIM), &
        rhoR    => U_R(DENS_PRIM), &
        pL      => U_L(PRES_PRIM), &
        pR      => U_R(PRES_PRIM), &
        gammaL  => U_L(GAMC_PRIM), &
        gammaR  => U_R(GAMC_PRIM), &
        psiL    => U_L(GLMP_PRIM), &
        psiR    => U_R(GLMP_PRIM), &
        uL      => U_L(VELX_PRIM), &
        uR      => U_R(VELX_PRIM), &
        vL      => U_L(VELY_PRIM), &
        vR      => U_R(VELY_PRIM), &
        wL      => U_L(VELZ_PRIM), &
        wR      => U_R(VELZ_PRIM), &
        B1L     => U_L(MAGX_PRIM), &
        B1R     => U_R(MAGX_PRIM), &
        B2L     => U_L(MAGY_PRIM), &
        B2R     => U_R(MAGY_PRIM), &
        B3L     => U_L(MAGZ_PRIM), &
        B3R     => U_R(MAGZ_PRIM))

    gammaA = 0.5*(gammaL + gammaR)

    !! Inverse temperature.
    betaL = 0.5 * rhoL / pL
    betaR = 0.5 * rhoR / pR

    ! Compute discrete wave speeds (eq. C.17, C.18)
    bb1A = 0.25 * ABS((B1L + B1R) * (B1L/rhoL + B1R/rhoR))
    bb2A = 0.25 * ABS((B2L + B2R) * (B2L/rhoL + B2R/rhoR))
    bb3A = 0.25 * ABS((B3L + B3R) * (B3L/rhoL + B3R/rhoR))

    bbA = bb1A + bb2A + bb3A

    !! Alfven speed.
    caA = SQRT(bb1A)

    !! Sound speed.
    aaA = 0.5*gammaA*(pL + pR) * (0.5/rhoL + 0.5/rhoR)

    xx = max(0.0,aaA + bbA - 2.0*sqrt(aaA*bb1A))

    !! Fast magneto-acoustic wave (always positive).
    cfA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) + sqrt(xx))

    !! Slow magneto-acoustic wave (positive or zero).
    csA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) - sqrt(xx))

    !! Compute averages for fluxes.
    rhoA    = 0.5 * (rhoL + rhoR)
    rhoLN   = LN_MEAN(rhoL, rhoR)
    betaA   = 0.5 * (betaL + betaR)
    betaLN  = LN_MEAN(betaL, betaR)
    pA      = 0.5 * rhoA / betaA
    pLN     = 0.5 * rhoLN / betaLN
    u1A     = 0.5 * (uL + uR)
    v1A     = 0.5 * (vL + vR)
    w1A     = 0.5 * (wL + wR)
    B1A     = 0.5 * (B1L + B1R)
    B2A     = 0.5 * (B2L + B2R)
    B3A     = 0.5 * (B3L + B3R)
    u2A     = 0.5 * (uL * uL + uR * uR)
    v2A     = 0.5 * (vL * vL + vR * vR)
    w2A     = 0.5 * (wL * wL + wR * wR)
    B1B1A   = 0.5 * (B1L * B1L + B1R * B1R)
    B2B2A   = 0.5 * (B2L * B2L + B2R * B2R)
    B3B3A   = 0.5 * (B3L * B3L + B3R * B3R)
    B1B2A   = 0.5 * (B1L * B2L + B1R * B2R)
    B1B3A   = 0.5 * (B1L * B3L + B1R * B3R)
    betauA  = 0.5 * (betaL * uL + betaR * uR)
    betavA  = 0.5 * (betaL * vL + betaR * vR)
    betawA  = 0.5 * (betaL * wL + betaR * wR)
    uB1B1A  = 0.5 * (uL * B1L * B1L + uR * B1R * B1R)
    uB2B2A  = 0.5 * (uL * B2L * B2L + uR * B2R * B2R)
    uB3B3A  = 0.5 * (uL * B3L * B3L + uR * B3R * B3R)
    vB1B2A  = 0.5 * (vL * B1L * B2L + vR * B1R * B2R)
    wB1B3A  = 0.5 * (wL * B1L * B3L + wR * B1R * B3R)
    uB1A    = 0.5 * (uL * B1L + uR * B1R)
    vB2A    = 0.5 * (vL * B2L + vR * B2R)
    wB3A    = 0.5 * (wL * B3L + wR * B3R)
    B1psiA  = 0.5 * (psiL * B1L + psiR * B1R)
    psiA    = 0.5 * (psiL + psiR)

    f_EC = 0.0

    !! Entropy-conserving flux.
    F_EC(DENS_FLUX) = rhoLN*u1A
    F_EC(XMOM_FLUX) = u1A*F_EC(DENS_FLUX) - B1A*B1A + pA + 0.5*(B1B1A + B2B2A + B3B3A)
    F_EC(YMOM_FLUX) = v1A*F_EC(DENS_FLUX) - B1A*B2A
    F_EC(ZMOM_FLUX) = w1A*F_EC(DENS_FLUX) - B1A*B3A
    F_EC(MAGX_FLUX) = hy_ch*PsiA
    F_EC(MAGY_FLUX) = u1A*B2A - v1A*B1A
    F_EC(MAGZ_FLUX) = u1A*B3A - w1A*B1A
    F_EC(GLMP_FLUX) = hy_ch*B1A
    F_EC(ENER_FLUX) = (0.5/(gammaA-1.0)/betaLN - 0.5*(u2A+v2A+w2A) ) * F_EC(DENS_FLUX) &
               + u1A*F_EC(XMOM_FLUX) + v1A*F_EC(YMOM_FLUX) + w1A*F_EC(ZMOM_FLUX) &
               + B1A*F_EC(MAGX_FLUX) + B2A*F_EC(MAGY_FLUX) + B3A*F_EC(MAGZ_FLUX) &
               - 0.5*(uB1B1A+uB2B2A+uB3B3A) &
               + uB1A*B1A + vB2A*B1A + wB3A*B1A &
               + psiA*F_EC(GLMP_FLUX) - hy_ch*B1psiA

    f = F_EC

! # if NSPECIES > 0
!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = 0.5*(rhoL*uL*U_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END) + rhoR*uR*U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
! # endif
! # if NMASS_SCALARS > 0
!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = 0.5*(rhoL*uL*U_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END) + rhoR*uR*U_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
! # endif
    end associate

contains

    pure elemental function LN_MEAN(a,b) result(r)

        real, intent(in)    :: a,b
        real                :: x,u,r

        real, parameter :: eps = 1e-4
        real, parameter :: c1 = 1.0/6.0
        real, parameter :: c2 = 2.0/45.0
        real, parameter :: c3 = 22.0/945.0

        x = a/b

        if (abs(a-b) < eps) then
            u  = (x*(x-2.0)+1.0)/(x*(x+2.0)+1.0)
            r = (a+b)*(0.5 - u*(c1 - u*(c2 - c3*u)))
        else
            r = (a-b)/log(x)
        endif

    end function

end function

pure function Nine_wave_LLF_prim(U_L,U_R) result(F)

    implicit none

    real, intent(in)    :: U_L(NPRIM_VARS)
    real, intent(in)    :: U_R(NPRIM_VARS)
    real                :: F(NFLUXES)

    real :: LambdaMax

    !! Averages for wavespeeds.
    real :: bb1A, bb2A, bb3A, bbA, caA, csA, cfA, aaA
    real :: u1A, gammaA, xx

    real :: betaL, betaR

    !! Averages for flux computations.
    real :: rhoA, rhoLN, pA, pLN
    real :: B1A,B2A,B3A
    real :: u2A,v2A,w2A,betaLN,betaA,betauA,betavA,betawA
    real :: B1B1A,B2B2A,B3B3A,B1B2A,B1B3A
    real :: v1A, w1A, uB1A, vB2A, wB3A
    real :: uB1B1A, uB2B2A, uB3B3A, vB1B2A, wB1B3A
    real :: B1psiA, psiA

    real :: F_EC(NFLUXES)

    !! MHD waves selective dissipation.
    real :: bperpA, psiFplus, psiFminus, psiSplus, psiSminus
    real :: beta2A, beta3A, alphas, alphaf, sgnb1, abeta, aA, aLN
    real :: aabbA, ca, cs, cf

    !! LLF stabilization.
    real :: jump(NFLUXES)
    real :: EL, ER

    associate(&
        rhoL    => U_L(DENS_PRIM), &
        rhoR    => U_R(DENS_PRIM), &
        pL      => U_L(PRES_PRIM), &
        pR      => U_R(PRES_PRIM), &
        gammaL  => U_L(GAMC_PRIM), &
        gammaR  => U_R(GAMC_PRIM), &
        psiL    => U_L(GLMP_PRIM), &
        psiR    => U_R(GLMP_PRIM), &
        uL      => U_L(VELX_PRIM), &
        uR      => U_R(VELX_PRIM), &
        vL      => U_L(VELY_PRIM), &
        vR      => U_R(VELY_PRIM), &
        wL      => U_L(VELZ_PRIM), &
        wR      => U_R(VELZ_PRIM), &
        B1L     => U_L(MAGX_PRIM), &
        B1R     => U_R(MAGX_PRIM), &
        B2L     => U_L(MAGY_PRIM), &
        B2R     => U_R(MAGY_PRIM), &
        B3L     => U_L(MAGZ_PRIM), &
        B3R     => U_R(MAGZ_PRIM))

    gammaA = 0.5*(gammaL + gammaR)

    !! Inverse temperature.
    betaL = 0.5 * rhoL / pL
    betaR = 0.5 * rhoR / pR

    ! Compute discrete wave speeds (eq. C.17, C.18)
    bb1A = 0.25 * ABS((B1L + B1R) * (B1L/rhoL + B1R/rhoR))
    bb2A = 0.25 * ABS((B2L + B2R) * (B2L/rhoL + B2R/rhoR))
    bb3A = 0.25 * ABS((B3L + B3R) * (B3L/rhoL + B3R/rhoR))

    bbA = bb1A + bb2A + bb3A

    !! Alfven speed.
    caA = SQRT(bb1A)

    !! Sound speed.
    aaA = 0.5*gammaA*(pL + pR) * (0.5/rhoL + 0.5/rhoR)

    xx = max(0.0,aaA + bbA - 2.0*sqrt(aaA*bb1A))

    !! Fast magneto-acoustic wave (always positive).
    cfA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) + sqrt(xx))

    !! Slow magneto-acoustic wave (positive or zero).
    csA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) - sqrt(xx))

    !! Compute averages for fluxes.
    rhoA    = 0.5 * (rhoL + rhoR)
    rhoLN   = LN_MEAN(rhoL, rhoR)
    betaA   = 0.5 * (betaL + betaR)
    betaLN  = LN_MEAN(betaL, betaR)
    pA      = 0.5 * rhoA / betaA
    pLN     = 0.5 * rhoLN / betaLN
    u1A     = 0.5 * (uL + uR)
    v1A     = 0.5 * (vL + vR)
    w1A     = 0.5 * (wL + wR)
    B1A     = 0.5 * (B1L + B1R)
    B2A     = 0.5 * (B2L + B2R)
    B3A     = 0.5 * (B3L + B3R)
    u2A     = 0.5 * (uL * uL + uR * uR)
    v2A     = 0.5 * (vL * vL + vR * vR)
    w2A     = 0.5 * (wL * wL + wR * wR)
    B1B1A   = 0.5 * (B1L * B1L + B1R * B1R)
    B2B2A   = 0.5 * (B2L * B2L + B2R * B2R)
    B3B3A   = 0.5 * (B3L * B3L + B3R * B3R)
    B1B2A   = 0.5 * (B1L * B2L + B1R * B2R)
    B1B3A   = 0.5 * (B1L * B3L + B1R * B3R)
    betauA  = 0.5 * (betaL * uL + betaR * uR)
    betavA  = 0.5 * (betaL * vL + betaR * vR)
    betawA  = 0.5 * (betaL * wL + betaR * wR)
    uB1B1A  = 0.5 * (uL * B1L * B1L + uR * B1R * B1R)
    uB2B2A  = 0.5 * (uL * B2L * B2L + uR * B2R * B2R)
    uB3B3A  = 0.5 * (uL * B3L * B3L + uR * B3R * B3R)
    vB1B2A  = 0.5 * (vL * B1L * B2L + vR * B1R * B2R)
    wB1B3A  = 0.5 * (wL * B1L * B3L + wR * B1R * B3R)
    uB1A    = 0.5 * (uL * B1L + uR * B1R)
    vB2A    = 0.5 * (vL * B2L + vR * B2R)
    wB3A    = 0.5 * (wL * B3L + wR * B3R)
    B1psiA  = 0.5 * (psiL * B1L + psiR * B1R)
    psiA    = 0.5 * (psiL + psiR)

    F_EC = 0.0

    !! Entropy-conserving flux.
    F_EC(DENS_FLUX) = rhoLN*u1A
    F_EC(XMOM_FLUX) = u1A*F_EC(DENS_FLUX) - B1A*B1A + pA + 0.5*(B1B1A + B2B2A + B3B3A)
    F_EC(YMOM_FLUX) = v1A*F_EC(DENS_FLUX) - B1A*B2A
    F_EC(ZMOM_FLUX) = w1A*F_EC(DENS_FLUX) - B1A*B3A
    F_EC(MAGX_FLUX) = hy_ch*PsiA
    F_EC(MAGY_FLUX) = u1A*B2A - v1A*B1A
    F_EC(MAGZ_FLUX) = u1A*B3A - w1A*B1A
    F_EC(GLMP_FLUX) = hy_ch*B1A
    F_EC(ENER_FLUX) = (0.5/(gammaA-1.0)/betaLN - 0.5*(u2A+v2A+w2A) ) * F_EC(DENS_FLUX) &
               + u1A*F_EC(XMOM_FLUX) + v1A*F_EC(YMOM_FLUX) + w1A*F_EC(ZMOM_FLUX) &
               + B1A*F_EC(MAGX_FLUX) + B2A*F_EC(MAGY_FLUX) + B3A*F_EC(MAGZ_FLUX) &
               - 0.5*(uB1B1A+uB2B2A+uB3B3A) &
               + uB1A*B1A + vB2A*B1A + wB3A*B1A &
               + psiA*F_EC(GLMP_FLUX) - hy_ch*B1psiA

    !! LLF stabilization.
    EL = 0.5*rhoL*(uL*uL+vL*vL+wL*wL) + pL/(gammaL-1.) + &
         0.5*(B1L*B1L + B2L*B2L + B3L*B3L + psiL*psiL)
    ER = 0.5*rhoR*(uR*uR+vR*vR+wR*wR) + pR/(gammaR-1.) + &
         0.5*(B1R*B1R + B2R*B2R + B3R*B3R + psiR*psiR)

    jump = 0.0

    jump(DENS_FLUX) = rhoR - rhoL
    jump(XMOM_FLUX) = rhoR*uR - rhoL*uL
    jump(YMOM_FLUX) = rhoR*vR - rhoL*vL
    jump(ZMOM_FLUX) = rhoR*wR - rhoL*wL
    jump(ENER_FLUX) = ER - EL
    jump(MAGX_FLUX) = B1R - B1L
    jump(MAGY_FLUX) = B2R - B2L
    jump(MAGZ_FLUX) = B3R - B3L
    jump(GLMP_FLUX) = psiR - psiL

    !! Maximum eigenvalue is the fast magneto-acoustic speed.
    LambdaMax = 0.5*(ABS(u1A) + cfA)

    !! Compute entropy-stable flux.
    !F = F_EC
    !return

    !! Compute entropy-stable flux.
    F = F_EC - LambdaMax*jump

!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = &
!!         F(DENS_FLUX)*merge(U_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END),U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = &
!!         F(DENS_FLUX)*merge(U_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END),U_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif

!! ------------------------------------ !!
!! Multi-species variant via Rusanov flux.
!!
!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = 0.5*(rhoL*uL*U_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END) + rhoR*uR*U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = 0.5*(rhoL*uL*U_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END) + rhoR*uR*U_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
!! # endif
!! 
!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = F(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) &
!!         + LambdaMax*(u_L(DENS_PRIM)*u_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END) - u_R(DENS_PRIM)*u_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = F(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) &
!!         + LambdaMax*(u_L(DENS_PRIM)*u_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END) - u_R(DENS_PRIM)*u_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
!! # endif
    end associate

    contains

    pure elemental function LN_MEAN(a,b) result(r)

        real, intent(in)    :: a,b
        real                :: x,u,r

        real, parameter :: eps = 1e-4
        real, parameter :: c1 = 1.0/6.0
        real, parameter :: c2 = 2.0/45.0
        real, parameter :: c3 = 22.0/945.0

        x = a/b

        if (abs(a-b) < eps) then
            u  = (x*(x-2.0)+1.0)/(x*(x+2.0)+1.0)
            r = (a+b)*(0.5 - u*(c1 - u*(c2 - c3*u)))
        else
            r = (a-b)/log(x)
        endif

    end function

end function

pure function Nine_Wave_HYBRIDWAVE_prim(U_L,U_R) result(F)

    implicit none

    real, intent(in)    :: U_L(NPRIM_VARS)
    real, intent(in)    :: U_R(NPRIM_VARS)
    real                :: F(NFLUXES)

    real :: LambdaMax

    !! Averages for wavespeeds.
    real :: bb1A, bb2A, bb3A, bbA, caA, csA, cfA, aaA
    real :: u1A, gammaA, xx

    real :: betaL, betaR

    !! Averages for flux computations.
    real :: rhoA, rhoLN, pA, pLN
    real :: B1A,B2A,B3A
    real :: u2A,v2A,w2A,betaLN,betaA,betauA,betavA,betawA
    real :: B1B1A,B2B2A,B3B3A,B1B2A,B1B3A
    real :: v1A, w1A, uB1A, vB2A, wB3A
    real :: uB1B1A, uB2B2A, uB3B3A, vB1B2A, wB1B3A
    real :: B1psiA, psiA, u2avg

    real :: F_EC(NFLUXES)

    !! MHD waves selective dissipation.
    real :: bperpA, psiFplus, psiFminus, psiSplus, psiSminus
    real :: beta2A, beta3A, alphas, alphaf, sgnb1, abeta, aA, aLN
    real :: aabbA, ca, cs, cf
    real :: dgfv_phi, xxx

    !! LLF stabilization.
    real :: jump(9), diss(9)

    !! MHD waves selective dissipation
    real, dimension(9,9) :: Dmatrix, Tmatrix, Rmatrix, RT

    real, parameter :: tiny = 1.e-32

    associate(&
        rhoL    => U_L(DENS_PRIM), &
        rhoR    => U_R(DENS_PRIM), &
        pL      => U_L(PRES_PRIM), &
        pR      => U_R(PRES_PRIM), &
        gammaL  => U_L(GAMC_PRIM), &
        gammaR  => U_R(GAMC_PRIM), &
        psiL    => U_L(GLMP_PRIM), &
        psiR    => U_R(GLMP_PRIM), &
        uL      => U_L(VELX_PRIM), &
        uR      => U_R(VELX_PRIM), &
        vL      => U_L(VELY_PRIM), &
        vR      => U_R(VELY_PRIM), &
        wL      => U_L(VELZ_PRIM), &
        wR      => U_R(VELZ_PRIM), &
        B1L     => U_L(MAGX_PRIM), &
        B1R     => U_R(MAGX_PRIM), &
        B2L     => U_L(MAGY_PRIM), &
        B2R     => U_R(MAGY_PRIM), &
        B3L     => U_L(MAGZ_PRIM), &
        B3R     => U_R(MAGZ_PRIM))

    gammaA = 0.5*(gammaL + gammaR)

    !! Inverse temperature.
    betaL = 0.5 * rhoL / pL
    betaR = 0.5 * rhoR / pR

    ! Compute discrete wave speeds (eq. C.17, C.18)
    bb1A = 0.25 * ABS((B1L + B1R) * (B1L/rhoL + B1R/rhoR))
    bb2A = 0.25 * ABS((B2L + B2R) * (B2L/rhoL + B2R/rhoR))
    bb3A = 0.25 * ABS((B3L + B3R) * (B3L/rhoL + B3R/rhoR))

    bbA = bb1A + bb2A + bb3A

    !! Alfven speed.
    caA = SQRT(bb1A)

    !! Sound speed.
    aaA = 0.5*gammaA*(pL + pR) * (0.5/rhoL + 0.5/rhoR)

    xx = max(0.0,aaA + bbA - 2.0*sqrt(aaA*bb1A))

    !! Fast magneto-acoustic wave (always positive).
    cfA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) + sqrt(xx))

    !! Slow magneto-acoustic wave (positive or zero).
    csA = 0.5*(sqrt(aaA + bbA + 2.0*sqrt(aaA*bb1A)) - sqrt(xx))

    !! Compute averages for fluxes.
    rhoA    = 0.5 * (rhoL + rhoR)
    rhoLN   = LN_MEAN(rhoL, rhoR)
    betaA   = 0.5 * (betaL + betaR)
    betaLN  = LN_MEAN(betaL, betaR)
    pA      = 0.5 * rhoA / betaA
    pLN     = 0.5 * rhoLN / betaLN
    u1A     = 0.5 * (uL + uR)
    v1A     = 0.5 * (vL + vR)
    w1A     = 0.5 * (wL + wR)
    B1A     = 0.5 * (B1L + B1R)
    B2A     = 0.5 * (B2L + B2R)
    B3A     = 0.5 * (B3L + B3R)
    u2A     = 0.5 * (uL * uL + uR * uR)
    v2A     = 0.5 * (vL * vL + vR * vR)
    w2A     = 0.5 * (wL * wL + wR * wR)
    B1B1A   = 0.5 * (B1L * B1L + B1R * B1R)
    B2B2A   = 0.5 * (B2L * B2L + B2R * B2R)
    B3B3A   = 0.5 * (B3L * B3L + B3R * B3R)
    B1B2A   = 0.5 * (B1L * B2L + B1R * B2R)
    B1B3A   = 0.5 * (B1L * B3L + B1R * B3R)
    betauA  = 0.5 * (betaL * uL + betaR * uR)
    betavA  = 0.5 * (betaL * vL + betaR * vR)
    betawA  = 0.5 * (betaL * wL + betaR * wR)
    uB1B1A  = 0.5 * (uL * B1L * B1L + uR * B1R * B1R)
    uB2B2A  = 0.5 * (uL * B2L * B2L + uR * B2R * B2R)
    uB3B3A  = 0.5 * (uL * B3L * B3L + uR * B3R * B3R)
    vB1B2A  = 0.5 * (vL * B1L * B2L + vR * B1R * B2R)
    wB1B3A  = 0.5 * (wL * B1L * B3L + wR * B1R * B3R)
    uB1A    = 0.5 * (uL * B1L + uR * B1R)
    vB2A    = 0.5 * (vL * B2L + vR * B2R)
    wB3A    = 0.5 * (wL * B3L + wR * B3R)
    B1psiA  = 0.5 * (psiL * B1L + psiR * B1R)
    psiA    = 0.5 * (psiL + psiR)

    !! Entropy-conserving flux.
    F_EC(DENS_FLUX) = rhoLN*u1A
    F_EC(XMOM_FLUX) = u1A*F_EC(DENS_FLUX) - B1A*B1A + pA + 0.5*(B1B1A + B2B2A + B3B3A)
    F_EC(YMOM_FLUX) = v1A*F_EC(DENS_FLUX) - B1A*B2A
    F_EC(ZMOM_FLUX) = w1A*F_EC(DENS_FLUX) - B1A*B3A
    F_EC(MAGX_FLUX) = hy_ch*PsiA
    F_EC(MAGY_FLUX) = u1A*B2A - v1A*B1A
    F_EC(MAGZ_FLUX) = u1A*B3A - w1A*B1A
    F_EC(GLMP_FLUX) = hy_ch*B1A
    F_EC(ENER_FLUX) = (0.5/(gammaA-1.0)/betaLN - 0.5*(u2A+v2A+w2A) ) * F_EC(DENS_FLUX) &
               + u1A*F_EC(XMOM_FLUX) + v1A*F_EC(YMOM_FLUX) + w1A*F_EC(ZMOM_FLUX) &
               + B1A*F_EC(MAGX_FLUX) + B2A*F_EC(MAGY_FLUX) + B3A*F_EC(MAGZ_FLUX) &
               - 0.5*(uB1B1A+uB2B2A+uB3B3A) &
               + uB1A*B1A + vB2A*B1A + wB3A*B1A &
               + psiA*F_EC(GLMP_FLUX) - hy_ch*B1psiA

    !! Compute additional averages.
    u2avg = uL*uR + vL*vR + wL*wR

    aA = SQRT(gammaA*pA/rhoLN)
    aLN = SQRT(gammaA*pLN/rhoLN)
    abeta = SQRT(0.5*gammaA/betaA)
    bb1A = B1A/SQRT(rhoLN)
    bb2A = B2A/SQRT(rhoLN)
    bb3A = B3A/SQRT(rhoLN)
    bbA = SQRT(bb1A*bb1A + bb2A*bb2A + bb3A*bb3A)
    aabbA = aA*aA + bbA*bbA

    !! Alfven speed.
    ca = ABS(bb1A)

    !! Control round-off errors
    xx = max(0.0,aabbA*aabbA - 4.0*aA*aA*bb1A*bb1A)
    xxx = aabbA + SQRT(xx)

    !! Fast magnetoacoustic speed.
    cf = SQRT(0.5 * xxx)

    ! Control round-off errors.
    xxx = max(0.0,aabbA - SQRT(xx))

    !! Slow magnetoacoustic speed.
    cs = SQRT(0.5 * xxx)

    bperpA = SQRT(bb2A*bb2A + bb3A*bb3A)

    !! In case of very small bperpA, the values of betaA_{2,3}
    !! are indeterminable so we make them pairwise orthogonal.
    if (bperpA .gt. tiny) then
      beta2A = bb2A/bperpA
      beta3A = bb3A/bperpA
    else
      bperpA = 0.0
      beta2A = 1.0/sqrt(2.0)
      beta3A = 1.0/sqrt(2.0)
    endif

    !! Avoid negative round-off errors when computing alphaf (auxiliary variable).
    xx = 0.0
    if ((cf*cf-cs*cs) .gt. 0.0) then
        xx = (aA*aA-cs*cs)/(cf*cf-cs*cs)
    end if

    if (xx .gt. 0.0) then
        alphaf = SQRT(xx)
    else
        alphaf = 0.0
    end if

    !! Avoid negative round-off errors when computing alphas (auxiliary variable).
    xx = 0.0

    if ((cf*cf-cs*cs) .gt. 0.0) then
        xx = (cf*cf-aA*aA)/(cf*cf-cs*cs)
    end if

    if (xx .gt. 0.0) then
        alphas = SQRT(xx)
    else
        alphas = 0.0
    end if

    sgnb1 = sign(1.0,bb1A)

    !! Derigs et al. (2018), (4.63).
    psiSplus =  0.5*alphas*rhoLN*u2avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN/(gammaA-1.0) + alphas*cs*rhoLN*u1A + &
                alphaf*cf*rhoLN*sgnb1*(v1A*beta2A + w1A*beta3A)
    psiSminus = 0.5*alphas*rhoLN*u2avg - abeta*alphaf*rhoLN*bperpA + &
                alphas*rhoLN*aLN*aLN/(gammaA-1.0) - alphas*cs*rhoLN*u1A - &
                alphaf*cf*rhoLN*sgnb1*(v1A*beta2A + w1A*beta3A)

    psiFplus =  0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN/(gammaA-1.0) + alphaf*cf*rhoLN*u1A - &
                alphas*cs*rhoLN*sgnb1*(v1A*beta2A + w1A*beta3A)
    psiFminus = 0.5*alphaf*rhoLN*u2avg + abeta*alphas*rhoLN*bperpA + &
                alphaf*rhoLN*aLN*aLN/(gammaA-1.0) - alphaf*cf*rhoLN*u1A + &
                alphas*cs*rhoLN*sgnb1*(v1A*beta2A + w1A*beta3A)

    ! + fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,1) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(u1A + cf), &
                      rhoLN*(alphaf*v1A - alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*w1A - alphas*cs*beta3A*sgnb1), &
                      psiFplus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,2) = (/ 0.0, &
                      0.0,  &
                      rhoLN*SQRT(rhoA)*beta3A, &
                      -rhoLN*SQRT(rhoA)*beta2A, &
                      -rhoLN*SQRT(rhoA)*(beta2A*w1A - beta3A*v1A), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! + slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,3) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(u1A + cs), &
                      rhoLN*(alphas*v1A + alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*w1A + alphaf*cf*beta3A*sgnb1), &
                      psiSplus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! + GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,4) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      B1A + psiA, &
                      1.0, &
                      0.0, &
                      0.0, &
                      1.0 /)

    ! Entropy wave
    ! Derigs et al. (2018), (4.66)
    Rmatrix(:,5) = (/ 1.0, &
                      u1A, &
                      v1A, &
                      w1A, &
                      0.5*u2avg, &
                      0.0, &
                      0.0, &
                      0.0, &
                      0.0 /)

    ! - GLM wave
    ! Dergs et al. (2018), eq. (4.65)
    Rmatrix(:,6) = (/ 0.0, &
                      0.0, &
                      0.0, &
                      0.0, &
                      B1A - psiA, &
                      1.0, &
                      0.0, &
                      0.0, &
                      -1.0 /)

    ! - slow magnetoacoustic wave
    ! Derigs et al. (2018), (4.69)
    Rmatrix(:,7) = (/ alphas*rhoLN, &
                      alphas*rhoLN*(u1A - cs), &
                      rhoLN*(alphas*v1A - alphaf*cf*beta2A*sgnb1), &
                      rhoLN*(alphas*w1A - alphaf*cf*beta3A*sgnb1), &
                      psiSminus, &
                      0.0, &
                      -alphaf*abeta*beta2A*SQRT(rhoLN), &
                      -alphaf*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    ! - Alfven wave
    ! Derigs et al. (2018), (4.67)
    Rmatrix(:,8) = (/ 0.0, &
                      0.0, &
                      -rhoLN*SQRT(rhoA)*beta3A, &
                      rhoLN*SQRT(rhoA)*beta2A, &
                      rhoLN*SQRT(rhoA)*(beta2A*w1A - beta3A*v1A), &
                      0.0, &
                      -rhoLN*beta3A, &
                      rhoLN*beta2A, &
                      0.0 /)

    ! - fast magnetoacoustic wave
    ! Derigs et al. (2018), (4.68)
    Rmatrix(:,9) = (/ alphaf*rhoLN, &
                      alphaf*rhoLN*(u1A - cf), &
                      rhoLN*(alphaf*v1A + alphas*cs*beta2A*sgnb1), &
                      rhoLN*(alphaf*w1A + alphas*cs*beta3A*sgnb1), &
                      psiFminus, &
                      0.0, &
                      alphas*abeta*beta2A*SQRT(rhoLN), &
                      alphas*abeta*beta3A*SQRT(rhoLN), &
                      0.0 /)

    !! Blending factor.
    dgfv_phi = sqrt(ABS(1.0-(pR/pL))/(1.0+(pR/pL)))

    !! Obtain wavespeed from maximum of non-hybrid D matrix.
    LambdaMax = abs(u1A) + abs(cf)

    Dmatrix = 0.0
    Dmatrix(1,1) = (1.0-dgfv_phi)*ABS( u1A + cf     ) + dgfv_phi*LambdaMax
    Dmatrix(2,2) = (1.0-dgfv_phi)*ABS( u1A + ca     ) + dgfv_phi*LambdaMax
    Dmatrix(3,3) = (1.0-dgfv_phi)*ABS( u1A + cs     ) + dgfv_phi*LambdaMax
    Dmatrix(4,4) = (1.0-dgfv_phi)*ABS( u1A + hy_ch  ) + dgfv_phi*LambdaMax
    Dmatrix(5,5) = (1.0-dgfv_phi)*ABS( u1A          ) + dgfv_phi*LambdaMax
    Dmatrix(6,6) = (1.0-dgfv_phi)*ABS( u1A - hy_ch  ) + dgfv_phi*LambdaMax
    Dmatrix(7,7) = (1.0-dgfv_phi)*ABS( u1A - cs     ) + dgfv_phi*LambdaMax
    Dmatrix(8,8) = (1.0-dgfv_phi)*ABS( u1A - ca     ) + dgfv_phi*LambdaMax
    Dmatrix(9,9) = (1.0-dgfv_phi)*ABS( u1A - cf     ) + dgfv_phi*LambdaMax

    !! Diagonal scaling matrix as described in Winters et al., eq. (4.15).
    Tmatrix = 0.0
    Tmatrix(1,1) = 0.5/gammaA/rhoLN ! + f
    Tmatrix(2,2) = 0.25/betaA/rhoLN/rhoLN ! + a
    Tmatrix(3,3) = Tmatrix(1,1) ! + s
    Tmatrix(4,4) = 0.25 / betaA ! + GLM
    Tmatrix(5,5) = rhoLN*(gammaA-1.0)/gammaA ! E
    Tmatrix(6,6) = Tmatrix(4,4) ! - GLM
    Tmatrix(7,7) = Tmatrix(1,1) ! - s
    Tmatrix(8,8) = Tmatrix(2,2) ! - a
    Tmatrix(9,9) = Tmatrix(1,1) ! - f

    !! Scale D matrix.
    Dmatrix = MATMUL(Dmatrix,Tmatrix)

    !! Compute jump in entropy variables.
    jump = entropy_vector(rhoR, uR, vR, wR, pR, B1R, B2R, B3R, psiR, gammaR) - &
           entropy_vector(rhoL, uL, vL, wL, pL, B1L, B2L, B3L, psiL, gammaL)

    RT = TRANSPOSE(Rmatrix)

    diss = 0.5*MATMUL(MATMUL(Rmatrix,MATMUL(Dmatrix,RT)),jump)

    f = 0.0

    f(DENS_FLUX) = f_EC(DENS_FLUX) + diss(1)
    f(XMOM_FLUX) = f_EC(XMOM_FLUX) + diss(2)
    f(YMOM_FLUX) = f_EC(YMOM_FLUX) + diss(3)
    f(ZMOM_FLUX) = f_EC(ZMOM_FLUX) + diss(4)
    f(ENER_FLUX) = f_EC(ENER_FLUX) + diss(5)
    f(MAGX_FLUX) = f_EC(MAGX_FLUX) + diss(6)
    f(MAGY_FLUX) = f_EC(MAGY_FLUX) + diss(7)
    f(MAGZ_FLUX) = f_EC(MAGZ_FLUX) + diss(8)
    f(GLMP_FLUX) = f_EC(GLMP_FLUX) + diss(9)

!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = 0.5*(rhoL*uL*U_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END) + rhoR*uR*U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = 0.5*(rhoL*uL*U_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END) + rhoR*uR*U_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
!! # endif
!! 
!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = F(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) &
!!         + 0.5*LambdaMax*(u_L(DENS_PRIM)*u_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END) - u_R(DENS_PRIM)*u_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = F(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) &
!!         + 0.5*LambdaMax*(u_L(DENS_PRIM)*u_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END) - u_R(DENS_PRIM)*u_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
!! # endif
    end associate

    contains

    pure elemental function LN_MEAN(a,b) result(r)

        real, intent(in)    :: a,b
        real                :: x,u,r

        real, parameter :: eps = 1e-4
        real, parameter :: c1 = 1.0/6.0
        real, parameter :: c2 = 2.0/45.0
        real, parameter :: c3 = 22.0/945.0

        x = a/b

        if (abs(a-b) < eps) then
            u  = (x*(x-2.0)+1.0)/(x*(x+2.0)+1.0)
            r = (a+b)*(0.5 - u*(c1 - u*(c2 - c3*u)))
        else
            r = (a-b)/log(x)
        endif

    end function

    pure function entropy_vector(rho, u, v, w, p, B1, B2, B3, psi, gamma)

        !! Computes the entropy vector for given (primitive) quantities

        implicit none

        real                :: entropy_vector(9)
        real, intent(IN)    :: rho, u, v, w, p, B1, B2, B3, psi, gamma
        real :: s

        s = LOG(p) - gamma*LOG(rho)

        entropy_vector(1) = -((gamma-s)/(gamma - 1.0) - 0.5*rho*(u*u+v*v+w*w)/p)
        entropy_vector(2) = -rho*u/p
        entropy_vector(3) = -rho*v/p
        entropy_vector(4) = -rho*w/p
        entropy_vector(5) =  rho/p
        entropy_vector(6) = -rho*B1/p
        entropy_vector(7) = -rho*B2/p
        entropy_vector(8) = -rho*B3/p
        entropy_vector(9) = -rho*psi/p

    end function

end function

pure function xriem_fv_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)

    f = Nine_Wave_LLF_prim(uL,uR)
    !f = Nine_Wave_HYBRIDWAVE_prim(uL,uR)
    !f = Nine_Wave_none_prim(uL,uR)
    !f = hll_prim(uL,uR)
    !f = hllc_prim(uL,uR)
    !f = rusanov_prim(uL,uR)

end function

pure function xriem_dg_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)

    f = Nine_Wave_LLF_prim(uL,uR)
    !f = Nine_Wave_none_prim(uL,uR)
    !f = hll_prim(uL,uR)
    !f = hllc_prim(uL,uR)
    !f = rusanov_prim(uL,uR)

end function

pure function yriem_fv_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)
     
    f = yrotate_back_flux(xriem_fv_prim(yrotate_prim(uL),yrotate_prim(uR)))

end function

pure function zriem_fv_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)
     
    f = zrotate_back_flux(xriem_fv_prim(zrotate_prim(uL),zrotate_prim(uR)))

end function

pure function yriem_dg_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)
     
    f = yrotate_back_flux(xriem_dg_prim(yrotate_prim(uL),yrotate_prim(uR)))

end function

pure function zriem_dg_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)
     
    f = zrotate_back_flux(xriem_dg_prim(zrotate_prim(uL),zrotate_prim(uR)))

end function

!! Rusanov flux for primitive variables.
pure function rusanov_prim(uL,uR) result(f)

    real, intent(in) :: uL(NPRIM_VARS),uR(NPRIM_VARS)
    real             :: f(NFLUXES)
    
    real :: enerL,enerR,lmax

    f = 0.0

    !f = two_point_xflux(uL,uR)
    f = 0.5*(xflux_prim(uL) + xflux_prim(uR))

    lMax = 0.5*max(abs(uL(VELX_PRIM)) + Fastest_Signal_Speed_prim(uL), &
               abs(uR(VELX_PRIM)) + Fastest_Signal_Speed_prim(uR))

    enerL = uL(DENS_PRIM)*uL(ENER_PRIM) + 0.5*(sum(uL(MAGX_PRIM:MAGZ_PRIM)**2) + uL(GLMP_PRIM)**2)
    enerR = uR(DENS_PRIM)*uR(ENER_PRIM) + 0.5*(sum(uR(MAGX_PRIM:MAGZ_PRIM)**2) + uR(GLMP_PRIM)**2)

    f(DENS_FLUX) = F(DENS_FLUX) + lMax*(uL(DENS_PRIM)               - uR(DENS_PRIM))
    f(XMOM_FLUX) = F(XMOM_FLUX) + lMax*(uL(DENS_PRIM)*uL(VELX_PRIM) - uR(DENS_PRIM)*uR(VELX_PRIM))
    f(YMOM_FLUX) = F(YMOM_FLUX) + lMax*(uL(DENS_PRIM)*uL(VELY_PRIM) - uR(DENS_PRIM)*uR(VELY_PRIM))
    f(ZMOM_FLUX) = F(ZMOM_FLUX) + lMax*(uL(DENS_PRIM)*uL(VELZ_PRIM) - uR(DENS_PRIM)*uR(VELZ_PRIM))
    f(ENER_FLUX) = F(ENER_FLUX) + lMax*(enerL - enerR)
    f(MAGX_FLUX) = F(MAGX_FLUX) + lMax*(uL(MAGX_PRIM) - uR(MAGX_PRIM))
    f(MAGY_FLUX) = F(MAGY_FLUX) + lMax*(uL(MAGY_PRIM) - uR(MAGY_PRIM))
    f(MAGZ_FLUX) = F(MAGZ_FLUX) + lMax*(uL(MAGZ_PRIM) - uR(MAGZ_PRIM))
    f(GLMP_FLUX) = F(GLMP_FLUX) + lMax*(uL(GLMP_PRIM) - uR(GLMP_PRIM))

!! Not needed anymore.
!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = &
!!         F(DENS_FLUX)*merge(UL(SPEC_PRIM_BEGIN:SPEC_PRIM_END),UR(SPEC_PRIM_BEGIN:SPEC_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = &
!!         F(DENS_FLUX)*merge(UL(MSCL_PRIM_BEGIN:MSCL_PRIM_END),UR(MSCL_PRIM_BEGIN:MSCL_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif

!! Rusano variant.
!! # if NSPECIES > 0
!!     f(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = F(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) &
!!         + lmax*(uL(DENS_PRIM)*uL(SPEC_PRIM_BEGIN:SPEC_PRIM_END) - uR(DENS_PRIM)*uR(SPEC_PRIM_BEGIN:SPEC_PRIM_END))
!! # endif
!! # if NMASS_SCALARS > 0
!!     f(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = F(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) &
!!         + lmax*(uL(DENS_PRIM)*uL(MSCL_PRIM_BEGIN:MSCL_PRIM_END) - uR(DENS_PRIM)*uR(MSCL_PRIM_BEGIN:MSCL_PRIM_END))
!! # endif
end function

!! !! HLLE Riemann solver
!! pure function HLL_prim(U_L,U_R) result(F)
!! 
!!     implicit none
!! 
!!     real, intent(in) :: u_L(NPRIM_VARS),u_R(NPRIM_VARS)
!!     real             :: f(NFLUXES)
!!  
!!     real :: F_L(NFLUXES), F_R(NFLUXES)
!!     real :: cu_L(NFLUXES), cu_R(NFLUXES)
!! 
!!     real :: c_L,c_R !! sound speed
!!     real :: H_L,H_R !! enthalpy
!!     real :: p_L,p_R !! pressure
!!     real :: S_L,S_R !! signal wave speed
!!     real :: V_L(3),V_R(3) !! velocities
!! 
!!     real :: sqrtrho_L,sqrtrho_R,ssqrtrho,absvel
!!     real :: roeVel(3),roeH,roeC
!!     real :: beta, kappa
!! 
!!     kappa = 0.5*(u_L(GAMC_PRIM) + u_R(GAMC_PRIM))
!! 
!!     beta = SQRT(0.5*(kappa-1.0)/kappa)
!! 
!!     !! velocities
!!     V_L = U_L(VELX_PRIM:VELZ_PRIM)
!!     V_R = U_R(VELX_PRIM:VELZ_PRIM)
!! 
!!     !! pressures
!!     p_L = U_L(PRES_PRIM)
!!     p_R = U_R(PRES_PRIM)
!! 
!!     !! sound speeds
!!     c_L = SQRT(kappa*p_L/U_L(DENS_PRIM))
!!     c_R = SQRT(kappa*p_R/U_R(DENS_PRIM))
!! 
!!     !! enthalpy
!!     H_L = U_L(ENER_PRIM) + p_L/U_L(DENS_PRIM)
!!     H_R = U_R(ENER_PRIM) + p_R/U_R(DENS_PRIM)
!! 
!!     !! auxiliaries
!!     SqrtRho_L = SQRT(U_L(DENS_PRIM))
!!     SqrtRho_R = SQRT(U_R(DENS_PRIM))
!!     sSqrtRho  = 1.0/(SqrtRho_L+SqrtRho_R)
!! 
!!     !! Roe averages
!!     RoeVel    = (SqrtRho_R*V_R + SqrtRho_L*V_L) * sSqrtRho
!!     RoeH      = (SqrtRho_R*H_R + SqrtRho_L*H_L) * sSqrtRho
!!     absVel    = DOT_PRODUCT(RoeVel,RoeVel)
!!     RoeC      = SQRT((kappa-1.0)*(RoeH-0.5*absVel))
!! 
!!     !! signal velocities
!!     S_L = MIN(RoeVel(1) - RoeC, V_L(1) - beta*c_L, 0.0)
!!     S_R = MAX(RoeVel(1) + RoeC, V_R(1) + beta*c_R, 0.0)
!! 
!!     if (S_L .GE. 0.0 .and. S_R .GT. 0.0) then
!! 
!!         F = xflux_prim(U_L)
!! 
!!     elseif (S_R .LE. 0.0 .and. S_L .LT. 0.0) then
!! 
!!         F = xflux_prim(U_R)
!! 
!!     else
!! 
!!         F_L = xflux_prim(U_L)
!!         F_R = xflux_prim(U_R)
!! 
!!         cu_L = 0.0
!!         cu_L(DENS_FLUX) = u_L(DENS_PRIM)
!!         cu_L(XMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELX_PRIM)
!!         cu_L(YMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELY_PRIM)
!!         cu_L(ZMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELZ_PRIM)
!!         cu_L(ENER_FLUX) = u_L(DENS_PRIM)*u_L(ENER_PRIM)
!! 
!! # if NSPECIES > 0
!!         cu_L(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = u_L(DENS_PRIM)*u_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!         cu_L(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = u_L(DENS_PRIM)*u_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!! 
!!         cu_R = 0.0
!!         cu_R(DENS_FLUX) = u_R(DENS_PRIM)
!!         cu_R(XMOM_FLUX) = u_R(DENS_PRIM)*u_R(VELX_PRIM)
!!         cu_R(YMOM_FLUX) = u_R(DENS_PRIM)*u_R(VELY_PRIM)
!!         cu_R(ZMOM_FLUX) = u_R(DENS_PRIM)*u_R(VELZ_PRIM)
!!         cu_R(ENER_FLUX) = u_R(DENS_PRIM)*u_R(ENER_PRIM)
!! 
!! # if NSPECIES > 0
!!         cu_R(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = u_R(DENS_PRIM)*u_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!         cu_R(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = u_R(DENS_PRIM)*u_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!! 
!!         F = (S_R*F_L - S_L*F_R + S_L*S_R*(cU_R - cU_L)) / (S_R - S_L)
!! 
!!     end if 
!! 
!! end function

!! HLLC Riemann solver
!! pure function HLLC_prim(U_L,U_R) result(F)
!! 
!!     implicit none
!! 
!!     real, intent(in) :: u_L(NPRIM_VARS),u_R(NPRIM_VARS)
!!     real             :: f(NFLUXES)
!!  
!!     real :: F_L(NFLUXES), F_R(NFLUXES)
!! 
!!     real :: c_L,c_R !! sound speed
!!     real :: H_L,H_R !! enthalpy
!!     real :: p_L,p_R !! pressure
!!     real :: S_L,S_R !! signal wave speed
!!     real :: V_L,V_R !! velocities
!! 
!!     real :: HTilde, vTilde, aTilde
!!     real :: SStar, QStar, UStar(NFLUXES), cons(NFLUXES)
!! 
!!     ! real :: smooth
!!     ! real :: F_HLLE(N_VARS)
!! 
!!     real :: beta, kappa
!! 
!!     kappa = 0.5*(u_L(GAMC_PRIM) + u_R(GAMC_PRIM))
!!     beta = SQRT(0.5*(kappa-1)/kappa)
!! 
!!     !! velocities
!!     v_L = U_L(VELX_PRIM)
!!     v_R = U_R(VELX_PRIM)
!! 
!!     !! pressures
!!     p_L = U_L(PRES_PRIM)
!!     p_R = U_R(PRES_PRIM)
!! 
!!     !! sound speed
!!     c_L = SQRT(kappa*p_L/U_L(DENS_PRIM))
!!     c_R = SQRT(kappa*p_R/U_R(DENS_PRIM))
!! 
!!     !! Signal velocity estimates by Davis (1988) and Einfeldt (1988)
!! 
!!     !! enthalpy
!!     H_L = U_L(ENER_PRIM) + p_L/U_L(DENS_PRIM)
!!     H_R = U_R(ENER_PRIM) + p_R/U_R(DENS_PRIM)
!! 
!!     !! Roe averages
!!     hTilde = (SQRT(U_L(DENS_PRIM))*H_L + SQRT(U_R(DENS_PRIM))*H_R)/(SQRT(U_L(DENS_PRIM)) + SQRT(U_R(DENS_PRIM)))
!!     vTilde = (SQRT(U_L(DENS_PRIM))*v_L + SQRT(U_R(DENS_PRIM))*v_R)/(SQRT(U_L(DENS_PRIM)) + SQRT(U_R(DENS_PRIM)))
!!     aTilde = SQRT((kappa-1.0)*(HTilde - 0.5*vTilde**2))
!! 
!!     !! signal velocities
!!     S_L = MIN(vTilde - aTilde, V_L - beta*c_L, 0.0)
!!     S_R = MAX(vTilde + aTilde, V_R + beta*c_R, 0.0)
!! 
!!     if (S_L .GE. 0.0) then
!!         F = xflux_prim(U_L)
!! 
!!     elseif (S_R .LE. 0.0) then
!!         F = xflux_prim(U_R)
!! 
!!     else
!!         F_L = xflux_prim(U_L)
!!         F_R = xflux_prim(U_R)
!! 
!!         ! smooth = SQRT(abs(p_L - p_R)/(p_L + p_R))
!!         ! F_HLLE = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)
!! 
!!         SStar = (p_R - p_L + U_L(DENS_PRIM)*U_L(VELX_PRIM)*(S_L - V_L) - U_R(DENS_PRIM)*U_R(VELX_PRIM)*(S_R - V_R)) / (U_L(DENS_PRIM)*(S_L - V_L) - U_R(DENS_PRIM)*(S_R - V_R))
!! 
!!         UStar = 0.0
!! 
!!         if (S_L .le. 0.0 .and. 0.0 .le. SStar) then
!!             QStar = U_L(DENS_PRIM)*(S_L - V_L)/(S_L-SStar)
!! 
!!             UStar(DENS_FLUX) = QStar
!!             UStar(XMOM_FLUX) = QStar * SStar
!!             UStar(YMOM_FLUX) = QStar * U_L(VELY_PRIM)
!!             UStar(ZMOM_FLUX) = QStar * U_L(VELZ_PRIM)
!!             UStar(ENER_FLUX) = QStar * (U_L(ENER_PRIM) + (SStar-V_L)*(SStar + p_L/U_L(DENS_PRIM)/(S_L - V_L)))
!! 
!! # if NSPECIES > 0
!!             UStar(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = QStar*u_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!             UStar(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = QStar*u_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!!             !! Convert to conservative state variables.
!!             cons = 0.0
!!             cons(DENS_FLUX) = u_L(DENS_PRIM)
!!             cons(XMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELX_PRIM)
!!             cons(YMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELY_PRIM)
!!             cons(ZMOM_FLUX) = u_L(DENS_PRIM)*u_L(VELZ_PRIM)
!!             cons(ENER_FLUX) = u_L(DENS_PRIM)*u_L(ENER_PRIM)
!! 
!! # if NSPECIES > 0
!!             cons(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = u_L(DENS_PRIM)*u_L(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!             cons(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = u_L(DENS_PRIM)*u_L(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!! 
!!             !F = smooth*F_HLLE + (1.0_dp-smooth)*(F_L + S_L*(UStar - U_L))
!!             F = F_L + S_L*(UStar - cons)
!! 
!!         else
!!             QStar = U_R(DENS_PRIM)*(S_R - V_R)/(S_R-SStar)
!! 
!!             UStar(DENS_FLUX) = QStar
!!             UStar(XMOM_FLUX) = QStar * SStar
!!             UStar(YMOM_FLUX) = QStar * U_R(VELY_PRIM)
!!             UStar(ZMOM_FLUX) = QStar * U_R(VELZ_PRIM)
!!             UStar(ENER_FLUX) = QStar * (U_R(ENER_PRIM) + (SStar-V_R)*(SStar + p_R/U_R(DENS_PRIM)/(S_R - V_R)))
!! 
!! # if NSPECIES > 0
!!             UStar(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = QStar*U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!             UStar(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = QStar*u_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!!             !! Convert to conservative state variables.
!!             cons = 0.0
!!             cons(DENS_FLUX) = U_R(DENS_PRIM)
!!             cons(XMOM_FLUX) = U_R(DENS_PRIM)*U_R(VELX_PRIM)
!!             cons(YMOM_FLUX) = U_R(DENS_PRIM)*U_R(VELY_PRIM)
!!             cons(ZMOM_FLUX) = U_R(DENS_PRIM)*U_R(VELZ_PRIM)
!!             cons(ENER_FLUX) = U_R(DENS_PRIM)*U_R(ENER_PRIM)
!! 
!! # if NSPECIES > 0
!!             cons(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = U_R(DENS_PRIM)*U_R(SPEC_PRIM_BEGIN:SPEC_PRIM_END)
!! # endif
!! # if NMASS_SCALARS > 0
!!             cons(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = U_R(DENS_PRIM)*U_R(MSCL_PRIM_BEGIN:MSCL_PRIM_END)
!! # endif
!! 
!!             !F = smooth*F_HLLE + (1.0_dp-smooth)*(F_R + S_R*(UStar - U_R))
!!             F = F_R + S_R*(UStar - cons)
!! 
!!         end if 
!!     end if 
!! 
!! end function

pure function TWO_POINT_XFLUX_CONS(uL,uR) result(f)

    implicit none

    real, intent(in)    :: uL(NCONS_VARS)
    real, intent(in)    :: uR(NCONS_VARS)
    real                :: f(NFLUXES)

    ! f = 0.5*(xflux_prim(cons2prim(uL)) + xflux_prim(cons2prim(uR)))

    f = Nine_wave_none_prim(cons2prim(uL),cons2prim(uR))

end function

!! pure function TWO_POINT_XFLUX_CONS(uL,uR) result(fstar)
!! 
!!     implicit none
!! 
!!     real, intent(in)    :: uL(NCONS_VARS)
!!     real, intent(in)    :: uR(NCONS_VARS)
!!     real                :: fstar(NFLUXES)
!! 
!!     real, parameter :: smu_0 = 1.0
!! 
!!     real :: betaLN,beta_R,beta_L
!!     real :: rhoLN,B2_L,B2_R,v2_L,v2_R
!!     real :: pTilde,p_L,p_R
!!     real :: v_L(3),v_R(3)
!!     real :: BAvg(3),vAvg(3)
!!     real :: v1_B2Avg
!!     real :: vB_Avg
!!     real :: psiAvg
!! 
!!     real :: kappa
!!     real :: kappaM1
!!     real :: skappaM1
!! 
!!     kappa    = 0.5*(UL(GAMC_CONS) + UR(GAMC_CONS))
!!     kappaM1  = kappa - 1.0
!!     skappaM1 = 1.0/kappaM1
!! 
!!     associate(&
!!         rho_L   => UL(DENS_CONS), &
!!         rho_R   => UR(DENS_CONS), &
!!         rhoV_L  => UL(MOMX_CONS:MOMZ_CONS), &
!!         rhoV_R  => UR(MOMX_CONS:MOMZ_CONS), &
!!         E_L     => UL(ENER_CONS) - 0.5*smu_0*UL(GLMP_CONS)**2, &
!!         E_R     => UR(ENER_CONS) - 0.5*smu_0*UR(GLMP_CONS)**2, &
!!         psi_L   => UL(GLMP_CONS), &
!!         psi_R   => UR(GLMP_CONS), &
!!         B_L     => UL(MAGX_CONS:MAGZ_CONS), &
!!         B_R     => UR(MAGX_CONS:MAGZ_CONS))
!! 
!!     !! Get the inverse density, velocity, and pressure on left and right.
!!     v_L = rhoV_L(:)/rho_L
!!     v_R = rhoV_R(:)/rho_R
!! 
!!     v2_L = SUM(v_L(:)*v_L(:))
!!     v2_R = SUM(v_R(:)*v_R(:))
!!     B2_L = SUM(B_L(:)*B_L(:))
!!     B2_R = SUM(B_R(:)*B_R(:))
!! 
!!     p_L    = kappaM1*(E_L - 0.5*(rho_L*v2_L+smu_0*B2_L))
!!     p_R    = kappaM1*(E_R - 0.5*(rho_R*v2_R+smu_0*B2_R))
!!     beta_L = 0.5*rho_L/p_L
!!     beta_R = 0.5*rho_R/P_R
!! 
!!     !! Get the averages for the numerical flux.
!!     rhoLN      = logmean(rho_L, rho_R)
!!     betaLN     = logmean(beta_L,beta_R)
!!     vAvg       = 0.5 * (v_L + v_R)
!!     BAvg       = 0.5 * (B_L + B_R)
!!     v1_B2Avg   = 0.5 * (v_L(1)*B2_L + v_R(1)*B2_R)
!!     vB_Avg     = 0.5 * (SUM(V_L(:)*B_L(:)) + SUM(V_R(:)*B_R(:)))
!!     pTilde     = 0.5 * ((rho_L+rho_R)/(beta_L+beta_R)+smu_0*0.5*(B2_L+B2_R))
!!     psiAvg     = 0.5 * (psi_L+psi_R)
!! 
!!     !! Entropy conserving and kinetic energy conserving flux.
!!     Fstar(DENS_FLUX) = rhoLN*vAvg(1)
!!     Fstar(XMOM_FLUX) = Fstar(DENS_FLUX)*vAvg(1) - smu_0*BAvg(1)*BAvg(1) + pTilde
!!     Fstar(YMOM_FLUX) = Fstar(DENS_FLUX)*vAvg(2) - smu_0*BAvg(1)*BAvg(2)
!!     Fstar(ZMOM_FLUX) = Fstar(DENS_FLUX)*vAvg(3) - smu_0*BAvg(1)*BAvg(3)
!!     Fstar(MAGX_FLUX) = vAvg(1)*Bavg(2) - BAvg(1)*vAvg(2)
!!     Fstar(MAGY_FLUX) = vAvg(1)*Bavg(3) - BAvg(1)*vAvg(3)
!!     Fstar(MAGZ_FLUX) = hy_ch*psiAvg
!!     Fstar(GLMP_FLUX) = hy_ch*BAvg(1)
!! 
!!     !! Total energy flux.
!!     Fstar(ENER_FLUX) = Fstar(DENS_FLUX)*0.5*(skappaM1/betaLN - 0.5*(v2_L+v2_R))  &
!!                + SUM(vAvg(:)*Fstar(XMOM_FLUX:ZMOM_FLUX)) &
!!                + smu_0*(SUM(BAvg(:)*Fstar(MAGX_FLUX:MAGZ_FLUX)) &
!!                        - 0.5*v1_B2Avg + BAvg(1)*vB_Avg &
!!                        + Fstar(GLMP_FLUX)*psiAvg-hy_ch*0.5*(psi_L*B_L(1)+psi_R*B_R(1)))
!!     end associate 
!! 
!!     contains
!! 
!!     pure elemental function logmean(a,b) result(r)
!! 
!!         implicit none
!! 
!!         real, intent(in)    :: a,b
!!         real                :: x,u,r
!! 
!!         real, parameter  :: eps = 1e-4
!!         real, parameter :: c1 = 1.0/6.0
!!         real, parameter :: c2 = 2.0/45.0
!!         real, parameter :: c3 = 22.0/945.0
!! 
!!         x = a/b
!! 
!!         if (abs(a-b) < eps) then
!!             u  = (x*(x-2.0)+1.0)/(x*(x+2.0)+1.0)
!!             r = (a+b)*(0.5 - u*(c1 - u*(c2 - c3*u)))
!!         else
!!             r = (a-b)/log(x)
!!         endif
!! 
!!     end function
!! 
!! end function

pure function TWO_POINT_YFLUX_CONS(uL,uR) result(f)

    implicit none

    real, intent(in) :: uL(NCONS_VARS), uR(NCONS_VARS)
    real             :: f(NFLUXES)

    f = yrotate_back_flux(TWO_POINT_XFLUX_CONS(yrotate_cons(uL),yrotate_cons(uR)))

end function

pure function TWO_POINT_ZFLUX_CONS(uL,uR) result(f)

    implicit none

    real, intent(in) :: uL(NCONS_VARS), uR(NCONS_VARS)
    real             :: f(NFLUXES)

    f = zrotate_back_flux(TWO_POINT_XFLUX_CONS(zrotate_cons(uL),zrotate_cons(uR)))

end function

pure function yrotate_prim(u) result(ru)

    !! Rotate state (primitive) vector to x-direction.

    implicit none

    real, intent(in)    :: u(NPRIM_VARS)
    real                :: ru(NPRIM_VARS)

    ru = u

    ru(VELX_PRIM) = u(VELY_PRIM)
    ru(VELY_PRIM) = u(VELZ_PRIM)
    ru(VELZ_PRIM) = u(VELX_PRIM)

    ru(MAGX_PRIM) = u(MAGY_PRIM)
    ru(MAGY_PRIM) = u(MAGZ_PRIM)
    ru(MAGZ_PRIM) = u(MAGX_PRIM)

end function

pure function zrotate_prim(u) result(ru)

    !! Rotate state vector to x-direction.

    implicit none

    real, intent(in)    :: u(NPRIM_VARS)
    real                :: ru(NPRIM_VARS)

    ru = u

    ru(VELX_PRIM) = u(VELZ_PRIM)
    ru(VELY_PRIM) = u(VELX_PRIM)
    ru(VELZ_PRIM) = u(VELY_PRIM)

    ru(MAGX_PRIM) = u(MAGZ_PRIM)
    ru(MAGY_PRIM) = u(MAGX_PRIM)
    ru(MAGZ_PRIM) = u(MAGY_PRIM)

end function

pure function yrotate_cons(u) result(ru)

    !! Rotate state vector to x-direction.

    implicit none

    real, intent(in)    :: u(NCONS_VARS)
    real                :: ru(NCONS_VARS)

    ru = u

    ru(MOMX_CONS) = u(MOMY_CONS)
    ru(MOMY_CONS) = u(MOMZ_CONS)
    ru(MOMZ_CONS) = u(MOMX_CONS)

    ru(MAGX_CONS) = u(MAGY_CONS)
    ru(MAGY_CONS) = u(MAGZ_CONS)
    ru(MAGZ_CONS) = u(MAGX_CONS)

end function

pure function zrotate_cons(u) result(ru)

    !! Rotate state vector to x-direction.

    implicit none

    real, intent(in)    :: u(NCONS_VARS)
    real                :: ru(NCONS_VARS)

    ru = u

    ru(MOMX_CONS) = u(MOMZ_CONS)
    ru(MOMY_CONS) = u(MOMX_CONS)
    ru(MOMZ_CONS) = u(MOMY_CONS)

    ru(MAGX_CONS) = u(MAGZ_CONS)
    ru(MAGY_CONS) = u(MAGX_CONS)
    ru(MAGZ_CONS) = u(MAGY_CONS)

end function

pure function yrotate_flux(f) result(rf)

    implicit none

    real, intent(in)    ::  f(NFLUXES)
    real                :: rf(NFLUXES)

    rf = f

    rf(XMOM_FLUX) = f(YMOM_FLUX)
    rf(YMOM_FLUX) = f(ZMOM_FLUX)
    rf(ZMOM_FLUX) = f(XMOM_FLUX)

    rf(MAGX_FLUX) = f(MAGY_FLUX)
    rf(MAGY_FLUX) = f(MAGZ_FLUX)
    rf(MAGZ_FLUX) = f(MAGX_FLUX)

end function

pure function zrotate_flux(f) result(rf)

    implicit none

    real, intent(in)    ::  f(NFLUXES)
    real                :: rf(NFLUXES)

    rf = f

    rf(XMOM_FLUX) = f(ZMOM_FLUX)
    rf(YMOM_FLUX) = f(XMOM_FLUX)
    rf(ZMOM_FLUX) = f(YMOM_FLUX)

    rf(MAGX_FLUX) = f(MAGZ_FLUX)
    rf(MAGY_FLUX) = f(MAGX_FLUX)
    rf(MAGZ_FLUX) = f(MAGY_FLUX)

end function

pure function yrotate_back_flux(rf) result(f)

    implicit none

    real, intent(in)    :: rf(NFLUXES)
    real                ::  f(NFLUXES)

    f = rf

    f(XMOM_FLUX) = rf(ZMOM_FLUX)
    f(YMOM_FLUX) = rf(XMOM_FLUX)
    f(ZMOM_FLUX) = rf(YMOM_FLUX)

    f(MAGX_FLUX) = rf(MAGZ_FLUX)
    f(MAGY_FLUX) = rf(MAGX_FLUX)
    f(MAGZ_FLUX) = rf(MAGY_FLUX)

end function

pure function zrotate_back_flux(rf) result(f)

    implicit none

    real, intent(in)    :: rf(NFLUXES)
    real                ::  f(NFLUXES)

    f = rf

    f(XMOM_FLUX) = rf(YMOM_FLUX)
    f(YMOM_FLUX) = rf(ZMOM_FLUX)
    f(ZMOM_FLUX) = rf(XMOM_FLUX)

    f(MAGX_FLUX) = rf(MAGY_FLUX)
    f(MAGY_FLUX) = rf(MAGZ_FLUX)
    f(MAGZ_FLUX) = rf(MAGX_FLUX)

end function

!! ========================================================================== !!

!! pure function xriemann_species_surface(f,uL,uR) result(sf)
!! 
!!     real, intent(in) :: f(NFLUXES),uL(NPRIM_VARS),uR(NPRIM_VARS)
!!     real             :: sf(NFLUXES)
!! 
!!     sf = f
!! 
!! # if NSPECIES > 0
!!     sf(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(SPEC_PRIM_BEGIN:SPEC_PRIM_END),UR(SPEC_PRIM_BEGIN:SPEC_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! # if NMASS_SCALARS > 0
!!     sf(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(MSCL_PRIM_BEGIN:MSCL_PRIM_END),UR(MSCL_PRIM_BEGIN:MSCL_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! 
!! end function
!! 
!! pure function yriemann_species_surface(f,uL,uR) result(sf)
!! 
!!     real, intent(in) :: f(NFLUXES),uL(NPRIM_VARS),uR(NPRIM_VARS)
!!     real             :: sf(NFLUXES)
!!      
!!     sf = f
!! 
!! # if NSPECIES > 0
!!     sf(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(SPEC_PRIM_BEGIN:SPEC_PRIM_END),UR(SPEC_PRIM_BEGIN:SPEC_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! # if NMASS_SCALARS > 0
!!     sf(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(MSCL_PRIM_BEGIN:MSCL_PRIM_END),UR(MSCL_PRIM_BEGIN:MSCL_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! 
!! end function
!! 
!! pure function zriemann_species_surface(f,uL,uR) result(sf)
!! 
!!     real, intent(in) :: f(NFLUXES),uL(NPRIM_VARS),uR(NPRIM_VARS)
!!     real             :: sf(NFLUXES)
!!      
!!     sf = f
!! 
!! # if NSPECIES > 0
!!     sf(SPECIES_FLUX_BEGIN:SPECIES_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(SPEC_PRIM_BEGIN:SPEC_PRIM_END),UR(SPEC_PRIM_BEGIN:SPEC_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! # if NMASS_SCALARS > 0
!!     sf(MASS_SCALARS_FLUX_BEGIN:MASS_SCALARS_FLUX_END) = &
!!         f(DENS_FLUX)*merge(UL(MSCL_PRIM_BEGIN:MSCL_PRIM_END),UR(MSCL_PRIM_BEGIN:MSCL_PRIM_END),F(DENS_FLUX) > 0.0)
!! # endif
!! 
!! end function

end module
