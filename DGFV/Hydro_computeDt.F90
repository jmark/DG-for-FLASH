!!****if* source/physics/Hydro/HydroMain/split/ES/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN) :: blockID,
!!                  real(IN)    :: x(:),
!!                  real(IN)    :: dx(:), 
!!                  real(IN)    :: uxgrid(:),
!!                  real(IN)    :: y(:), 
!!                  real(IN)    :: dy(:), 
!!                  real(IN)    :: uygrid(:), 
!!                  real(IN)    :: z(:), 
!!                  real(IN)    :: dz(:), 
!!                  real(IN)    :: uzgrid(:), 
!!                  integer(IN) :: blkLimits(2,MDIM)
!!                  integer(IN) :: blkLimitsGC(2,MDIM)
!!                  real,pointer   :: solnData(:,:,:,:),
!!                  real,(INOUT)   :: dtCheck, 
!!                  integer(INOUT) :: dtMinLoc(:),
!!                  real(INOUT), optional :: extraInfo)
!!
!!
!! DESCRIPTION
!!
!!  Computes the timestep DGFV solver. The Courant-Fredrichs-Lewy 
!!  criterion is used.  The maximum wave speed is computed and together with the velocities, 
!!  is used to constrain the timestep such that no information can propagate more than 
!!  one zone per timestep.
!!
!!
!! ARGUMENTS
!!
!!  blockID -       local blockID
!!  x, y, z -       coordinates
!!  dx, dy, dz -    deltas in each {x, y z} directions
!!  uxgrid, uygrid, uzgrid - velocity of grid expansion in {x, y z} directions
!!  blkLimits -     the indices for the interior endpoints of the block
!!  blkLimitsGC -   the indices for endpoints including the guardcells
!!  solnData -      the physical, solution data from grid
!!
!!  ALL ABOVE VARIABLES are dummy for this implementation of MHD.
!!  ONLY kept to conform to API standard
!! 
!!  dtCheck -      variable to hold timestep constraint
!!  dtMinLoc(5) -  array to hold location of cell responsible for minimum dt:
!!                 dtMinLoc(1) = i index
!!                 dtMinLoc(2) = j index
!!                 dtMinLoc(3) = k index
!!                 dtMinLoc(4) = blockID
!!                 dtMinLoc(5) = hy_meshMe
!! extraInfo    -  Driver_computeDt can provide extra info to the caller
!!                 using this argument.
!!
!! NOTES
!!
!!  For this implementation of MHD -- most of the arguments are dummy variables.
!!  More effient to just store the minimum timestep when you are doing the sweeps.
!!  This routine only serves to return the minimum values that have already been computed
!!  minimum hy_dt is actually computed in hy_es_sweep, where all required
!!  information is evaluated anyway. hy_es_sweep calls hy_es_setTimestep 
!!
!!
!!***

subroutine Hydro_computeDt ( blockID,  &
                             x, dx, uxgrid, &
                             y, dy, uygrid, &
                             z, dz, uzgrid, &
                             blkLimits,     &
                             blkLimitsGC,   &
                             solnData,      &
                             dtCheck, dtMinLoc,&
                             extraInfo )



  use Hydro_data,     ONLY : hy_dtmin, hy_dtminloc, hy_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks

  implicit none

#include "Flash.h"
#include "constants.h"

  !! Argument list ---------------------------------------------
  integer, intent(IN)       :: blockID 
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)        :: dtCheck
  integer,INTENT(INOUT)     :: dtMinLoc(5)
  real, pointer             :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif
  real,OPTIONAL,intent(INOUT) :: extraInfo
  !! -----------------------------------------------------------

  integer, dimension(MAXBLOCKS) :: blockList
  integer :: localNumBlocks, i

  if ( hy_dtmin < dtCheck ) then
     dtCheck  = hy_dtmin
     dtMinLoc = hy_dtminloc
  endif

  if(dtCheck <= 0.0) call Driver_abortFlash("[Hydro]: Computed dt is not > 0! Aborting!")

  !! Set default value to be null in not needed.
  if (present(extraInfo)) extraInfo = 0.

end subroutine Hydro_computeDt
