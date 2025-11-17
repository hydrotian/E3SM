module rof2iacMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module deals with arrays for exchanging data from MOSART
  ! to IAC (GCAM). Provides annual water availability metrics.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use RtmVar          , only : iulog, wrmflag
  use rof_cpl_indices , only : nt_nliq
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! rof -> iac variables structure
  ! Fields are dimensioned (ngrid)
  type, public :: rof2iac_type
     real(r8), pointer :: wr_avail(:) => null()       ! Main channel water availability (m3)
     real(r8), pointer :: wt_avail(:) => null()       ! Tributary water availability (m3)
     real(r8), pointer :: wtot_avail(:) => null()     ! Total surface water availability (m3)
     real(r8), pointer :: reservoir_stor(:) => null() ! Reservoir storage capacity (m3)
     real(r8), pointer :: streamflow(:) => null()     ! Annual mean streamflow (m3/s)

   contains
     procedure, public  :: Init
     procedure, public  :: update_rof2iac
  end type rof2iac_type

  ! !PUBLIC MEMBER FUNCTIONS:

contains

  !-------------------------------
  subroutine Init(this, begr, endr)

    ! !DESCRIPTION:
    ! Initialize MOSART variables required by IAC
    !
    ! !ARGUMENTS:
    class(rof2iac_type) :: this
    integer, intent(in) :: begr, endr  ! MOSART grid bounds

    ! Allocate arrays
    allocate(this%wr_avail(begr:endr))
    allocate(this%wt_avail(begr:endr))
    allocate(this%wtot_avail(begr:endr))
    allocate(this%reservoir_stor(begr:endr))
    allocate(this%streamflow(begr:endr))

    ! Initialize to zero
    this%wr_avail(:) = 0.0_r8
    this%wt_avail(:) = 0.0_r8
    this%wtot_avail(:) = 0.0_r8
    this%reservoir_stor(:) = 0.0_r8
    this%streamflow(:) = 0.0_r8

  end subroutine Init

  !------------------------------------------------------
  subroutine update_rof2iac(this, begr, endr, wr, wt, wh, erout, &
                             reservoir_stor_in, accum_time)
    !
    ! !DESCRIPTION:
    ! Calculate water availability metrics from MOSART state
    ! Uses accumulated/averaged MOSART data over the year
    !
    ! !ARGUMENTS:
    class(rof2iac_type), intent(inout) :: this
    integer, intent(in) :: begr, endr
    real(r8), intent(in) :: wr(:,:)              ! Main channel storage (m3)
    real(r8), intent(in) :: wt(:,:)              ! Tributary storage (m3)
    real(r8), intent(in) :: wh(:,:)              ! Hillslope storage (m)
    real(r8), intent(in) :: erout(:,:)           ! Channel outflow (m3/s)
    real(r8), intent(in) :: reservoir_stor_in(:) ! Reservoir storage (m3) [optional]
    real(r8), intent(in) :: accum_time           ! Total accumulation time (s)

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'update_rof2iac'
    integer :: n
    real(r8) :: env_flow_req
    real(r8), parameter :: ENV_FLOW_FRAC = 0.30_r8  ! 30% for environmental flows
    real(r8), parameter :: SEC_PER_YEAR = 31536000.0_r8

    ! Calculate water availability for each grid cell
    do n = begr, endr

       ! Annual mean streamflow (already in m3/s from accumulation)
       this%streamflow(n) = erout(n, nt_nliq)

       ! Total water passed through system (m3/year)
       ! This is the integral of flow over the year
       ! Using mean flow * seconds per year
       env_flow_req = ENV_FLOW_FRAC * this%streamflow(n) * SEC_PER_YEAR

       ! Extractable water = current storage (which represents annual average)
       ! We use the storage as a proxy for available water
       ! Main channel water availability
       this%wr_avail(n) = wr(n, nt_nliq)

       ! Tributary water availability
       this%wt_avail(n) = wt(n, nt_nliq)

       ! Total surface water availability
       ! Note: wh is in meters, need to convert to m3 using area
       ! For now, we'll use the sum of wr and wt as they're already in m3
       this%wtot_avail(n) = this%wr_avail(n) + this%wt_avail(n)

       ! Ensure non-negative
       if (this%wr_avail(n) < 0.0_r8) this%wr_avail(n) = 0.0_r8
       if (this%wt_avail(n) < 0.0_r8) this%wt_avail(n) = 0.0_r8
       if (this%wtot_avail(n) < 0.0_r8) this%wtot_avail(n) = 0.0_r8

       ! Reservoir storage (if WRM is enabled)
       if (wrmflag) then
          this%reservoir_stor(n) = reservoir_stor_in(n)
       else
          this%reservoir_stor(n) = 0.0_r8
       endif
    end do

  end subroutine update_rof2iac

end module rof2iacMod
