module iac2rofMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle coupled data from IAC (GCAM) for use in MOSART
  ! IAC is on a different grid (economic regions) than MOSART (river network)
  !
  ! !USES:
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use RtmVar         , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! iac -> rof structure
  ! Dimensioned by (ngrid)
  type, public :: iac2rof_type
     real(r8), pointer :: demand_irrig(:) => null()   ! Irrigation water demand (m3/year)
     real(r8), pointer :: demand_indust(:) => null()  ! Industrial water demand (m3/year)
     real(r8), pointer :: demand_munic(:) => null()   ! Municipal water demand (m3/year)
     real(r8), pointer :: demand_energy(:) => null()  ! Energy/cooling water demand (m3/year)
     real(r8), pointer :: demand_total(:) => null()   ! Total water demand (m3/year)
     real(r8), pointer :: consump_frac(:) => null()   ! Consumptive fraction (0-1)

   contains
     procedure, public :: Init
     procedure, public :: update_iac2rof
  end type iac2rof_type

  ! Number of water demand sectors
  integer, parameter, public :: nsectors_gcam = 4

contains

  !-------------------------------
  subroutine Init(this, begr, endr)
    !
    ! !DESCRIPTION
    ! Allocate and initialize iac variables used by MOSART
    !
    ! !ARGUMENTS:
    class(iac2rof_type) :: this
    integer, intent(in) :: begr, endr  ! MOSART grid bounds

    ! Allocate arrays
    allocate(this%demand_irrig(begr:endr))
    allocate(this%demand_indust(begr:endr))
    allocate(this%demand_munic(begr:endr))
    allocate(this%demand_energy(begr:endr))
    allocate(this%demand_total(begr:endr))
    allocate(this%consump_frac(begr:endr))

    ! Initialize to zero
    this%demand_irrig(:) = 0.0_r8
    this%demand_indust(:) = 0.0_r8
    this%demand_munic(:) = 0.0_r8
    this%demand_energy(:) = 0.0_r8
    this%demand_total(:) = 0.0_r8
    this%consump_frac(:) = 0.0_r8

  end subroutine Init

  !------------------------------------------------------
  subroutine update_iac2rof(this, begr, endr, qdem_gcam, qdem_total, consump_frac)
    !
    ! !DESCRIPTION:
    ! Convert IAC water demand from annual (m3/year) to rate (m3/s)
    ! and apply to MOSART arrays
    !
    ! !ARGUMENTS:
    class(iac2rof_type), intent(in) :: this
    integer, intent(in) :: begr, endr
    real(r8), intent(inout) :: qdem_gcam(:,:)    ! GCAM demand by sector (m3/s)
    real(r8), intent(inout) :: qdem_total(:)     ! Total GCAM demand (m3/s)
    real(r8), intent(inout) :: consump_frac(:)   ! Consumptive fraction

    ! !LOCAL VARIABLES:
    integer :: n
    real(r8), parameter :: SEC_PER_DAY = 86400.0_r8
    real(r8), parameter :: DAYS_PER_YEAR = 365.0_r8
    real(r8) :: seconds_per_year
    character(len=*), parameter :: subname = 'update_iac2rof'

    seconds_per_year = SEC_PER_DAY * DAYS_PER_YEAR

    ! Convert annual demand (m3/year) to rate (m3/s) and apply
    do n = begr, endr
       ! Convert each sector's demand from m3/year to m3/s
       qdem_gcam(n, 1) = this%demand_irrig(n) / seconds_per_year
       qdem_gcam(n, 2) = this%demand_indust(n) / seconds_per_year
       qdem_gcam(n, 3) = this%demand_munic(n) / seconds_per_year
       qdem_gcam(n, 4) = this%demand_energy(n) / seconds_per_year

       ! Total demand
       qdem_total(n) = this%demand_total(n) / seconds_per_year

       ! Consumptive fraction (no conversion needed)
       consump_frac(n) = this%consump_frac(n)
    end do

  end subroutine update_iac2rof

end module iac2rofMod
