# MOSART-GCAM Water Management Coupling Design

## Executive Summary

This document presents a detailed design for coupling the water management components between MOSART (Model for Scale Adaptive River Transport) and GCAM (Global Change Analysis Model) within the E3SM framework. The design leverages the existing E3SM-GCAM coupling infrastructure and extends it to enable bidirectional water information exchange: MOSART provides surface water availability to GCAM, and GCAM returns water demand feedback to MOSART.

**Design Date:** 2025-11-17
**Target Branch:** `claude/review-water-demand-01GsYnd5pLYDLfhhLhCQuXJr`
**Framework:** Extends existing IAC coupling architecture

---

## 1. Motivation and Objectives

### 1.1 Scientific Motivation

Current E3SM-GCAM coupling addresses:
- ✅ Terrestrial productivity (ELM → GCAM)
- ✅ Land use changes (GCAM → ELM)
- ✅ CO2 emissions (GCAM → EAM)
- ❌ **Water resource management** (NOT YET COUPLED)

**Gap:** GCAM makes water allocation decisions without knowledge of actual surface water availability, and MOSART does not account for human water demands from GCAM's economic modeling.

### 1.2 Objectives

**Primary Objectives:**
1. Enable MOSART to send surface water availability to GCAM annually
2. Enable GCAM to send water demand back to MOSART for hydrologic feedback
3. Leverage existing IAC coupling infrastructure (minimal code changes)
4. Maintain consistency with current coupling patterns

**Secondary Objectives:**
1. Support multiple water demand sectors (irrigation, industrial, municipal, energy)
2. Enable spatial resolution matching between MOSART grid and GCAM regions
3. Provide temporal resolution appropriate for water management (annual/seasonal)
4. Ensure restart capability and conservation

---

## 2. Architecture Overview

### 2.1 Coupling Pattern

**Leverage Existing Infrastructure:**
```
┌─────────────────────────────────────────────────────────────────┐
│                    E3SM Coupling Infrastructure                  │
│                                                                  │
│  ┌──────────┐         ┌──────────┐         ┌──────────┐        │
│  │   ELM    │────────▶│          │────────▶│   ELM    │        │
│  │          │◀────────│   GCAM   │◀────────│          │        │
│  └──────────┘  NPP,   │   (IAC)  │  Land   └──────────┘        │
│                frac    │          │  Use                        │
│                        │          │                             │
│  ┌──────────┐         │          │         ┌──────────┐        │
│  │   EAM    │◀────────│          │         │   EAM    │        │
│  │          │  CO2    │          │         │          │        │
│  └──────────┘         │          │         └──────────┘        │
│                        │          │                             │
│  ┌──────────┐         │          │         ┌──────────┐        │
│  │  MOSART  │────────▶│          │────────▶│  MOSART  │        │
│  │          │  Water  │          │  Water  │          │        │
│  │          │  Avail  │          │  Demand │          │        │
│  └──────────┘         └──────────┘         └──────────┘        │
│                            NEW                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Key Design Principle:** Parallel structure to ELM-GCAM coupling:
- Use same driver infrastructure (`prep_iac_mod.F90`)
- Use same MCT attribute vectors
- Use same annual coupling frequency
- Use same grid mapping mechanisms

### 2.2 Data Flow Diagram

```
Annual Cycle:
├─ Throughout Year (Days 1-365):
│  ├─ MOSART runs every hour:
│  │  ├─ Calculates surface water storage (wr, wt, wh)
│  │  ├─ Tracks channel flow (erout)
│  │  ├─ Monitors reservoir storage (if wrmflag)
│  │  └─ Accumulates water availability metrics
│  │
│  └─ Accumulation (similar to ELM NPP accumulation):
│     └─ prep_iac_accum_rof(): Accumulate MOSART water data
│
├─ Year Boundary (Day 1, 00:30 UTC):
│  ├─ Finalize MOSART accumulation:
│  │  ├─ prep_iac_accum_avg_rof(): Annual average
│  │  ├─ Calculate total water availability per grid cell
│  │  └─ Map MOSART grid → IAC grid
│  │
│  ├─ GCAM runs:
│  │  ├─ Receives: Surface water availability by region
│  │  ├─ Computes: Water allocation across sectors
│  │  │  ├─ Irrigation demand
│  │  │  ├─ Industrial demand
│  │  │  ├─ Municipal demand
│  │  │  └─ Energy/cooling demand
│  │  └─ Sends: Water demand by sector and region
│  │
│  └─ Apply to MOSART:
│     ├─ Map IAC grid → MOSART grid
│     ├─ Update MOSART water demand (rtmCTL%qdem_gcam)
│     └─ Temporal interpolation (annual → monthly/daily)
│
└─ Post-processing:
   └─ Reset accumulators for next year
```

---

## 3. Coupling Fields Definition

### 3.1 MOSART → GCAM (Annual Water Availability)

**New Fields in `seq_flds_mod.F90`:**

```fortran
! River/runoff to IAC (r2z) - Annual water availability
character(CXX) :: seq_flds_r2z_states

! Field list (to be added to seq_flds_r2z_states):
"Sr_wr_avail"        ! Main channel water availability (m³/year or mm/year)
"Sr_wt_avail"        ! Tributary/subnetwork water availability (m³/year)
"Sr_wtot_avail"      ! Total surface water availability (m³/year)
"Sr_reservoir_stor"  ! Reservoir storage capacity (m³) [if wrmflag]
"Sr_streamflow"      ! Annual mean streamflow (m³/s)
```

**Rationale for Each Field:**

| Field | Purpose | GCAM Usage |
|-------|---------|------------|
| `Sr_wr_avail` | Main channel water | Direct withdrawals from rivers |
| `Sr_wt_avail` | Tributary water | Local/distributed withdrawals |
| `Sr_wtot_avail` | Total availability | Overall water constraint |
| `Sr_reservoir_stor` | Storage capacity | Reservoir-based water supply |
| `Sr_streamflow` | Mean flow rate | Environmental flow requirements |

**Units:**
- Storage: m³ or mm equivalent over gridcell area
- Flow rate: m³/s
- Availability: m³ over the year (integral of flow)

### 3.2 GCAM → MOSART (Annual Water Demand)

**New Fields in `seq_flds_mod.F90`:**

```fortran
! IAC to river/runoff (z2r) - Water demand by sector
character(CXX) :: seq_flds_z2r_fluxes

! Field list (to be added to seq_flds_z2r_fluxes):
"Sz_demand_irrig"    ! Irrigation water demand (m³/year)
"Sz_demand_indust"   ! Industrial water demand (m³/year)
"Sz_demand_munic"    ! Municipal water demand (m³/year)
"Sz_demand_energy"   ! Energy/cooling water demand (m³/year)
"Sz_demand_total"    ! Total water demand (m³/year)
"Sz_consump_frac"    ! Consumptive fraction (dimensionless, 0-1)
```

**Rationale:**

| Field | Purpose | MOSART Usage |
|-------|---------|------------|
| `Sz_demand_irrig` | Crop irrigation | Irrigation water withdrawal from rivers/reservoirs |
| `Sz_demand_indust` | Industrial use | Industrial water allocation |
| `Sz_demand_munic` | Municipal use | Drinking water, urban supply |
| `Sz_demand_energy` | Power plant cooling | Thermal power generation |
| `Sz_demand_total` | Sum of all sectors | Total withdrawal constraint |
| `Sz_consump_frac` | Consumptive vs. return | How much water returns to stream |

**Units:**
- Demand: m³/year (total annual demand)
- Fraction: 0-1 (dimensionless)

**Conversion to MOSART internal units:**
```fortran
! m³/year → m³/s
rtmCTL%qdem_gcam(n, sector) = Sz_demand(n) / seconds_per_year
```

---

## 4. Implementation Details

### 4.1 New Modules and Files

#### 4.1.1 MOSART Export Module

**New File:** `components/mosart/src/cpl/rof2iacMod.F90`

**Purpose:** Package MOSART water data for export to IAC (parallel to `lnd2iacMod.F90`)

**Data Structure:**
```fortran
module rof2iacMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use RunoffMod, only : runoff_flow

  implicit none
  private

  type, public :: rof2iac_type
     ! Water availability (m³ or mm)
     real(r8), pointer :: wr_avail(:)       ! Main channel water availability
     real(r8), pointer :: wt_avail(:)       ! Tributary water availability
     real(r8), pointer :: wtot_avail(:)     ! Total surface water availability
     real(r8), pointer :: reservoir_stor(:) ! Reservoir storage [if wrmflag]
     real(r8), pointer :: streamflow(:)     ! Annual mean streamflow (m³/s)

   contains
     procedure, public :: Init
     procedure, public :: update_rof2iac
  end type rof2iac_type

contains

  subroutine Init(this, begr, endr)
    class(rof2iac_type) :: this
    integer, intent(in) :: begr, endr  ! MOSART grid bounds

    allocate(this%wr_avail(begr:endr))
    allocate(this%wt_avail(begr:endr))
    allocate(this%wtot_avail(begr:endr))
    allocate(this%reservoir_stor(begr:endr))
    allocate(this%streamflow(begr:endr))

    this%wr_avail(:) = 0.0_r8
    this%wt_avail(:) = 0.0_r8
    this%wtot_avail(:) = 0.0_r8
    this%reservoir_stor(:) = 0.0_r8
    this%streamflow(:) = 0.0_r8
  end subroutine Init

  subroutine update_rof2iac(this, begr, endr, rtmCTL, StorWater, &
                             accumulation_time)
    use RunoffMod, only : Trunoff
    use WRM_type_mod, only : WRMwater
    use RtmVar, only : wrmflag

    class(rof2iac_type), intent(inout) :: this
    integer, intent(in) :: begr, endr
    type(runoff_flow), intent(in) :: rtmCTL
    type(WRMwater), intent(in) :: StorWater
    real(r8), intent(in) :: accumulation_time  ! Total accumulation time (s)

    integer :: n
    real(r8) :: dt_year

    dt_year = accumulation_time  ! seconds

    do n = begr, endr
       ! Main channel water availability (integral over year)
       ! wr is in m³, integrate over timesteps
       this%wr_avail(n) = rtmCTL%wr(n, LIQ)

       ! Tributary water availability
       this%wt_avail(n) = rtmCTL%wt(n, LIQ)

       ! Total surface water availability
       this%wtot_avail(n) = rtmCTL%wr(n, LIQ) + rtmCTL%wt(n, LIQ) + &
                            rtmCTL%wh(n, LIQ)

       ! Mean streamflow (accumulated flow / time)
       this%streamflow(n) = rtmCTL%erout(n, LIQ)  ! Already in m³/s

       ! Reservoir storage (if WRM is enabled)
       if (wrmflag) then
          this%reservoir_stor(n) = StorWater%storage(n)
       else
          this%reservoir_stor(n) = 0.0_r8
       endif
    end do
  end subroutine update_rof2iac

end module rof2iacMod
```

#### 4.1.2 MOSART Import Module

**New File:** `components/mosart/src/cpl/iac2rofMod.F90`

**Purpose:** Apply GCAM water demand to MOSART (parallel to `iac2lndMod.F90`)

**Data Structure:**
```fortran
module iac2rofMod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  type, public :: iac2rof_type
     ! Water demand by sector (m³/year)
     real(r8), pointer :: demand_irrig(:)
     real(r8), pointer :: demand_indust(:)
     real(r8), pointer :: demand_munic(:)
     real(r8), pointer :: demand_energy(:)
     real(r8), pointer :: demand_total(:)
     real(r8), pointer :: consump_frac(:)  ! Consumptive fraction (0-1)

   contains
     procedure, public :: Init
     procedure, public :: update_iac2rof
  end type iac2rof_type

contains

  subroutine Init(this, begr, endr)
    class(iac2rof_type) :: this
    integer, intent(in) :: begr, endr

    allocate(this%demand_irrig(begr:endr))
    allocate(this%demand_indust(begr:endr))
    allocate(this%demand_munic(begr:endr))
    allocate(this%demand_energy(begr:endr))
    allocate(this%demand_total(begr:endr))
    allocate(this%consump_frac(begr:endr))

    this%demand_irrig(:) = 0.0_r8
    this%demand_indust(:) = 0.0_r8
    this%demand_munic(:) = 0.0_r8
    this%demand_energy(:) = 0.0_r8
    this%demand_total(:) = 0.0_r8
    this%consump_frac(:) = 0.0_r8
  end subroutine Init

  subroutine update_iac2rof(this, begr, endr, rtmCTL)
    use RunoffMod, only : Trunoff
    use RtmTimeManager, only : get_step_size

    class(iac2rof_type), intent(in) :: this
    integer, intent(in) :: begr, endr
    type(runoff_flow), intent(inout) :: rtmCTL

    integer :: n
    real(r8) :: seconds_per_year
    real(r8), parameter :: SEC_PER_DAY = 86400.0_r8
    real(r8), parameter :: DAYS_PER_YEAR = 365.0_r8

    seconds_per_year = SEC_PER_DAY * DAYS_PER_YEAR

    do n = begr, endr
       ! Convert annual demand (m³/year) to rate (m³/s)
       ! Store in new array rtmCTL%qdem_gcam (to be added to RunoffMod)
       rtmCTL%qdem_gcam(n, 1) = this%demand_irrig(n) / seconds_per_year
       rtmCTL%qdem_gcam(n, 2) = this%demand_indust(n) / seconds_per_year
       rtmCTL%qdem_gcam(n, 3) = this%demand_munic(n) / seconds_per_year
       rtmCTL%qdem_gcam(n, 4) = this%demand_energy(n) / seconds_per_year

       ! Total demand
       rtmCTL%qdem_total(n) = this%demand_total(n) / seconds_per_year

       ! Consumptive fraction (no conversion)
       rtmCTL%consump_frac(n) = this%consump_frac(n)
    end do
  end subroutine update_iac2rof

end module iac2rofMod
```

### 4.2 Modifications to Existing Files

#### 4.2.1 RunoffMod.F90 - Add GCAM Demand Fields

**File:** `components/mosart/src/riverroute/RunoffMod.F90`

**Additions to `runoff_flow` type:**
```fortran
type runoff_flow
   ! ... existing fields ...

   ! GCAM water demand (m³/s)
   real(r8), pointer :: qdem_gcam(:,:) => null()  ! (nr, nsectors)
   real(r8), pointer :: qdem_total(:) => null()    ! Total GCAM demand (m³/s)
   real(r8), pointer :: consump_frac(:) => null()  ! Consumptive fraction

   ! GCAM water supply (for feedback to GCAM)
   real(r8), pointer :: supply_gcam(:,:) => null()  ! Supply by sector (m³/s)
   real(r8), pointer :: deficit_gcam(:,:) => null() ! Deficit by sector (m³/s)

end type runoff_flow
```

**Allocation in initialization:**
```fortran
subroutine RunoffInit()
   ! ... existing code ...

   ! GCAM demand fields
   allocate(rtmCTL%qdem_gcam(begr:endr, 1:4))  ! 4 sectors
   allocate(rtmCTL%qdem_total(begr:endr))
   allocate(rtmCTL%consump_frac(begr:endr))
   allocate(rtmCTL%supply_gcam(begr:endr, 1:4))
   allocate(rtmCTL%deficit_gcam(begr:endr, 1:4))

   rtmCTL%qdem_gcam(:,:) = 0.0_r8
   rtmCTL%qdem_total(:) = 0.0_r8
   rtmCTL%consump_frac(:) = 0.0_r8
   rtmCTL%supply_gcam(:,:) = 0.0_r8
   rtmCTL%deficit_gcam(:,:) = 0.0_r8
end subroutine RunoffInit
```

#### 4.2.2 rof_cpl_indices.F90 - Add Coupling Indices

**File:** `components/mosart/src/cpl/rof_cpl_indices.F90`

**Add index variables:**
```fortran
! r2x (MOSART to coupler) - Water availability to IAC
integer, public :: index_r2x_Sr_wr_avail      = 0
integer, public :: index_r2x_Sr_wt_avail      = 0
integer, public :: index_r2x_Sr_wtot_avail    = 0
integer, public :: index_r2x_Sr_reservoir_stor = 0
integer, public :: index_r2x_Sr_streamflow    = 0

! x2r (coupler to MOSART) - Water demand from IAC
integer, public :: index_x2r_Sz_demand_irrig  = 0
integer, public :: index_x2r_Sz_demand_indust = 0
integer, public :: index_x2r_Sz_demand_munic  = 0
integer, public :: index_x2r_Sz_demand_energy = 0
integer, public :: index_x2r_Sz_demand_total  = 0
integer, public :: index_x2r_Sz_consump_frac  = 0
```

**Set indices in `rof_cpl_indices_set()`:**
```fortran
subroutine rof_cpl_indices_set()
   ! ... existing code ...

   ! IAC coupling (if iac_present)
   if (iac_present) then
      ! MOSART → IAC (water availability)
      index_r2x_Sr_wr_avail      = mct_aVect_indexRA(r2x, 'Sr_wr_avail')
      index_r2x_Sr_wt_avail      = mct_aVect_indexRA(r2x, 'Sr_wt_avail')
      index_r2x_Sr_wtot_avail    = mct_aVect_indexRA(r2x, 'Sr_wtot_avail')
      index_r2x_Sr_reservoir_stor = mct_aVect_indexRA(r2x, 'Sr_reservoir_stor')
      index_r2x_Sr_streamflow    = mct_aVect_indexRA(r2x, 'Sr_streamflow')

      ! IAC → MOSART (water demand)
      index_x2r_Sz_demand_irrig  = mct_aVect_indexRA(x2r, 'Sz_demand_irrig')
      index_x2r_Sz_demand_indust = mct_aVect_indexRA(x2r, 'Sz_demand_indust')
      index_x2r_Sz_demand_munic  = mct_aVect_indexRA(x2r, 'Sz_demand_munic')
      index_x2r_Sz_demand_energy = mct_aVect_indexRA(x2r, 'Sz_demand_energy')
      index_x2r_Sz_demand_total  = mct_aVect_indexRA(x2r, 'Sz_demand_total')
      index_x2r_Sz_consump_frac  = mct_aVect_indexRA(x2r, 'Sz_consump_frac')
   endif
end subroutine rof_cpl_indices_set
```

#### 4.2.3 rof_comp_mct.F90 - Export/Import

**File:** `components/mosart/src/cpl/rof_comp_mct.F90`

**Import from IAC (in `rof_import_mct()`):**
```fortran
subroutine rof_import_mct(x2r_r)
   ! ... existing imports ...

   ! Import GCAM water demand (if iac_present and iac_prognostic)
   if (iac_present .and. iac_prognostic) then
      do n = begr, endr
         ni = n - begr + 1

         iac2rof%demand_irrig(n)  = x2r_r%rAttr(index_x2r_Sz_demand_irrig, ni)
         iac2rof%demand_indust(n) = x2r_r%rAttr(index_x2r_Sz_demand_indust, ni)
         iac2rof%demand_munic(n)  = x2r_r%rAttr(index_x2r_Sz_demand_munic, ni)
         iac2rof%demand_energy(n) = x2r_r%rAttr(index_x2r_Sz_demand_energy, ni)
         iac2rof%demand_total(n)  = x2r_r%rAttr(index_x2r_Sz_demand_total, ni)
         iac2rof%consump_frac(n)  = x2r_r%rAttr(index_x2r_Sz_consump_frac, ni)
      end do

      ! Apply to MOSART internal arrays
      call iac2rof%update_iac2rof(begr, endr, rtmCTL)
   endif
end subroutine rof_import_mct
```

**Export to IAC (in `rof_export_mct()`):**
```fortran
subroutine rof_export_mct(r2x_r)
   ! ... existing exports ...

   ! Export MOSART water availability to GCAM (if iac_present)
   if (iac_present .and. iac_prognostic) then
      do n = begr, endr
         ni = n - begr + 1

         r2x_r%rAttr(index_r2x_Sr_wr_avail, ni) = rof2iac%wr_avail(n)
         r2x_r%rAttr(index_r2x_Sr_wt_avail, ni) = rof2iac%wt_avail(n)
         r2x_r%rAttr(index_r2x_Sr_wtot_avail, ni) = rof2iac%wtot_avail(n)
         r2x_r%rAttr(index_r2x_Sr_reservoir_stor, ni) = rof2iac%reservoir_stor(n)
         r2x_r%rAttr(index_r2x_Sr_streamflow, ni) = rof2iac%streamflow(n)
      end do
   endif
end subroutine rof_export_mct
```

#### 4.2.4 prep_iac_mod.F90 - Add MOSART Accumulation

**File:** `driver-mct/main/prep_iac_mod.F90`

**Add MOSART attribute vectors:**
```fortran
! MOSART export to IAC (on IAC grid, cpl PEs)
type(mct_aVect), pointer :: r2x_zx(:)

! MOSART accumulation (on MOSART grid, cpl PEs)
type(mct_aVect), pointer :: r2zacc_rx(:)   ! Accumulated MOSART export
integer, target          :: r2zacc_rx_cnt  ! Accumulation counter

! Mapper for MOSART to IAC
type(seq_map), pointer :: mapper_Sr2z
```

**Add accumulation functions:**
```fortran
subroutine prep_iac_accum_rof(timer)
   ! Similar to prep_iac_accum() but for MOSART
   integer :: eri
   type(mct_aVect), pointer :: r2x_rx

   call t_drvstartf(trim(timer), barrier=mpicom_CPLID)
   do eri = 1, num_inst_rof
      r2x_rx => component_get_c2x_cx(rof(eri))
      if (r2zacc_rx_cnt == 0) then
         call mct_avect_copy(r2x_rx, r2zacc_rx(eri))
      else
         call mct_avect_accum(r2x_rx, r2zacc_rx(eri))
      endif
   end do
   r2zacc_rx_cnt = r2zacc_rx_cnt + 1
   call t_drvstopf(trim(timer))
end subroutine prep_iac_accum_rof

subroutine prep_iac_accum_avg_rof(timer)
   ! Annual averaging for MOSART
   integer :: eri

   call t_drvstartf(trim(timer), barrier=mpicom_CPLID)
   if (r2zacc_rx_cnt > 1) then
      do eri = 1, num_inst_rof
         call mct_avect_avg(r2zacc_rx(eri), r2zacc_rx_cnt)
      end do
   endif
   r2zacc_rx_cnt = 0
   call t_drvstopf(trim(timer))
end subroutine prep_iac_accum_avg_rof

subroutine prep_iac_calc_r2x_zx(timer)
   ! Map MOSART accumulated data to IAC grid
   integer :: eri, ezi

   call t_drvstartf(trim(timer), barrier=mpicom_CPLID)
   do eri = 1, num_inst_rof
      ezi = mod((eri-1), num_inst_iac) + 1
      call seq_map_map(mapper_Sr2z, r2zacc_rx(eri), r2x_zx(ezi), &
                       fldlist=seq_flds_r2x_states, norm=.true.)
   end do
   call t_drvstopf(trim(timer))
end subroutine prep_iac_calc_r2x_zx
```

#### 4.2.5 seq_flds_mod.F90 - Field Definitions

**File:** `driver-mct/shr/seq_flds_mod.F90`

**Add field list variables:**
```fortran
! ROF to IAC fields
character(CXX), public :: seq_flds_r2z_states = ''

! IAC to ROF fields
character(CXX), public :: seq_flds_z2r_fluxes = ''
```

**Build field lists (in module initialization):**
```fortran
! MOSART to IAC (water availability)
seq_flds_r2z_states = 'Sr_wr_avail:Sr_wt_avail:Sr_wtot_avail:' // &
                      'Sr_reservoir_stor:Sr_streamflow'

! IAC to MOSART (water demand)
seq_flds_z2r_fluxes = 'Sz_demand_irrig:Sz_demand_indust:' // &
                      'Sz_demand_munic:Sz_demand_energy:' // &
                      'Sz_demand_total:Sz_consump_frac'
```

#### 4.2.6 cime_comp_mod.F90 - Driver Integration

**File:** `driver-mct/main/cime_comp_mod.F90`

**Add flags and variables:**
```fortran
logical :: rof_c2_iac    ! .true. => ROF to IAC coupling on
logical :: iac_c2_rof    ! .true. => IAC to ROF coupling on
```

**Add to coupling sequence in `cime_run_iac_setup_send()`:**
```fortran
subroutine cime_run_iac_setup_send()
   ! ... existing land accumulation ...

   ! MOSART accumulation and averaging
   if (rof_present .and. rof_c2_iac .and. do_iacrun_send) then
      call t_drvstartf('DRIVER_ROF_IAC_ACCUM')
      call prep_iac_accum_avg_rof('prep_iac_accum_avg_rof')
      call t_drvstopf('DRIVER_ROF_IAC_ACCUM')

      call t_drvstartf('DRIVER_R2Z_MAP')
      call prep_iac_calc_r2x_zx('prep_iac_calc_r2x_zx')
      call t_drvstopf('DRIVER_R2Z_MAP')

      call t_drvstartf('DRIVER_R2Z_MRG')
      call prep_iac_mrg_rof(infodata, 'prep_iac_mrg_rof')
      call t_drvstopf('DRIVER_R2Z_MRG')
   endif

   ! ... existing IAC exchange ...
end subroutine cime_run_iac_setup_send
```

---

## 5. Water Availability Calculation

### 5.1 Metrics for Water Availability

**Key Concept:** Water availability to GCAM should represent extractable surface water, not just instantaneous storage.

**Options for Calculation:**

**Option 1: Annual Mean Storage**
```fortran
! Simple average of water storage over year
wr_avail(n) = SUM(wr(n,LIQ) * dt) / total_time
```

**Option 2: Annual Integrated Flow**
```fortran
! Total water volume that flowed through
wr_avail(n) = SUM(erout(n,LIQ) * dt)  ! m³
```

**Option 3: Renewable Water Resource**
```fortran
! Flow minus environmental requirements
env_flow_req = 0.3 * mean_annual_flow  ! 30% for ecosystems
wr_avail(n) = SUM(erout(n,LIQ) * dt) - env_flow_req * seconds_per_year
```

**Recommendation:** Use **Option 3** with configurable environmental flow fraction.

### 5.2 Implementation in `rof2iacMod.F90`

```fortran
subroutine calculate_water_availability(this, begr, endr, rtmCTL, &
                                        annual_flow_accum, dt_year)
   class(rof2iac_type), intent(inout) :: this
   integer, intent(in) :: begr, endr
   type(runoff_flow), intent(in) :: rtmCTL
   real(r8), intent(in) :: annual_flow_accum(:)  ! Accumulated flow (m³)
   real(r8), intent(in) :: dt_year               ! Total time (s)

   integer :: n
   real(r8) :: mean_flow, env_flow_req
   real(r8), parameter :: ENV_FLOW_FRAC = 0.30_r8  ! 30% environmental flow

   do n = begr, endr
      ! Mean annual streamflow (m³/s)
      mean_flow = annual_flow_accum(n) / dt_year
      this%streamflow(n) = mean_flow

      ! Environmental flow requirement
      env_flow_req = ENV_FLOW_FRAC * mean_flow * dt_year  ! m³/year

      ! Extractable water = total flow - environmental flow
      this%wtot_avail(n) = annual_flow_accum(n) - env_flow_req

      ! Ensure non-negative
      if (this%wtot_avail(n) < 0.0_r8) this%wtot_avail(n) = 0.0_r8

      ! Main channel availability (assume 70% of total in main channel)
      this%wr_avail(n) = 0.70_r8 * this%wtot_avail(n)

      ! Tributary availability (remaining 30%)
      this%wt_avail(n) = 0.30_r8 * this%wtot_avail(n)
   end do
end subroutine calculate_water_availability
```

---

## 6. Water Demand Application in MOSART

### 6.1 Integration with WRM Module

**Current WRM Structure:**
- MOSART already has Water Resources Management (WRM) module
- WRM manages reservoirs and water allocation
- Current demand comes from ELM via `Flrl_demand`

**New GCAM Demand Integration:**
1. Add GCAM demand as separate source
2. WRM prioritizes demand sources
3. Feedback deficit to GCAM for next year

**Modified WRM Allocation Logic:**
```fortran
! In WRM_modules.F90
subroutine WRM_water_allocation(n)
   integer, intent(in) :: n  ! Grid cell index

   real(r8) :: total_demand, available_supply
   real(r8) :: demand_elim, demand_gcam_total
   real(r8) :: allocation_frac

   ! Total demand from both sources
   demand_elim = rtmCTL%qdem(n, LIQ)  ! From ELM
   demand_gcam_total = rtmCTL%qdem_total(n)  ! From GCAM

   total_demand = demand_elim + demand_gcam_total

   ! Available supply (from reservoirs + streamflow)
   available_supply = calculate_available_supply(n)

   ! Allocation fraction
   if (total_demand > 0.0_r8) then
      allocation_frac = min(1.0_r8, available_supply / total_demand)
   else
      allocation_frac = 1.0_r8
   endif

   ! Allocate to sectors proportionally
   do sector = 1, 4
      rtmCTL%supply_gcam(n, sector) = rtmCTL%qdem_gcam(n, sector) * allocation_frac
      rtmCTL%deficit_gcam(n, sector) = rtmCTL%qdem_gcam(n, sector) * (1.0_r8 - allocation_frac)
   end do

   ! Update reservoir storage after withdrawals
   call update_reservoir_storage(n, available_supply * allocation_frac)
end subroutine WRM_water_allocation
```

### 6.2 Temporal Downscaling

**Challenge:** GCAM provides annual demand, but MOSART runs hourly.

**Solution:** Temporal downscaling with seasonality

**Option 1: Uniform Distribution**
```fortran
! Divide annual demand equally across year
demand_hourly = demand_annual / hours_per_year
```

**Option 2: Monthly Distribution**
```fortran
! Use monthly fractions (provided by GCAM or climatology)
month_frac(1:12) = [0.08, 0.07, 0.08, 0.09, 0.10, 0.11, &
                     0.12, 0.11, 0.09, 0.08, 0.07, 0.06]  ! Summer peak
demand_monthly(m) = demand_annual * month_frac(m)
demand_hourly = demand_monthly(current_month) / hours_in_month
```

**Option 3: Daily Interpolation (like CO2)**
```fortran
! Linear interpolation between monthly mid-points
! Reuse iac_coupled_timeinterp() pattern
call gcam_demand_timeinterp(year_day, month, lower_bound, upper_bound, tfrac)
demand_daily = demand_monthly(lower_bound) * (1.0 - tfrac) + &
               demand_monthly(upper_bound) * tfrac
```

**Recommendation:** Use **Option 2** with monthly fractions provided by GCAM.

---

## 7. Grid Mapping

### 7.1 MOSART Grid to IAC Grid

**Challenge:** MOSART uses a river network grid (typically 1/8° or 1/4°), while GCAM uses economic regions (32-235 regions globally).

**Mapping Requirements:**
1. Conservative mapping for water volumes
2. Aggregation from fine MOSART grid to coarse GCAM regions
3. Disaggregation from GCAM regions back to MOSART grid

### 7.2 Mapping Files

**Generate Mapping Files (offline preprocessing):**

```bash
# SCRIP regridding weight generation
# MOSART grid → GCAM regions
ESMF_RegridWeightGen \
  --source mosart_grid_descriptor.nc \
  --destination gcam_regions.nc \
  --weight mosart_to_gcam_conserv.nc \
  --method conserve \
  --netcdf4

# GCAM regions → MOSART grid
ESMF_RegridWeightGen \
  --source gcam_regions.nc \
  --destination mosart_grid_descriptor.nc \
  --weight gcam_to_mosart_conserv.nc \
  --method conserve \
  --netcdf4
```

**Configuration in `seq_maps.rc`:**
```
rof2iac_smapname: mosart_to_gcam_conserv
rof2iac_smaptype: conserv
iac2rof_smapname: gcam_to_mosart_conserv
iac2rof_smaptype: conserv
```

### 7.3 Aggregation Strategy

**MOSART → GCAM (Upscaling):**
- Sum water volumes within each GCAM region
- Area-weighted averaging for intensive quantities

**GCAM → MOSART (Downscaling):**
- Distribute regional demand to MOSART cells
- Weight by cell area, population, or irrigated area
- Options:
  - Uniform: Equal demand per unit area
  - Population-weighted: Proportional to population density
  - Irrigation-weighted: Proportional to irrigated fraction

---

## 8. Namelist Configuration

### 8.1 MOSART Namelist

**File:** `components/mosart/bld/namelist_files/namelist_definition_mosart.xml`

**New Parameters:**
```xml
<entry id="gcam_coupling" type="logical" category="coupling">
  <values>
    <value>.false.</value>
  </values>
  <desc>Enable MOSART-GCAM water coupling</desc>
</entry>

<entry id="gcam_demand_downscale_method" type="char" category="coupling">
  <values>
    <value>monthly</value>
  </values>
  <valid_values>uniform,monthly,daily_interp</valid_values>
  <desc>Method for temporal downscaling of GCAM annual demand</desc>
</entry>

<entry id="gcam_env_flow_fraction" type="real" category="coupling">
  <values>
    <value>0.30</value>
  </values>
  <desc>Fraction of streamflow reserved for environmental flows (0-1)</desc>
</entry>

<entry id="gcam_spatial_disagg_method" type="char" category="coupling">
  <values>
    <value>area_weighted</value>
  </values>
  <valid_values>uniform,area_weighted,population_weighted,irrigation_weighted</valid_values>
  <desc>Method for spatial disaggregation of GCAM regional demand to MOSART grid</desc>
</entry>
```

### 8.2 Driver Namelist

**File:** `driver-mct/cime_config/namelist_definition_drv.xml`

**New Flags:**
```xml
<entry id="rof_c2_iac" type="logical">
  <values>
    <value>.false.</value>
  </values>
  <desc>ROF to IAC coupling flag</desc>
</entry>

<entry id="iac_c2_rof" type="logical">
  <values>
    <value>.false.</value>
  </values>
  <desc>IAC to ROF coupling flag</desc>
</entry>
```

---

## 9. Restart and History

### 9.1 Restart Variables

**MOSART Restart (RtmRestFile.F90):**
```fortran
! Add to restart file
call ncd_io('qdem_gcam', rtmCTL%qdem_gcam, flag='write', dim1name='gridcell', dim2name='sector')
call ncd_io('supply_gcam', rtmCTL%supply_gcam, flag='write', dim1name='gridcell', dim2name='sector')
call ncd_io('deficit_gcam', rtmCTL%deficit_gcam, flag='write', dim1name='gridcell', dim2name='sector')
call ncd_io('consump_frac', rtmCTL%consump_frac, flag='write', dim1name='gridcell')
```

**Driver Restart (seq_rest_mod.F90):**
```fortran
! Add MOSART accumulation state
call ncd_io('r2zacc_rx', r2zacc_rx, flag='write')
call ncd_io('r2zacc_rx_cnt', r2zacc_rx_cnt, flag='write')
```

### 9.2 History Output

**New MOSART History Fields (RtmHistFlds.F90):**
```fortran
! GCAM demand
call RtmHistAddfld('QDEM_GCAM_IRRIG', units='m3/s', avgflag='A', &
                   longname='GCAM irrigation water demand')
call RtmHistAddfld('QDEM_GCAM_INDUST', units='m3/s', avgflag='A', &
                   longname='GCAM industrial water demand')
call RtmHistAddfld('QDEM_GCAM_MUNIC', units='m3/s', avgflag='A', &
                   longname='GCAM municipal water demand')
call RtmHistAddfld('QDEM_GCAM_ENERGY', units='m3/s', avgflag='A', &
                   longname='GCAM energy/cooling water demand')

! Supply and deficit
call RtmHistAddfld('SUPPLY_GCAM_TOTAL', units='m3/s', avgflag='A', &
                   longname='Total GCAM water supply')
call RtmHistAddfld('DEFICIT_GCAM_TOTAL', units='m3/s', avgflag='A', &
                   longname='Total GCAM water deficit')

! Water availability sent to GCAM
call RtmHistAddfld('WATER_AVAIL_GCAM', units='m3', avgflag='I', &
                   longname='Annual water availability sent to GCAM')
```

---

## 10. Testing and Validation

### 10.1 Unit Tests

**Test 1: Field Index Consistency**
- Verify all coupling indices are properly set
- Check that optional fields are handled correctly

**Test 2: Unit Conversions**
- Verify m³/year ↔ m³/s conversions
- Check area normalizations

**Test 3: Temporal Interpolation**
- Test monthly downscaling at month boundaries
- Verify conservation of annual totals

### 10.2 Integration Tests

**Test 1: Water Balance**
```
Annual water balance closure:
  Input runoff = Output to ocean + Evaporation + Storage change + GCAM withdrawals
```

**Test 2: Demand-Supply Consistency**
```
For each sector:
  Supply + Deficit = Demand
  0 ≤ Supply ≤ Demand
  Supply ≤ Available water
```

**Test 3: Restart Reproducibility**
- Run for 2 years continuously
- Run for 1 year, restart, run 1 more year
- Verify identical results

### 10.3 Scientific Validation

**Comparison Datasets:**
1. FAO AQUASTAT water withdrawal data
2. USGS water use estimates (for US regions)
3. Global reservoir operations data (GRanD)
4. Environmental flow assessments

**Validation Metrics:**
1. Regional water demand magnitudes
2. Seasonal distribution of withdrawals
3. Water stress indicators (demand/availability ratio)
4. Reservoir storage dynamics

---

## 11. Phased Implementation Plan

### Phase 1: Infrastructure (Months 1-2)
- [ ] Create `rof2iacMod.F90` and `iac2rofMod.F90`
- [ ] Modify `RunoffMod.F90` to add GCAM demand fields
- [ ] Update `rof_cpl_indices.F90` with new indices
- [ ] Add field definitions to `seq_flds_mod.F90`
- [ ] Unit test each module independently

### Phase 2: Driver Integration (Months 3-4)
- [ ] Modify `prep_iac_mod.F90` for MOSART accumulation
- [ ] Update `cime_comp_mod.F90` coupling sequence
- [ ] Implement grid mapping (offline weight generation)
- [ ] Add to `seq_maps.rc` configuration
- [ ] Integration test: compile and run with stub IAC

### Phase 3: MOSART-Side Implementation (Months 5-6)
- [ ] Implement water availability calculation
- [ ] Modify `rof_comp_mct.F90` export/import
- [ ] Add temporal downscaling (monthly)
- [ ] Implement WRM integration
- [ ] Test with prescribed GCAM demand

### Phase 4: GCAM-Side Implementation (Months 7-8)
- [ ] Modify GCAM water module to read MOSART availability
- [ ] Update GCAM water demand calculation
- [ ] Return demand by sector to MOSART
- [ ] Two-way coupling test

### Phase 5: Validation and Tuning (Months 9-12)
- [ ] Scientific validation against observations
- [ ] Parameter tuning (env flow fraction, downscaling)
- [ ] Performance optimization
- [ ] Documentation and examples
- [ ] User guide and tutorials

---

## 12. Performance Considerations

### 12.1 Computational Cost

**Added Operations:**
1. MOSART accumulation every timestep: **O(ngrid)** → negligible
2. Annual averaging: **O(ngrid)** → once per year
3. Grid mapping: **O(nnz)** where nnz = nonzeros in weight matrix → once per year
4. Temporal downscaling: **O(ngrid)** → once per month

**Estimated Overhead:** < 1% of total runtime (dominated by GCAM execution)

### 12.2 Memory Footprint

**New Arrays:**
- `qdem_gcam(ngrid, 4)` - 4 sectors × 8 bytes × ngrid
- `supply_gcam(ngrid, 4)`
- `deficit_gcam(ngrid, 4)`
- `r2zacc_rx` - MCT aVect with 5 fields × ngrid

**For global 1/8° MOSART grid (~300,000 cells):**
- Memory increase: ~50 MB (negligible)

### 12.3 I/O Considerations

**Restart Files:**
- Additional ~50 MB per restart file
- Negligible compared to full model restart

**History Output:**
- Optional diagnostic fields
- Users can configure based on needs

---

## 13. Future Enhancements

### 13.1 Short-Term (1-2 Years)
1. **Groundwater Coupling:**
   - Add groundwater availability from ELM to GCAM
   - GCAM provides groundwater demand back to ELM

2. **Sub-Annual Coupling:**
   - Seasonal or monthly coupling instead of annual
   - Capture monsoon, snowmelt, and irrigation seasons

3. **Water Quality:**
   - Extend to nutrients, salinity, temperature
   - Link to BGC tracers in MOSART

### 13.2 Long-Term (3-5 Years)
1. **Multi-Model Ensemble:**
   - Couple multiple water models (MOSART, MRTM, etc.)
   - Provide uncertainty estimates to GCAM

2. **Floodplain Inundation:**
   - Use MOSART inundation model outputs
   - GCAM accounts for flood risk in economic decisions

3. **Managed Aquifer Recharge:**
   - GCAM provides recharge infrastructure decisions
   - MOSART/ELM simulate managed recharge

4. **Water Rights and Institutions:**
   - Encode water rights in GCAM-MOSART coupling
   - Respect priority systems, trans-boundary agreements

---

## 14. Summary of File Changes

### New Files:
1. `/components/mosart/src/cpl/rof2iacMod.F90` (~200 lines)
2. `/components/mosart/src/cpl/iac2rofMod.F90` (~150 lines)

### Modified Files:
1. `/components/mosart/src/riverroute/RunoffMod.F90` (~50 lines added)
2. `/components/mosart/src/cpl/rof_cpl_indices.F90` (~60 lines added)
3. `/components/mosart/src/cpl/rof_comp_mct.F90` (~100 lines added)
4. `/driver-mct/main/prep_iac_mod.F90` (~150 lines added)
5. `/driver-mct/shr/seq_flds_mod.F90` (~20 lines added)
6. `/driver-mct/main/cime_comp_mod.F90` (~80 lines added)
7. `/components/mosart/bld/namelist_files/namelist_definition_mosart.xml` (4 parameters)
8. `/driver-mct/cime_config/namelist_definition_drv.xml` (2 parameters)

**Total New Code:** ~800-1000 lines

---

## 15. References

1. **MOSART Documentation:**
   - Li, H., et al. (2013). "A physically based runoff routing model for land surface and Earth system models." J. Hydrometeor.

2. **GCAM Water Module:**
   - Kim, S., et al. (2016). "Balancing global water availability and use at basin scale in an integrated assessment model." Climatic Change.

3. **E3SM Coupling:**
   - Golaz, J., et al. (2019). "The DOE E3SM coupled model version 1." JAMES.

4. **Environmental Flows:**
   - Pastor, A., et al. (2014). "Accounting for environmental flow requirements in global water assessments." Hydrol. Earth Syst. Sci.

---

**End of Design Document**
