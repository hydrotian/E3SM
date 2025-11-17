# E3SM-GCAM Coupling Implementation Documentation

## Executive Summary

This document provides comprehensive documentation of the coupling scheme between E3SM (Energy Exascale Earth System Model) and GCAM (Global Change Analysis Model) implemented on the `tizhou/gcam/waterdemand/` branch. The coupling enables bidirectional information exchange between climate and land surface components (ELM, EAM) and the integrated assessment component (GCAM) for land use, carbon emissions, and economic modeling.

**Date:** 2025-11-17
**Branch:** `tizhou/gcam/waterdemand/`
**Base Model:** E3SM Version 3

---

## 1. Architecture Overview

### 1.1 Components Involved

The GCAM coupling integrates three main component types:

1. **IAC (Integrated Assessment Component)** - GCAM model wrapped as an E3SM component
2. **ELM (Energy Exascale Land Model)** - Land surface and vegetation model
3. **EAM (Energy Exascale Atmosphere Model)** - Atmospheric model

### 1.2 Coupling Framework

The coupling uses the **Model Coupling Toolkit (MCT)** framework with annual synchronization for land use and monthly interpolation for atmospheric CO2 emissions.

**Key Framework Components:**
- MCT (Model Coupling Toolkit) for data exchange
- Attribute vectors (mct_aVect) for field storage
- Conservative grid mapping between component grids
- Annual and monthly temporal resolution coupling

---

## 2. Component Locations and File Structure

### 2.1 GCAM Component

**Root Directory:** `/home/user/E3SM/components/gcam/`

**Structure:**
```
components/gcam/
├── src/                    # GCAM source code (submodule: E3SM-Project/giac)
├── bld/                    # Build configuration
│   └── namelist_files/
│       ├── namelist_definition_gcam.xml
│       └── namelist_defaults_gcam.xml
├── cime_config/            # CIME integration
│   ├── config_component.xml
│   ├── config_compsets.xml
│   └── user_nl_gcam
└── doc/                    # Documentation
    ├── NOTES.design
    └── NOTES.e3sm
```

### 2.2 Coupling Infrastructure

**Driver-Level Integration:**
- `/home/user/E3SM/driver-mct/main/cime_comp_mod.F90` - Main driver with IAC integration
- `/home/user/E3SM/driver-mct/main/prep_iac_mod.F90` - IAC data preparation and accumulation
- `/home/user/E3SM/driver-mct/shr/seq_flds_mod.F90` - Coupling field definitions
- `/home/user/E3SM/driver-mct/main/seq_rest_mod.F90` - Restart I/O
- `/home/user/E3SM/driver-mct/main/seq_hist_mod.F90` - History output

**Component-Level Coupling:**
- `/home/user/E3SM/components/elm/src/main/lnd2iacMod.F90` - ELM → IAC data export
- `/home/user/E3SM/components/elm/src/main/iac2lndMod.F90` - IAC → ELM data import
- `/home/user/E3SM/components/eam/src/chemistry/utils/iac_coupled_fields.F90` - IAC → EAM CO2 coupling

**MCT Interface Stub:**
- `/home/user/E3SM/components/stub_comps/siac/src/iac_comp_mct.F90` - IAC MCT interface stub

---

## 3. Coupling Data Exchange

### 3.1 ELM → GCAM (Annual Terrestrial Productivity)

**Purpose:** Provide GCAM with annual average terrestrial productivity data for economic land use decisions.

**Implementation Location:** `components/elm/src/main/lnd2iacMod.F90`

**Data Structure:**
```fortran
type lnd2iac_type
   real(r8), pointer :: hr(:,:)      ! Heterotrophic respiration per PFT (gC/m²/s)
   real(r8), pointer :: npp(:,:)     ! Net primary production per PFT (gC/m²/s)
   real(r8), pointer :: pftwgt(:,:)  ! PFT weight relative to gridcell (fraction)
end type
```

**Dimensions:**
- Grid cells × (numpft + 1)
- 17 PFTs total (0 = bare ground, 1-16 = vegetation types)
- Results in 51 coupled fields (17 PFTs × 3 variables)

**Coupling Fields in MCT (seq_flds_mod.F90:2587-2605):**
```fortran
Sl_hr_pft[0:16]      ! Total heterotrophic respiration per PFT
Sl_npp_pft[0:16]     ! Net primary production per PFT
Sl_pftwgt_pft[0:16]  ! PFT weight relative to gridcell
```

**Data Flow:**
1. ELM calculates NPP and respiration at every timestep (~30 min)
2. Driver accumulates these values via `prep_iac_accum()` (prep_iac_mod.F90:197-231)
3. At year boundary, driver averages via `prep_iac_accum_avg()` (prep_iac_mod.F90:235-268)
4. Averaged data mapped to IAC grid via `prep_iac_calc_l2x_zx()` (prep_iac_mod.F90:418-452)
5. Data merged and sent to GCAM

**Key Functions:**
- `update_lnd2iac()` (lnd2iacMod.F90:73-118) - Populates coupling arrays from ELM
- `prep_iac_accum()` - Annual accumulation of land data
- `prep_iac_accum_avg()` - Finalizes annual averaging

---

### 3.2 GCAM → ELM (Annual Land Use/Cover Updates)

**Purpose:** Apply GCAM-determined land use changes back to ELM for updated vegetation distribution and harvest rates.

**Implementation Location:** `components/elm/src/main/iac2lndMod.F90`

**Data Structure:**
```fortran
type iac2lnd_type
   real(r8), pointer :: frac_pft(:,:)      ! PFT fraction of gridcell (current year)
   real(r8), pointer :: frac_pft_prev(:,:) ! PFT fraction of gridcell (previous year)
   real(r8), pointer :: harvest_frac(:,:)  ! Harvest fraction by category
end type
```

**Dimensions:**
- Grid cells × (numpft + 1) for PFT fractions
- Grid cells × numharvest for harvest (5 categories)

**Coupling Fields in MCT (seq_flds_mod.F90:2617-2640):**
```fortran
Sz_pct_pft[0:16]      ! Percent PFT of vegetated land unit (current year)
Sz_pct_pft_prev[0:16] ! Percent PFT of vegetated land unit (previous year)
Sz_harvest_frac[0:4]  ! Harvest fraction by category
```

**Harvest Categories (numharvest = 5):**
0. Wood harvest
1. Crop harvest
2. Livestock grazing
3. Biofuel harvest
4. Other harvest

**Data Flow:**
1. GCAM runs annually and determines optimal land use
2. GCAM outputs PFT fractions and harvest rates
3. Driver maps IAC grid data to land grid via `prep_lnd_calc_z2x_lx()`
4. `update_iac2lnd()` (iac2lndMod.F90:104-267) applies changes to ELM
5. ELM performs temporal interpolation between current and previous year
6. Harvest rates applied to vegetation carbon pools

**Key Functions:**
- `update_iac2lnd()` - Converts (ngrid, pft) to patch-level data with temporal interpolation
- Time interpolation: `wt1 = 1.0 - get_curr_yearfrac()` for smooth transitions

**Unit Conversions:**
- IAC provides fraction of actual gridcell
- ELM requires fraction of land (must divide by `ldomain%frac` and `ldomain%mask`)
- Harvest rates normalized by vegetated column weight

---

### 3.3 GCAM → EAM (Monthly CO2 Emissions)

**Purpose:** Provide atmospheric model with anthropogenic CO2 emissions from GCAM for climate forcing.

**Implementation Location:** `components/eam/src/chemistry/utils/iac_coupled_fields.F90`

**Data Structure:**
```fortran
type iac_vertical_emiss_t
   integer  :: lchnk                        ! chunk index
   integer  :: ncol                         ! number of active columns
   real(r8), allocatable :: fco2_low_height(:)   ! CO2 at low altitude (kg/m²/s)
   real(r8), allocatable :: fco2_high_height(:)  ! CO2 at high altitude (kg/m²/s)
end type
```

**Coupling Fields in MCT (seq_flds_mod.F90:2660-2676):**
```fortran
! For each month (1-12):
Fazz_co2sfc_mon[1:12]    ! Surface flux of CO2 (moles/m²/s)
Fazz_co2airlo_mon[1:12]  ! Low altitude CO2 flux (moles/m²/s)
Fazz_co2airhi_mon[1:12]  ! High altitude CO2 flux (moles/m²/s)
```

**Total Coupled Fields:** 36 (12 months × 3 altitude levels)

**Vertical Distribution:**
- **Surface (co2sfc):** Ground level emissions
- **Low altitude (co2airlo):** Distributed below HIGH_LAYER (11 km)
- **High altitude (co2airhi):** Injected at HIGH_LAYER = 11,000 m

**Temporal Interpolation:**
- GCAM provides monthly values (mid-month)
- EAM performs daily interpolation via `iac_coupled_timeinterp()` (iac_coupled_fields.F90:232-359)
- Linear interpolation between monthly mid-points (day 15, 45, 74, ...)
- Special handling for Dec 16 - Jan 15 wraparound

**Data Flow:**
1. GCAM runs annually and calculates monthly CO2 emissions
2. Data distributed to atmosphere grid
3. `iac_coupled_fields_init()` allocates physics buffer storage
4. `iac_coupled_fields_adv()` (lines 139-230) updates physics buffer each timestep:
   - Finds vertical level for HIGH_LAYER (11 km)
   - Assigns high altitude flux to that level
   - Distributes low altitude flux among levels below HIGH_LAYER
5. Time interpolation provides daily values from monthly data

**Vertical Distribution Options (iac_low_height_option):**
- Option 1 (default): Equal distribution among all levels below HIGH_LAYER
- Additional options can be added for more sophisticated distributions

**Key Functions:**
- `iac_coupled_fields_register()` - Register CO2 field in physics buffer
- `iac_coupled_fields_init()` - Initialize data structures
- `iac_coupled_fields_adv()` - Advance CO2 fields each timestep
- `iac_coupled_timeinterp()` - Compute time interpolation fraction

---

## 4. Coupling Timing and Frequency

### 4.1 Temporal Resolution

**Configuration Location:** `components/gcam/cime_config/config_compsets.xml:43-80`

```xml
<entry id="NCPL_BASE_PERIOD">year</entry>
<entry id="IAC_NCPL">1</entry>        <!-- IAC runs once per year -->
<entry id="LND_NCPL">17520</entry>    <!-- ~30-minute timesteps -->
<entry id="ATM_NCPL">17520</entry>    <!-- ~30-minute timesteps -->
```

**Timing Summary:**
- **IAC/GCAM execution:** Once per year
- **Land accumulation:** Every land timestep (~30 min = 17,520 steps/year)
- **Annual averaging:** At year boundary (Day 1, 00:30 UTC)
- **CO2 interpolation:** Monthly to daily

### 4.2 Coupling Sequence

**Annual Cycle (from driver-mct/main/cime_comp_mod.F90):**

```
├─ Throughout Year (Days 1-365):
│  ├─ Every land timestep (~30 min):
│  │  ├─ ELM produces hr, npp, pftwgt
│  │  └─ prep_iac_accum(): Accumulate values in l2zacc_lx
│  │
│  └─ Every atmosphere timestep:
│     └─ EAM uses interpolated monthly CO2 values
│
├─ Year Boundary (Day 1, TOD=1800 seconds = 00:30 UTC):
│  ├─ cime_run_iac_setup_send() (cime_comp_mod.F90:4033-4084):
│  │  ├─ prep_iac_accum_avg(): Finalize annual averaging
│  │  ├─ prep_iac_calc_l2x_zx(): Map accumulated data to IAC grid
│  │  ├─ prep_iac_mrg(): Merge inputs
│  │  └─ component_exch(): Send to IAC (x2z flow)
│  │
│  ├─ IAC/GCAM runs:
│  │  ├─ Receives: Annual average NPP, respiration, PFT weights
│  │  ├─ Computes: Economic optimization, land use decisions, CO2 emissions
│  │  └─ Sends: PFT fractions, harvest rates, monthly CO2 fluxes
│  │
│  ├─ component_exch(): Receive from IAC (z2x flow)
│  │
│  ├─ cime_run_iac_recv_post() (cime_comp_mod.F90:4087-4121):
│  │  └─ Post-processing and diagnostics
│  │
│  ├─ For ELM:
│  │  ├─ prep_lnd_calc_z2x_lx(): Map Z2X fields to land grid
│  │  └─ ELM updates: frac_pft, harvest_frac
│  │
│  └─ For EAM:
│     ├─ prep_atm_calc_z2x_ax(): Map Z2X fluxes to atm grid
│     └─ EAM stores: 12 monthly CO2 flux fields × 3 altitudes
│
└─ Post-processing:
   ├─ prep_iac_zero_max(): Reset accumulators for next year
   └─ Restart/History I/O
```

### 4.3 Alarm Configuration

**IAC Alarms (cime_comp_mod.F90:277-498):**
- `iacrun_alarm`: Triggers IAC execution (annually)
- `iacrun_avg_alarm`: Triggers averaging before IAC run

---

## 5. Grid Mapping and Domain

### 5.1 Grid Configuration

**Components use different grids:**
- **Land (ELM):** Land grid (lnd_gnam)
- **Atmosphere (EAM):** Atmosphere grid (atm_gnam)
- **IAC (GCAM):** Economic regions grid (iac_gnam)

**Mapping Required:** Conservative remapping between:
- Land grid ↔ IAC grid
- Atmosphere grid ↔ IAC grid

### 5.2 Mapper Infrastructure

**Location:** `driver-mct/main/prep_iac_mod.F90:56-58`

```fortran
type(seq_map), pointer :: mapper_Sl2z  ! Land to IAC mapper
type(seq_map), pointer :: mapper_Sa2z  ! Atmosphere to IAC mapper
```

**Initialization (prep_iac_mod.F90:165-168):**
```fortran
call seq_map_init_rcfile(mapper_Sl2z, lnd(1), iac(1), &
     'seq_maps.rc','lnd2iac_smapname:','lnd2iac_smaptype:',samegrid_lz, &
     string='mapper_Sl2z initialization',esmf_map=esmf_map_flag)
```

**Mapping Configuration File:** `seq_maps.rc`
- `lnd2iac_smapname`: Mapping weight file
- `lnd2iac_smaptype`: Mapping algorithm (conservative, bilinear, etc.)
- `atm2iac_smapname`: Atmosphere to IAC mapping
- `atm2iac_smaptype`: Mapping type

**Usage:**
```fortran
! Map land data to IAC grid
call seq_map_map(mapper_Sl2z, l2zacc_lx(eli), l2x_zx(ezi), &
                 fldlist=seq_flds_l2x_states, norm=.true.)
```

---

## 6. Data Structures and Attribute Vectors

### 6.1 MCT Attribute Vectors

**Location:** `driver-mct/main/prep_iac_mod.F90:61-69`

```fortran
! Export to IAC (on IAC grid, coupler PEs)
type(mct_aVect), pointer :: l2x_zx(:)

! Accumulation (on land grid, coupler PEs)
type(mct_aVect), pointer :: l2zacc_lx(:)   ! Accumulated land export
integer, target          :: l2zacc_lx_cnt  ! Accumulation counter

! Maximum monthly values (currently unused)
type(mct_aVect), pointer :: l2zmax_lx(:)
```

**Dimensioning:**
- Arrays dimensioned by number of instances: `num_inst_lnd`, `num_inst_iac`
- Each aVect contains all coupling fields for grid points on local PE

### 6.2 Field Definitions

**Location:** `driver-mct/shr/seq_flds_mod.F90:233-236`

```fortran
character(CXX) :: seq_flds_x2z_states   ! States sent from coupler to IAC
character(CXX) :: seq_flds_z2x_states   ! States received from IAC (land use)
character(CXX) :: seq_flds_z2x_fluxes   ! Fluxes received from IAC (CO2)
character(CXX) :: seq_flds_x2z_fluxes   ! Fluxes sent to IAC
```

**Field Lists Built Dynamically:**
- Based on which components are active
- Based on namelist configuration (e.g., PFT count)

---

## 7. Restart and History

### 7.1 Restart Files

**IAC Restart Data (seq_rest_mod.F90):**
- `l2zacc_lx`: Accumulated land data (for mid-year restart)
- `l2zacc_lx_cnt`: Accumulation count
- `z2x_zx`: IAC export data

**GCAM-Specific Restart Files:**
- `gcam2glm_restart.r.YYYY.nc`: GCAM-to-GLM interface restart
- `output.glm.restart.state.YYYY.nc`: GLM (land use module) state
- Written by `iac_rpointer_write()` (iac2lndMod.F90:271-319)

**Restart Pointers:**
- `rpointer.gcam2glm`: Points to gcam2glm restart file
- `rpointer.glm`: Points to GLM restart file

### 7.2 History Output

**IAC History Fields (seq_hist_mod.F90):**
- `z2x_zx_avg`: Time-averaged IAC outputs
- Both states (PFT fractions) and fluxes (CO2)

**Component-Level History:**
- ELM writes lnd2iac fields
- EAM writes CO2 fields from IAC
- GCAM component writes internal diagnostics

---

## 8. Compsets and Configuration

### 8.1 Available Compsets

**Location:** `components/gcam/cime_config/config_compsets.xml:11-41`

| Compset | Description | Components |
|---------|-------------|------------|
| **Z** | IAC stub only | Data models + stub IAC |
| **ZLND** | IAC with data atmosphere | Data ATM + ELM + GCAM |
| **ZATM** | Fully coupled IAC | EAM + ELM + GCAM |
| **SSP245_ZATM** | SSP2-4.5 scenario | EAM + ELM + GCAM (RCP 4.5 forcing) |
| **SSP245_ZATM_BGC** | SSP2-4.5 with BGC | EAM + ELM + GCAM + BGC |
| **SSP370_ZATM_BGC** | SSP3-7.0 high emissions | EAM + ELM + GCAM + BGC |

### 8.2 Configuration Variables

**Component Configuration (config_component.xml):**
```xml
<entry id="COMP_IAC" valid_values="gcam">
  <desc>IAC component: GCAM integrated assessment model</desc>
</entry>

<entry id="USE_EHC" value="TRUE">
  <desc>Enable CXX library for GCAM compilation</desc>
</entry>
```

**Runtime Configuration:**
- Namelist files in `components/gcam/bld/namelist_files/`
- User customization via `user_nl_gcam`

---

## 9. Implementation Details

### 9.1 Accumulation Algorithm

**Purpose:** Convert 30-minute land model output to annual average for GCAM.

**Algorithm (prep_iac_mod.F90:197-268):**

```fortran
subroutine prep_iac_accum()
   ! Called every land timestep
   if (l2zacc_lx_cnt == 0) then
      ! First timestep: copy
      call mct_avect_copy(l2x_lx, l2zacc_lx(eli))
   else
      ! Subsequent timesteps: accumulate (sum)
      call mct_avect_accum(l2x_lx, l2zacc_lx(eli))
   endif
   l2zacc_lx_cnt = l2zacc_lx_cnt + 1
end subroutine

subroutine prep_iac_accum_avg()
   ! Called at year boundary
   if (l2zacc_lx_cnt > 1) then
      ! Average = sum / count
      call mct_avect_avg(l2zacc_lx(eli), l2zacc_lx_cnt)
   endif
   l2zacc_lx_cnt = 0  ! Reset for next year
end subroutine
```

**Notes:**
- Original design included monthly max option (commented out)
- Current implementation uses annual average only
- `l2zmax_lx` allocated but unused

### 9.2 PFT Indexing

**Convention:**
- PFT 0 = bare ground
- PFT 1-16 = vegetation types
- Arrays allocated `0:numpft` where `numpft = 16`

**In ELM:**
```fortran
do p = begp, endp
   g = veg_pp%gridcell(p)
   pft = veg_pp%itype(p)  ! PFT type (0-16)

   this%hr(g,pft) = col_cf%hr(c)
   this%npp(g,pft) = veg_cf%npp(p)
   this%pftwgt(g,pft) = veg_pp%wtgcell(p) * ldomain%frac(g) * ldomain%mask(g)
end do
```

**In GCAM:**
- Receives 17 separate fields per variable
- Economic optimization by PFT type
- Returns 17 separate fraction fields

### 9.3 Unit Conversions

**Land → IAC:**
- NPP: gC/m²/s (no conversion)
- Respiration: gC/m²/s (no conversion)
- PFT weight: fraction of gridcell (0-1)

**IAC → Land:**
- PFT fractions: percent (%) → fraction by dividing by 100
- Harvest: fraction of vegetated land unit

**IAC → Atmosphere:**
- CO2 flux: moles/m²/s
- Vertical distribution: surface, low (<11 km), high (11 km)

---

## 10. Key Design Patterns

### 10.1 Type-Bound Procedures

**ELM uses object-oriented design:**

```fortran
type, public :: lnd2iac_type
   real(r8), pointer :: hr(:,:)
   real(r8), pointer :: npp(:,:)
   real(r8), pointer :: pftwgt(:,:)
contains
   procedure, public :: Init
   procedure, public :: update_lnd2iac
end type

! Usage:
type(lnd2iac_type) :: lnd2iac
call lnd2iac%Init(bounds)
call lnd2iac%update_lnd2iac(bounds)
```

### 10.2 Grid-to-Patch Conversion

**IAC → ELM requires converting (grid, PFT) to patch:**

```fortran
! IAC data: iac2lnd%frac_pft(gridcell, pft)
! Target: veg_pp%wtgcell_iac(patch)

do p = begp, endp
   g = veg_pp%gridcell(p)
   pft = veg_pp%itype(p)

   ! Extract from 2D array and convert units
   temp = this%frac_pft(g, pft) / (ldomain%frac(g) * ldomain%mask(g))
   temp_prev = this%frac_pft_prev(g, pft) / (ldomain%frac(g) * ldomain%mask(g))

   ! Temporal interpolation
   wt1 = 1.0 - get_curr_yearfrac()
   veg_pp%wtgcell_iac(p) = temp + wt1 * (temp_prev - temp)
end do
```

### 10.3 Temporal Interpolation

**Monthly CO2 to Daily:**

```fortran
! Mid-month days: [15, 45, 74, 105, ..., 349]
! Current day in year: cyclic_curr_model_time

! Find bracketing months
lower_bound = findplb(mid_mon_num_days, tot_mon_in_a_year, cyclic_curr_model_time)
upper_bound = lower_bound + 1

! Special case: Dec 16 - Jan 15
if (lower_bound == 12) then
   upper_bound = 1
   day_at_upper_bnd = mid_mon_num_days(1) + num_days_in_current_year
endif

! Linear interpolation fraction
time_interp_frac = (model_time - day_at_lower_bnd) / (day_at_upper_bnd - day_at_lower_bnd)
```

---

## 11. Error Handling and Diagnostics

### 11.1 Validation Checks

**Time Interpolation Bounds:**
```fortran
if (time_interp_frac > 1.0 .or. time_interp_frac < 0.0) then
   call endrun('ERROR: time_interp_frac out of bounds')
endif
```

**Mid-Month Boundary Check:**
```fortran
if (day == mid_mon_num_days(month) + 1 .and. secs == 0) then
   if (time_interp_frac > tiny(time_interp_frac)) then
      call endrun('ERROR: time_interp_frac should be 0.0 at mid month days')
   endif
endif
```

### 11.2 Diagnostic Output

**MCT AVect Info:**
```fortran
call mct_avect_info(4, l2zacc_lx(eli), istr='TRS l2zacc')
```

**Optional NetCDF Diagnostic File:**
```fortran
! iac2lndMod.F90:229-264 (disabled by default)
if (.false.) then
   write(hfile, '(a)') './iac2lnd_update.nc'
   ! Write frac_pft, frac_pft_prev, harvest_frac to NetCDF
endif
```

---

## 12. Performance Considerations

### 12.1 Communication Patterns

**Annual Coupling Minimizes Overhead:**
- IAC runs once per year → minimal synchronization cost
- Land accumulation is local (no communication)
- Mapping operations only at year boundary

**Data Volume:**
- L2Z: ~51 fields × grid size
- Z2L: ~39 fields × grid size
- Z2A: 36 fields × grid size

### 12.2 Memory Footprint

**Persistent Attribute Vectors:**
- `l2x_zx(num_inst_lnd)`: IAC grid size
- `l2zacc_lx(num_inst_lnd)`: Land grid size
- `l2zmax_lx(num_inst_lnd)`: Land grid size (unused)

**EAM CO2 Storage:**
- `iac_vertical_emiss(begchunk:endchunk)`
- Per chunk: `fco2_low_height(pcols)` + `fco2_high_height(pcols)`

### 12.3 Computational Cost

**IAC Execution:**
- Runs once per year
- GCAM economic model solves for equilibrium
- Most expensive operation in coupling

**Accumulation:**
- O(1) cost per timestep (vector addition)
- 17,520 accumulations per year

**Mapping:**
- Conservative remapping once per year
- Cost proportional to overlap matrix size

---

## 13. Limitations and Assumptions

### 13.1 Current Limitations

1. **Spatial Resolution:**
   - IAC grid typically much coarser than land/atmosphere grids
   - Conservative mapping may smooth fine-scale heterogeneity

2. **Temporal Resolution:**
   - Annual coupling may miss sub-annual dynamics
   - No feedback from intra-annual climate variability to GCAM

3. **One-Way Heat/Moisture:**
   - Land use changes affect carbon/vegetation only
   - No feedback to energy/water budgets during year

4. **PFT Representation:**
   - Fixed to 17 PFT types
   - GCAM and ELM must use consistent PFT definitions

### 13.2 Key Assumptions

1. **Annual Averaging:**
   - Annual average NPP is sufficient for GCAM
   - Monthly CO2 is sufficient for EAM chemistry

2. **Grid Consistency:**
   - Land fraction (`ldomain%frac`) consistent across grids
   - Mapping preserves mass/energy

3. **Temporal Synchronization:**
   - Year boundary aligns across all components
   - IAC runs at beginning of year (TOD=1800)

4. **Restart Consistency:**
   - GCAM internal state captured in restart files
   - GLM land use state synchronized with ELM

---

## 14. References and Related Documentation

### 14.1 Source Documentation

- `components/gcam/doc/NOTES.design` - GCAM design notes
- `components/gcam/doc/NOTES.e3sm` - E3SM integration notes
- `driver-mct/main/cime_comp_mod.F90` - Inline comments on IAC integration

### 14.2 External References

- **GCAM Documentation:** http://jgcri.github.io/gcam-doc/
- **E3SM Documentation:** https://e3sm.org/model/
- **MCT User's Guide:** http://www.mcs.anl.gov/research/projects/mct/
- **CIME Documentation:** http://esmci.github.io/cime/

---

## 15. Glossary

| Term | Definition |
|------|------------|
| **aVect** | MCT attribute vector - container for coupling fields |
| **EAM** | Energy Exascale Atmosphere Model |
| **ELM** | Energy Exascale Land Model |
| **GCAM** | Global Change Analysis Model |
| **GLM** | GCAM Land Model (land use component) |
| **IAC** | Integrated Assessment Component |
| **MCT** | Model Coupling Toolkit |
| **PFT** | Plant Functional Type |
| **TOD** | Time of Day (seconds past midnight UTC) |
| **Mapper** | Grid mapping/regridding object |

---

**End of Documentation**
