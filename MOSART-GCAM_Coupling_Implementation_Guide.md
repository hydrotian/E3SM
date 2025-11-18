# MOSART-GCAM Water Management Coupling Implementation Guide

This document provides detailed instructions for implementing the GCAM-side coupling infrastructure to complement the E3SM-side implementation that has been completed in this repository.

## Overview

The E3SM repository now contains complete infrastructure for MOSART-GCAM water management coupling. To enable functional bidirectional water exchange, corresponding changes are required in:

1. **giac** repository (E3SM-Project/giac)
2. **gcam-core** repository (JGCRI/gcam-core)

## Base Code Versions

This implementation guide is based on:
- **giac**: https://github.com/E3SM-Project/giac/tree/0ffce23d17c01c0f4e598891243abd4bc37fe900
- **gcam-core**: https://github.com/JGCRI/gcam-core/tree/2901ad96aadf0a2ba272ef7718419dd8453af022

---

# Part 1: Changes Required in `giac` Repository

The IAC component wrapper (`giac`) serves as the bridge between E3SM's MCT coupling infrastructure and GCAM's internal representation.

## 1.1 File: `src/iac_comp_mct.F90`

### Add Water Coupling Field Indices

Add these indices to the module-level declarations (similar to existing land coupling indices):

```fortran
! In the index declarations section
integer, public :: index_r2z_Sr_wr_avail      = 0  ! MOSART main channel water availability
integer, public :: index_r2z_Sr_wt_avail      = 0  ! MOSART tributary water availability
integer, public :: index_r2z_Sr_wtot_avail    = 0  ! MOSART total water availability
integer, public :: index_r2z_Sr_reservoir_stor = 0 ! MOSART reservoir storage
integer, public :: index_r2z_Sr_streamflow    = 0  ! MOSART streamflow

integer, public :: index_z2r_Sz_demand_irrig  = 0  ! IAC irrigation water demand
integer, public :: index_z2r_Sz_demand_indust = 0  ! IAC industrial water demand
integer, public :: index_z2r_Sz_demand_munic  = 0  ! IAC municipal water demand
integer, public :: index_z2r_Sz_demand_energy = 0  ! IAC energy water demand
integer, public :: index_z2r_Sz_demand_total  = 0  ! IAC total water demand
integer, public :: index_z2r_Sz_consump_frac  = 0  ! IAC consumptive use fraction
```

### Initialize Field Indices

In the index initialization routine (e.g., `iac_cpl_indices_set`):

```fortran
subroutine iac_cpl_indices_set()
  ! ... existing code ...

  ! Set MOSART->IAC water availability indices
  if (rof_present .and. rof_c2_iac) then
     call mct_aVect_indexRA(r2z, 'Sr_wr_avail', index_r2z_Sr_wr_avail, perrWith='quiet')
     call mct_aVect_indexRA(r2z, 'Sr_wt_avail', index_r2z_Sr_wt_avail, perrWith='quiet')
     call mct_aVect_indexRA(r2z, 'Sr_wtot_avail', index_r2z_Sr_wtot_avail, perrWith='quiet')
     call mct_aVect_indexRA(r2z, 'Sr_reservoir_stor', index_r2z_Sr_reservoir_stor, perrWith='quiet')
     call mct_aVect_indexRA(r2z, 'Sr_streamflow', index_r2z_Sr_streamflow, perrWith='quiet')
  endif

  ! Set IAC->MOSART water demand indices
  if (rof_present .and. iac_c2_rof) then
     call mct_aVect_indexRA(z2r, 'Sz_demand_irrig', index_z2r_Sz_demand_irrig, perrWith='quiet')
     call mct_aVect_indexRA(z2r, 'Sz_demand_indust', index_z2r_Sz_demand_indust, perrWith='quiet')
     call mct_aVect_indexRA(z2r, 'Sz_demand_munic', index_z2r_Sz_demand_munic, perrWith='quiet')
     call mct_aVect_indexRA(z2r, 'Sz_demand_energy', index_z2r_Sz_demand_energy, perrWith='quiet')
     call mct_aVect_indexRA(z2r, 'Sz_demand_total', index_z2r_Sz_demand_total, perrWith='quiet')
     call mct_aVect_indexRA(z2r, 'Sz_consump_frac', index_z2r_Sz_consump_frac, perrWith='quiet')
  endif

end subroutine iac_cpl_indices_set
```

### Import Water Availability from MOSART

Modify `iac_import_mct` subroutine to extract water availability data:

```fortran
subroutine iac_import_mct(r2z_z)
  ! ... existing code ...

  ! Import water availability from MOSART if coupling is enabled
  if (index_r2z_Sr_wtot_avail > 0) then

     ! Allocate arrays if not already done
     if (.not. allocated(water_avail_main)) then
        allocate(water_avail_main(lsize))
        allocate(water_avail_trib(lsize))
        allocate(water_avail_total(lsize))
        allocate(reservoir_storage(lsize))
        allocate(mean_streamflow(lsize))
     endif

     ! Extract water availability fields from coupling attribute vector
     do n = 1, lsize
        water_avail_main(n)  = r2z_z%rattr(index_r2z_Sr_wr_avail, n)
        water_avail_trib(n)  = r2z_z%rattr(index_r2z_Sr_wt_avail, n)
        water_avail_total(n) = r2z_z%rattr(index_r2z_Sr_wtot_avail, n)
        reservoir_storage(n) = r2z_z%rattr(index_r2z_Sr_reservoir_stor, n)
        mean_streamflow(n)   = r2z_z%rattr(index_r2z_Sr_streamflow, n)
     enddo

     ! Map from IAC grid to GCAM regions and pass to GCAM
     call map_water_availability_to_gcam(water_avail_main, water_avail_trib, &
                                         water_avail_total, reservoir_storage, &
                                         mean_streamflow)

     ! Call GCAM interface to set water availability
     call gcam_set_water_supply(water_avail_total, lsize)

  endif

end subroutine iac_import_mct
```

### Export Water Demand to MOSART

Modify `iac_export_mct` subroutine to provide water demand data:

```fortran
subroutine iac_export_mct(z2r_z)
  ! ... existing code ...

  ! Export water demand to MOSART if coupling is enabled
  if (index_z2r_Sz_demand_total > 0) then

     ! Get water demand from GCAM
     call gcam_get_water_demand(demand_irrig, demand_indust, demand_munic, &
                                demand_energy, demand_total, consump_frac, lsize)

     ! Map from GCAM regions to IAC grid (inverse of import mapping)
     call map_water_demand_from_gcam(demand_irrig, demand_indust, demand_munic, &
                                     demand_energy, demand_total, consump_frac)

     ! Export to coupling attribute vector
     do n = 1, lsize
        z2r_z%rattr(index_z2r_Sz_demand_irrig, n)  = demand_irrig(n)
        z2r_z%rattr(index_z2r_Sz_demand_indust, n) = demand_indust(n)
        z2r_z%rattr(index_z2r_Sz_demand_munic, n)  = demand_munic(n)
        z2r_z%rattr(index_z2r_Sz_demand_energy, n) = demand_energy(n)
        z2r_z%rattr(index_z2r_Sz_demand_total, n)  = demand_total(n)
        z2r_z%rattr(index_z2r_Sz_consump_frac, n)  = consump_frac(n)
     enddo

  endif

end subroutine iac_export_mct
```

## 1.2 Grid Mapping Implementation

Create a new module for handling water-specific grid mapping:

### File: `src/water_grid_mapping.F90`

```fortran
module water_grid_mapping

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private

  public :: init_water_grid_mapping
  public :: map_water_availability_to_gcam
  public :: map_water_demand_from_gcam

  ! Mapping weights and indices
  integer, allocatable :: mosart_to_gcam_map(:,:)  ! (n_mosart_cells, max_gcam_regions)
  real(r8), allocatable :: mosart_to_gcam_weights(:,:)
  integer, allocatable :: gcam_to_mosart_map(:,:)  ! (n_gcam_regions, max_mosart_cells)
  real(r8), allocatable :: gcam_to_mosart_weights(:,:)

contains

  subroutine init_water_grid_mapping(iac_grid, mosart_grid)
    !---------------------------------------------------------------
    ! Initialize mapping between MOSART river network and GCAM regions
    ! This typically involves:
    ! 1. Reading basin-to-region mapping file
    ! 2. Computing area-weighted aggregation factors
    ! 3. Handling overlapping basins/regions
    !---------------------------------------------------------------

    ! Read mapping weights from file or compute on-the-fly
    ! For MOSART->GCAM: aggregate river cells to economic regions
    ! For GCAM->MOSART: disaggregate regional demands to river cells

    ! Implementation depends on GCAM region definition
    ! May require external mapping file relating MOSART basins to GCAM regions

  end subroutine init_water_grid_mapping

  subroutine map_water_availability_to_gcam(wr_avail, wt_avail, wtot_avail, &
                                            reservoir, streamflow)
    !---------------------------------------------------------------
    ! Aggregate MOSART water availability to GCAM regions
    ! Use conservative area-weighted aggregation
    !---------------------------------------------------------------
    real(r8), intent(in) :: wr_avail(:), wt_avail(:), wtot_avail(:)
    real(r8), intent(in) :: reservoir(:), streamflow(:)

    integer :: n_mosart, n_gcam, i, j
    real(r8) :: regional_avail(n_gcam_regions)

    ! Aggregate using mapping weights
    regional_avail(:) = 0.0_r8
    do i = 1, n_mosart_cells
       do j = 1, mosart_to_gcam_map(i,0)  ! number of regions for this cell
          regional_avail(mosart_to_gcam_map(i,j)) = &
               regional_avail(mosart_to_gcam_map(i,j)) + &
               wtot_avail(i) * mosart_to_gcam_weights(i,j)
       enddo
    enddo

    ! Store for GCAM interface

  end subroutine map_water_availability_to_gcam

  subroutine map_water_demand_from_gcam(demand_irrig, demand_indust, &
                                        demand_munic, demand_energy, &
                                        demand_total, consump_frac)
    !---------------------------------------------------------------
    ! Disaggregate GCAM water demand to MOSART grid
    ! Use spatial downscaling based on basin characteristics
    !---------------------------------------------------------------
    real(r8), intent(out) :: demand_irrig(:), demand_indust(:)
    real(r8), intent(out) :: demand_munic(:), demand_energy(:)
    real(r8), intent(out) :: demand_total(:), consump_frac(:)

    ! Get regional demands from GCAM
    ! Disaggregate to MOSART cells using appropriate weights
    ! (e.g., cropland area for irrigation, population for municipal)

  end subroutine map_water_demand_from_gcam

end module water_grid_mapping
```

## 1.3 GCAM Interface Module

Create C++ callable Fortran interface for GCAM:

### File: `src/gcam_water_interface.F90`

```fortran
module gcam_water_interface

  use iso_c_binding
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  ! Water data storage for GCAM exchange
  real(r8), allocatable, save :: gcam_water_supply(:)
  real(r8), allocatable, save :: gcam_demand_irrig(:)
  real(r8), allocatable, save :: gcam_demand_indust(:)
  real(r8), allocatable, save :: gcam_demand_munic(:)
  real(r8), allocatable, save :: gcam_demand_energy(:)
  real(r8), allocatable, save :: gcam_demand_total(:)
  real(r8), allocatable, save :: gcam_consump_frac(:)

contains

  subroutine gcam_set_water_supply(water_supply, n) bind(C, name="gcam_set_water_supply_f")
    !---------------------------------------------------------------
    ! Called from Fortran to store water supply for GCAM access
    !---------------------------------------------------------------
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: water_supply(n)

    if (.not. allocated(gcam_water_supply)) allocate(gcam_water_supply(n))
    gcam_water_supply(1:n) = water_supply(1:n)

  end subroutine gcam_set_water_supply

  subroutine gcam_get_water_demand(demand_i, demand_m, demand_u, demand_e, &
                                   demand_t, consump, n) bind(C, name="gcam_get_water_demand_f")
    !---------------------------------------------------------------
    ! Called from Fortran to retrieve water demand from GCAM
    !---------------------------------------------------------------
    integer(c_int), intent(in), value :: n
    real(c_double), intent(out) :: demand_i(n), demand_m(n), demand_u(n)
    real(c_double), intent(out) :: demand_e(n), demand_t(n), consump(n)

    if (allocated(gcam_demand_irrig)) then
       demand_i(1:n) = gcam_demand_irrig(1:n)
       demand_m(1:n) = gcam_demand_indust(1:n)
       demand_u(1:n) = gcam_demand_munic(1:n)
       demand_e(1:n) = gcam_demand_energy(1:n)
       demand_t(1:n) = gcam_demand_total(1:n)
       consump(1:n)  = gcam_consump_frac(1:n)
    else
       ! Return zeros if not initialized
       demand_i(:) = 0.0_r8
       demand_m(:) = 0.0_r8
       demand_u(:) = 0.0_r8
       demand_e(:) = 0.0_r8
       demand_t(:) = 0.0_r8
       consump(:)  = 0.5_r8  ! default consumptive fraction
    endif

  end subroutine gcam_get_water_demand

end module gcam_water_interface
```

---

# Part 2: Changes Required in `gcam-core` Repository

GCAM needs to implement the water resource management logic that uses MOSART water availability and generates sector-specific water demands.

## 2.1 Water Supply Module

### File: `cvs/objects/resources/include/water_supply.h`

```cpp
#ifndef _WATER_SUPPLY_H_
#define _WATER_SUPPLY_H_

#include <vector>
#include <map>
#include <string>
#include "util/base/include/inamed.h"

class WaterSupply : public INamed {
public:
    WaterSupply();
    ~WaterSupply();

    // Set water availability from E3SM/MOSART
    void setWaterAvailability(const std::string& aRegionName,
                             const double aTotalAvailability,
                             const double aReservoirStorage,
                             const double aStreamflow);

    // Get available water for extraction (after environmental flows)
    double getExtractableWater(const std::string& aRegionName,
                              const int aPeriod) const;

    // Calculate environmental flow requirements
    double calculateEnvironmentalFlows(const double aTotalAvailability) const;

    // Apply water extraction (update available volume)
    void applyExtraction(const std::string& aRegionName,
                        const double aExtraction,
                        const int aPeriod);

private:
    // Regional water availability (m³/year)
    std::map<std::string, double> mRegionalAvailability;

    // Reservoir storage capacity (m³)
    std::map<std::string, double> mReservoirStorage;

    // Mean streamflow (m³/s)
    std::map<std::string, double> mStreamflow;

    // Environmental flow requirement fraction (default 0.3)
    double mEnvironmentalFlowFraction;

    // Remaining extractable water after allocations
    std::map<std::string, double> mRemainingWater;
};

#endif // _WATER_SUPPLY_H_
```

### File: `cvs/objects/resources/source/water_supply.cpp`

```cpp
#include "resources/include/water_supply.h"
#include "util/base/include/xml_helper.h"

WaterSupply::WaterSupply()
: mEnvironmentalFlowFraction(0.3)  // Reserve 30% for environmental flows
{
}

WaterSupply::~WaterSupply() {
}

void WaterSupply::setWaterAvailability(const std::string& aRegionName,
                                       const double aTotalAvailability,
                                       const double aReservoirStorage,
                                       const double aStreamflow) {
    mRegionalAvailability[aRegionName] = aTotalAvailability;
    mReservoirStorage[aRegionName] = aReservoirStorage;
    mStreamflow[aRegionName] = aStreamflow;

    // Initialize remaining water as total minus environmental flows
    double envFlows = calculateEnvironmentalFlows(aTotalAvailability);
    mRemainingWater[aRegionName] = aTotalAvailability - envFlows;
}

double WaterSupply::getExtractableWater(const std::string& aRegionName,
                                       const int aPeriod) const {
    auto it = mRemainingWater.find(aRegionName);
    if (it != mRemainingWater.end()) {
        return it->second;
    }
    return 0.0;
}

double WaterSupply::calculateEnvironmentalFlows(const double aTotalAvailability) const {
    return aTotalAvailability * mEnvironmentalFlowFraction;
}

void WaterSupply::applyExtraction(const std::string& aRegionName,
                                 const double aExtraction,
                                 const int aPeriod) {
    mRemainingWater[aRegionName] -= aExtraction;

    // Ensure non-negative
    if (mRemainingWater[aRegionName] < 0.0) {
        mRemainingWater[aRegionName] = 0.0;
    }
}
```

## 2.2 Water Demand Calculator

### File: `cvs/objects/sectors/include/water_demand_sector.h`

```cpp
#ifndef _WATER_DEMAND_SECTOR_H_
#define _WATER_DEMAND_SECTOR_H_

#include <string>
#include <map>
#include "sectors/include/sector.h"

// Forward declarations
class Scenario;

class WaterDemandSector : public Sector {
public:
    WaterDemandSector();
    virtual ~WaterDemandSector();

    // Calculate water demand for this sector
    virtual void calcWaterDemand(const Scenario* aScenario,
                                const std::string& aRegionName,
                                const int aPeriod);

    // Get calculated demand by type
    double getIrrigationDemand(const std::string& aRegionName,
                              const int aPeriod) const;
    double getIndustrialDemand(const std::string& aRegionName,
                              const int aPeriod) const;
    double getMunicipalDemand(const std::string& aRegionName,
                             const int aPeriod) const;
    double getEnergyDemand(const std::string& aRegionName,
                          const int aPeriod) const;
    double getTotalDemand(const std::string& aRegionName,
                         const int aPeriod) const;

    // Get consumptive use fraction
    double getConsumptiveFraction(const std::string& aRegionName,
                                 const int aPeriod) const;

protected:
    // Sector-specific demand calculations
    double calcIrrigationDemand(const Scenario* aScenario,
                               const std::string& aRegionName,
                               const int aPeriod);
    double calcIndustrialDemand(const Scenario* aScenario,
                               const std::string& aRegionName,
                               const int aPeriod);
    double calcMunicipalDemand(const Scenario* aScenario,
                              const std::string& aRegionName,
                              const int aPeriod);
    double calcEnergyDemand(const Scenario* aScenario,
                           const std::string& aRegionName,
                           const int aPeriod);

private:
    // Demand by region and period (m³/year)
    std::map<std::pair<std::string, int>, double> mIrrigationDemand;
    std::map<std::pair<std::string, int>, double> mIndustrialDemand;
    std::map<std::pair<std::string, int>, double> mMunicipalDemand;
    std::map<std::pair<std::string, int>, double> mEnergyDemand;

    // Consumptive use fractions by sector
    std::map<std::pair<std::string, int>, double> mConsumptiveFraction;
};

#endif // _WATER_DEMAND_SECTOR_H_
```

### File: `cvs/objects/sectors/source/water_demand_sector.cpp`

```cpp
#include "sectors/include/water_demand_sector.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/model_time.h"

WaterDemandSector::WaterDemandSector() {
}

WaterDemandSector::~WaterDemandSector() {
}

void WaterDemandSector::calcWaterDemand(const Scenario* aScenario,
                                       const std::string& aRegionName,
                                       const int aPeriod) {
    // Calculate demand for each sector
    double irrigDemand = calcIrrigationDemand(aScenario, aRegionName, aPeriod);
    double industDemand = calcIndustrialDemand(aScenario, aRegionName, aPeriod);
    double municDemand = calcMunicipalDemand(aScenario, aRegionName, aPeriod);
    double energyDemand = calcEnergyDemand(aScenario, aRegionName, aPeriod);

    // Store demands
    auto key = std::make_pair(aRegionName, aPeriod);
    mIrrigationDemand[key] = irrigDemand;
    mIndustrialDemand[key] = industDemand;
    mMunicipalDemand[key] = municDemand;
    mEnergyDemand[key] = energyDemand;

    // Calculate weighted consumptive fraction
    double totalDemand = irrigDemand + industDemand + municDemand + energyDemand;
    if (totalDemand > 0) {
        // Irrigation: ~60% consumptive
        // Industrial: ~10% consumptive
        // Municipal: ~20% consumptive
        // Energy: ~5% consumptive
        double weightedConsump = (irrigDemand * 0.6 + industDemand * 0.1 +
                                 municDemand * 0.2 + energyDemand * 0.05) / totalDemand;
        mConsumptiveFraction[key] = weightedConsump;
    } else {
        mConsumptiveFraction[key] = 0.3;  // default
    }
}

double WaterDemandSector::calcIrrigationDemand(const Scenario* aScenario,
                                              const std::string& aRegionName,
                                              const int aPeriod) {
    // Calculate based on:
    // - Irrigated cropland area
    // - Crop water requirements
    // - Irrigation efficiency
    // - Climate/precipitation

    // Placeholder - implement based on GCAM's agriculture module
    return 0.0;
}

double WaterDemandSector::calcIndustrialDemand(const Scenario* aScenario,
                                              const std::string& aRegionName,
                                              const int aPeriod) {
    // Calculate based on:
    // - Industrial output/GDP
    // - Water intensity factors
    // - Technology improvements

    // Placeholder - implement based on GCAM's industry module
    return 0.0;
}

double WaterDemandSector::calcMunicipalDemand(const Scenario* aScenario,
                                             const std::string& aRegionName,
                                             const int aPeriod) {
    // Calculate based on:
    // - Population
    // - Per-capita water use
    // - Income/development level

    // Placeholder - implement based on GCAM's demographics
    return 0.0;
}

double WaterDemandSector::calcEnergyDemand(const Scenario* aScenario,
                                          const std::string& aRegionName,
                                          const int aPeriod) {
    // Calculate based on:
    // - Electricity generation by technology
    // - Cooling system types
    // - Water withdrawal/consumption factors

    // Placeholder - implement based on GCAM's energy module
    return 0.0;
}

double WaterDemandSector::getTotalDemand(const std::string& aRegionName,
                                        const int aPeriod) const {
    auto key = std::make_pair(aRegionName, aPeriod);
    double total = 0.0;

    auto it_i = mIrrigationDemand.find(key);
    if (it_i != mIrrigationDemand.end()) total += it_i->second;

    auto it_m = mIndustrialDemand.find(key);
    if (it_m != mIndustrialDemand.end()) total += it_m->second;

    auto it_u = mMunicipalDemand.find(key);
    if (it_u != mMunicipalDemand.end()) total += it_u->second;

    auto it_e = mEnergyDemand.find(key);
    if (it_e != mEnergyDemand.end()) total += it_e->second;

    return total;
}
```

## 2.3 Coupling Interface

### File: `cvs/objects/climate/include/e3sm_coupling_interface.h`

```cpp
#ifndef _E3SM_COUPLING_INTERFACE_H_
#define _E3SM_COUPLING_INTERFACE_H_

#include <vector>
#include <string>

// External C interface for Fortran interoperability
extern "C" {
    // Called from giac to provide water supply to GCAM
    void gcam_receive_water_supply(const double* water_availability,
                                   const int* region_ids,
                                   const int n_regions);

    // Called by giac to retrieve water demand from GCAM
    void gcam_provide_water_demand(double* demand_irrigation,
                                   double* demand_industrial,
                                   double* demand_municipal,
                                   double* demand_energy,
                                   double* demand_total,
                                   double* consumptive_fraction,
                                   const int* region_ids,
                                   const int n_regions);
}

class E3SMCouplingInterface {
public:
    static E3SMCouplingInterface* getInstance();

    void receiveWaterSupply(const std::vector<double>& aWaterAvailability,
                           const std::vector<int>& aRegionIDs);

    void provideWaterDemand(std::vector<double>& aDemandIrrig,
                          std::vector<double>& aDemandIndust,
                          std::vector<double>& aDemandMunic,
                          std::vector<double>& aDemandEnergy,
                          std::vector<double>& aDemandTotal,
                          std::vector<double>& aConsumpFrac,
                          const std::vector<int>& aRegionIDs);

private:
    E3SMCouplingInterface();
    static E3SMCouplingInterface* sInstance;

    std::map<int, double> mWaterSupplyByRegion;
    std::map<int, double> mDemandByRegion;
};

#endif // _E3SM_COUPLING_INTERFACE_H_
```

### File: `cvs/objects/climate/source/e3sm_coupling_interface.cpp`

```cpp
#include "climate/include/e3sm_coupling_interface.h"
#include "containers/include/world.h"
#include "sectors/include/water_demand_sector.h"

E3SMCouplingInterface* E3SMCouplingInterface::sInstance = nullptr;

E3SMCouplingInterface* E3SMCouplingInterface::getInstance() {
    if (!sInstance) {
        sInstance = new E3SMCouplingInterface();
    }
    return sInstance;
}

E3SMCouplingInterface::E3SMCouplingInterface() {
}

void E3SMCouplingInterface::receiveWaterSupply(
    const std::vector<double>& aWaterAvailability,
    const std::vector<int>& aRegionIDs) {

    for (size_t i = 0; i < aRegionIDs.size(); ++i) {
        mWaterSupplyByRegion[aRegionIDs[i]] = aWaterAvailability[i];
    }

    // Update GCAM's water supply objects
    // This requires integration with GCAM's region/resource framework
}

void E3SMCouplingInterface::provideWaterDemand(
    std::vector<double>& aDemandIrrig,
    std::vector<double>& aDemandIndust,
    std::vector<double>& aDemandMunic,
    std::vector<double>& aDemandEnergy,
    std::vector<double>& aDemandTotal,
    std::vector<double>& aConsumpFrac,
    const std::vector<int>& aRegionIDs) {

    // Retrieve calculated demands from GCAM
    // This requires accessing GCAM's water demand objects by region

    for (size_t i = 0; i < aRegionIDs.size(); ++i) {
        int regionID = aRegionIDs[i];

        // Get demands from water demand sector
        // (Requires integration with GCAM scenario/region structure)
        aDemandIrrig[i] = 0.0;   // Placeholder
        aDemandIndust[i] = 0.0;  // Placeholder
        aDemandMunic[i] = 0.0;   // Placeholder
        aDemandEnergy[i] = 0.0;  // Placeholder
        aDemandTotal[i] = 0.0;   // Placeholder
        aConsumpFrac[i] = 0.3;   // Placeholder
    }
}

// C interface implementations
extern "C" {
    void gcam_receive_water_supply(const double* water_availability,
                                   const int* region_ids,
                                   const int n_regions) {
        std::vector<double> avail(water_availability, water_availability + n_regions);
        std::vector<int> ids(region_ids, region_ids + n_regions);

        E3SMCouplingInterface::getInstance()->receiveWaterSupply(avail, ids);
    }

    void gcam_provide_water_demand(double* demand_irrigation,
                                   double* demand_industrial,
                                   double* demand_municipal,
                                   double* demand_energy,
                                   double* demand_total,
                                   double* consumptive_fraction,
                                   const int* region_ids,
                                   const int n_regions) {
        std::vector<double> d_irrig(n_regions);
        std::vector<double> d_indust(n_regions);
        std::vector<double> d_munic(n_regions);
        std::vector<double> d_energy(n_regions);
        std::vector<double> d_total(n_regions);
        std::vector<double> consump(n_regions);
        std::vector<int> ids(region_ids, region_ids + n_regions);

        E3SMCouplingInterface::getInstance()->provideWaterDemand(
            d_irrig, d_indust, d_munic, d_energy, d_total, consump, ids);

        // Copy results to output arrays
        std::copy(d_irrig.begin(), d_irrig.end(), demand_irrigation);
        std::copy(d_indust.begin(), d_indust.end(), demand_industrial);
        std::copy(d_munic.begin(), d_munic.end(), demand_municipal);
        std::copy(d_energy.begin(), d_energy.end(), demand_energy);
        std::copy(d_total.begin(), d_total.end(), demand_total);
        std::copy(consump.begin(), consump.end(), consumptive_fraction);
    }
}
```

---

# Part 3: Integration Steps

## 3.1 Build System Updates

### For `giac`:

Update `Makefile` or CMake configuration to include new water coupling modules:

```makefile
# Add to OBJS
OBJS += water_grid_mapping.o \
        gcam_water_interface.o

# Ensure linking with GCAM libraries
LDFLAGS += -L$(GCAM_LIB_DIR) -lgcam
```

### For `gcam-core`:

Update build system to compile new water modules:

```cmake
# In CMakeLists.txt or equivalent
set(WATER_SOURCES
    cvs/objects/resources/source/water_supply.cpp
    cvs/objects/sectors/source/water_demand_sector.cpp
    cvs/objects/climate/source/e3sm_coupling_interface.cpp
)

add_library(gcam ${GCAM_SOURCES} ${WATER_SOURCES})
```

## 3.2 Initialization Sequence

The coupling must be initialized in the correct order:

1. **E3SM Driver**: Initialize coupling flags (`rof_c2_iac`, `iac_c2_rof`)
2. **giac**: Initialize grid mappings between MOSART and GCAM
3. **GCAM**: Initialize water supply and demand objects
4. **First Exchange**: MOSART provides initial water availability
5. **GCAM Solve**: GCAM solves for water demands given availability
6. **Return Demands**: GCAM returns sector-specific demands to MOSART

## 3.3 Runtime Data Exchange

Annual coupling cycle:

```
Year N:
  1. MOSART accumulates water availability over 365 days
  2. At coupling time (Jan 1, 00:30 UTC):
     - MOSART averages accumulated data
     - E3SM driver maps to IAC grid
     - giac passes to GCAM via C interface
  3. GCAM solves for year N:
     - Uses water availability as constraint
     - Calculates sectoral demands
     - Computes allocation priorities
  4. GCAM returns demands:
     - giac receives via C interface
     - E3SM driver maps to MOSART grid
     - MOSART uses for year N+1 simulation
```

---

# Part 4: Testing and Validation

## 4.1 Unit Tests

### For `giac`:
- Test grid mapping weight calculations
- Verify conservation of water mass during aggregation/disaggregation
- Test Fortran-C interface boundary conditions

### For `gcam-core`:
- Test water demand calculations for each sector
- Verify allocation logic with various scarcity scenarios
- Test consumptive fraction calculations

## 4.2 Integration Tests

1. **Mass Balance Test**: Verify total water extracted ≤ total water available
2. **Coupling Frequency Test**: Ensure annual exchange timing is correct
3. **Multi-Year Simulation**: Test stability over multiple coupling cycles
4. **Extreme Scenarios**: Test with very low/high water availability

## 4.3 Validation Data

Compare model output against:
- Historical water withdrawal data by sector
- Regional water scarcity indicators
- Agricultural water use from FAO AQUASTAT
- Energy sector water use from IEA/EIA data

---

# Part 5: Configuration Files

## 5.1 Mapping Weight Files

Create basin-to-region mapping file (NetCDF or text format):

```
# Format: MOSART_cell_id  GCAM_region_id  weight
# weight = fraction of cell area in region
100001  1  0.80
100001  2  0.20
100002  1  1.00
...
```

## 5.2 GCAM Configuration

Update GCAM XML configuration to enable water coupling:

```xml
<scenario>
  <water-coupling enabled="1">
    <coupling-frequency>annual</coupling-frequency>
    <environmental-flow-fraction>0.30</environmental-flow-fraction>
    <sectors>
      <sector name="irrigation" consumptive-fraction="0.60"/>
      <sector name="industrial" consumptive-fraction="0.10"/>
      <sector name="municipal" consumptive-fraction="0.20"/>
      <sector name="energy" consumptive-fraction="0.05"/>
    </sectors>
  </water-coupling>
</scenario>
```

---

# Part 6: Known Issues and Future Work

## Current Limitations

1. **Spatial Resolution Mismatch**: MOSART river network vs. GCAM economic regions
2. **Temporal Aggregation**: Annual coupling may miss sub-annual variability
3. **Groundwater**: Current design focuses on surface water only
4. **Water Quality**: No tracking of water quality/pollution
5. **Irrigation Efficiency**: Fixed assumptions may not reflect technology changes

## Recommended Enhancements

1. Implement groundwater coupling
2. Add water quality tracking
3. Include reservoir operation optimization
4. Dynamic irrigation efficiency based on economic conditions
5. Sub-annual coupling for water-stressed regions

---

# Part 7: Contact and Support

For questions about implementation:

- **E3SM coupling infrastructure**: E3SM Water Group
- **MOSART model**: MOSART development team
- **GCAM model**: GCAM development team (JGCRI)
- **giac wrapper**: E3SM-GCAM coupling team

---

# Appendices

## Appendix A: Field Definitions

### MOSART → GCAM Fields

| Field Name | Units | Description |
|-----------|-------|-------------|
| `Sr_wr_avail` | m³ | Main channel water storage available for extraction |
| `Sr_wt_avail` | m³ | Tributary water storage available for extraction |
| `Sr_wtot_avail` | m³ | Total surface water available (main + tributary) |
| `Sr_reservoir_stor` | m³ | Reservoir storage capacity |
| `Sr_streamflow` | m³/s | Annual mean streamflow |

### GCAM → MOSART Fields

| Field Name | Units | Description |
|-----------|-------|-------------|
| `Sz_demand_irrig` | m³/year | Water demand for irrigation |
| `Sz_demand_indust` | m³/year | Water demand for industrial use |
| `Sz_demand_munic` | m³/year | Water demand for municipal use |
| `Sz_demand_energy` | m³/year | Water demand for energy production |
| `Sz_demand_total` | m³/year | Total water demand (sum of all sectors) |
| `Sz_consump_frac` | - | Fraction of withdrawn water consumed (0-1) |

## Appendix B: Coupling Timing

The coupling uses E3SM's IAC alarm system:
- **Frequency**: Annual
- **Timing**: First day of each year at 00:30 UTC
- **Accumulation**: MOSART accumulates daily, averages annually
- **Application**: GCAM demands applied to following year

## Appendix C: Unit Conversions

Common conversions needed:
- MOSART internal: m³/s
- GCAM water demands: m³/year
- Conversion: multiply m³/s by 31,536,000 (seconds/year)
- Area-specific: mm/year ↔ m³/year requires grid cell area

---

**Document Version**: 1.0
**Date**: 2025-11-18
**Based on E3SM Commit**: f0438dfdcc
**giac Base**: 0ffce23d17c01c0f4e598891243abd4bc37fe900
**gcam-core Base**: 2901ad96aadf0a2ba272ef7718419dd8453af022
