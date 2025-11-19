module MOSART_HeatBudgets_mod
! Description: MOSART global heat and energy budget diagnostics
!
! Created for river-atmosphere coupling with MOSART-heat
!-----------------------------------------------------------------------
  ! USES:
  use rof_cpl_indices, only : nt_rtm, nt_nliq, nt_nice
  use RtmVar         , only : iulog
  use RtmSpmd        , only : masterproc
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_sys_mod    , only : shr_sys_abort
  use shr_const_mod  , only : SHR_CONST_REARTH ! earth radius in m
  use shr_const_mod  , only : shr_const_pi

  implicit none
  private

  public MOSART_HeatBudget_Reset
  public MOSART_HeatBudget_Accumulate
  public MOSART_HeatBudget_Print

  !--- F for flux (W/m2 or energy flux density) ---

  integer, parameter :: f_swabs  = 1  ! Shortwave radiation absorbed
  integer, parameter :: f_lwnet  = 2  ! Net longwave radiation
  integer, parameter :: f_sensible = 3  ! Sensible heat flux
  integer, parameter :: f_latent = 4  ! Latent heat flux (evaporation)
  integer, parameter :: f_conductive = 5  ! Conductive heat exchange
  integer, parameter :: f_advective = 6  ! Advective heat transport

  integer, parameter, public :: f_size = f_advective

  character(len=20),parameter :: fname(f_size) = &
       (/&
       '   SW absorbed (in)', &
       '   LW net (in-out) ', &
       '  Sensible (out)   ', &
       '   Latent (out)    ', &
       ' Conductive (bed)  ', &
       '  Advective (net)  '  &
       /)

  !--- O for "other" term aggregated from flux to state ---

  integer, parameter :: o_evap = 1  ! Evaporated water mass (for water budget check)

  integer, parameter, public :: o_size = o_evap

  !--- S for state ---
  integer, parameter :: s_hcontent_beg = 1  ! River heat content beginning (J)
  integer, parameter :: s_hcontent_end = 2  ! River heat content end (J)
  integer, parameter :: s_hmain_beg    = 3  ! Main channel heat beginning (J)
  integer, parameter :: s_hmain_end    = 4  ! Main channel heat end (J)
  integer, parameter :: s_htrib_beg    = 5  ! Tributary heat beginning (J)
  integer, parameter :: s_htrib_end    = 6  ! Tributary heat end (J)

  integer, parameter, public :: s_size = s_htrib_end

  !--- P for period ---

  integer, parameter :: p_inst = 1
  integer, parameter :: p_day  = 2
  integer, parameter :: p_mon  = 3
  integer, parameter :: p_ann  = 4
  integer, parameter :: p_inf  = 5

  integer, parameter, public :: p_size = p_inf

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  real(r8), public :: rof_budg_hfluxG  (f_size, p_size) = 0.0_r8 ! global heat flux sum (W/m2)
  real(r8), public :: rof_budg_hother  (o_size, p_size) = 0.0_r8 ! global "other" term
  real(r8), public :: rof_budg_hstateG (s_size, p_size) = 0.0_r8 ! global heat state sum (J)
  real(r8), public :: rof_budg_hfluxN  (f_size, p_size) = 0.0_r8 ! counter, valid only on root pe

  !----- formats -----
  character(*),parameter :: FA0= "('    ',12x,(3x,a18,2x),' | ',(3x,a18,2x))"
  character(*),parameter :: FF = "('',a20,e18.8,' | ',e21.2)"
  character(*),parameter :: FS0= "(' ',12x,3(a22),' | ',(a22))"
  character(*),parameter :: FS2= "(' ',a12,22x,e22.8,22x,' | ',e22.8)"
  character(*),parameter :: FS3= "(' ',a12,3(e22.8),' | ',(e22.8))"

  !----- other variables -----
  real(r8)               :: unit_conversion

contains

  subroutine MOSART_HeatBudget_Reset(mode)
   !
   use RtmTimeManager, only : get_curr_date, get_prev_date
   !
   implicit none
   !
   character(len=*), intent(in),optional :: mode
   !
   integer :: year, mon, day, sec
   integer :: ip
   character(*),parameter :: subName = '(MOSART_HeatBudget_Reset) '

   if (.not.present(mode)) then
      call get_curr_date(year, mon, day, sec)

      do ip = 1,p_size
         if (ip == p_inst) then
            rof_budg_hfluxG(:,ip)  = 0.0_r8
            rof_budg_hfluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_day .and. sec==0) then
            rof_budg_hfluxG(:,ip)  = 0.0_r8
            rof_budg_hfluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_mon .and. day==1 .and. sec==0) then
            rof_budg_hfluxG(:,ip)  = 0.0_r8
            rof_budg_hfluxN(:,ip)  = 0.0_r8
         endif
         if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
            rof_budg_hfluxG(:,ip)  = 0.0_r8
            rof_budg_hfluxN(:,ip)  = 0.0_r8
         endif
      enddo

   else

      if (trim(mode) == 'inst') then
         rof_budg_hfluxG  (:,p_inst)   = 0.0_r8
         rof_budg_hother  (:,p_inst)   = 0.0_r8
         rof_budg_hstateG (:,p_inst)   = 0.0_r8
         rof_budg_hfluxN  (:,p_inst)   = 0.0_r8
      elseif (trim(mode) == 'day') then
         rof_budg_hfluxG  (:,p_day)    = 0.0_r8
         rof_budg_hother  (:,p_day)    = 0.0_r8
         rof_budg_hstateG (:,p_day)    = 0.0_r8
         rof_budg_hfluxN  (:,p_day)    = 0.0_r8
      elseif (trim(mode) == 'mon') then
         rof_budg_hfluxG  (:,p_mon)    = 0.0_r8
         rof_budg_hother  (:,p_mon)    = 0.0_r8
         rof_budg_hstateG (:,p_mon)    = 0.0_r8
         rof_budg_hfluxN  (:,p_mon)    = 0.0_r8
      elseif (trim(mode) == 'ann') then
         rof_budg_hfluxG  (:,p_ann)    = 0.0_r8
         rof_budg_hother  (:,p_ann)    = 0.0_r8
         rof_budg_hstateG (:,p_ann)    = 0.0_r8
         rof_budg_hfluxN  (:,p_ann)    = 0.0_r8
      elseif (trim(mode) == 'inf') then
         rof_budg_hfluxG  (:,p_inf)    = 0.0_r8
         rof_budg_hother  (:,p_inf)    = 0.0_r8
         rof_budg_hstateG (:,p_inf)    = 0.0_r8
         rof_budg_hfluxN  (:,p_inf)    = 0.0_r8
      elseif (trim(mode) == 'all') then
         rof_budg_hfluxG  (:,:)        = 0.0_r8
         rof_budg_hother  (:,:)        = 0.0_r8
         rof_budg_hstateG (:,:)        = 0.0_r8
         rof_budg_hfluxN  (:,:)        = 0.0_r8
      else
         call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
      endif
   endif

  end subroutine MOSART_HeatBudget_Reset

!--------------------------------------------------------------------
  subroutine MOSART_HeatBudget_Accumulate(sw_flux, lw_flux, sens_flux, lat_flux, &
                                          cond_flux, adv_flux, evap_mass, &
                                          hcontent_beg, hcontent_end, &
                                          hmain_beg, hmain_end, htrib_beg, htrib_end)
   !
   ! Accumulate heat budget terms for a single gridcell
   ! Called within the routing loop for each active gridcell
   !
   use RtmTimeManager, only : get_curr_date, get_prev_date, get_nstep
   implicit none
   !
   real(r8), intent(in) :: sw_flux       ! Shortwave absorbed (W)
   real(r8), intent(in) :: lw_flux       ! Net longwave (W)
   real(r8), intent(in) :: sens_flux     ! Sensible heat (W)
   real(r8), intent(in) :: lat_flux      ! Latent heat (W)
   real(r8), intent(in) :: cond_flux     ! Conductive heat (W)
   real(r8), intent(in) :: adv_flux      ! Advective heat (W)
   real(r8), intent(in) :: evap_mass     ! Evaporated mass (kg/s)
   real(r8), intent(in) :: hcontent_beg  ! Total heat content beginning (J)
   real(r8), intent(in) :: hcontent_end  ! Total heat content end (J)
   real(r8), intent(in) :: hmain_beg     ! Main channel heat beginning (J)
   real(r8), intent(in) :: hmain_end     ! Main channel heat end (J)
   real(r8), intent(in) :: htrib_beg     ! Tributary heat beginning (J)
   real(r8), intent(in) :: htrib_end     ! Tributary heat end (J)
   !
   integer                :: ip
   integer                :: year_prev, month_prev, day_prev, sec_prev
   integer                :: year_curr, month_curr, day_curr, sec_curr
   logical                :: update_state_beg, update_state_end
   !
   character(*),parameter :: subName = '(MOSART_HeatBudget_Accumulate) '

   ! Convert to W/m2 using Earth's surface area
   ! unit_conversion converts from total watts to W/m2
   unit_conversion = 1.d0/(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)

   ! Accumulate fluxes (instantaneous values)
   rof_budg_hfluxG(f_swabs,     p_inst) = rof_budg_hfluxG(f_swabs,     p_inst) + sw_flux   * unit_conversion
   rof_budg_hfluxG(f_lwnet,     p_inst) = rof_budg_hfluxG(f_lwnet,     p_inst) + lw_flux   * unit_conversion
   rof_budg_hfluxG(f_sensible,  p_inst) = rof_budg_hfluxG(f_sensible,  p_inst) + sens_flux * unit_conversion
   rof_budg_hfluxG(f_latent,    p_inst) = rof_budg_hfluxG(f_latent,    p_inst) + lat_flux  * unit_conversion
   rof_budg_hfluxG(f_conductive,p_inst) = rof_budg_hfluxG(f_conductive,p_inst) + cond_flux * unit_conversion
   rof_budg_hfluxG(f_advective, p_inst) = rof_budg_hfluxG(f_advective, p_inst) + adv_flux  * unit_conversion

   ! Accumulate other terms
   rof_budg_hother(o_evap, p_inst) = rof_budg_hother(o_evap, p_inst) + evap_mass

   ! Accumulate states
   rof_budg_hstateG(s_hcontent_beg, p_inst) = rof_budg_hstateG(s_hcontent_beg, p_inst) + hcontent_beg
   rof_budg_hstateG(s_hcontent_end, p_inst) = rof_budg_hstateG(s_hcontent_end, p_inst) + hcontent_end
   rof_budg_hstateG(s_hmain_beg,    p_inst) = rof_budg_hstateG(s_hmain_beg,    p_inst) + hmain_beg
   rof_budg_hstateG(s_hmain_end,    p_inst) = rof_budg_hstateG(s_hmain_end,    p_inst) + hmain_end
   rof_budg_hstateG(s_htrib_beg,    p_inst) = rof_budg_hstateG(s_htrib_beg,    p_inst) + htrib_beg
   rof_budg_hstateG(s_htrib_end,    p_inst) = rof_budg_hstateG(s_htrib_end,    p_inst) + htrib_end

   ! Propagate to other time periods
   call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
   call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

   do ip = p_inst+1, p_size
      rof_budg_hfluxG(:,ip) = rof_budg_hfluxG(:,ip) + rof_budg_hfluxG(:,p_inst)
      rof_budg_hother(:,ip) = rof_budg_hother(:,ip) + rof_budg_hother(:,p_inst)

      update_state_beg = .false.
      update_state_end = .false.

      select case (ip)
      case (p_day)
         if (sec_prev == 0) update_state_beg = .true.
         if (sec_curr == 0) update_state_end = .true.
      case (p_mon)
         if (sec_prev == 0 .and. day_prev == 1) update_state_beg = .true.
         if (sec_curr == 0 .and. day_curr == 1) update_state_end = .true.
      case (p_ann)
         if (sec_prev == 0 .and. day_prev == 1 .and. month_prev == 1) update_state_beg = .true.
         if (sec_curr == 0 .and. day_curr == 1 .and. month_curr == 1) update_state_end = .true.
      case (p_inf)
         if (get_nstep() == 1) update_state_beg = .true.
         update_state_end = .true.
      end select

      if (update_state_beg) then
         rof_budg_hstateG(:,ip) = rof_budg_hstateG(:, p_inst)
      endif

      if (update_state_end) then
         rof_budg_hstateG(:,ip) = rof_budg_hstateG(:, p_inst)
      endif

   end do

   rof_budg_hfluxN(:,:) = rof_budg_hfluxN(:,:) + 1._r8

  end subroutine MOSART_HeatBudget_Accumulate

!-----------------------------------------------------------------------
  subroutine MOSART_HeatBudget_Print()
   !
   use RtmTimeManager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
   !
   implicit none
   !
   integer :: budg_print_inst  = 0
   integer :: budg_print_daily = 0
   integer :: budg_print_month = 1
   integer :: budg_print_ann   = 1
   integer :: budg_print_ltann = 1
   integer :: budg_print_ltend = 0
   !
   ! !LOCAL VARIABLES:
   integer :: f,ip
   integer :: plev
   integer :: year, mon, day, sec
   integer :: cdate
   logical :: sumdone
   real(r8) :: budg_hfluxGpr (f_size,p_size)
   real(r8) :: net_heat_flux, dH_dt, residual

   character(*),parameter :: subName = '(MOSART_HeatBudget_Print) '

   sumdone = .false.

   if (get_nstep() <= 1) then
      call get_prev_date(year, mon, day, sec);
   else
      call get_curr_date(year, mon, day, sec);
   end if

   cdate = year*10000 + mon*100 + day

   do ip = 1,p_size
      plev = 0
      if (ip == p_inst) then
         plev = max(plev,budg_print_inst)
      endif
      if (ip==p_day .and. sec==0) then
         plev = max(plev,budg_print_daily)
      endif
      if (ip==p_mon .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_month)
      endif
      if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ann)
      endif
      if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
         plev = max(plev,budg_print_ltann)
      endif

      if (plev > 0) then
         if (.not.sumdone) then
            sumdone = .true.
            budg_hfluxGpr = rof_budg_hfluxG
            budg_hfluxGpr = budg_hfluxGpr/rof_budg_hfluxN
         end if

         if (ip == p_day .and. get_nstep() == 1) cycle
         if (ip == p_mon .and. get_nstep() == 1) cycle
         if (ip == p_ann .and. get_nstep() == 1) cycle
         if (ip == p_inf .and. get_nstep() == 1) cycle

         if (masterproc) then
            write(iulog,*)''
            write(iulog,*)'=========================================='
            write(iulog,*)'RIVER HEAT FLUXES: period ',trim(pname(ip)),': date = ',cdate,sec
            write(iulog,*)'=========================================='
            write(iulog,FA0)'  Time averaged','  Global total '
            write(iulog,FA0)'     (W/m2)    ','      (PW)     '
            write(iulog,'(42("-"),"|",25("-"))')
            do f = 1, f_size
               write(iulog,FF)fname(f),budg_hfluxGpr(f,ip), &
                    rof_budg_hfluxG(f,ip)/(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)*1.0e-15_r8
            end do
            write(iulog,'(42("-"),"|",25("-"))')

            ! Calculate net heat flux (positive = heating)
            net_heat_flux = budg_hfluxGpr(f_swabs,ip) + budg_hfluxGpr(f_lwnet,ip) &
                          - budg_hfluxGpr(f_sensible,ip) - budg_hfluxGpr(f_latent,ip) &
                          + budg_hfluxGpr(f_conductive,ip) + budg_hfluxGpr(f_advective,ip)

            write(iulog,FF)'   *NET FLUX*', net_heat_flux, &
                 net_heat_flux*(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2)*1.0e-15_r8
            write(iulog,'(42("-"),"|",25("-"))')
            write(iulog,*)''

            write(iulog,*)'RIVER HEAT STATES (PJ): period ',trim(pname(ip)),': date = ',cdate,sec
            write(iulog,FS0) &
                  '    Main Channel ', &
                  '       Tributary ', &
                  '      TOTAL HEAT '
            write(iulog,'(78("-"),"|",25("-"))')
            write(iulog,FS3) '         beg', &
                  rof_budg_hstateG(s_hmain_beg, ip)*1.0e-15_r8, &
                  rof_budg_hstateG(s_htrib_beg, ip)*1.0e-15_r8, &
                  rof_budg_hstateG(s_hcontent_beg, ip)*1.0e-15_r8
            write(iulog,FS3) '         end', &
                  rof_budg_hstateG(s_hmain_end, ip)*1.0e-15_r8, &
                  rof_budg_hstateG(s_htrib_end, ip)*1.0e-15_r8, &
                  rof_budg_hstateG(s_hcontent_end, ip)*1.0e-15_r8
            write(iulog,FS3) ' *CHANGE(dH)*', &
                  (rof_budg_hstateG(s_hmain_end,ip) - rof_budg_hstateG(s_hmain_beg,ip))*1.0e-15_r8, &
                  (rof_budg_hstateG(s_htrib_end,ip) - rof_budg_hstateG(s_htrib_beg,ip))*1.0e-15_r8, &
                  (rof_budg_hstateG(s_hcontent_end,ip) - rof_budg_hstateG(s_hcontent_beg,ip))*1.0e-15_r8
            write(iulog,'(78("-"),"|",25("-"))')

            ! Energy balance check: dH/dt should equal net flux
            dH_dt = (rof_budg_hstateG(s_hcontent_end,ip) - rof_budg_hstateG(s_hcontent_beg,ip))
            residual = net_heat_flux*(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2) * get_step_size() - dH_dt

            write(iulog,*) ''
            write(iulog,*) 'ENERGY CONSERVATION CHECK:'
            write(iulog,*) '  Net flux * dt (PJ):  ', net_heat_flux*(4.0_r8*shr_const_pi*SHR_CONST_REARTH**2) * get_step_size() * 1.0e-15_r8
            write(iulog,*) '  Heat change dH (PJ): ', dH_dt * 1.0e-15_r8
            write(iulog,*) '  Residual (PJ):       ', residual * 1.0e-15_r8
            if (abs(residual) > 1.0e9_r8) then  ! 1 GJ threshold
               write(iulog,*) '  ***** WARNING: Energy budget residual exceeds threshold *****'
            endif
            write(iulog,*) ''
            write(iulog,*) 'EVAPORATION:'
            write(iulog,*) '  Total evap (kg/s):   ', rof_budg_hother(o_evap, p_inst)
            write(iulog,*) '  Total evap (Sv):     ', rof_budg_hother(o_evap, p_inst) * 1.0e-9_r8
            write(iulog,*) '=========================================='
            write(iulog,*)''
         end if
      end if
   end do

  end subroutine MOSART_HeatBudget_Print

end module MOSART_HeatBudgets_mod
