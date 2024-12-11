
MODULE MOSARTinund_Core_MOD

!--------------------------------------------------------------------------------------
! DESCRIPTION: Core simulation of MOSART-Inundation.
! 
! HISTORY:
! 2011: Most subroutines/functions were initially developed.
! 2014-2016: The subroutines/functions were revised or created in offline MOSART-Inundation.
! 2017: Integrated with ACME.
! ... ... ...
!
!--------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_sys_mod, only: shr_sys_abort
  use RunoffMod, only: rtmCTL, Tctl, TUnit, TRunoff, &
      SMatP_upstrm, avsrc_upstrm, avdst_upstrm, SMatP_dnstrm
  use MOSARTinund_PreProcs_MOD, only: con1Em3
  use RtmVar, only: barrier_timers, iulog, inundflag
  use RtmSpmd, only: mpicom_rof, masterproc
  use mct_mod
  
  implicit none
  private
  
  public ManningEq, ChnlFPexchg
                   
  contains
  
  !function CRVRMAN_imp(slp_, n_, rr_) result(v_)
  real( r8 ) function ManningEq( slp_, n_, rr_ )
    ! DESCRIPTION: Estimate flow velocity with Manning equation.
    
    ! HISTORY:
    ! 2011: Initially developed (H.-Y. Li).
    ! 2014-2016: Revised in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none
    real( r8 ), intent(in) :: slp_    ! Slope (dimensionless).
    real( r8 ), intent(in) :: n_      ! Manning roughness coefficient ( s * m^(-1/3) ).
    real( r8 ), intent(in) :: rr_     ! Hydraulic radius (m).
    !real( r8 ), intent(out) :: v_    ! Flow velocity (m/s).
    real( r8 ) :: ftemp
    character( len = * ), parameter :: subname = '(ManningEq)'
    
    if ( rr_ .gt. 0._r8 ) then
    
      if ( slp_ .gt. 0._r8 ) then
        ftemp = 2._r8 / 3._r8
        ManningEq = ( rr_ ** ftemp ) * sqrt( slp_ ) / n_     ! Manning equation.
      elseif ( slp_ .lt. 0._r8 ) then
        ftemp = 2._r8 / 3._r8
        ManningEq = - ( rr_ ** ftemp ) * sqrt( - slp_ ) / n_ ! Manning equation.
      else      ! slp_ = 0
        ManningEq = 0._r8
      end if
      
    elseif ( rr_ .eq. 0._r8 ) then
      ManningEq = 0._r8
    else      ! rr_ < 0
      write( iulog, * ) trim( subname ) // ' ERROR: Hydraulic radius is negative !'
      call shr_sys_abort( trim( subname ) // ' ERROR: Hydraulic radius is negative !' )
    end if
    
    return
  end function ManningEq

  subroutine ChnlFPexchg ( )
    ! DESCRIPTION: Calculate water exchange between the main channel and floodplains (assuming instantaneous exchange).
    
    ! HISTORY:
    ! 2014-2016: Created and improved in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...

    ! Input:  
    !   + Channel storage and floodplain storage.
    ! Output: 
    !   + New channel storage and floodplain storage after exchange;
    !   + Flooded fraction and flooded area for the floodplains;
    !   + New channel water depth and floodplain max water depth (i.e., elevation difference between final water level (after exchange) and banktop);
    !   + Exchange amount.
    
    implicit none
    integer :: iu
    real( r8 ) :: wr_rcd    ! For recording channel storage before exchange (m^3).    
    real( r8 ) :: w_over    ! = channel storage + floodplain storage - channel storage capacity (m^3).
    integer :: j            ! Index.
    real( r8 ) :: d_s       ! The storage between the current water level and the below water level through the point " j " in the elevation profile (m^3).
    real( r8 ) :: d_e       ! Elevation difference between the current water level and the point " j " in the elevation profile (m).
    real( r8 ) :: hf        ! Floodplain max water depth, namely the elevation difference between the final water level and the banktop (i.e., channel bankfull water level) (m).
    real( r8 ) :: ff_unit   ! Flooded fraction in the computation unit (including channel area) (dimensionless).
    real( r8 ) :: wr_over   ! Channel water storage between the final water level and the banktop ( = channel storage - channel storage capacity ) (m^3).
    real( r8 ) :: Vcri      ! Critial volume to trigger inundation
    character( len = * ), parameter :: subname = '(ChnlFPexchg)'
    
    !$OMP PARALLEL DO PRIVATE(wr_rcd, w_over, j, d_s, d_e, hf, ff_unit, wr_over, Vcri) SCHEDULE(GUIDED)
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

        ! If the channel water level is higher than the floodplain water level, or on the contrary :
        if (TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) .gt. TRunoff%hf_ini(iu) + con1Em3 .or. &
            (TRunoff%hf_ini(iu) .gt. con1Em3 .and. TRunoff%hf_ini(iu) .gt. TRunoff%yr( iu, 1 ) - TUnit%rdepth(iu) + con1Em3)) then

          wr_rcd = TRunoff%wr( iu, 1 )
          Vcri = TUnit%wr_bf(iu)*TUnit%linear_a(iu)
          ! Calculate w_over ( = channel storage + floodplain storage - channel storage capacity ) :
          w_over = TRunoff%wr( iu, 1 ) + TRunoff%wf_ini( iu ) - Vcri
        
          ! if ( wr + wf > Vcri ) :
          if (w_over .gt. 0._r8) then
            ! ---------------------------------  
            ! Calculate: (1) Flooded fraction in the computation unit (including channel area); (2) The difference between final water level and banktop elevation :
            ! --------------------------------- 
            do j = 2, TUnit%npt_eprof3(iu) - 1    ! Note: j starts from 2 .
              if (TUnit%s_eprof3(iu, j) < w_over .and. w_over <= TUnit%s_eprof3(iu, j+1)) then
                ! Water volume above the level of point " j " :
                d_s = w_over - TUnit%s_eprof3(iu, j)

                ! --------------------------------- 
                ! The relationship between d_s and d_e is expressed as a quadratic equation: p3*(d_e)^2 + q3*(d_e) - (d_s) = 0.  
                ! TUnit%p3( :, : ) is the coefficient "p3" of this quadratic equation (unit: m). 
                ! --------------------------------- 
                if (TUnit%p3(iu, j) > 0._r8) then
                  
                  ! Calculate the elevation difference between the current water level and the point " j " ( i.e., the solution of the above quadratic equation ) :
                  d_e = ( - TUnit%q3(iu, j) + sqrt(TUnit%q3(iu, j)**2._r8 + 4._r8*TUnit%p3(iu, j)*d_s)) / (2._r8*TUnit%p3(iu, j))

                  ! Calculate the flooded fraction in the computation unit (including channel area) :
                  ff_unit = d_e * TUnit%alfa3(iu, j) + TUnit%a_eprof3(iu, j)

                ! if ( TUnit%p3(iu, j) = 0 )
                else

                  ! Calculate the elevation difference "d_e" ( = Storage difference / Flooded area ) :
                  !d_e = d_s / TUnit%area(iu) / TUnit%a_eprof3(iu, j)
                  d_e = d_s / ( TUnit%area(iu) * TUnit%frac(iu) * TUnit%a_eprof3(iu, j) )

                  ! Calculate the flooded fraction in the computation unit (including channel area) :
                  ff_unit = TUnit%a_eprof3(iu, j)
                endif

                ! The elevation difference between final water level and banktop :
                hf = TUnit%e_eprof3(iu, j) + d_e
                exit
              endif
            enddo
          
            ! Channel water storage between final water level and banktop :
            wr_over = hf * TUnit%rwidth( iu ) * TUnit%rlen( iu )            

            ! Channel storage after exchange :
            TRunoff%wr_exchg( iu ) = Vcri + wr_over        

            ! Floodplain storage after exchange :
            TRunoff%wf_exchg(iu) = w_over - wr_over              
          
            ! Ratio of flooded  area to computation unit area :
            TRunoff%ff_unit(iu) = ff_unit                                                                  
            ! Ratio of flooded floodplain area to computation unit area :
            TRunoff%ff_fp(iu) = ff_unit - TUnit%a_chnl( iu )              

            ! Area of flooded floodplain :
            !TRunoff%fa_fp(iu) = TUnit%area(iu) * TRunoff%ff_fp(iu)      
            TRunoff%fa_fp(iu) = TUnit%area(iu) * TUnit%frac(iu) * TRunoff%ff_fp(iu)

            ! Floodplain max water depth after exchange :
            TRunoff%hf_exchg(iu) = hf                             

            ! Channel water depth after exchange :
            TRunoff%yr_exchg( iu ) = TUnit%rdepth( iu ) + hf            
          
            if ( - con1Em3 < TRunoff%wf_exchg(iu) .and. TRunoff%wf_exchg(iu) < 0._r8) then  
              TRunoff%wf_exchg(iu) = 0._r8            
            elseif ( TRunoff%wf_exchg(iu) <= - con1Em3 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: Floodplain water volume is negative !'
              call shr_sys_abort( trim( subname ) // ' ERROR: Floodplain water volume is negative !' )
            endif   
          
          ! if ( wr + wf <= Vcri )
          else    
        
            ! The channel water volume after exchange ( all the water remains in the channel ) :
            TRunoff%wr_exchg(iu) = TRunoff%wr(iu,1) + TRunoff%wf_ini(iu)      

            ! The floodplain water volume after exchange is zero :
            TRunoff%wf_exchg(iu) = 0._r8                                  
          
            TRunoff%ff_fp(iu) = 0._r8
            TRunoff%fa_fp(iu) = 0._r8
            TRunoff%hf_exchg(iu) = 0._r8
            TRunoff%ff_unit(iu) = TUnit%a_chnl( iu )

            ! Channel water depth after exchange ( = water volume / channel area ) :
            TRunoff%yr_exchg( iu ) = TRunoff%wr_exchg( iu ) / TUnit%rwidth( iu ) / TUnit%rlen( iu )
          
          endif      ! if ( wr + wf > wr_bf )
        
          ! Channel--floodplain exchange amount (Positive: flow from channel to floodplain; vice versa) (m^3) :
          TRunoff%netchange(iu) = wr_rcd - TRunoff%wr_exchg(iu)
        
        ! ---------------------------------  
        ! No channel--floodplain exchange for two situations: 
        ! (1) Floodplain is inundated: channel water level equals floodplain water level; 
        ! (2) Floodplain is not inundated: channel water level is below banktop.
        ! ---------------------------------  
        else
          TRunoff%wr_exchg( iu ) = TRunoff%wr( iu, 1 )
          TRunoff%yr_exchg( iu ) = TRunoff%yr( iu, 1 )
          TRunoff%wf_exchg( iu ) = TRunoff%wf_ini( iu )
          TRunoff%hf_exchg( iu ) = TRunoff%hf_ini( iu )
          TRunoff%ff_fp( iu ) = TRunoff%ff_ini ( iu )
          TRunoff%fa_fp( iu ) = TUnit%area(iu) * TUnit%frac(iu) * TRunoff%ff_fp(iu)
          TRunoff%ff_unit( iu ) = TRunoff%ffunit_ini( iu )
          TRunoff%netchange( iu ) = 0._r8
        end if    ! if ( the channel water level is higher than the floodplain water level, or on the contrary )
      end if      ! if ( TUnit%mask( iu ) .gt. 0 )
    end do
    !$OMP END PARALLEL DO
    
  end subroutine ChnlFPexchg
  
!#endif
  
end MODULE MOSARTinund_Core_MOD
