MODULE MOSARTinund_data_driven_MOD
!------------------------------------------------------------------------------!
! DESCRIPTION: MOSART inundation with linear log relationship.
!              FF = a * log(Vtotal) + b
! 
!------------------------------------------------------------------------------!

    use shr_kind_mod, only: r8 => shr_kind_r8
    use shr_sys_mod, only: shr_sys_abort
    use RunoffMod, only: rtmCTL, Tctl, TUnit, TRunoff

    implicit none
    private

    public inundation_run

    contains

    subroutine inundation_run(iunit)

        implicit none
        integer, intent(in)  :: iunit
        real(r8)             :: Vcri  ! Critial volume calculated based on a and b
        real(r8)             :: Vtotal

        if (TUnit%linear_a(iunit) > -999._r8) then
            Vcri = exp((TUnit%a_chnl(iunit) - TUnit%linear_b(iunit))/TUnit%linear_a(iunit))
            if (Vcri > TUnit%wr_bf(iunit)) then
                Vcri = TUnit%wr_bf(iunit)
            endif
            Vtotal = TRunoff%wr(iunit,1) + TRunoff%wf_ini(iunit)
            TRunoff%ff_unit(iunit) = calculate_flooded_fraction(Vtotal,TUnit%linear_a(iunit),TUnit%linear_b(iunit),TUnit%a_chnl(iunit))
            if (TRunoff%ff_unit(iunit) > TUnit%a_chnl(iunit) ) then
                TRunoff%ff_fp(iunit) = TRunoff%ff_unit(iunit) - TUnit%a_chnl(iunit) 
            else
                TRunoff%ff_fp(iunit) = 0._r8
            endif
        else
            TRunoff%ff_fp(iunit) = 0._r8
            TRunoff%ff_fp(iunit) = 0._r8
            Vcri = TUnit%wr_bf(iunit)
        endif

        if (TRunoff%ff_fp(iunit) > 0._r8) then
            call calculate_volume_exchange(TRunoff%wr(iunit,1),TRunoff%wf_ini(iunit), &
                                           TUnit%wr_bf(iunit),TRunoff%ff_unit(iunit), &
                                           TRunoff%ff_fp(iunit),TUnit%a_chnl(iunit), Vcri)
            TRunoff%yr(iunit,1) = TRunoff%wr(iunit,1) / TUnit%rwidth(iunit) / TUnit%rlen(iunit) ! water depth is updated here
        endif

        TRunoff%ff_ini(iunit)     = TRunoff%ff_fp(iunit)
        TRunoff%ffunit_ini(iunit) = TRunoff%ff_unit(iunit)

    end subroutine inundation_run

    real (r8) function calculate_flooded_fraction(Vtotal,a,b,f_channel)

        implicit none
        real(r8), intent(in) :: Vtotal ! total volume [m^{3}], Vchannel + Vfloodplain
        real(r8), intent(in) :: a, b, f_channel   ! FF = a * log(Vtotal) + b
        character( len = * ), parameter :: subname = '(calculate_flooded_fraction)'

        if (Vtotal < 1._r8) then
            calculate_flooded_fraction = 0._r8
        else
            calculate_flooded_fraction = a*log(Vtotal) + b
            if (calculate_flooded_fraction > 1._r8) then
                calculate_flooded_fraction= 1._r8
            elseif (calculate_flooded_fraction < f_channel) then
                calculate_flooded_fraction = f_channel ! minimum flooded fraction of a gridcell should be the channel fraction, not zero.
            endif
        endif

        return

    end function calculate_flooded_fraction

    subroutine calculate_volume_exchange(wr,wf,wr_bf,ff_unit,ff_fp,f_channel,Vcri)

        implicit none
        real(r8), intent(inout) :: wr, wf, ff_fp, ff_unit
        real(r8), intent(in)    :: wr_bf, f_channel
        real(r8)                :: w_over     ! channel storage + floodplain storage - channel storage capacity (m^3).
        real(r8)                :: hbf_excess ! bankfull excess water depth
        real(r8), intent(in)    :: Vcri   ! critical volume, FF = 0 when Vtotal <= Vcri

        character( len = * ), parameter :: subname = '(calculate_volume_exchange)'

        if (wr + wf > Vcri) then
            w_over = wr + wf - Vcri
            hbf_excess = w_over / ff_unit
            wr         = Vcri + hbf_excess * (ff_unit - ff_fp) ! It's possible that ff_fp = 0 from `calculate_flooded_fraction`, so w_over will stay in the channel
            wf         = hbf_excess * ff_fp
        else
            wr = wr + wf
            wf = 0._r8
            ff_fp = 0._r8 ! It's also possible that ff_fp > 0 from `calculate_flooded_fraction`, so we make them consistent
            ff_unit = f_channel
        endif

    end subroutine

end MODULE MOSARTinund_data_driven_MOD