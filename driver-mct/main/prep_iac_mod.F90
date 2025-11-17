module prep_iac_mod

#include "shr_assert.h"
  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: cxx => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_iac, num_inst_frc, num_inst_rof
  use seq_comm_mct,     only: CPLID, IACID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use shr_log_mod     , only: errMsg => shr_log_errMsg
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: iac, lnd, atm, rof

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_iac_init
  public :: prep_iac_mrg

  public :: prep_iac_accum
  public :: prep_iac_accum_avg

  public :: prep_iac_zero_max

  public :: prep_iac_calc_l2x_zx

  public :: prep_iac_get_l2x_zx
  public :: prep_iac_get_l2zacc_lx
  public :: prep_iac_get_l2zacc_lx_cnt
  public :: prep_iac_get_mapper_Sl2z
  public :: prep_iac_get_mapper_Sa2z

  ! MOSART-IAC coupling
  public :: prep_iac_accum_rof
  public :: prep_iac_accum_avg_rof
  public :: prep_iac_calc_r2x_zx
  public :: prep_iac_get_r2x_zx
  public :: prep_iac_get_r2zacc_rx
  public :: prep_iac_get_r2zacc_rx_cnt
  public :: prep_iac_get_mapper_Sr2z

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sl2z
  type(seq_map), pointer :: mapper_Sa2z
  type(seq_map), pointer :: mapper_Sr2z  ! MOSART to IAC mapper

  ! attribute vectors - land
  type(mct_aVect), pointer :: l2x_zx(:) ! Lnd export, iac grid, cpl pes

  ! accumulation variables - land
  type(mct_aVect), pointer :: l2zacc_lx(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: l2zacc_lx_cnt  ! l2zacc_lx: number of time samples accumulated

  ! This holds our max of monthly averages
  type(mct_aVect), pointer :: l2zmax_lx(:)   ! lnd export, lnd grid, cpl pes

  ! attribute vectors - MOSART
  type(mct_aVect), pointer :: r2x_zx(:) ! MOSART export, iac grid, cpl pes

  ! accumulation variables - MOSART
  type(mct_aVect), pointer :: r2zacc_rx(:)   ! MOSART export, rof grid, cpl pes
  integer        , target  :: r2zacc_rx_cnt  ! r2zacc_rx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID                            ! MPI cpl communicator

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_iac_init(infodata, lnd_c2_iac, rof_c2_iac)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: lnd_c2_iac ! .true.  => lnd to iac coupling on
    logical                 , intent(in)    :: rof_c2_iac ! .true.  => rof to iac coupling on
    !
    ! Local Variables
    integer                  :: lsize_z, lsize_l, lsize_r
    integer                  :: eli, erl, eri
    logical                  :: samegrid_lz   ! samegrid land and iac
    logical                  :: samegrid_az   ! samegrid atm and iac
    logical                  :: samegrid_rz   ! samegrid rof and iac
    logical                  :: lnd_present   ! .true. => land is present
    logical                  :: iac_present   ! .true. => iac is present
    logical                  :: atm_present   ! .true. => atm is present
    logical                  :: rof_present   ! .true. => rof is present
    logical                  :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                  :: esmf_map_flag ! .true. => use esmf for mapping
    character(CL)            :: lnd_gnam      ! lnd grid
    character(CL)            :: iac_gnam      ! iac grid
    character(CL)            :: atm_gnam      ! atm grid
    character(CL)            :: rof_gnam      ! rof grid
    type(mct_avect), pointer :: z2x_zx, x2z_zx, l2x_lx, r2x_rx
    character(*), parameter  :: subname = '(prep_iac_init)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------
    call seq_infodata_getData(infodata, &
         lnd_present=lnd_present,       &
         iac_present=iac_present,       &
         atm_present=atm_present,       &
         rof_present=rof_present,       &
         lnd_gnam=lnd_gnam,             &
         iac_gnam=iac_gnam,             &
         atm_gnam=atm_gnam,             &
         rof_gnam=rof_gnam,             &
         esmf_map_flag=esmf_map_flag)


    if (iac_present) then
       allocate(mapper_Sl2z)
       allocate(mapper_Sa2z)
       allocate(mapper_Sr2z)

       z2x_zx => component_get_c2x_cx(iac(1))
       lsize_z = mct_aVect_lsize(z2x_zx)

       allocate(l2x_zx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_avect_init(l2x_zx(eli), rList=seq_flds_l2x_fields, lsize=lsize_z)
          call mct_avect_zero(l2x_zx(eli))
       end do

       allocate(r2x_zx(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_avect_init(r2x_zx(eri), rList=seq_flds_r2x_fields, lsize=lsize_z)
          call mct_avect_zero(r2x_zx(eri))
       end do
    end if

    if (iac_present .and. lnd_present) then
       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       ! iac imports
       x2z_zx => component_get_x2c_cx(iac(1))
       lsize_z = mct_aVect_lsize(x2z_zx)

       ! lnd exports
       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       ! Now create into our accumulator avect (on the land grid,
       ! still), using shared fields in both import and export avects
       allocate(l2zacc_lx(num_inst_lnd))
       allocate(l2zmax_lx(num_inst_lnd))

       do erl = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2z_zx, l2zacc_lx(erl), lsize=lsize_l)
          call mct_aVect_zero(l2zacc_lx(erl))       
          call mct_aVect_initSharedFields(l2x_lx, x2z_zx, l2zmax_lx(erl), lsize=lsize_l)
          call mct_aVect_zero(l2zmax_lx(erl))       
       end do

       samegrid_lz = .true.
       if (trim(lnd_gnam) /= trim(iac_gnam)) samegrid_lz = .false.       

       if (lnd_c2_iac) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sl2z'
          end if

          call seq_map_init_rcfile(mapper_Sl2z, lnd(1), iac(1), &
               'seq_maps.rc','lnd2iac_smapname:','lnd2iac_smaptype:',samegrid_lz, &
               string='mapper_Sl2z initialization',esmf_map=esmf_map_flag)
       end if
       call shr_sys_flush(logunit)
    end if

    ! currently lnd and atm must be present (maybe not on same grid)
    !   if iac is present
    ! the atm2iac mapper is used for the iac domain settings
    ! there currently are no variables passed from atm2iac
    if (iac_present .and. atm_present) then
       samegrid_az = .true.
       if (trim(atm_gnam) /= trim(iac_gnam)) samegrid_az = .false.

       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Sa2z'
       end if
       
       call seq_map_init_rcfile(mapper_Sa2z, atm(1), iac(1), &
          'seq_maps.rc','atm2iac_smapname:','atm2iac_smaptype:',samegrid_az, &
           string='mapper_Sa2z initialization',esmf_map=esmf_map_flag)
       call shr_sys_flush(logunit)
    end if

    ! MOSART-IAC coupling initialization
    if (iac_present .and. rof_present) then
       if (.not. allocated(mpicom_CPLID)) then
          call seq_comm_getData(CPLID, &
               mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)
       end if

       ! iac imports
       x2z_zx => component_get_x2c_cx(iac(1))
       lsize_z = mct_aVect_lsize(x2z_zx)

       ! rof exports
       r2x_rx => component_get_c2x_cx(rof(1))
       lsize_r = mct_aVect_lsize(r2x_rx)

       ! Create accumulator avect on the rof grid
       allocate(r2zacc_rx(num_inst_rof))

       do eri = 1,num_inst_rof
          call mct_aVect_initSharedFields(r2x_rx, x2z_zx, r2zacc_rx(eri), lsize=lsize_r)
          call mct_aVect_zero(r2zacc_rx(eri))
       end do

       r2zacc_rx_cnt = 0

       samegrid_rz = .true.
       if (trim(rof_gnam) /= trim(iac_gnam)) samegrid_rz = .false.

       if (rof_c2_iac) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sr2z'
          end if

          call seq_map_init_rcfile(mapper_Sr2z, rof(1), iac(1), &
               'seq_maps.rc','rof2iac_smapname:','rof2iac_smaptype:',samegrid_rz, &
               string='mapper_Sr2z initialization',esmf_map=esmf_map_flag)
       end if
       call shr_sys_flush(logunit)
    end if

  end subroutine prep_iac_init

  !================================================================================================

  subroutine prep_iac_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate land input to iac
    ! This will lead to a yearly average of lnd values
    ! To take a max of monthly means, then we'll have to add in a check
    ! for end of month, then run the average, then finally compare with
    ! previous max for this year and store in yet another avect
    !
    ! Actually just use the annual average because max months don't
    !   correspond to each other
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    type(mct_aVect), pointer :: l2x_lx
    character(*), parameter  :: subname = '(prep_iac_accum)'
    !---------------------------------------------------------------
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       l2x_lx => component_get_c2x_cx(lnd(eli))
       if (l2zacc_lx_cnt == 0) then
          call mct_avect_copy(l2x_lx, l2zacc_lx(eli))
       else
          ! Add the current field values to the iac accumulator avect
          call mct_avect_accum(l2x_lx, l2zacc_lx(eli))
       endif
    end do
    l2zacc_lx_cnt = l2zacc_lx_cnt + 1
    call t_drvstopf (trim(timer))

  end subroutine prep_iac_accum

  !================================================================================================

  subroutine prep_iac_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of land input to iac component, and build
    ! the max of monthly means
    !
    ! Actually just use the annual average
    !    the timer has been changed to do this average yearly
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: ezi, eli
    character(*), parameter :: subname = '(prep_iac_accum_avg)'
    !---------------------------------------------------------------
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2zacc_lx_cnt > 1 ) then 
       do eli = 1,num_inst_lnd
          write(logunit,*) 'TRS eli = ', eli, num_inst_lnd

          call mct_avect_avg(l2zacc_lx(eli),l2zacc_lx_cnt)
          
          call mct_avect_info(4,l2zacc_lx(eli),istr='TRS l2zacc')
       end do
    endif
    l2zacc_lx_cnt = 0
    call t_drvstopf (trim(timer))
    

  end subroutine prep_iac_accum_avg

  !================================================================================================
 subroutine prep_iac_zero_max()
   !---------------------------------------------------------------
   ! Description
   ! Sets l2zmax_lx to zero, to be called at the start of the year
   ! after iac_run
   !
   ! Arguments
   !
   ! Local Variables
   integer :: erl
   character(*), parameter :: subname = '(prep_iac_zero_max)'
   !---------------------------------------------------------------
   do erl = 1,num_inst_lnd
      call mct_aVect_zero(l2zmax_lx(erl))       
   end do
 end subroutine prep_iac_zero_max

  !================================================================================================

  subroutine prep_iac_mrg(infodata, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge iac inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: ezi
    type(mct_aVect), pointer :: x2z_zx
    character(*), parameter  :: subname = '(prep_iac_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)
    do ezi = 1,num_inst_iac
       x2z_zx => component_get_x2c_cx(iac(ezi))  ! This is actually modifying x2z_zx
       call prep_iac_merge(l2x_zx(ezi), x2z_zx)
    end do
    call t_drvstopf (trim(timer_mrg))

  end subroutine prep_iac_mrg

  !================================================================================================

  subroutine prep_iac_merge(l2x_z, x2z_z)
    !============================================================================
    ! Description
    ! Merge lnd for iac input
    !
    ! Arguments
    type(mct_aVect),intent(in)    :: l2x_z
    type(mct_aVect),intent(inout) :: x2z_z
    !
    ! Local variables
    integer       :: i, n, mdx
    integer       :: nflds,lsize
    integer       :: index_l2x, index_x2z
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL) :: field        ! field string
    character(CL),allocatable :: mrgstr(:)   ! temporary string, for logging
    character(*), parameter   :: subname = '(prep_iac_merge) '

    !-------------------------------------------------
    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2z_z)

    ! Loop over field indices
    nflds = shr_string_listGetNum(trim(seq_flds_x2z_states))

    allocate(mrgstr(nflds))
    mdx=1

    do i = 1, nflds
       ! Field i name
       call seq_flds_getField(field,i,seq_flds_x2z_states)
       ! Field i indices in input and output avects
       index_l2x = mct_aVect_indexRA(l2x_z, trim(field))
       index_x2z = mct_aVect_indexRA(x2z_z, trim(field))
       if (first_time) then
          mrgstr(mdx) = subname//'x2z%'//trim(field)//' =' //&
               ' = l2z%'//trim(field)
       endif

       ! l2z fields are all state fields, so it's a straight copy.
       ! It seems like mct_aVect_copy would work here...
       do n = 1, lsize
          x2z_z%rAttr(index_x2z,n) = l2x_z%rAttr(index_l2x,n)
       end do
       mdx = mdx+1
    end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,nflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.
  end subroutine prep_iac_merge

  !================================================================================================

  subroutine prep_iac_calc_l2x_zx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_zx (note that l2x_zx is a local module variable)
    !
    ! Arguments
    ! Don't know if we need these fractions just yet
    ! type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables

    integer :: eli, ezi
    type(mct_aVect), pointer :: l2x_lx
    character(*), parameter :: subname = '(prep_lnd_calc_l2x_zx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       ezi = mod((eli-1), num_inst_iac) + 1
       !l2x_lx => component_get_c2x_cx(lnd(eli))
       !call seq_map_map(mapper_Sl2z, l2_lx, l2x_zx(eli), fldlist=seq_flds_l2x_states, norm=.true.)
       !Need to map the max vector l2zmax_lx, which is local
       !call seq_map_map(mapper_Sl2z, l2zmax_lx(eli), l2x_zx(ezi), fldlist=seq_flds_l2x_states, norm=.true.)
       !Actually just use the annual average
       call seq_map_map(mapper_Sl2z, l2zacc_lx(eli), l2x_zx(ezi),&
          fldlist=seq_flds_l2x_states, norm=.true.)
    enddo

    ! Now zero out max avects for the next year of lnd runs
    call prep_iac_zero_max()

    call t_drvstopf  (trim(timer))

  end subroutine prep_iac_calc_l2x_zx

  !================================================================================================

  function prep_iac_get_l2x_zx()
    type(mct_aVect), pointer :: prep_iac_get_l2x_zx(:)
    prep_iac_get_l2x_zx => l2x_zx(:)
  end function prep_iac_get_l2x_zx

  function prep_iac_get_l2zacc_lx()
    type(mct_aVect), pointer :: prep_iac_get_l2zacc_lx(:)
    prep_iac_get_l2zacc_lx => l2zacc_lx(:)
  end function prep_iac_get_l2zacc_lx

  function prep_iac_get_l2zacc_lx_cnt()
    integer, pointer :: prep_iac_get_l2zacc_lx_cnt
    prep_iac_get_l2zacc_lx_cnt => l2zacc_lx_cnt
  end function prep_iac_get_l2zacc_lx_cnt

  function prep_iac_get_mapper_Sl2z()
    type(seq_map), pointer :: prep_iac_get_mapper_Sl2z
    prep_iac_get_mapper_Sl2z => mapper_Sl2z
  end function prep_iac_get_mapper_Sl2z

  function prep_iac_get_mapper_Sa2z()
    type(seq_map), pointer :: prep_iac_get_mapper_Sa2z
    prep_iac_get_mapper_Sa2z => mapper_Sa2z
  end function prep_iac_get_mapper_Sa2z

  !================================================================================================
  ! MOSART-IAC coupling functions
  !================================================================================================

  subroutine prep_iac_accum_rof(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate MOSART input to IAC
    ! This will lead to a yearly average of MOSART values
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_aVect), pointer :: r2x_rx
    character(*), parameter  :: subname = '(prep_iac_accum_rof)'
    !---------------------------------------------------------------
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))
       if (r2zacc_rx_cnt == 0) then
          call mct_avect_copy(r2x_rx, r2zacc_rx(eri))
       else
          ! Add the current field values to the iac accumulator avect
          call mct_avect_accum(r2x_rx, r2zacc_rx(eri))
       endif
    end do
    r2zacc_rx_cnt = r2zacc_rx_cnt + 1
    call t_drvstopf (trim(timer))

  end subroutine prep_iac_accum_rof

  !================================================================================================

  subroutine prep_iac_accum_avg_rof(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of MOSART input to IAC component
    ! Calculate annual average
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: ezi, eri
    character(*), parameter :: subname = '(prep_iac_accum_avg_rof)'
    !---------------------------------------------------------------
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (r2zacc_rx_cnt > 1 ) then
       do eri = 1,num_inst_rof
          call mct_avect_avg(r2zacc_rx(eri),r2zacc_rx_cnt)
          call mct_avect_info(4,r2zacc_rx(eri),istr='TRS r2zacc')
       end do
    endif
    r2zacc_rx_cnt = 0
    call t_drvstopf (trim(timer))

  end subroutine prep_iac_accum_avg_rof

  !================================================================================================

  subroutine prep_iac_calc_r2x_zx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_zx (note that r2x_zx is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri, ezi
    character(*), parameter :: subname = '(prep_iac_calc_r2x_zx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       ezi = mod((eri-1), num_inst_iac) + 1
       ! Map the accumulated MOSART data to IAC grid
       call seq_map_map(mapper_Sr2z, r2zacc_rx(eri), r2x_zx(ezi),&
          fldlist=seq_flds_r2z_states, norm=.true.)
    enddo

    call t_drvstopf  (trim(timer))

  end subroutine prep_iac_calc_r2x_zx

  !================================================================================================

  function prep_iac_get_r2x_zx()
    type(mct_aVect), pointer :: prep_iac_get_r2x_zx(:)
    prep_iac_get_r2x_zx => r2x_zx(:)
  end function prep_iac_get_r2x_zx

  function prep_iac_get_r2zacc_rx()
    type(mct_aVect), pointer :: prep_iac_get_r2zacc_rx(:)
    prep_iac_get_r2zacc_rx => r2zacc_rx(:)
  end function prep_iac_get_r2zacc_rx

  function prep_iac_get_r2zacc_rx_cnt()
    integer, pointer :: prep_iac_get_r2zacc_rx_cnt
    prep_iac_get_r2zacc_rx_cnt => r2zacc_rx_cnt
  end function prep_iac_get_r2zacc_rx_cnt

  function prep_iac_get_mapper_Sr2z()
    type(seq_map), pointer :: prep_iac_get_mapper_Sr2z
    prep_iac_get_mapper_Sr2z => mapper_Sr2z
  end function prep_iac_get_mapper_Sr2z

end module prep_iac_mod
