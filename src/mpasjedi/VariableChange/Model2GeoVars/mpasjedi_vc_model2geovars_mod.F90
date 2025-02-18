! (C) Copyright 2020-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_vc_model2geovars_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log

!oops
use kinds, only : kind_real
use oops_variables_mod, only: oops_variables

!ufo
use gnssro_mod_transform, only: geometric2geop
use ufo_vars_mod

!MPAS-Model
use mpas_derived_types
use mpas_field_routines
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
use mpas_geom_mod, only: mpas_geom
use mpas4da_mod, only: da_template_pool

!TODO: package the variable changes somewhere
use mpas2ufo_vars_mod
!use mpas_vc_utils

implicit none

private
public :: mpasjedi_vc_model2geovars

type :: mpasjedi_vc_model2geovars
  character(len=10) :: tropprs_method
 contains
  procedure, public :: create
  procedure, public :: delete
  procedure, public :: changevar
end type mpasjedi_vc_model2geovars

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(mpasjedi_vc_model2geovars), intent(inout) :: self
type(mpas_geom),                  intent(in)    :: geom
type(fckit_configuration),        intent(in)    :: conf

character(len=10) :: str

! Method to use for tropopause pressure ([gsi] or thompson)
!if (.not. conf%get("tropopause pressure method", str)) str = 'thompson'

self%tropprs_method = 'thompson'

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(mpasjedi_vc_model2geovars), intent(inout) :: self

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self, geom, xm, xg)

  class(mpasjedi_vc_model2geovars), intent(inout) :: self
  type(mpas_geom),                  intent(in)    :: geom !< mpas mesh descriptors
  class(mpas_fields),               intent(in)    :: xm   !< model state fields
  class(mpas_fields),               intent(inout) :: xg   !< state containing geovar fields

  ! pool-related pointers
  type(mpas_pool_type), pointer :: mFields
  type (mpas_pool_type), pointer :: sfc_input
  type(mpas_pool_data_type), pointer :: mdata, gdata

  ! reusable fields
  type(field2DReal), pointer :: fieldr2_a, fieldr2_b

  real(kind=RKIND), dimension(:), allocatable :: uu, vv, f10m

  ! reusable arrays
  real(kind=RKIND), dimension(:), pointer :: ptrr1_a, ptrr1_b
  real(kind=RKIND), dimension(:,:), pointer :: ptrr2_a, ptrr2_b
  real(kind=RKIND), dimension(:,:), allocatable :: r2_a, r2_b

  ! iteration-specific variables
  character(len=MAXVARLEN) :: geovar
  integer :: nCells, nVertLevels, nVertLevelsP1
  integer :: iVar, iCell, iLevel, i
  real (kind=RKIND) :: lat
  real (kind=kind_real) :: rz

  ! surface variables
  character(len=MAXVARLEN), parameter :: &
    MPASSfcNames(5) = &
      [ character(len=MAXVARLEN) :: 'ivgtyp', 'isltyp', 'landmask', 'seaice_fraction', 'snowc']
  character(len=MAXVARLEN), parameter :: &
    ValidCRTMSfcNames(8) = &
      [var_sfc_landtyp_usgs, var_sfc_landtyp_igbp, &
       var_sfc_vegtyp, var_sfc_soiltyp, &
       var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac]
  type(oops_variables) :: RequestedCRTMSfcNames
  type(mpas_pool_type), pointer :: RequestedCRTMSfcFields => null()
  character(len=MAXVARLEN), parameter :: &
    AllCRTMLandTypeNames(3) = &
      [var_sfc_landtyp_usgs, var_sfc_landtyp_igbp, var_sfc_landtyp_npoess]
  integer, dimension(:), pointer :: vegtyp, soiltyp
  real(kind=RKIND), dimension(:), pointer :: wfrac, lfrac, ifrac, sfrac

  ! air pressure on w levels
  real(kind=RKIND), allocatable :: plevels(:,:)

  integer, dimension(:), pointer :: domainMask

  ! config members
  character(len=StrKIND), pointer :: &
    config_microp_scheme, config_radt_cld_scheme, mminlu
  logical, pointer :: config_microp_re

  ! convenient local variables
  mFields => xm % subFields
  nCells = geom%nCellsSolve
  nVertLevels = geom%nVertLevels
  nVertLevelsP1 = geom%nVertLevelsP1

  ! CRTM surface type and fraction fields are interdependent
  ! pre-calculate all such fields if any are requested
  if ( any(xg%has(ValidCRTMSfcNames)) ) then

    if (.not. all(xm%has(MPASSfcNames))) then
      call abor1_ftn('mpasjedi_vc_model2geovars::changevar: xm must include MPASSfcNames!')
    end if

    ! land type check
    if (count(xg%has(AllCRTMLandTypeNames)) > 1) then
      call abor1_ftn('mpasjedi_vc_model2geovars::changevar: '&
                    &'can only produce one CRTM land surface classification at a time')
    end if

    ! mminlu: MPAS land surface classification system
    call mpas_pool_get_subpool(geom % domain % blocklist % structs, 'sfc_input', sfc_input)
    call mpas_pool_get_array(sfc_input, 'mminlu', mminlu)
    ! possible values (see VEGPARM.TBL): "USGS", "MODIFIED_IGBP_MODIS_NOAH", "NLCD40", "USGS-RUC", "MODI-RUC"

    RequestedCRTMSfcNames = oops_variables()
    do iVar = 1, size(ValidCRTMSfcNames)
      geovar = ValidCRTMSfcNames(ivar)
      if (xg%has(geovar)) then
        if (geovar == var_sfc_landtyp_usgs .and. trim(mminlu) /= 'USGS') then
          ! fail if MPAS ivgtyp is inconsistent with var_sfc_landtyp_usgs
          call abor1_ftn('mpasjedi_vc_model2geovars::changevar: '&
                        &'the "USGS" CRTM land surface classification '&
                        &'for CRTM ObsOperators associated with infrared instruments '&
                        &'(i.e., obs operator.options.IRVISlandCoeff) '&
                        &'can only be used with mminlu=USGS, not mminlu='//trim(mminlu))
        end if
        if (geovar == var_sfc_landtyp_igbp .and. trim(mminlu) /= 'MODIFIED_IGBP_MODIS_NOAH') then
          ! fail if MPAS ivgtyp is inconsistent with var_sfc_landtyp_igbp
          call abor1_ftn('mpasjedi_vc_model2geovars::changevar: '&
                        &'the "IGBP" CRTM land surface classification '&
                        &'for CRTM ObsOperators associated with infrared instruments '&
                        &'(i.e., obs operator.options.IRVISlandCoeff) '&
                        &'can only be used with mminlu=MODIFIED_IGBP_MODIS_NOAH, not mminlu='//trim(mminlu))
        end if
        call RequestedCRTMSfcNames%push_back(geovar)
      end if
    end do

    call da_template_pool(geom, &
                          RequestedCRTMSfcFields, &
                          RequestedCRTMSfcNames%nvars(), &
                          RequestedCRTMSfcNames%varlist())

    !! surface types
    ! land type
    if (RequestedCRTMSfcNames%has(var_sfc_landtyp_usgs)) then
      call xm%copy_to('ivgtyp', RequestedCRTMSfcFields, var_sfc_landtyp_usgs)
    elseif (RequestedCRTMSfcNames%has(var_sfc_landtyp_igbp)) then
      call xm%copy_to('ivgtyp', RequestedCRTMSfcFields, var_sfc_landtyp_igbp)
    end if

    ! veg type
    ! uses ivgtyp as input
    if (RequestedCRTMSfcNames%has(var_sfc_vegtyp)) then
      call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_vegtyp, vegtyp)
      call xm%get('ivgtyp', mdata)

      select case (trim(mminlu))
        case ('USGS')
          do iCell = 1, nCells
            vegtyp(iCell) = convert_type_veg_usgs(mdata%i1%array(iCell))
          end do
        case ('MODIFIED_IGBP_MODIS_NOAH')
          do iCell = 1, nCells
            vegtyp(iCell) = convert_type_veg_igbp(mdata%i1%array(iCell))
          end do
        case default
          call abor1_ftn('mpasjedi_vc_model2geovars::changevar: invalid mminlu, must be one of'&
                        &'[USGS, MODIFIED_IGBP_MODIS_NOAH]')
      end select
    end if

    ! soil type
    if (RequestedCRTMSfcNames%has(var_sfc_soiltyp)) then
      call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_soiltyp, soiltyp)
      call xm%get('isltyp', mdata)
      do iCell = 1, nCells
        soiltyp(iCell) = convert_type_soil(mdata%i1%array(iCell))
      end do
    end if

    ! For now, destruct RequestedCRTMSfcNames here, assuming that the ValidCRTMSfcNames referenced below are
    ! requested in "xg". RequestedCRTMSfcNames could be used below for more error/consistency checking if UFO
    ! GeoVar requests become more heterogeneous.
    call RequestedCRTMSfcNames%destruct()

    !! surface fractions
    ! land, will be adjusted later
    ! TODO: a binary landmask is a crude indicator of sub-grid land fraction
    call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_lfrac, lfrac)
    call xm%get('landmask', mdata) !'land-ocean mask (1=>land ; 0=>ocean)'
    lfrac(1:nCells) = real(mdata%i1%array(1:nCells), RKIND)

    ! ice
    call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_ifrac, ifrac)
    call xm%get('seaice_fraction', mdata)     !'fractional area coverage of sea-ice'
    ifrac(1:nCells) = mdata%r1%array(1:nCells)

    ! snow
    ! TODO: Investigate. snowc varies between 0. and 1., but Registry description indicates it is binary.
    !       Similar comment to landmask.
    call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_sfrac, sfrac)
    call xm%get('snowc', mdata)    !'flag for snow on ground (0=>no snow; 1=>otherwise)'
    sfrac(1:nCells) = mdata%r1%array(1:nCells)

    ! water+snow+land
    call mpas_pool_get_array(RequestedCRTMSfcFields, var_sfc_wfrac, wfrac)
    do iCell = 1, nCells
      if (ifrac(iCell) > MPAS_JEDI_ZERO_kr) then
        sfrac(iCell) = max(min(sfrac(iCell), MPAS_JEDI_ONE_kr - ifrac(iCell)), MPAS_JEDI_ZERO_kr)
        wfrac(iCell) = max(MPAS_JEDI_ONE_kr - ifrac(iCell) - sfrac(iCell), MPAS_JEDI_ZERO_kr)
      else if (sfrac(iCell) > MPAS_JEDI_ZERO_kr) then
        wfrac(iCell) = max(MPAS_JEDI_ONE_kr - sfrac(iCell), MPAS_JEDI_ZERO_kr)
      else
        wfrac(iCell) = max(MPAS_JEDI_ONE_kr - lfrac(iCell), MPAS_JEDI_ZERO_kr)
        ifrac(iCell) = MPAS_JEDI_ZERO_kr
        sfrac(iCell) = MPAS_JEDI_ZERO_kr
      end if
      lfrac(iCell) = max(MPAS_JEDI_ONE_kr - ifrac(iCell) - sfrac(iCell) - wfrac(iCell), MPAS_JEDI_ZERO_kr)
    end do
  end if


  ! pre-calculate pressure on w levels
  allocate(plevels(1:nVertLevelsP1,1:nCells))
  call xm%get('air_pressure', ptrr2_a)
  call xm%get('air_pressure_at_surface', ptrr1_a)
  call pressure_half_to_full(ptrr2_a(:,1:nCells), geom%zgrid(:,1:nCells), ptrr1_a(1:nCells), &
                             nCells, nVertLevels, plevels)

  ! populate xg
  do iVar = 1, xg % nf

    geovar = trim(xg % fldnames(iVar))

    if (geom%has_identity(geovar)) then
      ! "identity" variable changes are controlled in an external configuration file
      ! through the geom object. Refer to the mpas_geom type for more information.
      if (xm%has(geom%identity(geovar))) then
        call xm%copy_to(geom%identity(geovar), xg, geovar)
      else
        call abor1_ftn('mpasjedi_vc_model2geovars::changevar: '&
                      &'state missing identity field for geovar => '//trim(geovar))
      end if
    else

      call xg%get(geovar, gdata)

      select case (trim(geovar))

        case ( var_tv ) !-virtual_temperature
          call xm%get('air_temperature', ptrr2_a)
          call xm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_b)
          allocate(r2_a(1:nVertLevels, 1:nCells))
          call q_to_w( ptrr2_b(:,1:nCells), r2_a(:,1:nCells) )
          call tw_to_tv( ptrr2_a(:,1:nCells), r2_a(:,1:nCells), gdata%r2%array(:,1:nCells) )
          deallocate(r2_a)

        case ( var_mixr ) !-water_vapor_mixing_ratio_wrt_dry_air
          call xm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_a)
          call q_to_w(ptrr2_a(:,1:nCells), gdata%r2%array(:,1:nCells))

          ! Ensure positive-definite mixing ratios
          !  with respect to precision of crtm::CRTM_Parameters::ZERO.
          ! TODO: this should be moved to the mpas_fields%read step
          !       only other place it should be needed is add_incr (already there)
          where(gdata%r2%array(:,1:nCells) < MPAS_JEDI_ZERO_kr)
            gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr
          end where

        case ( var_airdens ) ! moist_air_density
          call xm%get('dry_air_density', ptrr2_a) ! dry air density, kg/m^3
          call xm%get('water_vapor_mixing_ratio_wrt_dry_air',  ptrr2_b) ! water vapor mixing ratio, kg/kg
          call dryrho_to_moistrho(ptrr2_a, ptrr2_b, nCells, nVertLevels)
          gdata%r2%array(:,1:nCells) = ptrr2_a(:,1:nCells)

        case ( var_prsi ) !-air_pressure_levels
          gdata%r2%array(1:nVertLevelsP1,1:nCells) = plevels(1:nVertLevelsP1,1:nCells)

        case ( var_w ) !-upward_air_velocity
          call xm%get('w', ptrr2_a)
          call geom%full_to_half(ptrr2_a(:,1:nCells), gdata%r2%array(:,1:nCells), nCells)

        case ( var_oz ) !-mole_fraction_of_ozone_in_air :TODO: not directly available from MPAS
          !call xm%get('mole_fraction_of_ozone_in_air', mdata)
          gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr !mdata%r2%array(:,1:nCells)

        case ( var_co2 ) !-mole_fraction_of_carbon_dioxide_in_air :TODO: not directly available from MPAS
          !call xm%get('mole_fraction_of_carbon_dioxide_in_air', mdata)
          gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr !mdata%r2%array(:,1:nCells)

        case ( var_clw_wp ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
          call q_fields_forward('cloud_liquid_water', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_cli_wp ) !-mass_content_of_cloud_ice_in_atmosphere_layer
          call q_fields_forward('cloud_liquid_ice', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clr_wp ) !-mass_content_of_rain_in_atmosphere_layer
          call q_fields_forward('rain_water', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_cls_wp ) !-mass_content_of_snow_in_atmosphere_layer
          call q_fields_forward('snow_water', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clg_wp ) !-mass_content_of_graupel_in_atmosphere_layer
          call q_fields_forward('graupel', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clh_wp ) !-mass_content_of_hail_in_atmosphere_layer
          call q_fields_forward('hail', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clwefr ) !-effective_radius_of_cloud water particle
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            call xm%get('re_cloud', ptrr2_a) !- [m]
            gdata%r2%array(:,1:nCells) = ptrr2_a(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
          else
            gdata%r2%array(:,1:nCells) = 10.0_kind_real ![micron]  ! Or, we can refer to rre_cloud which is default for MPAS radiance scheme
          end if

        case ( var_cliefr ) !-effective_radius_of_cloud ice particle
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            call xm%get('re_ice', ptrr2_a) !- [m]
            gdata%r2%array(:,1:nCells) = ptrr2_a(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
          else
            gdata%r2%array(:,1:nCells) = 30.0_kind_real ! [micron]
          end if

        case ( var_clrefr ) !-effective_radius_of_rain water_particle
          ! effective_radius_of_rain_water is not calculated in MPAS model physics
          call xm%get('rain_water', fieldr2_a)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            allocate(r2_a(1:nVertLevels, 1:nCells))
            allocate(r2_b(1:nVertLevels, 1:nCells))
            if (trim(config_microp_scheme) == 'mp_thompson') then
               call xm%get('rain_number_concentration', ptrr2_a) !- [nb kg^{-1}]: MPAS output for 2-moment MP scheme
               r2_b(:,1:nCells) = ptrr2_a(:,1:nCells)
            else
               r2_b = MPAS_JEDI_ONE_kr
            end if
            call xm%get('dry_air_density', fieldr2_b) !- [kg m^{-3}]: Dry air density

            call effectRad_rainwater(fieldr2_a%array(:,1:nCells), fieldr2_b%array(:,1:nCells),&
                                     r2_b(:,1:nCells), r2_a(:,1:nCells), config_microp_scheme, &
                                     nCells, nVertLevels)
            gdata%r2%array(:,1:nCells) = r2_a(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
            deallocate(r2_a)
            deallocate(r2_b)
          else
            gdata%r2%array(:,1:nCells) = 999.0_kind_real ! [micron]
          end if

        case ( var_clsefr ) !-effective_radius_of_snow particle
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            call xm%get('re_snow', mdata) !- [m]
            gdata%r2%array(:,1:nCells) = mdata%r2%array(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
          else
            gdata%r2%array(:,1:nCells) = 600.0_kind_real ! [micron]
          end if

        case ( var_clgefr ) !-effective_radius_of_graupel_particle
          !effective_radius_of_graupel is not calculated in MPAS model physics
          call xm%get('graupel', fieldr2_a)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            allocate(r2_a(1:nVertLevels, 1:nCells))
            call xm%get('dry_air_density', fieldr2_b) !- [kg m^{-3}]: Dry air density
            call effectRad_graupel(fieldr2_a%array(:,1:nCells), fieldr2_b%array(:,1:nCells), &
                                   r2_a(:,1:nCells), config_microp_scheme, &
                                   nCells, nVertLevels)
            gdata%r2%array(:,1:nCells) = r2_a(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
            deallocate(r2_a)
          else
            gdata%r2%array(:,1:nCells) = 600.0_kind_real ! [micron]
          end if

        case ( var_clhefr ) !-effective_radius_of_hail_particle :TODO: currently filled w/ default value
          !The current MP schemes do not include hail (wsm7 has)
          !call xm%get('re_hail', mdata)
          !gdata%r2%array(:,1:nCells) = mdata%r2%array(:,1:nCells) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
          gdata%r2%array(:,1:nCells) = 600.0_kind_real ! [micron]

        case ( var_cldfrac ) !-cloud_area_fraction_in_atmosphere_layer
          call mpas_pool_get_config(geom % domain % blocklist % configs, &
                                   'config_radt_cld_scheme', config_radt_cld_scheme)
          call mpas_pool_get_config(geom % domain % blocklist % configs, &
                                   'config_microp_scheme', config_microp_scheme)

          if ( trim(config_radt_cld_scheme) == MPAS_JEDI_OFF .or. &
              trim(config_microp_scheme) == MPAS_JEDI_OFF ) then
            gdata%r2%array(:,1:nCells) = MPAS_JEDI_LESSONE_kr
          else if (xm%has('cldfrac')) then
            call xm%copy_to('cldfrac', xg, geovar)
          else
            call abor1_ftn('mpasjedi_vc_model2geovars::changevar: cldfrac must be added to the state &
              & variables in order to populate the var_cldfrac geovar with the MPAS diagnostic cloud &
              & fraction')
          end if

        case ( var_z ) !-geopotential_height, geopotential heights at midpoint
          ! calculate midpoint geometricZ (unit: m):
          gdata%r2%array(:,1:nCells)=geom%height(:,1:nCells)

        case ( var_zi ) !-geopotential_height_levels, geopotential heights at w levels
          gdata%r2%array(:,1:nCells)=geom%zgrid(:,1:nCells)

        case ( var_geomz ) !-height_above_mean_sea_level
          ! calculate midpoint geometricZ (unit: m):
          gdata%r2%array(:,1:nCells) = geom%height(:,1:nCells)

        case ( var_tropprs ) !-tropopause pressure
          call xm%get('air_pressure',    ptrr2_a)
          call xm%get('air_temperature', ptrr2_b)
          if (trim(self%tropprs_method) == "thompson") then
            call tropopause_pressure_th(ptrr2_a(:,1:nCells), geom%zgrid(:,1:nCells), ptrr2_b(:,1:nCells), &
                                        nCells, nVertLevels, gdata%r1%array(1:nCells))
          !elseif (trim(self%tropprs_method) == "wmo") then
          !  call tropopause_pressure_wmo(ptrr2_a(:,1:nCells), geom%zgrid(:,1:nCells), ptrr2_b(:,1:nCells), &
          !                               nCells, nVertLevels, gdata%r1%array(1:nCells))
          else
            call abor1_ftn('mpasjedi_vc_model2geovars::changevar: invalid tropopause pressure determination method, &
                           & must be one of [thompson]')
          endif

!! begin surface variables
        case ( var_sfc_z ) !-geopotential_height_at_surface
          gdata%r1%array(1:nCells) = geom%zgrid(1,1:nCells)

        case ( var_sfc_geomz ) !-height_above_mean_sea_level_at_surface
          gdata%r1%array(1:nCells) = geom%zgrid(1,1:nCells)

        case ( var_sfc_sdepth ) !-surface_snow_thickness
          call xm%get('snowh', mdata)
          gdata%r1%array(1:nCells) = mdata%r1%array(1:nCells) * MPAS_JEDI_THOUSAND_kr ! [m] -> [mm]

!! NOTE: With introducing "mpas io scaling factor", this is handled with "identity" VC.
!  If users comment out "mpas identity field" and "mpas io scaling factor" 
!  of "field name: vegetation_area_fraction" in the geovars.yaml,
!  following case-statement will be activated.
!  TODO: This should be removed in future.
        case ( var_sfc_vegfrac ) !-vegetation_area_fraction
          call xm%get('vegetation_area_fraction', mdata)
          gdata%r1%array(1:nCells) = mdata%r1%array(1:nCells) / 100.0_kind_real ! [unitless, 0~100] = [%, 0~1]

        case ( var_sfc_soilm ) !-volume_fraction_of_condensed_water_in_soil
          ! NOTE: use 1st level
          ! NOTE: units are equivalent with water density 1 g/cm3
          call xm%get('smois', mdata)
          gdata%r1%array(1:nCells) = mdata%r2%array(1,1:nCells) ! [m3/m3] -> [g/cm3]

        case ( var_sfc_soilt ) !-soil_temperature
          ! NOTE: use 1st level
          call xm%get('tslb', mdata)
          gdata%r1%array(1:nCells) = mdata%r2%array(1,1:nCells) ! [K] -> [K]

        case ( var_sfc_landtyp_usgs, var_sfc_landtyp_igbp, &
               var_sfc_vegtyp, var_sfc_soiltyp, &
               var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac )
          ! RequestedCRTMSfcNames
          !-land_type_index, vegetation_type_index, soil_type
          !-water_area_fraction, land_area_fraction, ice_area_fraction, surface_snow_area_fraction
          call xg%copy_from(geovar, RequestedCRTMSfcFields)

        case ( var_sfc_wspeed ) !-wind_speed_at_surface
          call xm%get('eastward_wind_at_10m', ptrr1_a)
          call xm%get('northward_wind_at_10m', ptrr1_b)
          gdata%r1%array(1:nCells)=sqrt( ptrr1_a(1:nCells)**2 + ptrr1_b(1:nCells)**2 )

        case ( var_observable_domain_mask ) !-domain check
          call mpas_pool_get_array(geom%domain%blocklist%allFields, 'bdyMaskCell', domainMask)
! pass only the domain interior points.
          gdata%r1%array(1:nCells)= real(domainMask(1:nCells))

        case ( var_sfc_fact10 ) !-wind_reduction_factor_at_10m
          call xm%get('eastward_wind_at_10m', ptrr1_a)
          call xm%get('northward_wind_at_10m', ptrr1_b)

! get model wind at surface layer
          call xm%get('eastward_wind', mdata)
          allocate(uu(1:nCells))
          uu(1:nCells) = mdata%r2%array(1, 1:nCells)

          call xm%get('northward_wind', mdata)
          allocate(vv(1:nCells))
          vv(1:nCells) = mdata%r2%array(1, 1:nCells)
! get model wind at surface layer

          allocate(f10m(1:nCells))
          f10m(1:nCells) = sqrt( ptrr1_a(1:nCells)**2 + ptrr1_b(1:nCells)**2 )

          do i=1, nCells
           if (f10m(i) > 0.0_kind_real) then
               f10m(i) = f10m(i)/ sqrt( uu(i)**2 + vv(i)**2 )
           else
               f10m(i) = 1.0_kind_real
           end if
          end do

          gdata%r1%array(1:nCells) = f10m(1:nCells)

          deallocate(uu)
          deallocate(vv)

!! end surface variables
        case default
          call abor1_ftn('mpasjedi_vc_model2geovars::changevar: geovar not implemented => '//trim(geovar))

      end select

    end if

  end do !iVar

  if (associated(RequestedCRTMSfcFields)) then
    call mpas_pool_destroy_pool(RequestedCRTMSfcFields)
  end if

  deallocate(plevels)

end subroutine changevar

! --------------------------------------------------------------------------------------------------

end module mpasjedi_vc_model2geovars_mod

