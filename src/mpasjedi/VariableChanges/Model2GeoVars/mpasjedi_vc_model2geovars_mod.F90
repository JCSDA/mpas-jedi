! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_vc_model2geovars_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log

!oops
use kinds, only : kind_real

!ufo
use gnssro_mod_transform, only: geometric2geop
use ufo_vars_mod

!MPAS-Model
use mpas_derived_types
use mpas_field_routines
use mpas_kind_types, only: RKIND
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
  type(mpas_pool_data_type), pointer :: mdata, gdata

  ! reusable fields
  type(field2DReal), pointer :: fieldr2_a, fieldr2_b

  ! reusable arrays
  real(kind=kind_real), dimension(:), pointer :: ptrr1_a, ptrr1_b
  real(kind=kind_real), dimension(:,:), pointer :: ptrr2_a, ptrr2_b
  real(kind=kind_real), dimension(:,:), allocatable :: r2_a, r2_b

  ! iteration-specific variables
  character(len=MAXVARLEN) :: geovar
  integer :: nCells, nVertLevels, nVertLevelsP1
  integer :: iVar, iCell, iLevel
  real (kind=kind_real) :: lat

  ! surface variables
  character(len=MAXVARLEN), parameter :: &
    MPASSfcClassifyNames(5) = &
      [ character(len=MAXVARLEN) :: 'ivgtyp', 'isltyp', 'landmask', 'xice', 'snowc']
  character(len=MAXVARLEN), parameter :: &
    CRTMSfcClassifyNames(7) = &
      [var_sfc_vegtyp, var_sfc_landtyp, var_sfc_soiltyp, &
       var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac]
  type(mpas_pool_type), pointer :: CRTMSfcClassifyFields => null()
  integer, dimension(:), pointer :: vegtyp, soiltyp
  real(kind=kind_real), dimension(:), pointer :: wfrac, lfrac, ifrac, sfrac

  ! air pressure on w levels
  real(kind=kind_real), allocatable :: plevels(:,:)

  ! config members
  character(len=StrKIND), pointer :: &
    config_microp_scheme, config_radt_cld_scheme
  logical, pointer :: config_microp_re

  ! convenient local variables
  mFields => xm % subFields
  nCells = geom%nCellsSolve
  nVertLevels = geom%nVertLevels
  nVertLevelsP1 = geom%nVertLevelsP1

  ! CRTM surface type and fraction fields are interdependent
  ! pre-calculate all such fields if any are requested
  if ( any(xg%has(CRTMSfcClassifyNames)) ) then

    if (.not. all(xm%has(MPASSfcClassifyNames))) then
      call abor1_ftn('mpasjedi_vc_model2geovars::changevar: xm must include MPASSfcClassifyNames!')
    end if

    call da_template_pool(geom, CRTMSfcClassifyFields, size(CRTMSfcClassifyNames), CRTMSfcClassifyNames)

    !! surface types
    ! land type
    call xm%copy_to('ivgtyp', CRTMSfcClassifyFields, var_sfc_landtyp)

    ! veg type
    ! uses ivgtyp as input
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_vegtyp, vegtyp)
    call xm%get('ivgtyp', mdata)
    do iCell = 1, nCells
      vegtyp(iCell) = convert_type_veg(mdata%i1%array(iCell))
    end do

    ! soil type
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_soiltyp, soiltyp)
    call xm%get('isltyp', mdata)
    do iCell = 1, nCells
      soiltyp(iCell) = convert_type_soil(mdata%i1%array(iCell))
    end do

    !! surface fractions
    ! land, will be adjusted later
    ! TODO: a binary landmask is a crude indicator of sub-grid land fraction
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_lfrac, lfrac)
    call xm%get('landmask', mdata) !'land-ocean mask (1=>land ; 0=>ocean)'
    lfrac(1:nCells) = real(mdata%i1%array(1:nCells), RKIND)

    ! ice
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_ifrac, ifrac)
    call xm%get('xice', mdata)     !'fractional area coverage of sea-ice'
    ifrac(1:nCells) = mdata%r1%array(1:nCells)

    ! snow
    ! TODO: Investigate. snowc varies between 0. and 1., but Registry description indicates it is binary.
    !       Similar comment to landmask.
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_sfrac, sfrac)
    call xm%get('snowc', mdata)    !'flag for snow on ground (0=>no snow; 1=>otherwise)'
    sfrac(1:nCells) = mdata%r1%array(1:nCells)

    ! water+snow+land
    call mpas_pool_get_array(CRTMSfcClassifyFields, var_sfc_wfrac, wfrac)
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
  call xm%get('pressure', ptrr2_a)
  call xm%get('surface_pressure', ptrr1_a)
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
          call xm%get('temperature', ptrr2_a)
          call xm%get('spechum', ptrr2_b)
          allocate(r2_a(1:nVertLevels, 1:nCells))
          call q_to_w( ptrr2_b(:,1:nCells), r2_a(:,1:nCells) )
          call tw_to_tv( ptrr2_a(:,1:nCells), r2_a(:,1:nCells), gdata%r2%array(:,1:nCells) )
          deallocate(r2_a)

        case ( var_mixr ) !-humidity_mixing_ratio
          call xm%get('spechum', ptrr2_a)
          call q_to_w(ptrr2_a(:,1:nCells), gdata%r2%array(:,1:nCells))
          gdata%r2%array(:,1:nCells) = gdata%r2%array(:,1:nCells) * MPAS_JEDI_THOUSAND_kr ! [kg/kg] -> [g/kg]

          ! Ensure positive-definite mixing ratios
          !  with respect to precision of crtm::CRTM_Parameters::ZERO.
          ! TODO: this should be moved to the mpas_fields%read step
          !       only other place it should be needed is add_incr (already there)
          where(gdata%r2%array(:,1:nCells) < MPAS_JEDI_ZERO_kr)
            gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr
          end where

        case ( var_prsi ) !-air_pressure_levels
          gdata%r2%array(1:nVertLevelsP1,1:nCells) = plevels(1:nVertLevelsP1,1:nCells)

        case ( var_oz ) !-mole_fraction_of_ozone_in_air :TODO: not directly available from MPAS
          !call xm%get('o3', mdata)
          gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr !mdata%r2%array(:,1:nCells)

        case ( var_co2 ) !-mole_fraction_of_carbon_dioxide_in_air :TODO: not directly available from MPAS
          !call xm%get('co2', mdata)
          gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr !mdata%r2%array(:,1:nCells)

        case ( var_clw ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
          call q_fields_forward('qc', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_cli ) !-mass_content_of_cloud_ice_in_atmosphere_layer
          call q_fields_forward('qi', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clr ) !-mass_content_of_rain_in_atmosphere_layer
          call q_fields_forward('qr', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_cls ) !-mass_content_of_snow_in_atmosphere_layer
          call q_fields_forward('qs', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clg ) !-mass_content_of_graupel_in_atmosphere_layer
          call q_fields_forward('qg', mFields, gdata%r2, plevels, nCells, nVertLevels)

        case ( var_clh ) !-mass_content_of_hail_in_atmosphere_layer
          call q_fields_forward('qh', mFields, gdata%r2, plevels, nCells, nVertLevels)

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
          call xm%get('qr', fieldr2_a)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            allocate(r2_a(1:nVertLevels, 1:nCells))
            allocate(r2_b(1:nVertLevels, 1:nCells))
            if (trim(config_microp_scheme) == 'mp_thompson') then
               call xm%get('nr', ptrr2_a) !- [nb kg^{-1}]: MPAS output for 2-moment MP scheme
               r2_b(:,1:nCells) = ptrr2_a(:,1:nCells)
            else
               r2_b = MPAS_JEDI_ONE_kr
            end if
            call xm%get('rho', fieldr2_b) !- [kg m^{-3}]: Dry air density

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
          call xm%get('qg', fieldr2_a)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
          call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
          if (config_microp_re) then
            allocate(r2_a(1:nVertLevels, 1:nCells))
            call xm%get('rho', fieldr2_b) !- [kg m^{-3}]: Dry air density
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
          allocate(r2_a(1:nVertLevels,1:nCells))
          call geometricZ_full_to_half(geom%zgrid(:,1:nCells), nCells, &
                                       nVertLevels, r2_a(:,1:nCells))
          do iCell = 1, nCells
            lat = geom%latCell(iCell) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
            do iLevel = 1, nVertLevels
               call geometric2geop(lat, r2_a(iLevel,iCell), gdata%r2%array(iLevel,iCell))
            enddo
          enddo
          deallocate(r2_a)

        case ( var_geomz ) !-height
          ! calculate midpoint geometricZ (unit: m):
          call geometricZ_full_to_half(geom%zgrid(:,1:nCells), nCells, &
                                       nVertLevels, gdata%r2%array(:,1:nCells))

!! begin surface variables
        case ( var_sfc_z ) !-surface_geopotential_height
          do iCell=1,nCells
            lat = geom%latCell(iCell) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
            call geometric2geop(lat, geom%zgrid(1,iCell), gdata%r1%array(iCell))
          enddo

        case ( var_sfc_geomz ) !-surface_altitude
          gdata%r1%array(1:nCells) = geom%zgrid(1,1:nCells)

        case ( var_sfc_sdepth ) !-surface_snow_thickness
          call xm%get('snowh', mdata)
          gdata%r1%array(1:nCells) = mdata%r1%array(1:nCells) * MPAS_JEDI_THOUSAND_kr ! [m] -> [mm]

        case ( var_sfc_vegfrac ) !-vegetation_area_fraction
          call xm%get('vegfra', mdata)
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

        case ( var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp, &
               var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac )
          ! CRTMSfcClassifyNames
          !-land_type_index, vegetation_type_index, soil_type
          !-water_area_fraction, land_area_fraction, ice_area_fraction, surface_snow_area_fraction
          call xg%copy_from(geovar, CRTMSfcClassifyFields)

        case ( var_sfc_wspeed ) !-surface_wind_speed
          call xm%get('u10', ptrr1_a)
          call xm%get('v10', ptrr1_b)
          gdata%r1%array(1:nCells)=sqrt( ptrr1_a(1:nCells)**2 + ptrr1_b(1:nCells)**2 )

!! end surface variables

        case default
          call abor1_ftn('mpasjedi_vc_model2geovars::changevar: geovar not implemented => '//trim(geovar))

      end select

    end if

  end do !iVar

  if (associated(CRTMSfcClassifyFields)) then
    call mpas_pool_destroy_pool(CRTMSfcClassifyFields)
  end if

  deallocate(plevels)

end subroutine changevar

! --------------------------------------------------------------------------------------------------

end module mpasjedi_vc_model2geovars_mod

