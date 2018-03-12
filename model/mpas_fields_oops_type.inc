
type :: mpas_field
  type (domain_type), pointer :: domain
  type (core_type), pointer :: corelist
  !type (dm_info), pointer :: dminfo
  type (mpas_geom), pointer :: geom                                 !< Number of unstructured grid cells
  integer :: nf                           !< Number of variables in fld
  !integer :: sum_scalar                   !< Number of variables in fld
  !integer :: sum_aero                     !< Number of variables in fld
  !integer :: ns                           !< Number of surface fields (x1d [nCells])
  character(len=22), allocatable  :: fldnames(:)      !< Variable identifiers
  type (mpas_pool_type), pointer  :: subFields
end type mpas_field
