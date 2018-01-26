/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "FieldsMPAS.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "util/Logger.h"
#include "Fortran.h"
#include "GeometryMPAS.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "VariablesMPAS.h" // Adding
// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
FieldsMPAS::FieldsMPAS(const GeometryMPAS & geom, const oops::Variables & vars,
                         const util::DateTime & time):
  geom_(new GeometryMPAS(geom)), vars_(vars), time_(time)
{
  mpas_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
}
// -----------------------------------------------------------------------------
FieldsMPAS::FieldsMPAS(const FieldsMPAS & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  mpas_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  if (copy) {
    mpas_field_copy_f90(keyFlds_, other.keyFlds_);
  } else {
    mpas_field_zero_f90(keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsMPAS::FieldsMPAS(const FieldsMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  mpas_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  mpas_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsMPAS::FieldsMPAS(const FieldsMPAS & other, const GeometryMPAS & geom)
  : geom_(new GeometryMPAS(geom)), vars_(other.vars_), time_(other.time_)
{
  mpas_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  mpas_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsMPAS::FieldsMPAS(const FieldsMPAS & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
  mpas_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  mpas_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsMPAS::~FieldsMPAS() {
  mpas_field_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsMPAS & FieldsMPAS::operator=(const FieldsMPAS & rhs) {
  mpas_field_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
FieldsMPAS & FieldsMPAS::operator+=(const FieldsMPAS & rhs) {
  mpas_field_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsMPAS & FieldsMPAS::operator-=(const FieldsMPAS & rhs) {
  mpas_field_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsMPAS & FieldsMPAS::operator*=(const double & zz) {
  mpas_field_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void FieldsMPAS::zero() {
  mpas_field_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::zero(const util::DateTime & time) {
  mpas_field_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void FieldsMPAS::axpy(const double & zz, const FieldsMPAS & rhs) {
  mpas_field_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double FieldsMPAS::dot_product_with(const FieldsMPAS & fld2) const {
  double zz;
  mpas_field_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsMPAS::schur_product_with(const FieldsMPAS & dx) {
    mpas_field_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::random() {
  mpas_field_random_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::interpolate(const ufo::Locations & locs, const oops::Variables & vars,
                              ufo::GeoVaLs & gom) const {
  const VariablesMPAS varmpas(vars);
  mpas_field_interp_tl_f90(keyFlds_, locs.toFortran(), varmpas.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::interpolateTL(const ufo::Locations & locs, const oops::Variables & vars,
                                ufo::GeoVaLs & gom) const {
  const VariablesMPAS varmpas(vars);
  mpas_field_interp_tl_f90(keyFlds_, locs.toFortran(), varmpas.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::interpolateAD(const ufo::Locations & locs, const oops::Variables & vars,
                                const ufo::GeoVaLs & gom) {
  const VariablesMPAS varmpas(vars);
  mpas_field_interp_ad_f90(keyFlds_, locs.toFortran(), varmpas.toFortran(), gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::changeResolution(const FieldsMPAS & other) {
  mpas_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::add(const FieldsMPAS & rhs) {
  mpas_field_add_incr_f90(keyFlds_, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::diff(const FieldsMPAS & x1, const FieldsMPAS & x2) {
  mpas_field_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::define(oops::UnstructuredGrid & ug) const {
  mpas_field_define_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::convert_to(oops::UnstructuredGrid & ug) const {
  mpas_field_convert_to_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::convert_from(const oops::UnstructuredGrid & ug) {
  mpas_field_convert_from_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsMPAS::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  mpas_field_read_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsMPAS::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  mpas_field_write_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double FieldsMPAS::norm() const {
  double zz = 0.0;
  mpas_field_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsMPAS::print(std::ostream & os) const {
  int nx = -1;
  int ny = -1;
  int nf = -1;
  int nb = -1;
//  mpas_field_sizes_f90(keyFlds_, nx, ny, nf, nb);
  os << std::endl << "  Resolution = " << nx << ", " << ny
     << ", Fields = " << nf << ", " << nb;
  nf += nb;
  std::vector<double> zstat(3*nf);
  mpas_field_gpnorm_f90(keyFlds_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2];
  }
}
// -----------------------------------------------------------------------------
}  // namespace mpas
