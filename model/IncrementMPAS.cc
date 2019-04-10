/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/IncrementMPAS.h"

#include <algorithm>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "model/GeometryMPAS.h"
#include "model/GetValuesTrajMPAS.h"
#include "model/StateMPAS.h"

#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & geom,
                             const oops::Variables & vars,
                             const util::DateTime & time):
  geom_(new GeometryMPAS(geom)), vars_(vars), time_(time)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  mpas_increment_zero_f90(keyInc_);
  oops::Log::trace() << "IncrementMPAS constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & resol,
                             const IncrementMPAS & other)
  : geom_(new GeometryMPAS(resol)), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  mpas_increment_change_resol_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "IncrementMPAS constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  if (copy) {
    mpas_increment_copy_f90(keyInc_, other.keyInc_);
  } else {
    mpas_increment_zero_f90(keyInc_);
  }
  oops::Log::trace() << "IncrementMPAS copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  mpas_increment_copy_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "IncrementMPAS copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::~IncrementMPAS() {
  mpas_increment_delete_f90(keyInc_);
  oops::Log::trace() << "IncrementMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementMPAS::diff(const StateMPAS & x1, const StateMPAS & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  oops::Log::debug() << "IncrementMPAS:diff x1 " << x1.toFortran() << std::endl;
  oops::Log::debug() << "IncrementMPAS:diff x2 " << x2.toFortran() << std::endl;
  mpas_increment_diff_incr_f90(keyInc_, x1.toFortran(), x2.toFortran());
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator=(const IncrementMPAS & rhs) {
  mpas_increment_copy_f90(keyInc_, rhs.keyInc_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator+=(const IncrementMPAS & dx) {
  ASSERT(this->validTime() == dx.validTime());
  mpas_increment_self_add_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator-=(const IncrementMPAS & dx) {
  ASSERT(this->validTime() == dx.validTime());
  mpas_increment_self_sub_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator*=(const double & zz) {
  mpas_increment_self_mul_f90(keyInc_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::zero() {
  mpas_increment_zero_f90(keyInc_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::zero(const util::DateTime & vt) {
  mpas_increment_zero_f90(keyInc_);
  time_ = vt;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::axpy(const double & zz, const IncrementMPAS & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, dx.keyInc_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::axpy(const double & zz, const StateMPAS & xx,
                       const bool check) {
  ASSERT(!check || this->validTime() == xx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::accumul(const double & zz, const StateMPAS & xx) {
  mpas_increment_axpy_state_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::schur_product_with(const IncrementMPAS & dx) {
  mpas_increment_self_schur_f90(keyInc_, dx.keyInc_);
}
// -----------------------------------------------------------------------------
double IncrementMPAS::dot_product_with(const IncrementMPAS & other) const {
  double zz;
  mpas_increment_dot_prod_f90(keyInc_, other.keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::random() {
  mpas_increment_random_f90(keyInc_);
}
// -----------------------------------------------------------------------------
/// Get increment values at observation locations
// -----------------------------------------------------------------------------
void IncrementMPAS::getValuesTL(const ufo::Locations & locs,
                                const oops::Variables & vars,
                                ufo::GeoVaLs & gom,
                                const GetValuesTrajMPAS & traj) const {
  oops::Log::trace() << "IncrementMPAS::getValuesTL starting" << std::endl;
  const eckit::Configuration * conf = &vars.toFortran();
  mpas_increment_getvalues_tl_f90(keyInc_, locs.toFortran(), &conf,
                              gom.toFortran(), traj.toFortran());
  oops::Log::trace() << "IncrementMPAS::getValuesTL done" << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::getValuesAD(const ufo::Locations & locs,
                                const oops::Variables & vars,
                                const ufo::GeoVaLs & gom,
                                const GetValuesTrajMPAS & traj) {
  const eckit::Configuration * conf = &vars.toFortran();
  oops::Log::trace() << "IncrementMPAS::getValuesAD starting" << std::endl;
  mpas_increment_getvalues_ad_f90(keyInc_, locs.toFortran(), &conf,
                                  gom.toFortran(), traj.toFortran());
  oops::Log::trace() << "IncrementMPAS::getValuesAD done" << std::endl;
}
// -----------------------------------------------------------------------------
/// Unstructured grid
// -----------------------------------------------------------------------------
void IncrementMPAS::ug_coord(oops::UnstructuredGrid & ug) const {
  mpas_increment_ug_coord_f90(keyInc_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::field_to_ug(oops::UnstructuredGrid & ug,
                                const int & its) const {
  mpas_increment_increment_to_ug_f90(keyInc_, ug.toFortran(), its);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::field_from_ug(const oops::UnstructuredGrid & ug,
                                const int & its) {
  mpas_increment_increment_from_ug_f90(keyInc_, ug.toFortran(), its);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementMPAS::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  mpas_increment_read_file_f90(keyInc_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  mpas_increment_write_file_f90(keyInc_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double IncrementMPAS::norm() const {
  double zz = 0.0;
  mpas_increment_rms_f90(keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::print(std::ostream & os) const {
  int nc = 0;
  int nf = 0;
  os << std::endl << "  Valid time: " << validTime() << std::endl;
  mpas_increment_sizes_f90(keyInc_, nc, nf);
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<double> zstat(3*nf);
  mpas_increment_gpnorm_f90(keyInc_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "Fld=" << jj+1 << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2]
       << " : " << vars_[jj];
  }
}
// -----------------------------------------------------------------------------
void IncrementMPAS::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  mpas_increment_dirac_f90(keyInc_, &conf);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
