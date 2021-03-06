/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & geom,
                             const oops::Variables & vars,
                             const util::DateTime & time):
  geom_(new GeometryMPAS(geom)), vars_(vars), time_(time)
{
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  mpas_increment_zero_f90(keyInc_);
  oops::Log::trace() << "IncrementMPAS constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & resol,
                             const IncrementMPAS & other)
  : geom_(new GeometryMPAS(resol)), vars_(other.vars_), time_(other.time_)
{
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  mpas_increment_change_resol_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "IncrementMPAS constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
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
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
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

  // If the states x1, x2 have a different geometry than the increment, need to
  // convert them.
  std::shared_ptr<const GeometryMPAS> stateGeom = x1.geometry();
  ASSERT(stateGeom->isEqual(*(x2.geometry())));
  if (geom_->isEqual(*stateGeom)) {
    mpas_increment_diff_incr_f90(keyInc_, x1.toFortran(), x2.toFortran());
  } else {
  // Note: this is likely a high-to-low resolution interpolation that should
  //       probably not be done with barycentric?
    StateMPAS x1_ir(*geom_, x1);
    StateMPAS x2_ir(*geom_, x2);
    mpas_increment_diff_incr_f90(keyInc_, x1_ir.toFortran(), x2_ir.toFortran());
  }
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
// ------------------------------------------------------------------------------
void IncrementMPAS::ones() {
  mpas_increment_ones_f90(keyInc_);
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
/// ATLAS
// -----------------------------------------------------------------------------
void IncrementMPAS::setAtlas(atlas::FieldSet * afieldset) const {
  mpas_increment_set_atlas_f90(toFortran(), geom_->toFortran(), vars_,
    afieldset->get());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::toAtlas(atlas::FieldSet * afieldset) const {
  mpas_increment_to_atlas_f90(toFortran(), geom_->toFortran(), vars_,
    afieldset->get());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::fromAtlas(atlas::FieldSet * afieldset) {
  mpas_increment_from_atlas_f90(toFortran(), geom_->toFortran(), vars_,
    afieldset->get());
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementMPAS::read(const eckit::Configuration & config) {
  mpas_increment_read_file_f90(keyInc_, config, time_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::write(const eckit::Configuration & config) const {
  mpas_increment_write_file_f90(keyInc_, config, time_);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t IncrementMPAS::serialSize() const {
  // Field
  size_t nn;
  mpas_increment_serial_size_f90(keyInc_, nn);

  // Magic value
  nn += 1;

  // Date and time
  nn += time_.serialSize();
  return nn;
}

// -----------------------------------------------------------------------------
constexpr double SerializeCheckValue = -54321.98765;
void IncrementMPAS::serialize(std::vector<double> & vect) const {
  // Serialize the field
  size_t nn;
  mpas_increment_serial_size_f90(keyInc_, nn);
  std::vector<double> vect_field(nn, 0.0);
  mpas_increment_serialize_f90(keyInc_, nn, vect_field.data());
  vect.insert(vect.end(), vect_field.begin(), vect_field.end());

  // Magic value placed in serialization; used to validate deserialization
  vect.push_back(SerializeCheckValue);

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::deserialize(const std::vector<double> & vect,
  size_t & index) {
  mpas_increment_deserialize_f90(keyInc_, vect.size(), vect.data(), index);

  // Use magic value to validate deserialization
  ASSERT(vect.at(index) == SerializeCheckValue);
  ++index;

  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
double IncrementMPAS::norm() const {
  double zz = 0.0;
  mpas_increment_rms_f90(keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::print(std::ostream & os) const {
  // store os fmt state
  std::ios oldState(nullptr);
  oldState.copyfmt(os);

  int nc = 0;
  int nf = 0;
  mpas_increment_sizes_f90(keyInc_, nc, nf);

  os << std::endl << "  Valid time: " << validTime();
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<double> zstat(3*nf);
  mpas_increment_gpnorm_f90(keyInc_, nf, zstat[0]);
  os << std::setprecision(9);
  os << std::scientific;
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "Fld=" << jj+1 << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2]
       << " : " << vars_[jj];
  }

  // restore os fmt state
  os.copyfmt(oldState);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::dirac(const eckit::Configuration & config) {
  mpas_increment_dirac_f90(keyInc_, config);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
