/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

//?  #include "oops/generic/InterpolatorUnstructured.h"
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
  oops::Log::trace() << "Increment::Increment (from geom, vars and time) done" << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & resol,
                             const IncrementMPAS & other)
  : geom_(new GeometryMPAS(resol)), vars_(other.vars_), time_(other.time_)
{
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  mpas_increment_change_resol_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "Increment::Increment (from geom/resol and other) done" << std::endl;
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
  oops::Log::trace() << "Increment::Increment (from other and bool copy) done" << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  mpas_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  mpas_increment_copy_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "Increment::Increment (from other) done" << std::endl;
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
  vars_ = rhs.vars_;
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
IncrementMPAS & IncrementMPAS::operator*=(const real_type & zz) {
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
void IncrementMPAS::axpy(const real_type & zz, const IncrementMPAS & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, dx.keyInc_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::axpy(const real_type & zz, const StateMPAS & xx,
                       const bool check) {
  ASSERT(!check || this->validTime() == xx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::accumul(const real_type & zz, const StateMPAS & xx) {
  mpas_increment_axpy_state_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::schur_product_with(const IncrementMPAS & dx) {
  mpas_increment_self_schur_f90(keyInc_, dx.keyInc_);
}
// -----------------------------------------------------------------------------
real_type IncrementMPAS::dot_product_with(const IncrementMPAS & other) const {
  real_type zz;
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
// Here toAtlas and fromAtlas are used in saber (B^-1 term) without halo or vertical flip
void IncrementMPAS::setAtlas(atlas::FieldSet * afieldset) const {
  oops::Log::trace() << "mpasjedi::Increment::setAtlas starting" << std::endl;
  mpas_increment_set_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get(), false);
  oops::Log::trace() << "mpasjedi::Increment::setAtlas done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IncrementMPAS::toAtlas(atlas::FieldSet * afieldset) const {
  oops::Log::trace() << "mpasjedi::Increment::toAtlas starting" << std::endl;
  mpas_increment_to_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get(), false, false);
  oops::Log::trace() << "mpasjedi::Increment::toAtlas done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IncrementMPAS::fromAtlas(atlas::FieldSet * afieldset) {
  oops::Log::trace() << "mpasjedi::Increment::fromAtlas starting" << std::endl;
  mpas_increment_from_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get(), false, false);
  oops::Log::trace() << "mpasjedi::Increment::fromAtlas done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IncrementMPAS::getFieldSet(const oops::Variables & vars, atlas::FieldSet & fset) const {
  // Quenstionable : true for halo and flip_vert_lev
  const bool include_halo = true;
  const bool flip_vert_lev = true;
  oops::Log::trace() << "mpasjedi::Increment:::getFieldSet starting" << std::endl;
  mpas_increment_set_atlas_f90(keyInc_, geom_->toFortran(), vars, fset.get(), include_halo);
  mpas_increment_to_atlas_f90(keyInc_, geom_->toFortran(), vars, fset.get(), include_halo,
                              flip_vert_lev);
  oops::Log::trace() << "mpasjedi::Increment::getFieldSet done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IncrementMPAS::getFieldSetAD(const oops::Variables & vars, const atlas::FieldSet & fset) {
  // Quenstionable : true for halo and flip_vert_lev
  const bool include_halo = true;
  const bool flip_vert_lev = true;
  oops::Log::trace() << "mpasjedi::Increment:::getFieldSetAD starting" << std::endl;
  mpas_increment_to_atlas_ad_f90(keyInc_, geom_->toFortran(), vars, fset.get(), include_halo,
                                 flip_vert_lev);
  oops::Log::trace() << "mpasjedi::Increment::getFieldSetAD done" << std::endl;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementMPAS::read(const IncrementMPASReadParameters & params) {
  mpas_increment_read_file_f90(keyInc_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::write(const IncrementMPASWriteParameters & params) const {
  mpas_increment_write_file_f90(keyInc_, params.toConfiguration(), time_);
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
constexpr real_type SerializeCheckValue = -54321.98765;
void IncrementMPAS::serialize(std::vector<real_type> & vect) const {
  // Serialize the field
  size_t nn;
  mpas_increment_serial_size_f90(keyInc_, nn);
  std::vector<real_type> vect_field(nn, 0.0);
  mpas_increment_serialize_f90(keyInc_, nn, vect_field.data());
  vect.insert(vect.end(), vect_field.begin(), vect_field.end());

  // Magic value placed in serialization; used to validate deserialization
  vect.push_back(SerializeCheckValue);

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::deserialize(const std::vector<real_type> & vect,
  size_t & index) {
  mpas_increment_deserialize_f90(keyInc_, vect.size(), vect.data(), index);

  // Use magic value to validate deserialization
  ASSERT(vect.at(index) == SerializeCheckValue);
  ++index;

  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
real_type IncrementMPAS::norm() const {
  real_type zz = 0.0;
  mpas_increment_rms_f90(keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::print(std::ostream & os) const {
  int nc = 0;
  int nf = 0;
  mpas_increment_sizes_f90(keyInc_, nc, nf);

  os << std::endl << "  Valid time: " << validTime();
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<real_type> zstat(3*nf);
  mpas_increment_gpnorm_f90(keyInc_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "Fld=" << jj+1 << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2]
       << " : " << vars_[jj];
  }
}
// -----------------------------------------------------------------------------
void IncrementMPAS::dirac(const DiracParameters & params) {
  mpas_increment_dirac_f90(keyInc_, params.toConfiguration());
}
// -----------------------------------------------------------------------------

}  // namespace mpas
