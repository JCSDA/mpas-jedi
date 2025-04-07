/*
 * (C) Copyright 2017-2023 UCAR
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

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/State/State.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                             const oops::Variables & vars,
                             const util::DateTime & time):
  geom_(geom), vars_(vars), time_(time), iterLevs_(geom_.nIterLevs(vars_)),
  iterVals_(std::accumulate(iterLevs_.begin(), iterLevs_.end(), 0))
{
  mpas_increment_create_f90(keyInc_, geom_.toFortran(), vars_);
  mpas_increment_zero_f90(keyInc_);
  oops::Log::trace() << "Increment::Increment (from geom, vars and time) done" << std::endl;
}

// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & resol,
                             const Increment & other)
  : geom_(resol), vars_(other.vars_), time_(other.time_),
  iterLevs_(other.iterLevs_), iterVals_(other.iterVals_)
{
  mpas_increment_create_f90(keyInc_, geom_.toFortran(), vars_);

  // If both increments have same resolution, then copy instead of interpolating
  if (geom_.isEqual(other.geom_)) {
    mpas_increment_copy_f90(keyInc_, other.keyInc_);
    time_ = other.time_;
    return;
  }

  // Build oops interpolator -- this takes a few extra steps at this level
  const oops::GeometryData source_geom(other.geom_.functionSpace(),
                                       other.geom_.fields(),
                                       other.geom_.levelsAreTopDown(),
                                       other.geom_.getComm());
  const atlas::FunctionSpace target_fs = geom_.functionSpace();
  eckit::LocalConfiguration conf;
  // Use oops interpolator for consistency with State resolution change.
  conf.set("local interpolator type", "oops unstructured grid interpolator");
  oops::GlobalInterpolator interp(conf, source_geom, target_fs, geom_.getComm());

  atlas::FieldSet source{};
  atlas::FieldSet target{};

  // Interpolate atlas::FieldSet representation of mpas data
  other.toFieldSet(source);
  interp.apply(source, target);
  this->fromFieldSet(target);

  oops::Log::trace() << "Increment::Increment (from geom/resol and other) done" << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
  iterLevs_(other.iterLevs_), iterVals_(other.iterVals_)
{
  mpas_increment_create_f90(keyInc_, geom_.toFortran(), vars_);
  if (copy) {
    mpas_increment_copy_f90(keyInc_, other.keyInc_);
  } else {
    mpas_increment_zero_f90(keyInc_);
  }
  oops::Log::trace() << "Increment::Increment (from other and bool copy) done" << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
  iterLevs_(other.iterLevs_), iterVals_(other.iterVals_)
{
  mpas_increment_create_f90(keyInc_, geom_.toFortran(), vars_);
  mpas_increment_copy_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "Increment::Increment (from other) done" << std::endl;
}
// -----------------------------------------------------------------------------
Increment::~Increment() {
  mpas_increment_delete_f90(keyInc_);
  oops::Log::trace() << "Increment destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  oops::Log::debug() << "Increment:diff x1 " << x1.toFortran() << std::endl;
  oops::Log::debug() << "Increment:diff x2 " << x2.toFortran() << std::endl;
  // If the states x1, x2 have a different geometry than the increment, need to
  // convert them.
  const Geometry & stateGeom = x1.geometry();
  ASSERT(stateGeom.isEqual(x2.geometry()));
  if (geom_.isEqual(stateGeom)) {
    mpas_increment_diff_incr_f90(keyInc_, x1.toFortran(), x2.toFortran());
  } else {
  // Note: this is likely a high-to-low resolution interpolation that should
  //       probably not be done with barycentric?
    State x1_ir(geom_, x1);
    State x2_ir(geom_, x2);
    mpas_increment_diff_incr_f90(keyInc_, x1_ir.toFortran(), x2_ir.toFortran());
  }
}
// -----------------------------------------------------------------------------
Increment & Increment::operator=(const Increment & rhs) {
  mpas_increment_copy_f90(keyInc_, rhs.keyInc_);
  time_ = rhs.time_;
  vars_ = rhs.vars_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  mpas_increment_self_add_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  mpas_increment_self_sub_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator*=(const real_type & zz) {
  mpas_increment_self_mul_f90(keyInc_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void Increment::zero() {
  mpas_increment_zero_f90(keyInc_);
}
// -----------------------------------------------------------------------------
void Increment::zero(const util::DateTime & vt) {
  mpas_increment_zero_f90(keyInc_);
  time_ = vt;
}
// ------------------------------------------------------------------------------
void Increment::ones() {
  mpas_increment_ones_f90(keyInc_);
}
// ------------------------------------------------------------------------------
void Increment::sqrt() {
  atlas::FieldSet fset{};
  toFieldSet(fset);
  util::sqrtFieldSet(fset);
  fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
void Increment::axpy(const real_type & zz, const Increment & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, dx.keyInc_);
}
// -----------------------------------------------------------------------------
void Increment::axpy(const real_type & zz, const State & xx,
                       const bool check) {
  ASSERT(!check || this->validTime() == xx.validTime());
  mpas_increment_axpy_inc_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void Increment::accumul(const real_type & zz, const State & xx) {
  mpas_increment_axpy_state_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void Increment::schur_product_with(const Increment & dx) {
  mpas_increment_self_schur_f90(keyInc_, dx.keyInc_);
}
// -----------------------------------------------------------------------------
real_type Increment::dot_product_with(const Increment & other) const {
  real_type zz;
  mpas_increment_dot_prod_f90(keyInc_, other.keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void Increment::random() {
  mpas_increment_random_f90(keyInc_);
}
// -----------------------------------------------------------------------------
std::vector<double> Increment::rmsByLevel(const std::string & varname) const {
  throw eckit::NotImplemented("mpasjedi::Increment::rmsByLevel not implemented yet",
                              Here());
}
// -----------------------------------------------------------------------------
/// Getpoint/Setpoint
// -----------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) {
  mpas_increment_getpoint_f90(keyInc_, iter.toFortran(), iterVals_[0], iterVals_.size());
  return oops::LocalIncrement(vars_, iterVals_, iterLevs_);
}
// -------------------------------------------------------------------------------------------------
void Increment::setLocal(const oops::LocalIncrement & values, const GeometryIterator & iter) {
  const std::vector<double> vals = values.getVals();
  mpas_increment_setpoint_f90(keyInc_, iter.toFortran(), vals[0], vals.size());
}
// -----------------------------------------------------------------------------
/// ATLAS
// -----------------------------------------------------------------------------
void Increment::fromFieldSet(const atlas::FieldSet & fset) {
  const bool include_halo = false;  /* always false, only fill ceter of domain */
  const bool flip_vert_lev = true;
  oops::Log::trace() << "mpasjedi::Increment::fromFieldSet starting" << std::endl;
  mpas_increment_from_fieldset_f90(keyInc_, geom_.toFortran(), vars_, fset.get(), include_halo,
                                   flip_vert_lev);
  oops::Log::trace() << "mpasjedi::Increment::fromFieldSet done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Increment::toFieldSet(atlas::FieldSet & fset) const {
  const bool include_halo = false;
  const bool flip_vert_lev = true;
  oops::Log::trace() << "mpasjedi::Increment:::toFieldSet starting" << std::endl;
  mpas_increment_to_fieldset_f90(keyInc_, geom_.toFortran(), vars_, fset.get(), include_halo,
                              flip_vert_lev);
  oops::Log::trace() << "mpasjedi::Increment::toFieldSet done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void Increment::read(const eckit::Configuration & config) {
  ReadParameters_ params;
  params.deserialize(config);
  mpas_increment_read_file_f90(keyInc_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
void Increment::write(const eckit::Configuration & config) const {
  WriteParameters_ params;
  params.deserialize(config);
  mpas_increment_write_file_f90(keyInc_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t Increment::serialSize() const {
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
void Increment::serialize(std::vector<real_type> & vect) const {
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
void Increment::deserialize(const std::vector<real_type> & vect,
  size_t & index) {
  mpas_increment_deserialize_f90(keyInc_, vect.size(), vect.data(), index);

  // Use magic value to validate deserialization
  ASSERT(vect.at(index) == SerializeCheckValue);
  ++index;

  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
real_type Increment::norm() const {
  real_type zz = 0.0;
  mpas_increment_rms_f90(keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void Increment::print(std::ostream & os) const {
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
void Increment::dirac(const eckit::Configuration & config) {
  DiracParameters_ params;
  params.deserialize(config);
  mpas_increment_dirac_f90(keyInc_, params.toConfiguration());
}
// -----------------------------------------------------------------------------

}  // namespace mpas
