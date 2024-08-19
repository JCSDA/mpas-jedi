/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/State/State.h"

namespace mpas {


// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
State::State(const Geometry & geom,
                     const oops::Variables & vars,
                     const util::DateTime & time)
  : geom_(geom), vars_(vars), time_(time)
{
  oops::Log::trace() << "State::State create." << std::endl;
  mpas_state_create_f90(keyState_, geom_.toFortran(), stateVars(), vars_);
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & geom,
                     const eckit::Configuration & config)
  : geom_(geom), vars_(), time_(util::DateTime())
{
  oops::Log::trace() << "State::State create and read." << std::endl;
  StateParameters params;
  params.deserialize(config);

  // Set up vars
  vars_ = oops::Variables(params.state_variables.value());
  mpas_state_create_f90(keyState_, geom_.toFortran(), stateVars(), vars_);

  if (params.analytic_init.value() != boost::none) {
    analytic_init(params.toConfiguration());
  } else {
    // filename must be set if not analytic initial condition
    ASSERT(params.filename.value() != boost::none);
    read(params.toConfiguration());
  }

  oops::Log::trace() << "State::State created and read in."
                     << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & geom,
                     const State & other)
  : geom_(geom), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "State::State create by interpolation."
                     << std::endl;

  // create new state with geom
  mpas_state_create_f90(keyState_, geom_.toFortran(), stateVars(), vars_);

  // interpolate other to geom
  changeResolution(other);

  oops::Log::trace() << "State::State created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const oops::Variables & vars, const State & other)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
  oops::Log::trace() << "State::State create with variable change." << std::endl;

  mpas_state_create_f90(keyState_, geom_.toFortran(), stateVars(), vars_);
  mpas_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "State::State created with variablechange" << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const State & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "State::State before copied." << std::endl;

  mpas_state_create_f90(keyState_, geom_.toFortran(), stateVars(), vars_);
  mpas_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
State::~State() {
  mpas_state_delete_f90(keyState_);
  oops::Log::trace() << "State::State destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
State & State::operator=(const State & rhs) {
  mpas_state_copy_f90(keyState_, rhs.keyState_);
  time_ = rhs.time_;
  vars_ = rhs.vars_;
  return *this;
}
// -----------------------------------------------------------------------------
void State::toFieldSet(atlas::FieldSet & fset) const {
  const bool include_halo = false;
  const bool flip_vert_lev = true;
  mpas_state_to_fieldset_f90(keyState_, geom_.toFortran(), vars_, fset.get(), include_halo,
                          flip_vert_lev);
  oops::Log::trace() << "State toFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::fromFieldSet(const atlas::FieldSet & fset) {
  const bool include_halo = false;  /* always false, only fill ceter of domain */
  const bool flip_vert_lev = true;
  mpas_state_from_fieldset_f90(keyState_, geom_.toFortran(), vars_, fset.get(), include_halo,
                          flip_vert_lev);
  oops::Log::trace() << "State fromFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interpolate full state
// -----------------------------------------------------------------------------
void State::changeResolution(const State & other) {
  oops::Log::trace() << "State changed resolution starting" << std::endl;

  // If both states have same resolution, then copy instead of interpolating
  if (geom_.isEqual(other.geom_)) {
    mpas_state_copy_f90(keyState_, other.keyState_);
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
  // Use oops interpolator to handle integer/categorical fields correctly
  // Once the atlas interpolator gains support for this feature, we could make this configurable
  // from the user-facing yaml file; for now though, the atlas interpolator would be wrong for the
  // many integer fields of mpas-jedi.
  conf.set("local interpolator type", "oops unstructured grid interpolator");
  oops::GlobalInterpolator interp(conf, source_geom, target_fs, geom_.getComm());

  atlas::FieldSet source{};
  atlas::FieldSet target{};

  // Interpolate atlas::FieldSet representation of mpas data
  other.toFieldSet(source);
  interp.apply(source, target);
  this->fromFieldSet(target);

  oops::Log::trace() << "State changed resolution" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  oops::Log::trace() << "State add increment starting" << std::endl;
  ASSERT(validTime() == dx.validTime());
  // Interpolate increment to state resolution
  Increment dx_sr(geom_, dx);
  mpas_state_add_incr_f90(keyState_, dx_sr.toFortran());
  oops::Log::trace() << "State add increment done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t State::serialSize() const {
  // Field
  size_t nn;
  mpas_state_serial_size_f90(keyState_, nn);

  // Magic factor
  nn += 1;

  // Date and time
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
constexpr real_type SerializeCheckValue = -54321.98765;
void State::serialize(std::vector<real_type> & vect) const {
  // Serialize the field
  size_t nn;
  mpas_state_serial_size_f90(keyState_, nn);
  std::vector<real_type> vect_field(nn, 0);
  mpas_state_serialize_f90(keyState_, nn, vect_field.data());
  vect.insert(vect.end(), vect_field.begin(), vect_field.end());

  // Magic value placed in serialization; used to validate deserialization
  vect.push_back(SerializeCheckValue);

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void State::deserialize(const std::vector<real_type> & vect, size_t & index) {
  // Deserialize the field
  mpas_state_deserialize_f90(keyState_, vect.size(), vect.data(), index);

  // Use magic value to validate deserialization
  ASSERT(vect.at(index) == SerializeCheckValue);
  ++index;

  // Deserialize the date and time
  time_.deserialize(vect, index);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void State::analytic_init(const eckit::Configuration & config) {
  oops::Log::trace() << "State analytic init starting" << std::endl;
  mpas_state_analytic_init_f90(keyState_, config, time_);
  oops::Log::trace() << "State analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::read(const eckit::Configuration & config) {
  StateParameters params;
  params.deserialize(config);
  mpas_state_read_file_f90(keyState_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
void State::write(const eckit::Configuration & config) const {
  eckit::LocalConfiguration wconf(config);
  if (wconf.getBool("use_oops_filename", false)) {
    wconf.set("filename", wconf.getString("datadir") + "/" + wconf.getString("prefix") + ".nc");
  }
  mpas_state_write_file_f90(keyState_, wconf, time_);
}
// -----------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  int nc = 0;
  int nf = 0;
  mpas_state_sizes_f90(keyState_, nc, nf);

  os << std::endl << "  Valid time: " << validTime();
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<real_type> zstat(3*nf);
  mpas_state_gpnorm_f90(keyState_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "Fld=" << jj+1 << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2]
       << " : " << vars_[jj];
  }
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void State::zero() {
  mpas_state_zero_f90(keyState_);
}
// -----------------------------------------------------------------------------
void State::accumul(const real_type & zz, const State & xx) {
  mpas_state_axpy_f90(keyState_, zz, xx.keyState_);
}
// -----------------------------------------------------------------------------
real_type State::norm() const {
  real_type zz = 0.0;
  mpas_state_rms_f90(keyState_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
oops::Variables State::stateVars()
{
  // ---------------------------------------------------------------------------
  /// Temporary Auxilliary Variable Definitions
  // ---------------------------------------------------------------------------
  //--- TODO: nr and mpas_re_fields are still added in mpas_state_interface_mod
  //          they must be configured in the yaml and stream_list.atmosphere
  //          so that realistic values read from file.  Then we can remove this
  //          extra oops::Variables object.
  oops::Variables statevars(vars_);
  return statevars;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
