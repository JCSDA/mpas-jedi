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

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"


namespace mpas {


// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & geom,
                     const oops::Variables & vars,
                     const util::DateTime & time)
  : geom_(new GeometryMPAS(geom)), vars_(vars), time_(time)
{
  oops::Log::trace() << "StateMPAS::StateMPAS create." << std::endl;
  mpas_state_create_f90(keyState_, geom_->toFortran(), stateVars(), vars_);
  oops::Log::trace() << "StateMPAS::StateMPAS created." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & geom,
                     const StateMPASParameters & params)
  : geom_(new GeometryMPAS(geom)), time_(util::DateTime())
{
  oops::Log::trace() << "StateMPAS::StateMPAS create and read." << std::endl;

  // Set up vars
  vars_ = oops::Variables(params.state_variables.value());

  mpas_state_create_f90(keyState_, geom_->toFortran(), stateVars(), vars_);

  if (params.analytic_init.value() != boost::none) {
    analytic_init(params);
  } else {
    // filename must be set if not analytic initial condition
    ASSERT(params.filename.value() != boost::none);
    read(params);
  }

  oops::Log::trace() << "StateMPAS::StateMPAS created and read in."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & geom,
                     const StateMPAS & other)
  : geom_(new GeometryMPAS(geom)), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "StateMPAS::StateMPAS create by interpolation."
                     << std::endl;

  // create new state with geom
  mpas_state_create_f90(keyState_, geom_->toFortran(), stateVars(), vars_);

  // interpolate other to geom
  changeResolution(other);

  oops::Log::trace() << "StateMPAS::StateMPAS created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const StateMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "StateMPAS::StateMPAS before copied." << std::endl;

  mpas_state_create_f90(keyState_, geom_->toFortran(), stateVars(), vars_);
  mpas_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "StateMPAS::StateMPAS copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::~StateMPAS() {
  mpas_state_delete_f90(keyState_);
  oops::Log::trace() << "StateMPAS::StateMPAS destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateMPAS & StateMPAS::operator=(const StateMPAS & rhs) {
  mpas_state_copy_f90(keyState_, rhs.keyState_);
  time_ = rhs.time_;
  vars_ = rhs.vars_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate full state
// -----------------------------------------------------------------------------
void StateMPAS::changeResolution(const StateMPAS & other) {
  mpas_state_change_resol_f90(keyState_, other.toFortran());
  oops::Log::trace() << "StateMPAS changed resolution" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateMPAS & StateMPAS::operator+=(const IncrementMPAS & dx) {
  oops::Log::trace() << "StateMPAS add increment starting" << std::endl;
  ASSERT(validTime() == dx.validTime());
  // Interpolate increment to state resolution
  IncrementMPAS dx_sr(*geom_, dx);
  mpas_state_add_incr_f90(keyState_, dx_sr.toFortran());
  oops::Log::trace() << "StateMPAS add increment done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t StateMPAS::serialSize() const {
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
constexpr double SerializeCheckValue = -54321.98765;
void StateMPAS::serialize(std::vector<double> & vect) const {
  // Serialize the field
  size_t nn;
  mpas_state_serial_size_f90(keyState_, nn);
  std::vector<double> vect_field(nn, 0);
  mpas_state_serialize_f90(keyState_, nn, vect_field.data());
  vect.insert(vect.end(), vect_field.begin(), vect_field.end());

  // Magic value placed in serialization; used to validate deserialization
  vect.push_back(SerializeCheckValue);

  // Serialize the date and time
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void StateMPAS::deserialize(const std::vector<double> & vect, size_t & index) {
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
void StateMPAS::read(const StateMPASParameters & params) {
  mpas_state_read_file_f90(keyState_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
void StateMPAS::analytic_init(const StateMPASParameters & params) {
  oops::Log::trace() << "StateMPAS analytic init starting" << std::endl;
  mpas_state_analytic_init_f90(keyState_, params.toConfiguration(), time_);
  oops::Log::trace() << "StateMPAS analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateMPAS::write(const StateMPASWriteParameters & params) const {
  mpas_state_write_file_f90(keyState_, params.toConfiguration(), time_);
}
// -----------------------------------------------------------------------------
void StateMPAS::print(std::ostream & os) const {
  // store os fmt state
  std::ios oldState(nullptr);
  oldState.copyfmt(os);

  int nc = 0;
  int nf = 0;
  mpas_state_sizes_f90(keyState_, nc, nf);

  os << std::endl << "  Valid time: " << validTime();
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<double> zstat(3*nf);
  mpas_state_gpnorm_f90(keyState_, nf, zstat[0]);
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
/// For accumulator
// -----------------------------------------------------------------------------
void StateMPAS::zero() {
  mpas_state_zero_f90(keyState_);
}
// -----------------------------------------------------------------------------
void StateMPAS::accumul(const double & zz, const StateMPAS & xx) {
  mpas_state_axpy_f90(keyState_, zz, xx.keyState_);
}
// -----------------------------------------------------------------------------
double StateMPAS::norm() const {
  double zz = 0.0;
  mpas_state_rms_f90(keyState_, zz);
  return zz;
}
// -----------------------------------------------------------------------------

oops::Variables StateMPAS::stateVars()
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
