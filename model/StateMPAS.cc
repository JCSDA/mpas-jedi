/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/StateMPAS.h"

#include <algorithm>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "model/GeometryMPAS.h"
#include "model/IncrementMPAS.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Temporary Auxilliary Variable Definitions
// -----------------------------------------------------------------------------
//--- TODO: theta, rho, u should come from "State" variables in YAML
//          + var list for oops::State4D previously came from
//            cost_function/variables YAML, but needs to come from
//            cost_function/Jb/Background/state/variables
//          + vars_ initialized in ModelMPAS, TlmMPAS, and TlmIdMPAS
//            from YAML will need updating to reflect any changes
//          + Note: it may be possible to merge State and Increment back
//            into Fields class once YAML definitions used correctly
//--- TODO: all other aux variables should come from "VariableChange"
//          variables in YAML
const oops::Variables
    StateMPAS::auxvars_({ "theta", "rho", "u", "index_qv", "pressure",
        "landmask", "xice", "snowc", "skintemp", "ivgtyp", "isltyp",
        "snowh", "vegfra", "u10", "v10", "lai", "smois", "tslb", "w"});

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & geom,
                     const oops::Variables & incvars,
                     const util::DateTime & time)
  : geom_(new GeometryMPAS(geom)), vars_(incvars), time_(time)
{
  oops::Log::trace() << "StateMPAS::StateMPAS create." << std::endl;
  oops::Variables statevars(vars_);
  statevars += auxvars_;
  mpas_state_create_f90(keyState_, geom_->toFortran(), statevars, vars_);
  oops::Log::trace() << "StateMPAS::StateMPAS created." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol,
                     const oops::Variables & incvars,
                     const eckit::Configuration & config)
  : geom_(new GeometryMPAS(resol)), vars_(incvars), time_(util::DateTime())
{
  oops::Log::trace() << "StateMPAS::StateMPAS create and read." << std::endl;

  oops::Variables statevars(vars_);
  statevars += auxvars_;
  mpas_state_create_f90(keyState_, geom_->toFortran(), statevars, vars_);

  if (config.has("analytic_init")) {
    mpas_state_analytic_init_f90(keyState_, resol.toFortran(), config, time_);
  } else {
    mpas_state_read_file_f90(keyState_, config, time_);
  }

  oops::Log::trace() << "StateMPAS::StateMPAS created and read in."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol,
                     const StateMPAS & other)
  : geom_(new GeometryMPAS(resol)), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "StateMPAS::StateMPAS create by interpolation."
                     << std::endl;

  oops::Variables statevars(vars_);
  statevars += auxvars_;
  mpas_state_create_f90(keyState_, geom_->toFortran(), statevars, vars_);
  mpas_state_change_resol_f90(keyState_, other.keyState_);
  oops::Log::trace() << "StateMPAS::StateMPAS created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const StateMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "StateMPAS::StateMPAS before copied." << std::endl;

  oops::Variables statevars(vars_);
  statevars += auxvars_;
  mpas_state_create_f90(keyState_, geom_->toFortran(), statevars, vars_);
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
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate full state
// -----------------------------------------------------------------------------
void StateMPAS::changeResolution(const StateMPAS & other) {
  mpas_state_change_resol_f90(keyState_, other.keyState_);
  oops::Log::trace() << "StateMPAS changed resolution" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateMPAS & StateMPAS::operator+=(const IncrementMPAS & dx) {
  oops::Log::trace() << "StateMPAS add increment starting" << std::endl;
  ASSERT(this->validTime() == dx.validTime());
  mpas_state_add_incr_f90(keyState_, dx.toFortran());
  oops::Log::trace() << "StateMPAS add increment done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateMPAS::read(const eckit::Configuration & config) {
  mpas_state_read_file_f90(keyState_, config, time_);
}
// -----------------------------------------------------------------------------
void StateMPAS::analytic_init(const eckit::Configuration & config,
                              const GeometryMPAS & geom) {
  oops::Log::trace() << "StateMPAS analytic init starting" << std::endl;
  mpas_state_analytic_init_f90(keyState_, geom.toFortran(), config, time_);
  oops::Log::trace() << "StateMPAS analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateMPAS::write(const eckit::Configuration & config) const {
  mpas_state_write_file_f90(keyState_, config, time_);
}
// -----------------------------------------------------------------------------
void StateMPAS::print(std::ostream & os) const {
  int nc = 0;
  int nf = 0;
  os << std::endl << "  Valid time: " << validTime() << std::endl;
  mpas_state_sizes_f90(keyState_, nc, nf);
  os << std::endl << "  Resolution: nCellsGlobal = " << nc <<
     ", nFields = " << nf;
  std::vector<double> zstat(3*nf);
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

}  // namespace mpas
