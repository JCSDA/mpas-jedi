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

#include "model/GeometryMPAS.h"
#include "model/GetValuesTrajMPAS.h"
#include "model/IncrementMPAS.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

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
  const eckit::Configuration * conf = &vars_.toFortran();
  mpas_state_create_f90(keyState_, geom_->toFortran(), &conf);
  oops::Log::trace() << "StateMPAS::StateMPAS created." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol,
                     const oops::Variables & vars,
                     const eckit::Configuration & file)
  : geom_(new GeometryMPAS(resol)), vars_(vars), time_(util::DateTime())
{
  oops::Log::trace() << "StateMPAS::StateMPAS create and read." << std::endl;

  const eckit::Configuration * cvars = &vars_.toFortran();
  mpas_state_create_f90(keyState_, geom_->toFortran(), &cvars);

  const eckit::Configuration * conf = &file;
  util::DateTime * dtp = &time_;

  if (file.has("analytic_init")) {
    mpas_state_analytic_init_f90(keyState_, resol.toFortran(), &conf, &dtp);
  } else {
    mpas_state_read_file_f90(keyState_, &conf, &dtp);
  }

  oops::Log::trace() << "StateMPAS::StateMPAS created and read in."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol,
                     const StateMPAS & other)
  : geom_(new GeometryMPAS(resol)), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  oops::Log::trace() << "StateMPAS::StateMPAS create by interpolation."
                     << std::endl;
  mpas_state_create_f90(keyState_, geom_->toFortran(), &conf);
  mpas_state_change_resol_f90(keyState_, other.keyState_);
  oops::Log::trace() << "StateMPAS::StateMPAS created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const StateMPAS & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  oops::Log::trace() << "StateMPAS::StateMPAS before copied." << std::endl;
  mpas_state_create_f90(keyState_, geom_->toFortran(), &conf);
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
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void StateMPAS::getValues(const ufo::Locations & locs,
                          const oops::Variables & vars,
                          ufo::GeoVaLs & gom) const {
  oops::Log::trace() << "StateMPAS::getValues starting" << std::endl;
  const eckit::Configuration * conf = &vars.toFortran();
  mpas_state_getvalues_notraj_f90(keyState_, locs.toFortran(), &conf,
                                 gom.toFortran());
  oops::Log::trace() << "StateMPAS::getValues done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateMPAS::getValues(const ufo::Locations & locs,
                          const oops::Variables & vars,
                          ufo::GeoVaLs & gom,
                          const GetValuesTrajMPAS & traj) const {
  oops::Log::trace() << "StateMPAS::getValues traj starting" << std::endl;
  const eckit::Configuration * conf = &vars.toFortran();
  mpas_state_getvalues_f90(keyState_, locs.toFortran(), &conf,
                           gom.toFortran(), traj.toFortran());
  oops::Log::trace() << "StateMPAS::getValues traj done" << std::endl;
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
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  mpas_state_read_file_f90(keyState_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void StateMPAS::analytic_init(const eckit::Configuration & config,
                              const GeometryMPAS & geom) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  oops::Log::trace() << "StateMPAS analytic init starting" << std::endl;
  mpas_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
  oops::Log::trace() << "StateMPAS analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateMPAS::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  mpas_state_write_file_f90(keyState_, &conf, &dtp);
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
