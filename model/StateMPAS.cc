/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "StateMPAS.h"

#include <algorithm>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ModelBiasMPAS.h"
#include "FieldsMPAS.h"
#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "ModelMPAS.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol, const oops::Variables & vars,
                 const util::DateTime & vt)
  : fields_(new FieldsMPAS(resol, vars, vt)), stash_()
{
  oops::Log::trace() << "StateMPAS::StateMPAS created." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol, const eckit::Configuration & file)
  : fields_(), stash_()
{
// Should get variables from file. YT
//  eckit::LocalConfiguration modelvars;
//  modelvars.set("variables", "cv");
//  oops::Variables vars(modelvars);

//  eckit::LocalConfiguration modelvars;
//  const std::vector<std::string> vv{"cv"};
//  modelvars.set("variables", vv); 
//  oops::Variables vars(vv);

/*
// Should get variables from file. YT
  const std::vector<std::string> vv{"x","bc"};
  oops::Variables vars(vv);
  fields_.reset(new FieldsQG(resol, vars, util::DateTime()));
  fields_->read(file);

  ASSERT(fields_);
  Log::trace() << "StateQG::StateQG created and read in." << std::endl;

 */
// Should get variables from file. YT

// WORKING WITHOUT READING THE NAMELIST AS QG/
  oops::Log::trace() << "StateMPAS::GD0 enforcing to variable to cv" << std::endl;
//  --- For Interface ---
//  const std::vector<std::string> vv{"theta", "rho", "index_qv", "uReconstructZonal", "uReconstructMeridional"};
//  --- For HofX ---
//  const std::vector<std::string> vv{"theta", "index_qv", "pressure_base"};
  const std::vector<std::string> vv{"theta", "rho", "index_qv", "uReconstructZonal", "uReconstructMeridional", "pressure_base"};
//  --- For Dirac ---
//  const std::vector<std::string> vv{"theta"};
  oops::Log::trace() << "StateMPAS::GD1 enforcing to variable to cv" << std::endl;
  oops::Variables vars(vv);
  oops::Log::trace() << "StateMPAS::GD2" << std::endl;


  fields_.reset(new FieldsMPAS(resol, vars, util::DateTime()));
  oops::Log::trace() << "StateMPAS::GD3 before read" << std::endl;
  fields_->read(file);
  oops::Log::trace() << "StateMPAS::GD3 after read" << std::endl;

  ASSERT(fields_);

  oops::Log::trace() << "StateMPAS::StateMPAS created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const GeometryMPAS & resol, const StateMPAS & other)
  : fields_(new FieldsMPAS(*other.fields_, resol)), stash_()
{
  ASSERT(fields_);
  oops::Log::trace() << "StateMPAS::StateMPAS created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::StateMPAS(const StateMPAS & other)
  : fields_(new FieldsMPAS(*other.fields_)), stash_()
{
  
  oops::Log::trace() << "StateMPAS::StateMPAS before copied." << std::endl;
  ASSERT(fields_);
  oops::Log::trace() << "StateMPAS::StateMPAS copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateMPAS::~StateMPAS() {
  oops::Log::trace() << "StateMPAS::StateMPAS destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateMPAS & StateMPAS::operator=(const StateMPAS & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void StateMPAS::interpolate(const ufo::Locations & locs, const oops::Variables & vars, ufo::GeoVaLs & cols) const {
  fields_->interpolate(locs, vars, cols);
}
// -----------------------------------------------------------------------------
/// Interpolate full fields
// -----------------------------------------------------------------------------
void StateMPAS::changeResolution(const StateMPAS & other) {
  fields_->changeResolution(*other.fields_);
  oops::Log::trace() << "StateMPAS interpolated" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateMPAS & StateMPAS::operator+=(const IncrementMPAS & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  fields_->add(dx.fields());
  return *this;
}
// -----------------------------------------------------------------------------
/// Convert to/from unstructured grid
// -----------------------------------------------------------------------------
void StateMPAS::convert_to(oops::UnstructuredGrid & ug) const {
  fields_->convert_to(ug);
}
// -----------------------------------------------------------------------------
void StateMPAS::convert_from(const oops::UnstructuredGrid & ug) {
  fields_->convert_from(ug);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateMPAS::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void StateMPAS::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
void StateMPAS::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateMPAS::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void StateMPAS::accumul(const double & zz, const StateMPAS & xx) {
  fields_->axpy(zz, *xx.fields_);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
