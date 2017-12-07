/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "IncrementMPAS.h"

#include <algorithm>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ErrorCovarianceMPAS.h"
#include "FieldsMPAS.h"
#include "GeometryMPAS.h"
#include "StateMPAS.h"
#include "oops/base/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & resol, const oops::Variables & vars,
                         const util::DateTime & vt)
  : fields_(new FieldsMPAS(resol, vars, vt)), stash_()
{
  fields_->zero();
  oops::Log::trace() << "IncrementMPAS constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const GeometryMPAS & resol, const IncrementMPAS & other)
  : fields_(new FieldsMPAS(*other.fields_, resol)), stash_()
{
  oops::Log::trace() << "IncrementMPAS constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other, const bool copy)
  : fields_(new FieldsMPAS(*other.fields_, copy)), stash_()
{
  oops::Log::trace() << "IncrementMPAS copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::IncrementMPAS(const IncrementMPAS & other)
  : fields_(new FieldsMPAS(*other.fields_)), stash_()
{
  oops::Log::trace() << "IncrementMPAS copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementMPAS::~IncrementMPAS() {
  oops::Log::trace() << "IncrementMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementMPAS::diff(const StateMPAS & x1, const StateMPAS & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  oops::Log::debug() << "IncrementMPAS:diff incr " << *fields_ << std::endl;
  oops::Log::debug() << "IncrementMPAS:diff x1 " << x1.fields() << std::endl;
  oops::Log::debug() << "IncrementMPAS:diff x2 " << x2.fields() << std::endl;
  fields_->diff(x1.fields(), x2.fields());
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator=(const IncrementMPAS & rhs) {
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator+=(const IncrementMPAS & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator-=(const IncrementMPAS & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementMPAS & IncrementMPAS::operator*=(const double & zz) {
  *fields_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void IncrementMPAS::zero(const util::DateTime & vt) {
  fields_->zero(vt);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::axpy(const double & zz, const IncrementMPAS & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::accumul(const double & zz, const StateMPAS & xx) {
  fields_->axpy(zz, xx.fields());
}
// -----------------------------------------------------------------------------
void IncrementMPAS::schur_product_with(const IncrementMPAS & dx) {
  fields_->schur_product_with(*dx.fields_);
}
// -----------------------------------------------------------------------------
double IncrementMPAS::dot_product_with(const IncrementMPAS & other) const {
  return dot_product(*fields_, *other.fields_);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::random() {
  fields_->random();
}
// -----------------------------------------------------------------------------
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void IncrementMPAS::interpolateTL(const ufo::Locations & locs, const oops::Variables & vars, ufo::GeoVaLs & cols) const {
  oops::Log::debug() << "IncrementMPAS::interpolateTL fields in" << *fields_ << std::endl;
  fields_->interpolateTL(locs, vars, cols);
  oops::Log::debug() << "IncrementMPAS::interpolateTL gom " << cols << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementMPAS::interpolateAD(const ufo::Locations & locs, const oops::Variables & vars, const ufo::GeoVaLs & cols) {
  oops::Log::debug() << "IncrementMPAS::interpolateAD gom " << cols << std::endl;
  oops::Log::debug() << "IncrementMPAS::interpolateAD fields in" << *fields_ << std::endl;
  fields_->interpolateAD(locs, vars, cols);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementMPAS::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
void IncrementMPAS::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------

}  // namespace mpas
