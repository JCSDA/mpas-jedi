/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "VariablesMPAS.h"

#include<vector>

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

namespace mpas {

// -----------------------------------------------------------------------------

VariablesMPAS::VariablesMPAS(const oops::Variables & oopsvars) {
  oops::Log::debug() << "VariablesMPAS oopsvar:" << oopsvars.variables() << std::endl;
  this->setF90(oopsvars.variables());
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

VariablesMPAS::VariablesMPAS(const eckit::Configuration & config) {
  oops::Log::debug() << "VariablesMPAS config:" << config << std::endl;
  std::vector<std::string> vars;
  config.get("variables", vars);
  this->setF90(vars);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

void VariablesMPAS::setF90(const std::vector<std::string> vars) {
  size_t nv = vars.size();
  oops::Log::debug() << "setF90 " << nv << " vars = " << vars << std::endl;
  fvars_.resize(nv + 2);
  fvars_[0] = nv;
  for (size_t jj = 0; jj < nv; ++jj) {
     int ii = 0;
     if (vars[jj]=="x") ii = 1;
     if (vars[jj]=="q") ii = 2;
     if (vars[jj]=="u") ii = 3;
     if (vars[jj]=="v") ii = 4;
     if (vars[jj]=="bc") ii = 5;
     if (vars[jj]=="cv") ii = 1;
     ASSERT(ii > 0);
     fvars_[jj+1] = ii;
  }
  fvars_[nv+1] = 999;  // just for checking
  oops::Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

VariablesMPAS::~VariablesMPAS() {}

// -----------------------------------------------------------------------------

VariablesMPAS::VariablesMPAS(const VariablesMPAS & other): fvars_(other.fvars_) {}

// -----------------------------------------------------------------------------

void VariablesMPAS::print(std::ostream & os) const {
  os << "mpas::VariablesMPAS: vars = " << fvars_;
}

// -----------------------------------------------------------------------------

}  // namespace mpas
