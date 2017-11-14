/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_MPASVARIABLES_H_
#define MPAS_MODEL_MPASVARIABLES_H_

#include <ostream>
#include <string>

#include "util/Logger.h"
#include "model/Fortran.h"
#include "eckit/config/Configuration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Variables class to handle variables for MPAS model.

class Variables : public util::Printable,
              private util::ObjectCounter<Variables> {
 public:
  static const std::string classname() {return "mpas::Variables";}

  explicit Variables(const eckit::Configuration & config) {
    using oops::Log;
    Log::debug() << "Variables config:" << config << std::endl;
    const eckit::Configuration * conf = &config;
    mpas_var_create_f90(keyVar_, &conf);
  }
  explicit Variables(const int keyVar): keyVar_(keyVar) {}

  ~Variables() {mpas_var_delete_f90(keyVar_);}

  Variables(const Variables & other) {mpas_var_clone_f90(other.keyVar_, keyVar_);}

  int& toFortran() {return keyVar_;}
  const int& toFortran() const {return keyVar_;}

 private:
  void print(std::ostream &) const;
  int keyVar_;
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPAS_MODEL_MPASVARIABLES_H_
