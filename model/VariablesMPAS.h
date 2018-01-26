/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_VARIABLESMPAS_H_
#define MPAS_MODEL_VARIABLESMPAS_H_

#include <ostream>
#include <string>
#include <vector>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace mpas {

// -----------------------------------------------------------------------------

class VariablesMPAS : public util::Printable,
                    private util::ObjectCounter<VariablesMPAS> {
 public:
  static const std::string classname() {return "mpas::VariablesMPAS";}

  explicit VariablesMPAS(const oops::Variables &);
  explicit VariablesMPAS(const eckit::Configuration &);

  ~VariablesMPAS();

  VariablesMPAS(const VariablesMPAS &);

//model/FieldsMPAS.h:  int & toFortran() {return keyFlds_;}
//model/FieldsMPAS.h:  const int & toFortran() const {return keyFlds_;}
//  int * toFortran() {return &fvars_[0];}  //BJJ
  const int * toFortran() const {return &fvars_[0];}

 private:
  void print(std::ostream &) const;
  void setF90(const std::vector<std::string>);
  std::vector<int> fvars_;
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPAS_MODEL_VARIABLESMPAS_H_
