/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <memory>
#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

#include "mpasjedi/getvalues/GetValues.interface.h"

// -------------------------------------------------------------------------------------------------

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
}

namespace util {
  class DateTime;
}

namespace mpas {
  class GeometryMPAS;
  class StateMPAS;

// -------------------------------------------------------------------------------------------------

class GetValues : public util::Printable,
  private util::ObjectCounter<GetValues> {
 public:
  static const std::string classname() {return "mpas::GetValues";}

  GetValues(const GeometryMPAS &, const ufo::Locations &,
            const eckit::Configuration &);
  virtual ~GetValues();

  void fillGeoVaLs(const StateMPAS &, const util::DateTime &,
            const util::DateTime &, ufo::GeoVaLs &) const;

 private:
  void print(std::ostream &) const;
  F90getvalues keyGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const GeometryMPAS> geom_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas

