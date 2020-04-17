/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "model/GeometryMPAS.h"
#include "model/getvalues/LinearGetValues.interface.h"
#include "model/IncrementMPAS.h"
#include "model/StateMPAS.h"

// -------------------------------------------------------------------------------------------------

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace mpas {
  class StateMPAS;
  class GeometryMPAS;

// -------------------------------------------------------------------------------------------------

class LinearGetValues : public util::Printable, private util::ObjectCounter<LinearGetValues> {
 public:
  static const std::string classname() {return "mpas::LinearGetValues";}

  LinearGetValues(const GeometryMPAS &, const ufo::Locations &);
  virtual ~LinearGetValues();

  void setTrajectory(const StateMPAS & state, const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals);
  void fillGeoVaLsTL(const IncrementMPAS & inc, const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals) const;
  void fillGeoVaLsAD(IncrementMPAS & inc, const util::DateTime & t1, const util::DateTime & t2,
                     const ufo::GeoVaLs & geovals) const;

 private:
  void print(std::ostream &) const;
  F90lineargetvalues keyLinearGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const GeometryMPAS> geom_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
